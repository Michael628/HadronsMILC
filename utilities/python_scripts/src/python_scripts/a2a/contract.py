#! /usr/bin/env python3
import os
import sys
import logging
import itertools
import opt_einsum as oe
import h5py
import pickle
import typing as t

import pandas as pd

try:
    import cupy as xp
except ImportError:
    import numpy as xp

from time import perf_counter
from dataclasses import dataclass, field
from sympy.utilities.iterables import multiset_permutations
from mpi4py import MPI

import python_scripts
from python_scripts import utils
from python_scripts.a2a import config

COMM = MPI.COMM_WORLD


def convert_to_numpy(corr: xp.ndarray):
    """Converts a cupy array to a numpy array"""
    if hasattr(xp, 'asnumpy'):
        xp.asnumpy(corr)
    else:
        return corr

def make_contraction_key(contraction: t.Tuple[str], diagram_config: config.DiagramConfig):
    con_key = "".join(contraction)

    if diagram_config.has_high:
        con_key = con_key.replace(diagram_config.high_label, "")

    if diagram_config.has_low:
        con_key = con_key.replace(diagram_config.low_label, "e")

    return con_key

def time_average(cij: xp.ndarray, open_indices: t.Tuple = (0, -1)) -> xp.ndarray:
    """Takes an array with dim >= 2 and returns an array of (dim-1) where the
    i-th element in the last axis of the array is the sum of all input elements
    separated by i in the axes specified by `open_indices`

    Parameters
    ----------
    cij : ndarray
        A dim >= 2 array

    open_indices : tuple, optional
        axis indices to average over. Defaults to first and last index.
        Assumes both indices have same length

    Returns
    -------
    ndarray
        Array matching `cij` after averaging. Last dimension of output is
        average over separations in `open_indices` axes
        using periodic boundary.
    """

    cij = xp.asarray(cij)  # Remain cp/np agnostic for utility functions

    nt = cij.shape[open_indices[0]]
    dim = len(cij.shape)

    ones = xp.ones(cij.shape)
    t_range = xp.array(range(nt))

    t_start = [None] * dim
    t_start[open_indices[0]] = slice(None)
    t_start = tuple(t_start)

    t_end = [None] * dim
    t_end[open_indices[1]] = slice(None)
    t_end = tuple(t_end)

    t_mask = xp.mod(t_range[t_end] * ones
                      - t_range[t_start] * ones, xp.array([nt]))

    time_removed_indices = tuple(
        slice(None) if t_start[i] == t_end[i] else 0 for i in range(dim))

    corr = xp.zeros(cij[time_removed_indices].shape
                      + (nt,), dtype=xp.complex128)

    corr[:] = xp.array([cij[t_mask == t].sum() for t in range(nt)])

    return convert_to_numpy(corr / nt)


@dataclass
class MesonLoader:
    """Iterable object that loads meson fields for processing by a Contractor

    Parameters
    ----------
    mesonfiles : str
        File location of meson field to load.
    times : iterable
        An iterable object containing slices. The slices will be used to load
        the corresponding time slices from the meson hdf5 file in the order
        they appear in `times`
    shift_mass : bool, optional
        If true, it will be assumed that the meson field being loaded is
        constructed from eigenvectors with eigenvalues based on `oldmass`.
        The eigenvalues will then be replaced by corresponding values
        using `newmass` via the `meson_mass_alter` method.
        evalfile : str, optional
            File location of hdf5 file containing list of eigenvalues.
            Required if `shift_mass` is True.
    oldmass : str, optional
        Mass for evals in `evalfile`. Required if `shift_mass` is True.
    newmass : str, optional
        New mass to use with `evalfile`. Required if `shift_mass` is True.
    """
    mesonfiles: t.List[str]
    times: t.Tuple
    shift_mass: bool = False
    evalfile: str = ""
    oldmass: float = 0
    newmass: float = 0
    vmax_index: int = field(default_factory=lambda: slice(None))
    wmax_index: int = field(default_factory=lambda: slice(None))
    milc_mass: bool = True

    def meson_mass_alter(self, mat: xp.ndarray):
        """Shifts mass of `mat` according to `newmass` and `oldmass`
        properties. adjusts for MILC mass convention (factor of 2) if
        `milc_mass` == True.
        """

        with h5py.File(self.evalfile, 'r') as f:
            try:
                evals = xp.array(
                    f['/evals'][()].view(xp.float64), dtype=xp.float64)
            except KeyError:
                evals = xp.array(f['/EigenValueFile/evals']
                                   [()].view(xp.float64), dtype=xp.float64)

        evals = xp.sqrt(evals)

        mult_factor = 2. if self.milc_mass else 1.
        eval_scaling = xp.zeros((len(evals), 2), dtype=xp.complex128)
        eval_scaling[:, 0] = xp.divide(mult_factor * self.oldmass + 1.j * evals,
                                         mult_factor * self.newmass + 1.j * evals)
        eval_scaling[:, 1] = xp.conjugate(eval_scaling[:, 0])
        eval_scaling = eval_scaling.reshape((-1,))

        mat[:] = xp.multiply(
            mat, eval_scaling[xp.newaxis, xp.newaxis, :])

    def load_meson(self, file, time: slice = slice(None)):
        """Reads 3-dim array from hdf5 file.

        Parameters
        ----------
        time : slice, optional
            Designates time slice range to read from hdf5 file

        Returns
        -------
        ndarray
            The requested array from the hdf5 file

        Notes
        -----
        Assumes array is single precision complex. Promotes to double precision
        """
        t1 = perf_counter()

        with h5py.File(file, "r") as f:
            a_group_key = list(f.keys())[0]

            temp = f[a_group_key]['a2aMatrix']
            temp = xp.array(temp[
                                  time, self.wmax_index, self.vmax_index].view(xp.complex64),
                              dtype=xp.complex128)

        t2 = perf_counter()
        logging.debug(f"Loaded array {temp.shape} in {t2 - t1} sec")

        if self.shift_mass:
            fact = "2*" if self.milc_mass else ""
            logging.info(
                (f"Shifting mass from {fact}{self.oldmass:f} "
                 f"to {fact}{self.newmass:f}")
            )
            self.meson_mass_alter(temp)

        return temp

    def __iter__(self):
        self.mesonlist = [None for _ in range(len(self.mesonfiles))]
        self.iter_count = -1
        return self

    def __next__(self):
        """Each iteration returns

        """
        if self.iter_count < len(self.times[0]) - 1:

            self.iter_count += 1

            current_times = [self.times[i][self.iter_count]
                             for i in range(len(self.mesonfiles))]

            for i, (time, file) in \
                    enumerate(zip(current_times, self.mesonfiles)):
                try:
                    # Check if meson exists from last iter
                    if self.mesonlist[i] is not None:
                        if self.mesonlist[i][0] == time:
                            continue

                    # Check for matching time slice
                    matches = [j for j in range(
                        len(current_times[:i])) if current_times[j] == time]

                    # Check for matching file names
                    j = self.mesonfiles.index(file)

                    # Check that file matches desired time
                    if j not in matches:
                        raise ValueError

                    logging.debug(f"Found {time} at index {j}")
                    self.mesonlist[i] = (
                        time, self.mesonlist[j][1])  # Copy reference

                except ValueError:
                    logging.debug(f"Loading {time} from {file}")
                    self.mesonlist[i] = (time, self.load_meson(
                        file, time))  # Load new file

            return tuple(self.mesonlist)
        else:
            self.mesonlist = None
            raise StopIteration


def contract(m1: xp.ndarray, m2: xp.ndarray, m3: xp.ndarray = None,
             m4: xp.ndarray = None, open_indices: t.Tuple = (0, -1)):
    """Performs contraction of up to 4 3-dim arrays down to one 2-dim array

    Parameters
    ----------
    m1 : ndarray
    m2 : ndarray
    m3 : ndarray, optional
    m4 : ndarray, optional
    open_indices : tuple, default=(0,-1)
        A two-dimensional tuple containing the time indices that will
        not be contracted in the full product. The default=(0,-1) leaves
        indices of the first and last matrix open, summing over all others.

    Returns
    -------
    ndarray
        The resultant 2-dim array from the contraction
    """
    npoint = sum([1 for m in [m1, m2, m3, m4] if m is not None])

    if len(open_indices) > npoint:
        raise ValueError(
            (f"Length of open_indices must be "
             "<= diagram degree (maximum: 4)")
        )

    index_list = ['i', 'j', 'k', 'l'][:npoint]
    out_indices = "".join(index_list[i] for i in open_indices)

    if npoint == 2:  # two-point contractions
        cij = oe.contract(f'imn,jnm->{out_indices}', m1, m2)

    elif npoint == 3:  # three-point contractions
        cij = oe.contract(f'imn,jno,kom->{out_indices}', m1, m2, m3)

    else:  # four-point contractions
        cij = oe.contract(
            f'imn,jno,kop,lpm->{out_indices}', m1, m2, m3, m4)

    return cij


def execute(contraction: t.Tuple[str], diagram_config: config.DiagramConfig, run_config: config.RunContractConfig):
    def generate_time_sets():
        """Breaks meson field time extent into `comm_size` blocks and
        returns unique list of blocks for each `rank`.

        Returns
        -------
        tuple
            A tuple of length `diagram_config.npoint` containing lists of slices.
            One list for each meson field
        """

        workers = run_config.comm_size

        slice_indices = list(itertools.product(
            range(workers), repeat=diagram_config.npoint))

        if diagram_config.symmetric:  # filter for only upper-triangular slices
            slice_indices = list(
                filter(lambda x: list(x) == sorted(x), slice_indices))
            workers = int(
                (len(slice_indices) + workers - 1) / workers)

        offset = int(run_config.rank * workers)

        slice_indices = list(zip(*slice_indices[offset:offset + workers]))

        tspacing = int(run_config.time / run_config.comm_size)

        return tuple([slice(int(ti * tspacing),
                            int((ti + 1) * tspacing)) for ti in times]
                     for times in slice_indices)


    def conn_2pt():
        corr = {}

        times = generate_time_sets()

        for gamma in diagram_config.gammas:

            mesonfiles = tuple(diagram_config.mesonfile(w=contraction[i],
                                                        v=contraction[i + 1],
                                                        gamma=gamma) for i in [0, 2])

            mat_gen = MesonLoader(mesonfiles=mesonfiles,
                                  times=times, **diagram_config.meson_params)

            cij = xp.zeros((run_config.time, run_config.time), dtype=xp.complex128)

            for (t1, m1), (t2, m2) in mat_gen:

                logging.info(f"Contracting {gamma}: {t1},{t2}")

                cij[t1, t2] = contract(m1, m2)
                if diagram_config.symmetric and t1 != t2:
                    cij[t2, t1] = cij[t1, t2].T

            logging.debug("Contraction completed")

            if run_config.comm_size > 1:
                temp = None
                if run_config.rank == 0:
                    temp = xp.empty_like(cij)
                COMM.Barrier()
                COMM.Reduce(cij, temp, op=MPI.SUM, root=0)

                if run_config.rank == 0:
                    corr[gamma] = convert_to_numpy(time_average(temp))
            else:
                corr[gamma] = convert_to_numpy(time_average(cij))

            del m1, m2
        return corr


    def sib_conn_3pt():
        corr = {}

        times = generate_time_sets()

        for gamma in diagram_config.gammas:

            mesonfiles = tuple(
                diagram_config.mesonfile(w=contraction[i],
                                         v=contraction[i + 1],
                                         gamma=g)
                for i, g in zip([0, 2, 4], [gamma, "G1_G1", gamma])
            )

            mat_gen = MesonLoader(mesonfiles=mesonfiles,
                                  times=times, **diagram_config.meson_params)

            cij = xp.zeros((run_config.time, run_config.time, run_config.time), dtype=xp.complex128)

            for (t1, m1), (t2, m2), (t3, m3) in mat_gen:

                logging.info(f"Contracting {gamma}: {t1},{t2},{t3}")
                cij[t1, t2, t3] = contract(
                    m1, m2, m3, open_indices=[0, 1, 2])

                if diagram_config.symmetric:
                    raise Exception(
                        "Symmetric 3dim optimization not implemented.")

            logging.debug("Contraction completed.")

            if run_config.comm_size > 1:
                temp = None
                if run_config.rank == 0:
                    temp = xp.empty_like(cij)
                COMM.Barrier()
                COMM.Reduce(cij, temp, op=MPI.SUM, root=0)

                if run_config.rank == 0:
                    corr[gamma] = convert_to_numpy(temp)  # time_average(temp)
            else:
                corr[gamma] = convert_to_numpy(cij)  # time_average(cij)

        return corr


    def qed_conn_4pt(subdiagram: config.Diagrams) -> pd.DataFrame:
        corr = pd.DataFrame()

        times = generate_time_sets()

        for gamma in diagram_config.gammas:

            for i in range(diagram_config.n_em):
                emlabel = f"{diagram_config.emseedstring}_{i}"
                if subdiagram == config.Diagrams.photex:
                    ops = [gamma, emlabel, gamma, emlabel]
                elif subdiagram == config.Diagrams.selfen:
                    ops = [gamma, emlabel, emlabel, gamma]
                else:
                    raise ValueError("Invalid qed diagram.")
                mesonfiles = tuple(
                    diagram_config.mesonfile(w=contraction[i],
                                             v=contraction[i + 1],
                                             gamma=g)
                    for i, g in zip([0, 2, 4, 6], ops)
                )

                mat_gen = MesonLoader(mesonfiles=mesonfiles,
                                      times=times, **diagram_config.meson_params)

                cij = xp.zeros((run_config.time,)*4, dtype=xp.complex128)

                for (t1, m1), (t2, m2), (t3, m3), (t4, m4) in mat_gen:

                    logging.info(f"Contracting ({gamma},{emlabel}): {t1}, {t2}, {t3}, {t4}")
                    cij[t1, t2, t3, t4] = contract(
                        m1, m2, m3, m4, open_indices=[0, 1, 2, 3])

                    if diagram_config.symmetric:
                        raise Exception(
                            "Symmetric 4dim optimization not implemented.")

                logging.debug("Contraction completed.")

                if run_config.comm_size > 1:
                    temp = None
                    if run_config.rank == 0:
                        temp = xp.empty_like(cij)
                    COMM.Barrier()
                    COMM.Reduce(cij, temp, op=MPI.SUM, root=0)

                    if run_config.rank == 0:
                        corr[gamma] = convert_to_numpy(temp)  # time_average(temp)
                else:
                    corr[gamma] = convert_to_numpy(cij)  # time_average(cij)

        return corr


    # def qed_conn_4pt():
    #     cij = xp.zeros((run_config.time,)*4,
    #                    dtype=xp.complex128)
    #
    #     times = generate_time_sets()
    #
    #     seedkey = make_contraction_key(contraction, diagram_config)
    #
    #     gammas = [diagram_config.mesonKey.format(gamma=g) for g in ['X', 'Y', 'Z']]
    #
    #     corr = dict(
    #         zip(diagram_config.subdiagrams, [
    #             {seedkey: dict(zip(gammas, [{}] * len(gammas)))}
    #         ] * len(diagram_config.subdiagrams)))
    #
    #     matg1 = [diagram_config.mesonfile(w=contraction[0],
    #                                       v=contraction[1],
    #                                       gammam=gamma)
    #              for gamma in gammas]
    #     matg2_photex = []
    #     matg2_selfen = []
    #
    #     for gamma in gammas:
    #         if "photex" in diagram_config.subdiagrams:
    #             matg2_photex.append(diagram_config.mesonfile(w=contraction[4],
    #                                                          v=contraction[5],
    #                                                          gammam=gamma))
    #         if "selfen" in diagram_config.subdiagrams:
    #             matg2_selfen.append(diagram_config.mesonfile(w=contraction[6],
    #                                                          v=contraction[7],
    #                                                          gammam=gamma))
    #
    #     for i in range(diagram_config.n_em):
    #
    #         emlabel = f"{diagram_config.emseedstring}_{i}"
    #
    #         matp1 = diagram_config.mesonfile(w=contraction[2],
    #                                          v=contraction[3],
    #                                          gammam=emlabel)
    #         matp2_selfen = None
    #         matp2_photex = None
    #         if "selfen" in diagram_config.subdiagrams:
    #             matp2_selfen = diagram_config.mesonfile(w=contraction[4],
    #                                                     v=contraction[5],
    #                                                     gammam=emlabel)
    #         if "photex" in diagram_config.subdiagrams:
    #             matp2_photex = diagram_config.mesonfile(w=contraction[6],
    #                                                     v=contraction[7],
    #                                                     gammam=emlabel)
    #
    #         for j, gamma in enumerate(gammas):
    #                 if "photex" in diagram_config.subdiagrams:
    #                     cij[:] = contract(
    #                         m1=matg1[j],
    #                         m2=matp1,
    #                         m3=matg2_photex[j],
    #                         m4=matp2_photex,
    #                         open_indices=(0, 2)
    #                     )
    #                     corr['photex'][seedkey][gamma][emlabel] = time_average(cij)
    #                 if "selfen" in diagram_config.subdiagrams:
    #                     cij[:] = contract(
    #                         m1=matg1[j],
    #                         m2=matp1,
    #                         m3=matp2_selfen,
    #                         m4=matg2_selfen[j]
    #                     )
    #                     corr['selfen'][seedkey][gamma][emlabel] = time_average(cij)
    #
    #     return corr

    if hasattr(xp, 'cuda'):
        my_device = run_config.rank % xp.cuda.runtime.getDeviceCount()
        logging.debug(f"Rank {run_config.rank} is using gpu device {my_device}")
        xp.cuda.Device(my_device).use()

    logging.info(f"Processing mode: {', '.join(contraction)}")

    contraction_types = {
        'conn_2pt': conn_2pt,
        'sib_conn_3pt': sib_conn_3pt,
        'qed_conn_photex_4pt': lambda: qed_conn_4pt(config.Diagrams.photex),
        'qed_conn_selfen_4pt': lambda: qed_conn_4pt(config.Diagrams.selfen)
    }

    if diagram_config.contraction_type in contraction_types:
        run = contraction_types[diagram_config.contraction_type]
    else:
        raise ValueError(f"No contraction implementation for `{diagram_config.contraction_type}`.")

    return run()

def main(param_file: str):

    params = utils.load_param(param_file)

    run_config = config.get_contract_config(params)

    python_scripts.setup()

    if run_config.logging_level:
        logging.getLogger().setLevel(run_config.logging_level)

    if run_config.hardware == 'cpu':
        import numpy as xp
        globals()['xp'] = xp

    overwrite = run_config.overwrite_correlators

    diagrams = run_config.diagrams
    for diagram_config in diagrams:

        nmesons = diagram_config.npoint

        low_min = 0 if diagram_config.has_high else nmesons
        low_max = nmesons + 1 if diagram_config.has_low else 1

        perms = sum([
            list(multiset_permutations(['L'] * nlow + ['H'] * (nmesons - nlow)))
            for nlow in range(low_min, low_max)
        ], [])
        perms = list(map("".join, perms))
        # Overwrite chosen permutations with user input, if provided
        if diagram_config.perms:
            perms = diagram_config.perms

        logging.debug(f"Computing permutations: {perms}")

        for perm in perms:
            nlow = perm.count('L')

            permkey = "".join(
                sum(((perm[i], perm[(i + 1) % nmesons])
                     for i in range(nmesons)), ())
            )

            if diagram_config.has_high:
                # Build list of high source indices,
                # e.g. [[0,1], [0,2], ...]
                seeds = list(map(list, itertools.combinations(
                    list(range(diagram_config.high_count)), nmesons - nlow)))
            else:
                seeds = [[]]

            # Fill low-mode indices with None
            # e.g. [[None,0,1], [None,0,2], ...]
            _ = [
                seed.insert(i, None)
                for i in range(len(perm))
                if perm[i] == 'L'
                for seed in seeds
            ]

            # Double indices for <bra | ket> and cycle
            # e.g. [[None,0,0,1,1,None], [None,0,0,2,2,None], ...]
            seeds = [list(sum(zip(seed, seed), ())) for seed in seeds]
            seeds = [seed[1:] + seed[:1] for seed in seeds]

            outfile = diagram_config.outfile(permkey=permkey)

            if overwrite or not os.path.exists(outfile):
                logging.info(f'Contracting diagram: {diagram_config.diagram_label} ({permkey})')
            else:
                logging.info(f'Skipping write. File exists: {outfile}')
                continue

            contraction_list = [
                [
                    diagram_config.low_label
                    if seed[i] is None else
                    ('w' if i % 2 == 0 else 'v')
                    + diagram_config.high_label + s
                    for i, s in enumerate(map(str, seed))
                ]
                for seed in seeds
            ]

            start_time = perf_counter()

            corr = dict(zip(
                map(lambda x: make_contraction_key(x, diagram_config), contraction_list),
                map(lambda x: execute(x, diagram_config, run_config), contraction_list)
            ))

            stop_time = perf_counter()

            logging.debug('')
            logging.debug('    Total elapsed time for %s = %g seconds.' % (
                permkey, stop_time - start_time))
            logging.debug('')

            if run_config.rank < 1:
                os.makedirs(os.path.dirname(outfile), exist_ok=True)
                pickle.dump(corr, open(outfile, 'wb'))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise ValueError("Must provide input yaml file.")

    main(sys.argv[1])
