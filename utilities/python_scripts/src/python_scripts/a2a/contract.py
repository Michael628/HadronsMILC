#! /usr/bin/env python3
import os
import sys
import logging
import itertools
import numpy as np
import cupy as cp
import opt_einsum as oe
import h5py
import pickle
import re

from time import perf_counter
from dataclasses import dataclass
from sympy.utilities.iterables import multiset_permutations

from mpi4py import MPI
# from multiprocessing import Process, Pool, Lock, Manager

from python_scripts.processing.format import FilestemFormatBase as FSFormat
import python_scripts.utils as utils

cpnp = cp


def convert_to_numpy(corr):
    """Converts a cupy array to a numpy array"""
    return corr if type(corr) is np.ndarray else cp.asnumpy(corr)


def time_average(cij, open_indices=(0, -1)):
    """Takes an n,n array and returns a 1-dim array of length n, where the i-th
    output element is the sum of all input elements with indices separated by i

    Parameters
    ----------
    cij: ndarray
        A 2-dim square array

    Returns
    -------
    ndarray
        A 1-dim array with entries that are the average of input entries with
        equal index separations
    """

    cij = cpnp.asarray(cij)  # Remain cp/np agnostic for utility functions

    nt = cij.shape[open_indices[0]]
    dim = len(cij.shape)

    ones = cpnp.ones(cij.shape)
    t_range = cpnp.array(range(nt))

    t_start = [None] * dim
    t_start[open_indices[0]] = slice(None)
    t_start = tuple(t_start)

    t_end = [None] * dim
    t_end[open_indices[1]] = slice(None)
    t_end = tuple(t_end)

    t_mask = cpnp.mod(t_range[t_end] * ones -
                      t_range[t_start] * ones, cpnp.array([nt]))

    time_removed_indices = tuple(
        slice(None) if t_start[i] == t_end[i] else 0 for i in range(dim))

    corr = cpnp.zeros(cij[time_removed_indices].shape +
                      (nt,), dtype=np.complex128)

    # t1 = perf_counter()
    corr[:] = cpnp.array([cij[t_mask == t].sum() for t in range(nt)])
    # t2 = perf_counter()

    # print(f"claculation: {t2-t1}")

    return convert_to_numpy(corr / nt)


@dataclass
class MesonLoader:
    """Iterable object that loads meson fields for processing by a Contractor
    Parameters
    ----------
    file : str
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
    mesonfiles: list[str]
    times: tuple
    shift_mass: bool = False
    evalfile: str = ""
    oldmass: float = 0
    newmass: float = 0
    vmax_index: int = slice(None)
    wmax_index: int = slice(None)
    milc_mass: bool = True

    def meson_mass_alter(self, mat: cpnp.ndarray):

        with h5py.File(self.evalfile, 'r') as f:
            try:
                evals = cpnp.array(
                    f['/evals'][()].view(cpnp.float64), dtype=np.float64)
            except KeyError:
                evals = cpnp.array(f['/EigenValueFile/evals']
                                   [()].view(cpnp.float64), dtype=cpnp.float64)

        evals = cpnp.sqrt(evals)

        mult_factor = 2. if self.milc_mass else 1.
        eval_scaling = cpnp.zeros((len(evals), 2), dtype=cpnp.complex128)
        eval_scaling[:, 0] = cpnp.divide(mult_factor * self.oldmass + 1.j * evals,
                                         mult_factor * self.newmass + 1.j * evals)
        eval_scaling[:, 1] = cpnp.conjugate(eval_scaling[:, 0])
        eval_scaling = eval_scaling.reshape((-1,))

        mat[:] = cpnp.multiply(
            mat, eval_scaling[cpnp.newaxis, cpnp.newaxis, :])

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
            temp = cpnp.array(temp[
                time, self.wmax_index, self.vmax_index].view(np.complex64),
                dtype=np.complex128)

        t2 = perf_counter()
        logging.debug(f"Loaded array {temp.shape} in {t2-t1} sec")

        if self.shift_mass:
            fact = "2*" if self.milc_mass else ""
            logging.info(
                (f"Shifting mass from {fact}{self.oldmass:f} "
                    "to {fact}{self.newmass:f}")
            )
            self.meson_mass_alter(temp)

        return temp

    def __iter__(self):
        self.mesonlist = [None for _ in range(len(self.mesonfiles))]
        self.iter_count = -1
        return self

    def __next__(self):
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


class Contractor:
    """Performs contractions of MesonField
    """

    def __init__(self, series: str, cfg: str, **kwargs):
        self.low_max = None
        self.symmetric = False

        self.__dict__.update(kwargs)

        if 'comm' in kwargs:
            self.rank = self.comm.Get_rank()
            self.comm_size = self.comm.Get_size()

        self.series = series
        self.cfg = cfg
        self.has_high = False

        if 'high_label' in self.__dict__ and 'high_count' in self.__dict__:
            self.has_high = True

        if 'low_label' in self.__dict__:
            self.has_low = True

        self.meson_params = {
            "wmax_index": slice(self.low_max),
            "vmax_index": slice(self.low_max),
            "milc_mass": True
        }

        if 'newmass' in self.__dict__:
            self.meson_params['shift_mass'] = True
            self.meson_params['oldmass'] = float(f"0.{self.mass}")
            self.meson_params['newmass'] = float(f"0.{self.newmass}")
            self.meson_params['evalfile'] = self.evalfilestem.format(
                series=self.series, cfg=self.cfg)

        self.set_contraction_type()

    def set_contraction_type(self):
        self.run_contract = getattr(self, self.contraction_type)

        self.npoint = int(
            re.match(".*_([0-9])pt.*", self.contraction_type).group(1))

    def contract(self, m1: np.ndarray, m2: np.ndarray, m3: np.ndarray = None,
                 m4: np.ndarray = None, open_indices: tuple = (0, -1)):
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
        if len(open_indices) > self.npoint:
            raise Exception(
                (f"Length of open_indices must be "
                 "<= diagram degree ({self.npoint})")
            )

        index_list = ['i', 'j', 'k', 'l'][:self.npoint]
        out_indices = "".join(index_list[i] for i in open_indices)

        if self.npoint == 2:   # two-point contractions
            cij = oe.contract(f'imn,jnm->{out_indices}', m1, m2)

        elif self.npoint == 3:  # three-point contractions
            cij = oe.contract(f'imn,jno,kom->{out_indices}', m1, m2, m3)

        elif self.npoint == 4:  # four-point contractions
            cij = oe.contract(
                f'imn,jno,kop,lpm->{out_indices}', m1, m2, m3, m4)
        else:
            raise Exception("Expecting 2 <= self.npoint <= 4")

        return cij

    def generate_time_sets(self, symmetric: bool = False):
        """Breaks meson field time extent into `comm_size` blocks and
        returns unique list of blocks for each `rank`.

        Parameters
        ----------
        symmetric : bool, optional
            If true, time sets generated are upper triangular

        Returns
        -------
        tuple
            A tuple of length `npoint` containing lists of slices.
            One list for each meson field
        """

        workload = self.comm_size

        slice_indices = list(itertools.product(
            range(self.comm_size), repeat=self.npoint))

        if symmetric:  # filter for only upper-triangular slices
            slice_indices = list(
                filter(lambda x: list(x) == sorted(x), slice_indices))
            workload = int(
                (len(slice_indices) + self.comm_size - 1) / self.comm_size)

        offset = int(self.rank * workload)

        slice_indices = list(zip(*slice_indices[offset:offset + workload]))

        tspacing = int(self.nt / self.comm_size)

        return tuple([slice(int(ti * tspacing),
                            int((ti + 1) * tspacing)) for ti in times]
                     for times in slice_indices)

    def conn_2pt(self, contraction):

        corr = {}

        times = self.generate_time_sets(self.symmetric)

        mesonfile_replacements = FSFormat.formatdict(
            self.mesonfile, **self.__dict__)

        for gamma in self.gammas:

            mesonfiles = tuple(self.mesonfile.format(
                w=contraction[i],
                v=contraction[i + 1],
                gamma=gamma,
                **mesonfile_replacements) for i in [0, 2])

            mat_gen = MesonLoader(mesonfiles=mesonfiles,
                                  times=times, **self.meson_params)

            cij = cpnp.zeros((self.nt, self.nt), dtype=np.complex128)

            for (t1, m1), (t2, m2) in mat_gen:

                logging.info(f"Contracting {gamma}: {t1},{t2}")

                cij[t1, t2] = self.contract(m1, m2)
                if self.symmetric and t1 != t2:
                    cij[t2, t1] = cij[t1, t2].T

            logging.debug("Contraction completed")

            if 'comm' in self.__dict__ and self.comm_size > 1:
                temp = None
                if self.rank == 0:
                    temp = cpnp.empty_like(cij)
                self.comm.Barrier()
                self.comm.Reduce(cij, temp, op=MPI.SUM, root=0)

                if self.rank == 0:
                    corr[gamma] = convert_to_numpy(time_average(temp))
            else:
                corr[gamma] = convert_to_numpy(time_average(cij))

            del m1, m2
        return corr

    def sib_conn_3pt(self, contraction):

        corr = {}

        times = self.generate_time_sets(self.symmetric)

        mesonfile_replacements = FSFormat.formatdict(self.mesonfile,
                                                     **self.__dict__)

        for gamma in self.gammas:

            mesonfiles = tuple(
                self.mesonfile.format(
                    w=contraction[i],
                    v=contraction[i + 1],
                    gamma=g,
                    **mesonfile_replacements
                ) for i, g in zip([0, 2, 4], [gamma, "G1_G1", gamma])
            )

            mat_gen = MesonLoader(mesonfiles=mesonfiles,
                                  times=times, **self.meson_params)

            cij = cpnp.zeros((self.nt, self.nt, self.nt), dtype=np.complex128)

            for (t1, m1), (t2, m2), (t3, m3) in mat_gen:

                logging.info(f"Contracting {gamma}: {t1},{t2},{t3}")
                cij[t1, t2, t3] = self.contract(
                    m1, m2, m3, open_indices=[0, 1, 2])

                if self.symmetric:
                    raise Exception(
                        "Symmetric 3dim optimization not implemented.")

            logging.debug("Contraction completed.")

            if 'comm' in self.__dict__ and self.comm_size > 1:
                temp = None
                if self.rank == 0:
                    temp = cpnp.empty_like(cij)
                self.comm.Barrier()
                self.comm.Reduce(cij, temp, op=MPI.SUM, root=0)

                if self.rank == 0:
                    corr[gamma] = convert_to_numpy(temp)  # time_average(temp)
            else:
                corr[gamma] = convert_to_numpy(cij)  # time_average(cij)

        return corr

    def qed_conn_4pt(self, contraction):

        cij = np.zeros((self.nt, self.nt, self.nt, self.nt),
                       dtype=np.complex128)

        seedkey = "".join(contraction)

        gammas = [self.mesonKey.format(gamma=g) for g in ['X', 'Y', 'Z']]

        corr = dict(
            zip(self.subdiagrams, [
                {seedkey: dict(zip(gammas, [{}] * len(gammas)))}
            ] * len(self.subdiagrams)))

        matg1 = [self.load_meson(contraction[0], contraction[1], gamma)
                 for gamma in gammas]
        matg2_photex = []
        matg2_selfen = []

        for gamma in gammas:
            if "photex" in self.subdiagrams:
                matg2_photex.append(self.load_meson(
                    contraction[4], contraction[5], gamma))
            if "selfen" in self.subdiagrams:
                matg2_selfen.append(self.load_meson(
                    contraction[6], contraction[7], gamma))

        for i in range(self.Nem):

            emlabel = f"{self.emseedstring}_{i}"

            selfen_p2_key = contraction[4] + contraction[5] + emlabel
            photex_p2_key = contraction[6] + contraction[7] + emlabel

            matp1 = self.load_meson(contraction[2], contraction[3], emlabel)
            if "selfen" in self.subdiagrams:
                matp2_selfen = self.load_meson(
                    contraction[4], contraction[5], emlabel)
            if "photex" in self.subdiagrams:
                matp2_photex = self.load_meson(
                    contraction[6], contraction[7], emlabel)

            for j, gamma in enumerate(gammas):
                selfen_g2_key = contraction[6] + contraction[7] + gamma
                photex_g2_key = contraction[4] + contraction[5] + gamma

                for d in self.subdiagrams:

                    if d == "photex":
                        cij[:] = self.contract(
                            m1=np.array(matg1[j]),
                            m2=matp1,
                            m3=np.array(matg2_photex[j]),
                            m4=matp2_photex,
                            open_indices=(0, 2)
                        )
                    elif d == "selfen":
                        cij[:] = self.contract(
                            m1=np.array(matg1[j]),
                            m2=matp1,
                            m3=matp2_selfen,
                            m4=np.array(matg2_selfen[j])
                        )
                    else:
                        raise Exception(f"Unrecognized diagram label: {d}")

                    corr[d][seedkey][gamma][emlabel] = time_average(cij)

        return corr

    def execute(self, contraction):

        if cpnp.__name__ == 'cupy':
            my_device = self.rank % cp.cuda.runtime.getDeviceCount()
            logging.debug(f"Rank {self.rank} is using gpu device {my_device}")
            cp.cuda.Device(my_device).use()

        logging.info(f"Processing mode: {', '.join(contraction)}")

        return self.run_contract(contraction)

    def make_contraction_key(self, contraction):
        con_key = "".join(contraction)

        if self.has_high:
            con_key = con_key.replace(self.high_label, "")

        if self.has_low:
            con_key = con_key.replace(self.low_label, "e")

        return con_key


def main():

    comm = MPI.COMM_WORLD

    if len(sys.argv) != 2:
        logging.error(("Must provide sub-ensemble and "
                       "config data in 'series.ensemble' format"))
        exit()

    series, cfg = sys.argv[1].split('.')

    params = utils.load_param('params.yaml')

    if 'logging_level' in params['contract'] and comm.Get_rank() == 0:
        logging_level = params['contract']['logging_level']
    else:
        logging_level = logging.INFO

    logging.basicConfig(
        format="%(asctime)s - %(levelname)-5s - %(message)s",
        style="%",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging_level,
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )

    if 'hardware' in params['contract'] and \
            params['contract']['hardware'] == 'cpu':
        globals()['cpnp'] = np

    diagrams = params['contract']['diagrams']
    contractor_dict = dict(zip(diagrams, [
        Contractor(series=series,
                   cfg=cfg,
                   nt=int(params['lmi_param']['TIME']),
                   comm=comm,
                   **params['contract'][d]) for d in diagrams
    ]))

    for diagram, contractor in contractor_dict.items():

        outfile_replacements = FSFormat.formatdict(
            contractor.outfile, **contractor.__dict__)

        overwrite_outfile = getattr(contractor, 'overwrite', False)

        nmesons = contractor.npoint

        low_min = 0 if contractor.has_high else nmesons
        low_max = nmesons + 1 if contractor.has_low else 0

        perms = sum([
            list(multiset_permutations(['L'] * nlow + ['H'] * (nmesons - nlow)))
            for nlow in range(low_min, low_max)
        ], [])
        perms = list(map("".join, perms))
        # Overwrite chosen permutations with user input, if provided
        perms = getattr(contractor, 'perms', perms)

        logging.debug(f"Computing permutations: {perms}")

        for perm in perms:
            nlow = perm.count('L')

            permkey = "".join(
                sum(((perm[i], perm[(i + 1) % nmesons])
                     for i in range(nmesons)), ())
            )

            if contractor.has_high:
                # Build list of high source indices,
                # e.g. [[0,1], [0,2], ...]
                seeds = list(map(list, itertools.combinations(
                    list(range(contractor.high_count)), nmesons - nlow)))
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
            else:
                seeds = [[]]

            outfile = contractor.outfile.format(
                permkey=permkey, diagram=diagram, **outfile_replacements)

            if overwrite_outfile or not os.path.exists(outfile):
                logging.info(f'Contracting diagram: {diagram} ({permkey})')
            else:
                logging.info(f'Skipping write. File exists: {outfile}')
                continue

            contraction_list = [
                [
                    contractor.low_label
                    if seed[i] is None else
                    ('w' if i % 2 == 0 else 'v')
                    + contractor.high_label + s
                    for i, s in enumerate(map(str, seed))
                ]
                for seed in seeds
            ]

            start_time = perf_counter()

            corr = dict(zip(
                map(contractor.make_contraction_key, contraction_list),
                map(contractor.execute, contraction_list)
            ))

            stop_time = perf_counter()

            logging.debug('')
            logging.debug('    Total elapsed time for %s = %g seconds.' % (
                permkey, stop_time - start_time))
            logging.debug('')

            if ('comm' not in contractor.__dict__ or contractor.rank == 0):
                if not os.path.exists(os.path.dirname(outfile)):
                    os.makedirs(os.path.dirname(outfile))
                pickle.dump(corr, open(outfile, 'wb'))


if __name__ == '__main__':
    main()
