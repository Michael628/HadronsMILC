#! /usr/bin/env python3
import gc
import os
import sys
import itertools
import numpy as np
import cupy as cp
import cupy as cpnp
import h5py
import pickle
from time import perf_counter
from sympy.utilities.iterables import multiset_permutations
import opt_einsum as oe
import logging
from multiprocessing import Process, Pool, Lock, Manager
import concurrent.futures
import time as timer
import copy
from mpi4py import MPI
from dataclasses import dataclass
import yaml

def loadParam(file):
    """Read the YAML parameter file"""

    try:
        param = yaml.safe_load(open(file,'r'))
    except:
        print("WARNING: loadParam failed")
        sys.exit(1)

    return param

def build_subdictionary(key_list: list[str], source_dict: dict):
    """Builds a dictionary from elements of `source_dict`.

    Parameters
    ----------
    replacement_list : list[str]
        A list of keywords to find in `source_dict

    source_dict : dict
        The dictionary from which the returned dict is built

    Returns
    -------
    sub_dict : dict
        Dictionary for all elements in `source_dict` with values of type str
    """
    sub_dict = {}
    for key in key_list:
        if key in source_dict and type(source_dict[key]) is str:
            sub_dict[key] = source_dict[key]

    return sub_dict

@dataclass
class MesonLoader:
    """Iterable object that loads MesonFields for processing by a Contractor
    Parameters
    ----------
    file : str
        File location of meson field to load.
    times : iterable
        An iterable object containing slices. The slices will be used to load
        the corresponding time slices from the meson hdf5 file in the order
        they appear in `times`
    mass_shift : bool, optional
        If true, it will be assumed that the meson field being loaded is constructed
        from eigenvectors with eigenvalues based on `mass_old`. The eigenvalues will
        then be replaced by corresponding values using `mass_new` via the `shift_mass` method.
    evalfile : str, optional
        File location of hdf5 file containing list of eigenvalues. Required if `mass_shift` is True.
    oldmass : str, optional
        Mass for evals in `evalfile`. Required if `mass_shift` is True.
    newmass : str, optional
        New mass to use with `evalfile`. Required if `mass_shift` is True.
    """
    mesonfiles: list[str]
    times: tuple
    mass_shift: bool = False
    evalfile: str = ""
    oldmass: str = ""
    newmass: str = ""
    vmax_index: int = -1
    wmax_index: int = -1

    def shift_mass(self,mat: np.ndarray):
        raise Exception("Update this code")

        # with h5py.File(self.evalfilestem.format(series=self.series,cfg=self.cfg),'r') as f:
        #     try:
        #         evals = np.array(f['/evals'][()].view(np.float64),dtype=np.float64)
        #     except KeyError:
        #         evals = np.array(f['/EigenValueFile/evals'][()].view(np.float64),dtype=np.float64)

        # evals = np.sqrt(evals)
        # mass_old = float(f'0.{self.massold}')
        # mass_new = float(f'0.{self.mass}')

        # eval_scaling = np.zeros((len(evals),2),dtype=np.complex128)
        # eval_scaling[:,0] = np.divide(mass_old+1.j*evals,mass_new+1.j*evals)
        # eval_scaling[:,1] = np.conjugate(eval_scaling[:,0])
        # eval_scaling = eval_scaling.reshape((-1,))

        # mat[:] = np.multiply(mat,eval_scaling[np.newaxis,np.newaxis,:])

    def loadMeson(self, file, time: slice=slice(None)):
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
        Assumes array is single precision complex. Promotes to double precision.
        """
        t1 = timer.time()
        with h5py.File(file,"r") as f:
            a_group_key = list(f.keys())[0]

            temp = f[a_group_key]['a2aMatrix']
            temp = cpnp.array(temp[time,self.wmax_index,self.vmax_index].view(np.complex64),
                            dtype=np.complex128)
            t2 = round((timer.time() - t1),2)
            print(f"Loaded array {temp.shape} in {t2} sec")

            return temp

    def __iter__(self):
        self.mesonlist = [ None for _ in range(len(self.mesonfiles))]
        self.iter_count = 0
        return self

    def __next__(self):
        if self.iter_count < len(self.times[0]):

            current_times = [self.times[i][self.iter_count] for i in range(len(self.mesonfiles))]

            for i, (time,file) in enumerate(zip(current_times,self.mesonfiles)):
                try:
                    if self.mesonlist[i] is not None: # Check if meson exists from last iter
                        if self.mesonlist[i][0] == time:
                            continue

                    j = current_times[:i].index(time) # Check that time is the same
                    if file != self.mesonfiles[j]: # Check that file is the same
                        raise ValueError

                    print(f"Found {time} at index {j}")
                    self.mesonlist[i] = (time,self.mesonlist[j][1]) # Copy reference

                except ValueError:
                    print(f"Loading {time} from {file}")
                    self.mesonlist[i] = (time,self.loadMeson(file,time)) # Load new file

            self.iter_count += 1
            return tuple(self.mesonlist)
        else:
            self.mesonlist = None
            raise StopIteration


    
class Contractor:
    """Performs contractions of MesonField
    """
    def __init__(self, series: str, cfg: str, **kwargs):
        self.__dict__.update(kwargs)

        if 'comm' in kwargs:
            self.rank = self.comm.Get_rank()
            self.comm_size = self.comm.Get_size()

        # Set max indices of meson field v,w dimensions
        if 'wmax_index' not in kwargs:
            self.wmax_index = None
        if 'vmax_index' not in kwargs:
            self.vmax_index = None
        if 'symmetric' not in kwargs:
            self.symmetric = False

        self.series = series
        self.cfg    = cfg

    def printr0(self,text):
        if self.rank == 0:
            print(text)
        return

    def convert_to_numpy(self,corr):
        return corr if type(corr) is np.ndarray else cp.asnumpy(corr)

    def contract(self, m1: np.ndarray, m2: np.ndarray, m3: np.ndarray=None,
                    m4: np.ndarray=None, open_indices: tuple=(0,-1)):
        """Performs contraction of up to 4 3-dim arrays down to one 2-dim array

        Parameters
        ----------
        m1 : ndarray
        m2 : ndarray
        m3 : ndarray, optional
        m4 : ndarray, optional
        open_indices : tuple, default=(0,-1)
            A two-dimensional tuple containing the time indices that will not be 
            contracted in the full product. The default=(0,-1) leaves the indices
            of the first and last matrix open, summing over all others.

        Returns
        -------
        ndarray
            The resultant 2-dim array from the contraction
        """
        if len(open_indices) != 2:
            raise Exception("Length of open_indices must be 2")

        index_list = ['i','j','k','l'][:self.npoint]
        out_indices = "".join(index_list[i] for i in open_indices)

        if self.npoint == 2:   # two-point contractions
            cij = oe.contract(f'imn,jnm->{out_indices}',m1,m2)

        elif self.npoint == 3: # three-point contractions
            cij = oe.contract(f'imn,jno,kom->{out_indices}',m1,m2,m3)

        elif self.npoint == 4: # four-point contractions
            cij = oe.contract(f'imn,jno,kop,lpm->{out_indices}',m1,m2,m3,m4)
        else:
            raise Exception("Expecting 2 <= self.npoint <= 4")

        return cij

    def time_average(self,cij):
        """Takes an n,n array and returns a 1-dim array of length n, where the i-th
        output element is the sum of all input elements with indices separated by i. 

        Parameters
        ----------
        cij: ndarray
            A 2-dim square array

        Returns
        -------
        ndarray
            A 1-dim array with entries that are the average of input entries with equal
            index separations
        """

        C = cpnp.zeros((self.nt,),dtype=np.complex128)
        ones = cpnp.ones(cij.shape)
        t_range = cpnp.array(range(self.nt))

        t_start = (slice(None),None)
        t_end = (None,slice(None))
        t_mask = cpnp.mod(t_range[t_end]*ones-t_range[t_start]*ones,cpnp.array([self.nt]))

        for t in range(self.nt):
            C[t] = cij[t_mask==t].sum()

        #return cp.asnumpy(C/self.nt)
        return C/self.nt

    def generate_time_sets(self,symmetric: bool = False):
        """Breaks meson field time extent into `comm_size` blocks and returns unique
        list of blocks for each `rank`.

        Parameters
        ----------
        symmetric : bool, optional
            If true, time sets generated are upper triangular

        Returns
        -------
        tuple
            A tuple of length `npoint` containing lists of slices. One lsit for each meson field
        """

        workload = self.comm_size

        slice_indices = list(itertools.product(range(self.comm_size),repeat=self.npoint))

        if symmetric: # filter for only upper-triangular slices and split work evenly
            slice_indices = list(filter(lambda x: list(x) == sorted(x),slice_indices))
            workload = int((len(slice_indices)+self.comm_size-1)/self.comm_size)

        offset = int(self.rank*workload)

        slice_indices = list(zip(*slice_indices[offset:offset+workload]))

        tspacing = int(self.nt/self.comm_size)

        return tuple([slice(int(ti*tspacing),int((ti+1)*tspacing)) for ti in times]
                for times in slice_indices)


    def vec_conn_2pt(self,seedlist):

        corr = {}

        times = self.generate_time_sets(self.symmetric)

        mesonfile_replacements = {}
        if 'replacements' in self.mesonfile:
            mesonfile_replacements = build_subdictionary(self.mesonfile['replacements'],self.__dict__)

        for gamma in self.gammas:

            self.printr0(f"Processing {gamma}")

            mesonfiles = tuple(self.mesonfile['filestem'].format(
                                                w=seedlist[i],
                                                v=seedlist[i+1],
                                                gamma=gamma,
                                                **mesonfile_replacements) for i in [0,2])

            mat_gen = MesonLoader(mesonfiles=mesonfiles,
                                wmax_index=slice(self.wmax_index),
                                vmax_index=slice(self.vmax_index),
                                times=times)

            cij = cpnp.zeros((self.nt,self.nt),dtype=np.complex128)

            for (t1,m1),(t2,m2) in mat_gen:
                self.printr0(f"contracting {t1},{t1}")
                cij[t1,t2] = self.contract(m1,m2)
                if self.symmetric and t1 != t2:
                    cij[t2,t1] = cij[t1,t2].T

            self.printr0("contracted")
            if 'comm' in self.__dict__ and self.comm_size > 1:
                temp = None
                if self.rank == 0:
                    temp = cpnp.empty_like(cij)
                self.comm.Barrier()
                self.comm.Reduce(cij,temp,op=MPI.SUM,root=0)

                if self.rank == 0:
                    corr[gamma] = self.convert_to_numpy(self.time_average(temp))
            else:
                corr[gamma] = self.convert_to_numpy(self.time_average(cij))

            del m1,m2
        return corr

    def sib_conn_3pt(self,seedlist):

        corr = {}

        times1, times2, times3 = self.generate_time_sets(self.symmetric)

        mesonfile_replacements = {}
        if 'replacements' in self.mesonfile:
            mesonfile_replacements = build_subdictionary(self.mesonfile['replacements'],self.__dict__)

        mat1 = MesonLoader(mesonfile=self.mesonfile['filestem'].format(
                                            w=seedlist[2],
                                            v=seedlist[3],
                                            gamma="G1_G1",
                                            **mesonfile_replacements),
                            times=times2)

        for gamma in self.gammas:

            self.printr0(gamma)
            mat1 = MesonLoader(mesonfile=self.mesonfile['filestem'].format(
                                                w=seedlist[0],
                                                v=seedlist[1],
                                                gamma=gamma,
                                                **mesonfile_replacements),
                                times=times1)

            mat3 = MesonLoader(mesonfile=self.mesonfile['filestem'].format(
                                                w=seedlist[4],
                                                v=seedlist[5],
                                                gamma=gamma,
                                                **mesonfile_replacements),
                                times=times3)

            cij = cpnp.zeros((self.nt,self.nt,self.nt),dtype=np.complex128)

            for (t1,m1),(t2,m2),(t3,m3) in zip(mat1,mat2,mat3):
                cij[t1,t2,t3] = self.contract(m1,m2,m3)

            self.printr0("contracted")
            if 'comm' in self.__dict__ and self.comm_size > 1:
                temp = None
                if self.rank == 0:
                    temp = cpnp.empty_like(cij)
                self.comm.Barrier()
                self.comm.Reduce(cij,temp,op=MPI.SUM,root=0)

                if self.rank == 0:
                    corr[gamma] = temp#self.time_average(temp)
            else:
                corr[gamma] = cij#self.time_average(cij)

        return corr

    def qed_conn_4pt(self,seedlist):

        cij = np.zeros((self.nt,self.nt,self.nt,self.nt),dtype=np.complex128)
        
        seedkey = "".join(seedlist)

        gammas = [self.mesonKey.format(gamma=g) for g in ['X','Y','Z']]

        corr = dict(zip(self.subdiagrams,[{seedkey:dict(zip(gammas,
                                                         [{}]*len(gammas)))}]*len(self.subdiagrams)))

        #matg1 = [cp.asnumpy(self.loadMeson(seedlist[0],seedlist[1],gamma)) for gamma in gammas]
        matg1 = [self.loadMeson(seedlist[0],seedlist[1],gamma) for gamma in gammas]
        matg2_photex = []
        matg2_selfen = []

        for gamma in gammas:
            if "photex" in self.subdiagrams:
                #matg2_photex.append(cp.asnumpy(self.loadMeson(seedlist[4],seedlist[5],gamma)))
                matg2_photex.append(self.loadMeson(seedlist[4],seedlist[5],gamma))
            if "selfen" in self.subdiagrams:
                #matg2_selfen.append(cp.asnumpy(self.loadMeson(seedlist[6],seedlist[7],gamma)))
                matg2_selfen.append(self.loadMeson(seedlist[6],seedlist[7],gamma))

        for i in range(self.Nem):

            emlabel = f"{self.emseedstring}_{i}"

            selfen_p2_key = seedlist[4]+seedlist[5]+emlabel
            photex_p2_key = seedlist[6]+seedlist[7]+emlabel

            matp1 = self.loadMeson(seedlist[2],seedlist[3],emlabel)
            if "selfen" in self.subdiagrams:
                matp2_selfen = self.loadMeson(seedlist[4],seedlist[5],emlabel)
            if "photex" in self.subdiagrams:
                matp2_photex = self.loadMeson(seedlist[6],seedlist[7],emlabel)

            for j,gamma in enumerate(gammas):
                selfen_g2_key = seedlist[6]+seedlist[7]+gamma
                photex_g2_key = seedlist[4]+seedlist[5]+gamma

                for d in self.subdiagrams:

                    if d == "photex":
                        cij[:] = self.contract(m1=np.array(matg1[j]),m2=matp1,m3=np.array(matg2_photex[j]),m4=matp2_photex,
                                            open_indices=(0,2))
                    elif d == "selfen":
                        cij[:] = self.contract(m1=np.array(matg1[j]),m2=matp1,m3=matp2_selfen,m4=np.array(matg2_selfen[j]))
                    else:
                        raise Exception(f"Unrecognized diagram label: {d}")

                    corr[d][seedkey][gamma][emlabel] = self.time_average(cij)

        return corr

    def execute(self,diagram,seeds):

        permkey = "".join(sum(((self.perm[i], self.perm[(i+1) % self.npoint]) for i in range(self.npoint)),()))

        # os.system("nvidia-smi")
        with cp.cuda.Device(0):
            veclist = sum(tuple(zip([(mode,(self.low if mode == 'L' else "wseed%i" % seeds[i])) for mode, i in zip(self.perm,self.high_indices)],
                                    [(mode,(self.low if mode == 'L' else "vseed%i" % seeds[i])) for mode, i in zip(self.perm[1:]+self.perm[:1],self.high_indices[1:]+self.high_indices[:1])])),
                          ())

            seedlist = [v[1] for v in veclist]
            seedkey = "".join(seedlist)

            corr = {permkey:{diagram:{}}}

            # if seedkey not in C[permkey][self.diagrams[0]].keys():
            self.printr0(f"Processing modes: {veclist}")
            if diagram in ["vec_local","vec_onelink"]:
                corr[permkey][diagram][seedkey] = self.vec_conn_2pt(seedlist)
            if "sib" == diagram:
                corr[permkey]['sib'][seedkey] = self.sib_conn_3pt(seedlist)

            if "selfen" == diagram or "photex" == diagram:
                corr[permkey].update(self.qed_conn_4pt(seedlist))

            # cp.cuda.Device(index).synchronize()

            # i = barrier.wait()

            # if i == 0:
                    # Mesons.clear()
            # barrier.wait()

            return corr

def main():

    comm = MPI.COMM_WORLD

    if len(sys.argv) != 2:
        print("Must provide sub-ensemble and config data in 'series.ensemble' format")
        exit()

    series,cfg = sys.argv[1].split('.')
    processes = 1

    params = loadParam('params.yaml')

    if 'hardware' in params['contract'] and params['contract']['hardware'] == 'cpu':
        import numpy as cpnp

    diagrams = params['contract']['diagrams']
    contractor_dict = dict(zip(diagrams,
        [Contractor(series=series, 
                    cfg=cfg,
                    nt=int(params['LMIparam']['TIME']),
                    comm=comm,
                    **params['contract'][d]) for d in diagrams]))


    C = {}

    for diagram,contractor in contractor_dict.items():


        outfile_replacements = {}
        if 'replacements' in contractor.outfile:
            outfile_replacements = build_subdictionary(contractor.outfile['replacements'],contractor.__dict__)

        contractor.printr0(f'Contracting diagram: {diagram}')
        nmesons = contractor.npoint

        low_range = [nmesons] if 'high' not in contractor.__dict__ else range(contractor.low_start if 'low_start' in contractor.__dict__ else 0,nmesons+1)

        for nlow in low_range:
    
            perms = multiset_permutations(['L']*nlow+['H']*(nmesons-nlow))
            for perm in perms:

                contractor.perm = perm

                permkey = "".join(sum(((perm[i], perm[(i+1) % nmesons]) for i in range(nmesons)),()))

                if permkey not in C.keys():
                    C[permkey] = {}

                contractor.high_indices = [sum([0 if a =='L' else 1 for a in perm][:i]) for i in range(nmesons)]

                # Produces list of tuples like (('L','e1000'), ('H','wseed1'), ('H','vseed1'), ('L','e1000'))
                seeds = [()]
                if 'high' in contractor.__dict__:
                    seeds = itertools.combinations(list(range(contractor.high['size'])),nmesons-nlow)

                start_time = perf_counter( )

                #p = Pool(processes)

                #for result in p.starmap(execute,[(i % processes,item,diagram,contractor) for i,item in enumerate(items)])
                _ = [ 0 if not bool(result[permkey]) else 
                      C[permkey].update(result[permkey]) if diagram not in C[permkey] else 
                      C[permkey][diagram].update(v1) if seed not in C[permkey][diagram] else
                      C[permkey][diagram][seed].update(v2) if gamma not in C[permkey][diagram][seed] or type(C[permkey][diagram][seed][gamma]) is not dict else
                      C[permkey][diagram][seed][gamma].update(v3)
                      for result in [contractor.execute(diagram,seed) for i,seed in enumerate(seeds)]
                      for diagram,v1 in result[permkey].items()
                      for seed,v2 in v1.items()
                      for gamma, v3 in v2.items()]

                stop_time = perf_counter( )
                contractor.printr0('')
                contractor.printr0('    Elapsed wall clock time for IO+contraction = %g seconds.' % (stop_time - start_time) )
                contractor.printr0('')
                contractor.printr0(f"Finished {permkey}")
                
                if 'comm' not in contractor.__dict__ or contractor.rank == 0:
                    outfile = contractor.outfile['filestem'].format(permkey=permkey,diagram=diagram,**outfile_replacements)
                    d = C[permkey][diagram]
                    pickle.dump(C[permkey][diagram],open(outfile,'wb'))

if __name__ == '__main__':
    main()
