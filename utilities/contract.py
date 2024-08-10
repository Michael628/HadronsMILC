#! /usr/bin/env python3
import gc
import os
import sys
import itertools
import numpy as np
import cupy as cp
import h5py
import pickle
from time import perf_counter
from sympy.utilities.iterables import multiset_permutations
import opt_einsum as oe
import logging
from multiprocessing import Process, Pool, Lock, Manager
import concurrent.futures
import time
import copy
from mpi4py import MPI
from dataclasses import dataclass

sys.path.append("../scripts")

from TodoUtils import *

@dataclass
class MesonLoader:
    """Iterable object that loads MesonFields for processing by a Contractor
    Parameters
    ----------
    file : str
        File location of meson field to load.
    times: iterable
        An iterable object containing slices. The slices will be used to load
        the corresponding time slices from the meson hdf5 file in the order
        they appear in `times`
    mass_shift: bool, optional
        If true, it will be assumed that the meson field being loaded is constructed
        from eigenvectors with eigenvalues based on `mass_old`. The eigenvalues will
        then be replaced by corresponding values using `mass_new` via the `shift_mass` method.
    evalfile: str, optional
        File location of hdf5 file containing list of eigenvalues. Required if `mass_shift` is True.
    oldmass: str, optional
        Mass for evals in `evalfile`. Required if `mass_shift` is True.
    newmass: str, optional
        New mass to use with `evalfile`. Required if `mass_shift` is True.
    """
    mesonfile: str
    times: tuple
    mass_shift: bool = False
    evalfile: str = ""
    oldmass: str = ""
    newmass: str = ""

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

    def loadMeson(self,time: slice=slice(None)):
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
        with h5py.File(self.mesonfile,"r") as f:
            a_group_key = list(f.keys())[0]
            temp = f[a_group_key]['a2aMatrix']
            return np.array(temp[time].view(np.complex64),
                            dtype=np.complex128)

    def __iter__(self):
        self.iter_count = 0
        return self

    def __next__(self):
        if self.iter_count < len(self.times):
            self.iter_count += 1
            return (self.times[self.iter_count-1], self.loadMeson(time=self.times[self.iter_count-1]))
        else:
            raise StopIteration


    
class Contractor:
    """Performs contractions of MesonField
    """
    def __init__(self, series: str, cfg: str, **kwargs):
        self.__dict__.update(kwargs)

        if 'comm' in kwargs:
            self.rank = self.comm.Get_rank()
            self.comm_size = self.comm.Get_size()

        self.series = series
        self.cfg    = cfg

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

        index_list = ['i','j','k','l'][:self.Npoint]
        out_indices = "".join(index_list[i] for i in open_indices)

        if self.Npoint == 2:   # two-point contractions
            cij = oe.contract(f'imn,jnm->{out_indices}',m1,m2)

        elif self.Npoint == 3: # three-point contractions
            cij = oe.contract(f'imn,jno,kom->{out_indices}',m1,m2,m3)

        elif self.Npoint == 4: # four-point contractions
            cij = oe.contract(f'imn,jno,kop,lpm->{out_indices}',m1,m2,m3,m4)
        else:
            raise Exception("Expecting 2 <= self.Npoint <= 4")

        return cij

    def time_average(self,cij: np.ndarray):
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

        C = np.zeros((self.nt,),dtype=np.complex128)
        ones = np.ones(cij.shape)
        t_range = np.array(range(self.nt))

        t_start = (slice(None),None)
        t_end = (None,slice(None))
        t_mask = np.mod(t_range[t_end]*ones-t_range[t_start]*ones,np.array([self.nt]))

        for t in range(self.nt):
            C[t] = cij[t_mask==t].sum()

        #return cp.asnumpy(C/self.nt)
        return C/self.nt

    def generate_time_sets(self):
        """Breaks meson field time extent into `comm_size` blocks and returns unique
        list of blocks for each `rank`.

        Returns
        -------
        tuple
            A tuple of length `Npoint` containing lists of slices. One lsit for each meson field
        """

        spacing = self.nt/self.comm_size

        slice_indices = zip(*itertools.product(range(self.comm_size),repeat=self.Npoint))

        return tuple([slice(int(tindex*spacing),int((tindex+1)*spacing)) for tindex in tset][self.rank::self.comm_size] 
                for tset in slice_indices)


    def vec_conn_2pt(self,seedlist):

        corr = {}

        gammas = [self.mesonKey.format(gamma=g) for g in ['X','Y','Z']]

        times1, times2 = self.generate_time_sets()

        for gamma in gammas:

            print(gamma)
            mat1 = MesonLoader(mesonfile=self.mesonfilestem.format(
                                                w=seedlist[0],
                                                v=seedlist[1],
                                                g=gamma,
                                                series=self.series,
                                                cfg=self.cfg),
                                times=times1)

            mat2 = MesonLoader(mesonfile=self.mesonfilestem.format(
                                                w=seedlist[2],
                                                v=seedlist[3],
                                                g=gamma,
                                                series=self.series,
                                                cfg=self.cfg),
                                times=times2)

            cij = np.zeros((self.nt,self.nt),dtype=np.complex128)

            for (t1,m1),(t2,m2) in zip(mat1,mat2):
                cij[t1,t2] = self.contract(m1,m2)

            print("contracted")
            if 'comm' in self.__dict__ and self.comm_size > 1:
                temp = None
                if self.rank == 0:
                    temp = np.empty_like(cij)
                self.comm.Barrier()
                self.comm.Reduce(cij,temp,op=MPI.SUM,root=0)

                if self.rank == 0:
                    corr[gamma] = self.time_average(temp)
            else:
                corr[gamma] = self.time_average(cij)

        return corr

    def sib_conn_3pt(self,seedlist):

        corr = {}
        cij = np.zeros((self.nt,self.nt,self.nt),dtype=np.complex128)

        gammas = [self.mesonKey.format(gamma=g) for g in ['X','Y','Z']]

        mat2 = self.loadMeson(seedlist[2],seedlist[3],self.mesonKey.format(gamma='1'))
            
        for gamma in gammas:
            mat1 = self.loadMeson(seedlist[0],seedlist[1],gamma)
            mat3 = self.loadMeson(seedlist[4],seedlist[5],gamma)

            corr[gamma] = self.time_average(self.contract(m1=mat1,m2=mat2,m3=mat3))

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

def execute(index,seeds,diagram,contractor):

    permkey = "".join(sum(((contractor.perm[i], contractor.perm[(i+1) % contractor.Npoint]) for i in range(contractor.Npoint)),()))

    if True:
    #with cp.cuda.Device(index):
        veclist = sum(tuple(zip([(mode,(contractor.elabel if mode == 'L' else "wseed%i" % seeds[i])) for mode, i in zip(contractor.perm,contractor.highIndices)],
                                [(mode,(contractor.elabel if mode == 'L' else "vseed%i" % seeds[i])) for mode, i in zip(contractor.perm[1:]+contractor.perm[:1],contractor.highIndices[1:]+contractor.highIndices[:1])])),
                      ())

        seedlist = [v[1] for v in veclist]
        seedkey = "".join(seedlist)

        corr = {permkey:{diagram:{}}}

        # if seedkey not in C[permkey][contractor.diagrams[0]].keys():
        print(round(time.time()*1000),veclist)
        if diagram in ["vec_local","vec_onelink"]:
            corr[permkey][diagram][seedkey] = contractor.vec_conn_2pt(seedlist)
        if "sib" == diagram:
            corr[permkey]['sib'][seedkey] = contractor.sib_conn_3pt(seedlist)

        if "selfen" == diagram or "photex" == diagram:
            corr[permkey].update(contractor.qed_conn_4pt(seedlist))

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

    diagrams = params['contract_a2a']['diagrams']
    contractor_dict = dict(zip(diagrams,
        [Contractor(series=series, 
                    cfg=cfg,
                    nt=int(params['LMIparam']['TIME']),
                    comm=comm,
                    **params['contract_a2a'][d]) for d in diagrams]))


    C = {}

    for diagram,contractor in contractor_dict.items():

        print(f'Contracting diagram: {diagram}')
        Nmesons = contractor.Npoint

        for Nlow in range(contractor.lowStart,Nmesons+1):
    
            HLperms = multiset_permutations(['L']*Nlow+['H']*(Nmesons-Nlow))
            for perm in HLperms:

                contractor.perm = perm

                permkey = "".join(sum(((perm[i], perm[(i+1) % Nmesons]) for i in range(Nmesons)),()))

                #for s,d in zip(myContractor.subdirs,myContractor.diagrams):
                #    outfile = myContractor.outfilestem.format(permkey=permkey,diagram=d,
                #                                              subdir=s.format(**myContractor.__dict__),**myContractor.__dict__)
                #    if os.path.exists(outfile):
                #        diagrams.remove(d)

                if permkey not in C.keys():
                    C[permkey] = {}

                contractor.highIndices = [sum([0 if a =='L' else 1 for a in perm][:i]) for i in range(Nmesons)]

                # Produces list of tuples like (('L','e1000'), ('H','wseed1'), ('H','vseed1'), ('L','e1000'))
                items = itertools.combinations(list(range(contractor.Nseeds)),Nmesons-Nlow)

                start_time = perf_counter( )

                #p = Pool(processes)

                #for result in p.starmap(execute,[(i % processes,item,diagram,contractor) for i,item in enumerate(items)])
                _ = [ 0 if not bool(result[permkey]) else 
                      C[permkey].update(result[permkey]) if diagram not in C[permkey] else 
                      C[permkey][diagram].update(v1) if seed not in C[permkey][diagram] else
                      C[permkey][diagram][seed].update(v2) if gamma not in C[permkey][diagram][seed] or type(C[permkey][diagram][seed][gamma]) is not dict else
                      C[permkey][diagram][seed][gamma].update(v3)
                      for result in [execute(i % processes,item,diagram,contractor) for i,item in enumerate(items)]
                      for diagram,v1 in result[permkey].items()
                      for seed,v2 in v1.items()
                      for gamma, v3 in v2.items()]

                stop_time = perf_counter( )
                print('')
                print('    Elapsed wall clock time for IO+contraction = %g seconds.' % (stop_time - start_time) )
                print('')
                print(f"Finished {permkey}")
                
                if 'comm' not in contractor.__dict__ or contractor.rank == 0:
                    for s in contractor.subdirs:
                        outfile = contractor.outfilestem.format(permkey=permkey,diagram=diagram,
                                                                subdir=s.format(**contractor.__dict__),**contractor.__dict__)
                        pickle.dump(C[permkey][diagram],open(outfile,'wb'))

if __name__ == '__main__':
    main()
