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

if len(sys.argv) != 2:
    print("Must provide sub-ensemble and config data in 'series.ensemble' format")
    exit()

series,cfg = sys.argv[1].split('.')
processes = 1

manager = Manager()

class Contractor:
    def __init__(self,series,cfg):
        # self.Nt = 64
        self.Nt = 48
        self.evalfilestem = "eigs/evals/evalmassless3248f211b580m002426m06730m8447nv1000{series}.{cfg}.h5"
        self.mesonfilestem = "e1000n1dt1/mesons/m001555/mf_{series}.{cfg}/{g}_0_0_0.h5"
        #self.mesonfilestem = "e1000n1dt1/mesons/mf_{series}_{w}_{v}.{cfg}/{g}_0_0_0.h5"
        # self.mesonfilestem = "e2000n1dt1/mesons/m001907/mf_{series}.{cfg}/{g}_0_0_0.h5"
        # self.elabel = "e2000"
        self.elabel = "e1000"
        self.series = series
        self.cfg    = cfg
        self.Nseeds = 20
        self.Nem    = 10
        self.emseedstring  = 'em'
        self.seedstring    = 'seed'
        #self.diagrams      = ['selfen','photex']
        #self.subdirs       = ["em{Nem}","em{Nem}"]
        self.diagrams      = ['vec_local']
        self.massold = "001555"
        #self.massold = "002426"
        #self.massnew = "002426"
        #self.massnew = "001555"
        self.massnew = "003297"
        self.subdirs       = [f"m{self.massnew}"]
        #self.outfilestem   = 'a2a_corrs/{diagram}/{subdir}/corr_vec_local_conn_{diagram}_{permkey}_seed{Nseeds}_{series}.{cfg}.p'
        self.outfilestem   = 'a2a_corrs/{diagram}/{subdir}/corr_vec_local_conn_{diagram}_{permkey}_{series}.{cfg}.p'
        #self.Npoint        = 4
        # self.outfilestem   = 'a2a_corrs/{diagram}/{subdir}/corr_conn_{diagram}_{permkey}_{series}.{cfg}.p'
        # self.outfilestem   = 'a2a_corrs/test/corr_vec_local_conn_{diagram}_{permkey}_seed{Nseeds}_{series}.{cfg}.p'
        self.Npoint        = 2
        self.timingPrint = False
        # self.diagrams = ['vec_local']
        # self.subdirs = ['']
        # self.outfilestem ='a2a_corrs/test/corr_sib_vec_local_conn_{permkey}_seed{Nseeds}_{series}.{cfg}.p'
        # self.Npoint = 3

    def loadMeson(self,keyl,keyr,gamma,time=()):

        with h5py.File(self.mesonfilestem.format(w=keyl,v=keyr,g=gamma,series=self.series,cfg=self.cfg),"r") as f:
            a_group_key = list(f.keys())[0]
            return np.array(f[a_group_key]['a2aMatrix'][time].view(np.complex64),
                            dtype=np.complex128)
    
    def contract(self,m1,m2,m3=None,m4=None,open_indices=(0,-1)):
        # two-point functions
        if m3 is None:
            Cij = cp.array(oe.contract('imn,jnm->ij',m1,m2))
        # three-point functions
        elif m4 is None:
            index_list = ['i','j','k']
            out_indices = "".join(index_list[i] for i in open_indices)
            Cij = oe.contract(f'ilm,jmn,knl->{out_indices}',m1,m2,m3)
        else:
            index_list = ['i','j','k','l']
            out_indices = "".join(index_list[i] for i in open_indices)
            Cij = oe.contract(f'imn,jno,kop,lpm->{out_indices}',m1,m2,m3,m4)
        return Cij

    def time_average(self,Cij):
        dims = len(Cij.shape)
        if dims > 4:
            raise Exception("> 4-pt functions not yet implemented")

        C = cp.zeros((self.Nt,),dtype=np.complex128)
        ones = cp.ones(Cij.shape)
        trange = cp.array(range(self.Nt))
        # start_tuple = tuple(None if i != (corr_indices[0] % Nt) else slice(None) for i in range(dims))
        # end_tuple = tuple(None if i != (corr_indices[1] % Nt) else slice(None) for i in range(dims))
        start_tuple = (slice(None),None)
        end_tuple = (None,slice(None))
        tmask = cp.mod(trange[start_tuple]*ones-trange[end_tuple]*ones,cp.array([self.Nt]))

        # submask = tmask >= 0 #Always True

        # before_tuple = start_tuple
        # for n in range(1,dims-1):
            # current_tuple = tuple(None if i != n else slice(None) for i in range(dims))
            # after_tuple = tuple(None if i != (n+1) else slice(None) for i in range(dims))

            # Assigns True to normally ordered elements: t1<t2<t3 and False otherwise
            # if n == 1:
                # submask = cp.logical_xor(trange[before_tuple]*ones <= trange[current_tuple]*ones,
                            # cp.logical_xor(trange[before_tuple]*ones <= trange[after_tuple]*ones,
                                            # trange[current_tuple]*ones <= trange[after_tuple]*ones))
            # else:
                # submask = cp.logical_not(cp.logical_xor(submask,
                                            # cp.logical_xor(trange[before_tuple]*ones <= trange[current_tuple]*ones,
                                                # cp.logical_xor(trange[before_tuple]*ones <= trange[after_tuple]*ones,
                                                                # trange[current_tuple]*ones <= trange[after_tuple]*ones))))
            # before_tuple=current_tuple

        for t in range(self.Nt):
            C[t] = Cij[tmask==t].sum()
            # C[t] = Cij[cp.logical_and(tmask==t,submask)].sum() + Cij[cp.logical_and(tmask==t, cp.logical_not(submask))].sum()
            # for t1 in range(Nt):
            #     temp = Cij
            #     t2 = (t1+t) % Nt
            #     if len(Cij.shape) == 3:
            #         if t1 <= t2:
            #             temp = cp.subtract(cp.sum(temp[:,t1:t2+1],axis=1),cp.add(cp.sum(temp[:,t2:Nt],axis=1),cp.sum(temp[:,0:t1],axis=1)))
            #         else:
            #             temp = cp.subtract(cp.sum(temp[:,t2:t1],axis=1),cp.add(cp.sum(temp[:,t1:Nt],axis=1),cp.sum(temp[:,0:t2],axis=1)))
            #     C[t] += temp[t1,t2]    

        #return C/Nt
        return cp.asnumpy(C/self.Nt)

    def shiftMass(self,mat):
        with h5py.File(self.evalfilestem.format(series=self.series,cfg=self.cfg),'r') as f:
            try:
                evals = np.array(f['/evals'][()].view(np.float64),dtype=np.float64)
            except KeyError:
                evals = np.array(f['/EigenValueFile/evals'][()].view(np.float64),dtype=np.float64)

        evals = np.sqrt(evals)
        mass_old = float(f'0.{self.massold}')
        mass_new = float(f'0.{self.massnew}')

        eval_scaling = np.zeros((len(evals),2),dtype=np.complex128)
        eval_scaling[:,0] = np.divide(mass_old+1.j*evals,mass_new+1.j*evals)
        eval_scaling[:,1] = np.conjugate(eval_scaling[:,0])
        eval_scaling = eval_scaling.reshape((-1,))

        mat[:] = np.multiply(mat,eval_scaling[np.newaxis,np.newaxis,:])
        
    def vec_conn_2pt(self,seedlist,local):

        corr = {}

        shift_mass = False
        if self.massnew and self.massold and self.massnew != self.massold:
            shift_mass = True
            
        gammas = [f"G{g}_G{g}" if local else f"G{g}_G1" for g in ['X','Y','Z']]
        #gammas = [f"Gamma{g}" if local else f"G{g}_G1" for g in ['X','Y','Z']]

        for gamma in gammas:

            mat1 = self.loadMeson(seedlist[0],seedlist[1],gamma)
            if shift_mass and seedlist[1] == self.elabel:
                self.shiftMass(mat1)
            mat2 = self.loadMeson(seedlist[2],seedlist[3],gamma)
            if shift_mass and seedlist[3] == self.elabel:
                self.shiftMass(mat2)

            corr[gamma] = self.time_average(self.contract(mat1,mat2))

        return corr

    def local_sib_conn_3pt(self,seedlist):

        corr = {}
        Cij = cp.zeros((self.Nt,self.Nt,self.Nt),dtype=np.complex128)

        local_gammas = [f"G{g}_G{g}" for g in ['X','Y','Z']]

        mat2 = self.loadMeson(seedlist[2],seedlist[3],"G1_G1")
            
        for gamma in local_gammas:
            mat1 = self.loadMeson(seedlist[0],seedlist[1],gamma)
            mat3 = self.loadMeson(seedlist[4],seedlist[5],gamma)

            Cij[:] = contract(m1=mat1,m2=mat2,m3=mat3)

            corr[gammakey] = cp.asnumpy(self.time_average(Cij))

        return corr

    def local_qed_conn_4pt(self,seedlist):

        Cij = cp.zeros((self.Nt,self.Nt,self.Nt,self.Nt),dtype=np.complex128)
        
        seedkey = "".join(seedlist)

        local_gammas = [f"G{g}_G{g}" for g in ['X','Y','Z']]

        corr = dict(zip(self.diagrams,[{seedkey:dict(zip(local_gammas,
                                                         [{}]*len(local_gammas)))}]*len(self.diagrams)))

        matg1 = [cp.asnumpy(self.loadMeson(seedlist[0],seedlist[1],gamma)) for gamma in local_gammas]
        matg2_photex = []
        matg2_selfen = []

        for gamma in local_gammas:
            if "photex" in self.diagrams:
                matg2_photex.append(cp.asnumpy(self.loadMeson(seedlist[4],seedlist[5],gamma)))
            if "selfen" in self.diagrams:
                matg2_selfen.append(cp.asnumpy(self.loadMeson(seedlist[6],seedlist[7],gamma)))

        for i in range(self.Nem):

            emlabel = f"{self.emseedstring}_{i}"

            selfen_p2_key = seedlist[4]+seedlist[5]+emlabel
            photex_p2_key = seedlist[6]+seedlist[7]+emlabel

            matp1 = self.loadMeson(seedlist[2],seedlist[3],emlabel)
            if "selfen" in self.diagrams:
                matp2_selfen = self.loadMeson(seedlist[4],seedlist[5],emlabel)
            if "photex" in self.diagrams:
                matp2_photex = self.loadMeson(seedlist[6],seedlist[7],emlabel)

            for j,gamma in enumerate(local_gammas):
                selfen_g2_key = seedlist[6]+seedlist[7]+gamma
                photex_g2_key = seedlist[4]+seedlist[5]+gamma

                for d in self.diagrams:

                    if d == "photex":
                        Cij[:] = self.contract(m1=cp.array(matg1[j]),m2=matp1,m3=cp.array(matg2_photex[j]),m4=matp2_photex,
                                            open_indices=(0,2))
                    elif d == "selfen":
                        Cij[:] = self.contract(m1=cp.array(matg1[j]),m2=matp1,m3=matp2_selfen,m4=cp.array(matg2_selfen[j]))
                    else:
                        raise Exception(f"Unrecognized diagram label: {d}")

                    corr[d][seedkey][gamma][emlabel] = cp.asnumpy(self.time_average(Cij))

        return corr

myContractor = Contractor(series, cfg)

def execute(index,seeds,diagrams):

    permkey = "".join(sum(((myContractor.perm[i], myContractor.perm[(i+1) % myContractor.Npoint]) for i in range(myContractor.Npoint)),()))

    with cp.cuda.Device(index):
        veclist = sum(tuple(zip([(mode,(myContractor.elabel if mode == 'L' else "wseed%i" % seeds[i])) for mode, i in zip(myContractor.perm,myContractor.highIndices)],
                                [(mode,(myContractor.elabel if mode == 'L' else "vseed%i" % seeds[i])) for mode, i in zip(myContractor.perm[1:]+myContractor.perm[:1],myContractor.highIndices[1:]+myContractor.highIndices[:1])])),
                      ())

        seedlist = [v[1] for v in veclist]
        seedkey = "".join(seedlist)

        corr = {permkey:dict(zip(diagrams,[{}]*len(diagrams)))}

        # if seedkey not in C[permkey][myContractor.diagrams[0]].keys():
        print(round(time.time()*1000),veclist)
        if "vec_local" in diagrams:
            corr[permkey]["vec_local"][seedkey] = myContractor.vec_conn_2pt(seedlist,local=True)
        if "vec_onelink" in diagrams:
            corr[permkey]["vec_onelink"][seedkey] = myContractor.vec_conn_2pt(seedlist,local=False)
        if "sib" in diagrams:
            corr[permkey]['sib'][seedkey] = myContractor.local_sib_conn_3pt(seedlist)

        if "selfen" in diagrams or "photex" in diagrams:
            corr[permkey].update(myContractor.local_qed_conn_4pt(seedlist))

        # cp.cuda.Device(index).synchronize()

        # i = barrier.wait()

        # if i == 0:
                # Mesons.clear()
        # barrier.wait()

        return corr

# main

Nmesons = myContractor.Npoint
Nseeds = myContractor.Nseeds
Nem = myContractor.Nem
seedstring = myContractor.seedstring
C = {}
local_gammas = [f"G{g}_G{g}" for g in ['X','Y','Z']]


for Nlow in range(Nmesons+1):

    if Nlow < 2:
       continue
    
    HLperms = multiset_permutations(['L']*Nlow+['H']*(Nmesons-Nlow))
    for perm in HLperms:

        myContractor.perm = perm

        permkey = "".join(sum(((perm[i], perm[(i+1) % Nmesons]) for i in range(Nmesons)),()))

        diagrams = copy.deepcopy(myContractor.diagrams)
        print(diagrams)
        for s,d in zip(myContractor.subdirs,myContractor.diagrams):
            outfile = myContractor.outfilestem.format(permkey=permkey,Nseeds=Nseeds,series=series,cfg=cfg,
                                                diagram=d,subdir=s.format(Nem=Nem))

            if os.path.exists(outfile):
                diagrams.remove(d)

        if permkey not in C.keys():
            C[permkey] = {}
        print(diagrams)

        myContractor.highIndices = [sum([0 if a =='L' else 1 for a in perm][:i]) for i in range(Nmesons)]
        p = Pool(processes)
        # Produces list of tuples like (('L','e1000'), ('H','wseed1'), ('H','vseed1'), ('L','e1000'))
        items = itertools.combinations(list(range(Nseeds)),Nmesons-Nlow)

        start_time = perf_counter( )

        _ = [ 0 if not bool(result[permkey]) else 
             C[permkey].update(result[permkey]) if diagram not in C[permkey] else 
             C[permkey][diagram].update(v1) if seed not in C[permkey][diagram] else
             C[permkey][diagram][seed].update(v2) if gamma not in C[permkey][diagram][seed] or type(C[permkey][diagram][seed][gamma]) is not dict else
             C[permkey][diagram][seed][gamma].update(v3)
                for result in p.starmap(execute,[(i % processes,item,diagrams) for i,item in enumerate(items)])
                for diagram,v1 in result[permkey].items()
                for seed,v2 in v1.items()
                for gamma, v3 in v2.items()]

        stop_time = perf_counter( )
        print('')
        print('    Elapsed wall clock time for IO+contraction = %g seconds.' % (stop_time - start_time) )
        print('')
        print(f"Finished {permkey}")

        for s,d in zip(myContractor.subdirs,diagrams):
            outfile = myContractor.outfilestem.format(permkey=permkey,Nseeds=Nseeds,series=series,cfg=cfg,
                                                diagram=d,subdir=s.format(Nem=Nem))
            pickle.dump(C[permkey][d],open(outfile,'wb'))

