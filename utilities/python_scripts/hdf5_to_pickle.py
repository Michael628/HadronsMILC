#! /usr/bin/env python3

import os
import h5py
import shutil
import re
import pickle
from pathlib import Path
import numpy as np
import sys
import yaml
from processing_utils import *

params = yaml.safe_load(open('params.yaml','r'))['hdf5_to_pickle']

outfile_template = params['outfile']['dict']

purge_list = sys.argv[1:]
audit = False

if len(purge_list) > 0 and purge_list[0] == 'check':
    audit = True
    purge_list = purge_list[1:]
    
Nt = int(params['nt'])
dt = int(params['dt'])
Nslices = int(Nt/dt)

for infile in params['infiles']:

    split_time_files = infile['split_time_files'] if 'split_time_files' in infile else True
    print(split_time_files)
    
    for dset in infile['dsets']:
        for gamma in infile['gammas']:

            glabel,gindices = gamma['label'],gamma['indices']
            outfile_name = outfile_template.format(dset=dset,gamma=glabel)
            if not os.path.exists(os.path.dirname(outfile_name)):
                os.makedirs(os.path.dirname(outfile_name))

            # Determine appropriate output filename and load if already exists
            try:
                corr_data = pickle.load(open(outfile_name,"rb"))

                for purge_cfg in purge_list:
                    if purge_cfg in corr_data:
                        corr_data.pop(purge_cfg)
            except FileNotFoundError:
                corr_data = {}
            
            if dset == 'a2aLL':
                pattern = re.compile("^.*\.([0-9]+)\.h5$")
                configs = set()

                for series in [ s[-1] for s in os.listdir(infile['series_path'].format(dset=dset,gamma=glabel))]:

                    infile_path = infile['path'].format(dset=dset,gamma=glabel,series=series)
                    files = os.listdir(infile_path)

                    # Gather list of configs
                    for infile_base in files:

                        re_match = pattern.match(infile_base) # Grab series.cfg from filename
                        config = re_match[1]
                        seriescfg = f'{series}.{config}'
                        
                        configs.add(seriescfg)
                    
            else:
                if split_time_files:
                    pattern = re.compile("(.*)_t([0-9]+)_([a-z].[0-9]+)(.*)")
                    infile_path = infile['path'].format(dset=dset,gamma=glabel)
                    configs = set()
                    files = os.listdir(infile_path)

                    # Gather list of configs
                    for infile_base in filter(lambda x: True if "_t0" in x else False, files):

                        re_match = pattern.match(infile_base) # Grab series.cfg from filename
                        seriescfg = re_match[3]

                        configs.add(seriescfg)
                    print(configs)
                else:
                    pattern = re.compile(".*_([a-z].[0-9]+).h5$")
                    infile_path = infile['path'].format(dset=dset,gamma=glabel)
                    configs = set()
                    files = os.listdir(infile_path)

                    # Gather list of configs
                    for infile_base in files:
                        re_match = pattern.match(infile_base) # Grab series.cfg from filename
                        seriescfg = re_match[1]

                        configs.add(seriescfg)
                
            for seriescfg in configs:
                temp_corr = None

                if seriescfg in corr_data.keys() and not audit:
                    print(f"Complete: {seriescfg} already in output file '{outfile_name}'")
                    continue

                series,cfg = seriescfg.split(".")
                
                if dset == "a2aLL":

                    if len(gindices) == 0:
                        gindices = [0]
                        
                    for gindex,gdir in enumerate(gindices):

                        meson_label = gamma['meson_label'].format(gamma=gdir)
                        
                        infile_path = infile['path'].format(dset=dset,gamma=glabel,series=series)

                        infile_name = f"{infile_path}/" + infile['infile_label'].format(cfg=cfg,meson=meson_label)
            
                        if temp_corr is None:
                            temp_corr = np.zeros((len(gindices),Nt),dtype=np.complex128)

                        skip_cfg = False
                        try:
                            f = h5py.File(infile_name,"r")
                            temp_corr[gindex,:] = np.array(f[f'/{meson_label}/correlator'][()].view(np.complex128),dtype=np.complex128) 
                            f.close()
                        except FileNotFoundError:
                            print(f"Skipping {seriescfg} {glabel} {dset}: could not find file {infile_name}")
                            skip_cfg = True
                            break

                        except KeyError:
                            print(f"Failed to find meson with label '{meson_label}'. Trying 'vec5_0_vec5', assuming we're looking for pion data")
                            try:
                                temp_corr[gindex,:] = np.array(f[f'/vec5_0_vec5/correlator'][()].view(np.complex128),dtype=np.complex128)
                            except KeyError:
                                print(f"Skipping {seriescfg} {glabel} {dset}: could not get file contents {infile_name}")
                                skip_cfg = True
                                break
                            f.close()
                    
                    if skip_cfg:
                        continue
                else:
                    cfg_files = list(filter(lambda x: True if f"{seriescfg}.h5" in x else False, files))

                    if  split_time_files and len(cfg_files) != Nslices:
                        if not audit:
                            print(f"Skipping {seriescfg} {glabel} {dset}: expected {Nslices} time slices. Found {len(cfg_files)}")
                            continue
            
                    if not audit:
                        print(f"Processing: {seriescfg} {glabel} {dset}")

                    for infile_base in cfg_files:
                        if split_time_files:
                            tslice = int(pattern.search(infile_base)[2])
                            tindex = int(tslice/dt)

                            if temp_corr is None:
                                temp_corr = np.zeros((len(gindices),Nslices,Nt),dtype=np.complex128)

                            infile_name = f"{infile_path}/{infile_base}"

                            print(infile_name)
                            f = h5py.File(infile_name,"r")
                            temp_corr[...,tindex,:] = [np.array(grp.view(np.complex128),dtype=np.complex128) for grp in [f[f"/meson/meson_{k}/corr"][()] for k in gindices]]
                            f.close()
                        else:
                            if temp_corr is None:
                                temp_corr = np.zeros((len(gindices),Nslices,Nt),dtype=np.complex128)

                            infile_name = f"{infile_path}/{infile_base}"

                            print(infile_name)
                            f = h5py.File(infile_name,"r")
                            temp_corr[:] = [np.array(grp.view(np.complex128),dtype=np.complex128) for grp in [f[f"/meson/meson_{k}/srcCorrs"][()] for k in gindices]]
                            f.close()
                            

                if audit:
                    file_match = np.array_equal(corr_data[seriescfg],temp_corr if len(gindices) > 1 else temp_corr[0])
                    if not file_match:
                        print(f"Mismatch: {seriescfg}")
                    else:
                        print(f"Matching: {seriescfg}")
                else:
                    corr_data[seriescfg] = temp_corr if len(gindices) > 1 else temp_corr[0]

            if not audit:
                print(f"Writing: {outfile_name}")
                pickle.dump(corr_data,open(outfile_name,"wb"))
                if 'numpy' in params['outfile']:
                    pickle.dump(dictToCorr(corr_data),open(params['outfile']['numpy'].format(dset=dset,gamma=glabel),"wb"))
                    
