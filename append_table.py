#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 18:56:32 2016
combine all the existing tsv data
@author: sig
"""
# %% get actual dirs
import os 
import pandas as pd
os.chdir("/home/sig/project/nocoding")
dirs = list(filter(os.path.isdir, os.listdir()))
empty = [] # empty with a verbose path
for folder in dirs:
    folder = os.path.join("/home/sig/project/nocoding", folder) # add father directory
    if  not os.listdir(folder):
        empty.append(folder) 

empty_dirs = [one_dir.split("/")[-1] for one_dir in empty]
actual_dirs = list(set(dirs).difference(set(empty_dirs)))
# %%
final_big_tsv = pd.DataFrame()
small_dirs = actual_dirs[0:2]
# %%
for one_dir in actual_dirs:
    one_dir_raw = os.path.join("/home/sig/project/nocoding", one_dir)
    tmp_tsv = os.path.join(one_dir, os.listdir(one_dir_raw)[0])
    tmp_tsv = pd.read_csv(tmp_tsv, sep = "\t", index_col = 0)
    tmp_wgs = tmp_tsv[tmp_tsv["sequencing_strategy"] == "WGS"]
final_big_tsv = final_big_tsv.append(tmp_wgs)
