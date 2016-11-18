#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 18:41:09 2016
make dirs and move *.tsv.gz to dirs
@author: sig
"""
# %%
import os
file = open("./dirs.txt", "r")
dirs = []
while(1):
   line = file.readline()
   if line == "":
      break
   dirs.append(line.strip().split(" ")[0])
file.close()
# %%
os.chdir("/home/sig/project/nocoding") # change working directory
for one_dir in dirs:
   os.mkdir(one_dir) # make directories
   
# %% get gz files
os.chdir("/home/sig/project")
entities = os.listdir()
gzs = []
for entity in entities:
   if entity.endswith("tsv.gz"):
      gzs.append(entity)
# %% move files
import shutil
for gz in gzs:
   dist = gz.split(".")[2]
   shutil.move(gz, "./nocoding/" + dist)
