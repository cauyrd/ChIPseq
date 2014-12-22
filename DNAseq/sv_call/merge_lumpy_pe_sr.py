#!/usr/bin/python
#-*- coding: utf-8 -*-
#===============================================================================
#
#         FILE: test.py
#
#        USAGE: ./test.py  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Rendong Yang (cauyrd@gmail.com), 
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: Wed Apr 30 13:27:17 CDT 2014
#     REVISION: ---
#===============================================================================
import pandas as pd
import sys
import os
try:
	pesr1 = pd.read_table(sys.argv[1]+'.pesr1.bed',header=None)
except:
	sys.exit(0)
try:
	pesr2 = pd.read_table(sys.argv[1]+'.pesr2.bed',header=None)
except:
	sys.exit(0)
pesr1_id = set(map(tuple,pesr1.loc[:,[3,15]].values))
pesr2_id = set(map(tuple,pesr2.loc[:,[3,15]].values))
common = pesr1_id.intersection(pesr2_id)
if not common:
	sys.exit(0)
pe = pd.read_table(sys.argv[1]+'.lumpy.pe.bedpe',header=None)
sr = pd.read_table(sys.argv[1]+'.lumpy.sr.bedpe',header=None)
combine = pd.DataFrame()
for each in common:
	pe_values = {6:[each[0]]}
	sr_values = {6:[each[1]]}
	a = pe[pe.isin(pe_values).any(1)]
	b = sr[sr.isin(sr_values).any(1)]
	a.index = range(1)
	b.index = range(1)
	combine = combine.append(pd.concat([a, b],ignore_index=True,axis=1))
combine.to_csv(sys.argv[1]+'.lumpy.pesr.bedpe',sep='\t',header=False,index=False)
