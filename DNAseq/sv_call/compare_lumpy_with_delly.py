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
import sys
ifp1 = open('merge1.bed')
ifp2 = open('merge2.bed')
set1 = set()
set2 = set()
id_dict = {}
for line in ifp1:
	item = line.rstrip().split()
	set1.add((item[3],item[-2]))
for line in ifp2:
	item = line.rstrip().split()
	set2.add((item[3],item[-2]))
common = set1.intersection(set2)
ifp1.close()
ifp2.close()
for each in common:
	id_dict[each[0]] = each[1]
ifp = open(sys.argv[1]+'.tier2.lumpy.pesr.bedpe')
ofp = open(sys.argv[1]+'.tier1.lumpy_delly.bedpe','w')
for line in ifp:
	item = line.rstrip().split()
	if item[6] in id_dict:
		print >> ofp, line.rstrip()+'\t'+id_dict[item[6]]
