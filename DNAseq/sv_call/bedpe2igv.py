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
import os
import sys
color_map = {'TYPE:DELETION':('del','25,210,20'),'TYPE:DUPLICATION':('dup','110,0,110'), 'TYPE:INTERCHROM':('tra','255,0,0'), 'TYPE:INVERSION':('inv','0,0,255')}
if os.stat(sys.argv[1]+'.tier1.lumpy_delly.bedpe')[6] == 0:
	sys.exit(0)
else:
	ifp = open(sys.argv[1]+'.tier1.lumpy_delly.bedpe')
	ofp = open(sys.argv[1]+'.igv.bed','w')
	for line in ifp:
		items = line.rstrip().split()
		try:
			sv_type = color_map[items[10]][0]
			sv_color = color_map[items[10]][1]
		except KeyError:
			print items[10]+' not in color map'
			sys.exit(1)
		name1 = sv_type+','+items[3]+':'+items[4]+'-'+items[5]
		name2 = sv_type+','+items[0]+':'+items[1]+'-'+items[2]
		print >> ofp, '\t'.join(items[0:3])+'\t'+name1+'\t'+items[7]+'\t'+items[8]+'\t'+items[1]+'\t'+items[2]+'\t'+sv_color
		print >> ofp, '\t'.join(items[3:6])+'\t'+name2+'\t'+items[7]+'\t'+items[9]+'\t'+items[4]+'\t'+items[5]+'\t'+sv_color
	ifp.close()
	ofp.close()

