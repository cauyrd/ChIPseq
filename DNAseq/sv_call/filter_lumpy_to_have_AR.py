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
import HTSeq
ifp = open(sys.argv[1]+'.pesr.bedpe')
ofp1 = open(sys.argv[1]+'.tier3.lumpy.pe.bedpe','w')
ofp2 = open(sys.argv[1]+'.tier3.lumpy.sr.bedpe','w')
ofp3 = open(sys.argv[1]+'.tier2.lumpy.pesr.bedpe','w')
ar_iv = HTSeq.GenomicInterval( "chrX", 66763874, 66950461, "." )
for line in ifp:
	items = line.rstrip().split()
	pos1 = HTSeq.GenomicInterval(items[0], int(items[1]), int(items[2]), '.')
	pos2 = HTSeq.GenomicInterval(items[3], int(items[4]), int(items[5]), '.')
	if ar_iv.contains(pos1) or ar_iv.contains(pos2):
		id = items[11].split(':')[1].split(';')
		if len(id) > 1: # support by both pairend and split reads
			print >> ofp3, line.rstrip()
		elif id[0].split(',')[0] == '1':
			print >> ofp1, line.rstrip()
		else:
			print >> ofp2, line.rstrip()
ofp1.close()
ofp2.close()
ofp3.close()
ifp.close()
