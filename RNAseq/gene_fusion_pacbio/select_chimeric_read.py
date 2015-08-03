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
import pysam
import sys
import os
from itertools import combinations
import pybedtools
def select_pair(pair,ext_width):
	ar = pybedtools.BedTool('chrX 66763874 66950461',from_string=True)
	if ar.any_hits(pair[0]) and ar.any_hits(pair[1]):
		return ''
	if not ar.any_hits(pair[0]) and not ar.any_hits(pair[1]):
		return ''
	depth_cutoff = 5 
	left_read_names = set(pair[0].name.split(';'))
	right_read_names = set(pair[1].name.split(';'))
	common_reads = left_read_names & right_read_names
	if len(common_reads) < depth_cutoff:
		return ''
	left_start = pair[0].start + ext_width
	left_end = pair[0].end - ext_width
	left_chr = pair[0].chrom
	right_start = pair[1].start + ext_width
	right_end = pair[1].end - ext_width
	right_chr = pair[1].chrom
	left_pos = left_chr+'\t'+str(left_start)+'\t'+str(left_end)
	right_pos = right_chr+'\t'+str(right_start)+'\t'+str(right_end)
	return left_pos+'\t'+right_pos+'\t'+';'.join(list(common_reads))+'\t'+str(len(common_reads))

read_set = set()
mapq_cutoff = 40
ext_width = 1000
readfile = pysam.AlignmentFile(sys.argv[1],'r')
ofp = open('temp.bed','w')
for read in readfile.fetch():
	if read.query_name not in read_set and read.mapping_quality >= mapq_cutoff:
		try:
			alt_pos = read.opt('SA').split(';')[:-1]
			read_set.add(read.query_name)
			print >> ofp, readfile.getrname(read.tid)+'\t'+str(read.pos-ext_width)+'\t'+str(read.pos+ext_width)+'\t'+read.query_name
			for each in alt_pos:
				items = each.split(',')
				if items[4] >= mapq_cutoff:
					print >> ofp, items[0]+'\t'+str(int(items[1])-ext_width)+'\t'+str(int(items[1])+ext_width)+'\t'+read.query_name
		except KeyError:
			continue
ofp.close()
#sort and merge bed file
os.system('sortBed -i temp.bed >temp.sorted.bed')
os.system('mergeBed -i temp.sorted.bed -nms >temp.breakpoints.bed')
bedfile = pybedtools.BedTool('temp.breakpoints.bed')
all_pairs = combinations(bedfile,2)
ofp = open(sys.argv[2]+'.pairs.bedpe','w')
for pair in all_pairs:
	bedpe_string = select_pair(pair,ext_width)
	if bedpe_string:
		print >> ofp, bedpe_string
ofp.close()
