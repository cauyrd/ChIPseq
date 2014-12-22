#!/usr/bin/python
#-*- coding: utf-8 -*-
#===============================================================================
#
#         FILE: bwa2blat.py
#
#        USAGE: ./bwa2blat.py blatDir softclip_ratio input.bam output.bam
#
#  DESCRIPTION: the program is used to correct the CIGAR string with soft clippinb in BWA using BLAT mapping 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Rendong Yang (yang4414@umn.edu), 
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: Wed Nov 26 11:46:39 CST 2014
#===============================================================================
import sys
import os
import pysam
from Bio import SearchIO

def psl2sam(hsp,query_seq_len):
	"""psl2sam try to implement the psl2sam.pl script and return the cigar and mapping position estimated from psl file"""
	cigar = ''
	query_start = hsp.query_start
	query_end = hsp.query_end
	strand = hsp.query_strand_all[0] # may need replace by qery_strand
	soft_len = 0
	if strand == -1:
		query_start = query_seq_len - hsp.query_end
		query_end = query_seq_len - hsp.query_start
	if query_start:
		# 5'-end clipping
		soft_len = query_start
		cigar += str(query_start)+'S'
	x = hsp.query_span_all
	if strand == -1:
		y = [query_seq_len - item[1] for item in hsp.query_range_all] # may need r    eplace by query_start_all when the bug is fixed in Biopython
	else:
		y = [item[0] for item in hsp.query_range_all] # may need replace by query_start_all when the bug is fixed in Biopython
	z = hsp.hit_start_all
	y0, z0 = y[0], z[0]
	for i in range(1,len(hsp)):
		ly = y[i] - y[i-1] - x[i-1]
		lz = z[i] - z[i-1] - x[i-1]
		if ly < lz:
			# del: the reference gap is longer
			cigar += str(y[i] - y0)+'M'
			cigar += str(lz - ly)+'D'
			y0, z0 = y[i], z[i]
		elif lz < ly:
			# ins: the query gap is longer
			cigar += str(z[i] - z0)+'M'
			cigar += str(ly - lz)+'I'
			y0, z0 = y[i], z[i]
	cigar += str(query_end - y0)+'M'
	if (query_seq_len != query_end):
		# 3'-end clipping
		end3 = query_seq_len - query_end
		if end3 > soft_len:
			soft_len = end3
		cigar += str(end3)+'S'
	return cigar, soft_len

bwa_bam = pysam.Samfile(sys.argv[3],'rb')
blat_bam = pysam.Samfile(sys.argv[4],'wb', template=bwa_bam)
cutoff = float(sys.argv[2])

for read in bwa_bam.fetch():
	if read.is_secondary or read.is_unmapped:
		# secondary alignment or not mapped
		continue
	for each in read.cigar:
		if each[0] == 4 and each[1]/float(read.rlen) >= cutoff:
			# generating fasta file
			fa = open('temp.fa','w')
			print >> fa,'>'+read.qname
			print >> fa, read.seq
			fa.close()

			# run blat
			code = os.system('gfClient localhost 50000 '+sys.argv[1]+' temp.fa temp.psl >/dev/null 2>&1 ')
			if code != 0:
				print 'Execute gfClient failed!'
				sys.exit(1)
			try:
				blat = SearchIO.read('temp.psl','blat-psl')
			except:
				break
			hsps = blat.hsps
			hsps.sort(key=lambda k:k.score, reverse=True)
			if hsps[0].hit_id == bwa_bam.getrname(read.tid):
				# matching genomic coordinate
				if hsps[0].hit_start == read.pos or hsps[0].hit_end==read.aend:
					cigarstring, soft_len = psl2sam(hsps[0],blat.seq_len)
					if soft_len == 0:
						read.cigarstring, read.pos= cigarstring, hsps[0].hit_start
			break
	blat_bam.write(read)
bwa_bam.close()
blat_bam.close()
