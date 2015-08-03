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
def is_mapped_to_ar(read):
	ar_pos = ('chrX',66762874,66951461)
	if read.tid == 22: # tid 22 is chrX
		if read.pos < ar_pos[2] and read.pos > ar_pos[1]:
			return True
		else:
			return False
	return False
def soft_clip_ratio(read):
	if read.cigar[0][0] == 4:
		if read.cigar[-1][0] == 4:
			if read.cigar[0][1] > read.cigar[-1][1]:
				soft_len = read.cigar[0][1]
			else:
				soft_len = read.cigar[-1][1]
		else:
			soft_len = read.cigar[0][1]
	elif read.cigar[-1][0] == 4:
		soft_len = read.cigar[-1][1]
	else:
		soft_len = 0
	return soft_len/float(read.rlen)

mapq_cutoff = 40
soft_clip_ratio_cutoff = 0.1
# selecting candidate reads from GMAP mapping of pacbio reads
samfile = pysam.Samfile(sys.argv[1],'r')
outfile = pysam.Samfile('temp.bam','wb',template=samfile)
for read in samfile.fetch():
	if read.mapping_quality >= mapq_cutoff and is_mapped_to_ar(read):
		if soft_clip_ratio(read) > soft_clip_ratio_cutoff:
			outfile.write(read)
outfile.close()
# bam to fasta 
os.system('python bam_to_fasta.py temp.bam')
# mapping fasta to reference
reference = '/panfs/roc/rissdb/genomes/Homo_sapiens/hg19/bwa/hg19.fa'
os.system('bwa mem -x pacbio '+reference+' temp.fa >temp.sam')
os.system('samtools view -bS temp.sam >temp.bam')
os.system('samtools sort temp.bam temp.sorted')
os.system('samtools index temp.sorted.bam')
# predict fusion breakpoints
os.system('python select_chimeric_read.py temp.sorted.bam '+sys.argv[1])
#os.system('rm temp*')
