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
import os
ifp = open('bamfile.txt')
path = '/home/dehms/shared/chenzler/dnaseq/'
reference = '/panfs/roc/rissdb/genomes/Homo_sapiens/hg19/seq/hg19.fa'
for line in ifp:
	item = line.rstrip().split('/')
	sample = item[0]
	os.system('samtools mpileup -f '+reference+' -l AR.exon_100bp.bed -d 10000 -L 10000 -q 5 '+path+line.rstrip()+' >'+sample+'.mpileup')
	os.system('java -Xmx2g -jar ~/bin/VarScan.jar mpileup2snp '+sample+'.mpileup --output-vcf 1 --min-var-freq 0.01 | awk \'/^#/ || $7=="PASS"\' >'+sample+'.snp.vcf')
	os.system('java -Xmx2g -jar ~/bin/VarScan.jar mpileup2indel '+sample+'.mpileup --output-vcf 1 --min-var-freq 0.01 | awk \'/^#/ || $7=="PASS"\' >'+sample+'.indel.vcf')
