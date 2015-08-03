#!/usr/bin/python
#-*- coding: utf-8 -*-
#===============================================================================
#
#         FILE: test.py
#
#        USAGE: ./test.py vcf_file1 vcf_file2 vcf_file3 
#
#  DESCRIPTION: select the variants in vcf_file1 but not in vcf_file2 and vcf_file3
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
import vcf
import sys
myset = set()
vcf1 = vcf.Reader(open(sys.argv[1]))
vcf2 = vcf.Reader(open(sys.argv[2]))
uniq_vcf = vcf.Writer(open(sys.argv[1]+'.no_dbSNP.vcf','w'),vcf1)
for record in vcf2:
	if 'chr' not in record.CHROM:
		a = ['chr'+record.CHROM]
	else:
		a = [record.CHROM]
	a.append(record.POS)
	a.append(record.REF)
	a.extend(record.ALT)
	myset.add('_'.join(map(str,a)))
for  record in vcf1:
	if 'chr' not in record.CHROM:
		a = ['chr'+record.CHROM]
	else:
		a = [record.CHROM]
	a.append(record.POS)
	a.append(record.REF)
	a.extend(record.ALT)
	item = '_'.join(map(str,a))
	if item in myset:
		continue
	else:
		uniq_vcf.write_record(record)
