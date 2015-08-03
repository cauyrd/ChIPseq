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
import pybedtools
import os
max_num = 10
for dp in range(50,5000,50):
	print 'dp = ',dp,
	os.system('vcffilter -f " DP > '+str(dp)+' " wt_rep1.filtered.vcf.no_dbSNP.vcf > test1.vcf')
	os.system('vcffilter -f " DP > '+str(dp)+' " wt_rep2.filtered.vcf.no_dbSNP.vcf > test2.vcf')
	os.system('vcffilter -f " DP > '+str(dp)+' " wt_rep3.filtered.vcf.no_dbSNP.vcf > test3.vcf')
	wt1 = pybedtools.BedTool('test1.vcf')
	wt2 = pybedtools.BedTool('test2.vcf')
	wt3 = pybedtools.BedTool('test3.vcf')
	wt1_2 = len(wt1-wt2)
	wt1_3 = len(wt1-wt3)
	wt2_3 = len(wt2-wt3)
	print 'wt1-wt2 = ',wt1_2,';wt1-wt2 = ',wt1_3,';wt2-wt3 = ',wt2_3
	if wt1_2 < max_num and wt1_3 < max_num and wt2_3 < max_num:
		print 'found depth', dp
		break
print 'finished.'
