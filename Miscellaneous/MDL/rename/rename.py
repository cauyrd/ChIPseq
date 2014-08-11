#!/usr/bin/python
#-*- coding: utf-8 -*-
#===============================================================================
#
#         FILE: rename.py
#
#        USAGE: ./rename.py  
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
os.system('module load samtools')
if len(sys.argv) == 2:
	parent_path = '/home/yohes/shared/ngsdp_new/'
else:
	parent_path = '/home/thyagara/shared/ngsdp_new/projects/'
ifp = open('list.txt')
file_list = ['gene_panel_cvrg.tab','run_count','all_gene_cvrg.tab','full_bam.bam','gene_panel.vcf','plot.pdf','sample_sheet.tab']
for line in ifp:
	dirname = line.rstrip()
	child_path = parent_path+dirname+'/data/'
	os.chdir(child_path)
	for each in file_list:
		os.system('mv '+each+' '+dirname+'_'+each)
		if 'bam' in each:
			os.system('samtools index '+dirname+'_'+each)
ifp.close()
print "done for renameing."
