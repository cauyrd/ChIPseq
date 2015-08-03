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
id_table = {}
ifp = open('cuffmerge_gtf/merged.gtf')
for line in ifp:
	item = line.rstrip().split('\t')[-1].split(';')
	for each in item:
		if 'gene_id' in each:
			gene_id = each.split()[1].split('"')[1]
		if 'transcript_id' in each:
			transcript_id = each.split()[1].split('"')[1]
		if 'gene_name' in each:
			gene_name = each.split()[1].split('"')[1]
		if 'oId' in each:
			oId = each.split()[1].split('"')[1]
	if gene_id not in id_table:
		id_table[gene_id] = gene_name
	if transcript_id not in id_table:
		id_table[transcript_id] = oId
ifp.close()
ifp = open(sys.argv[1])
count_up = 0
count_down = 0
name = sys.argv[2]
header = ifp.readline().rstrip()
ofp1 = open(name+'.up.txt','w')
print >> ofp1, header
ofp2 = open(name+'.down.txt','w')
print >> ofp2, header

for line in ifp:
	item = line.rstrip().split()
	item[0] = id_table[item[0]]
	item[1] = id_table[item[1]]
	if item[-1] == 'yes' and item[6] == 'OK':
		if float(item[9]) < 0:
			count_down += 1
			print >> ofp1, '\t'.join(item)
		elif float(item[9]) > 0:
			count_up += 1
			print >> ofp2, '\t'.join(item)
print name+'\t'+str(count_up)+'\t'+str(count_down)
	
