#!/usr/bin/python
#-*- coding: utf-8 -*-
#===============================================================================
#
#         FILE: test.py
#
#        USAGE: ./test.py  
#
#  DESCRIPTION: input is the sample name
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
import vcf
sv_dict = {'DEL':'DELETION','DUP':'DUPLICATION','INV':'INVERSION','TRA':'INTERCHROM'}
#ifp = vcf.Reader(open(sys.argv[1]+'.del.vcf'))
#ofp = vcf.Writer(open(sys.argv[1]+'.treated.vcf','w'),ifp)
ofp = open(sys.argv[1]+'.delly.bedpe','w')
ar_iv = HTSeq.GenomicInterval( "chrX", 66763874, 66950461, "." )
for each in sv_dict.keys():
	ifp = vcf.Reader(open(sys.argv[1]+'.'+each.lower()+'.vcf'))
	for record in ifp:
		pos1 = HTSeq.GenomicPosition(record.CHROM,record.POS,'.')
		pos2 = HTSeq.GenomicPosition(record.INFO['CHR2'],record.INFO['END'],'.')
		if ar_iv.contains(pos1) or ar_iv.contains(pos2):
			#ofp.write_record(record)
			if record.var_subtype == 'TRA' and pos1.chrom > pos2.chrom:
				pos2,pos1 = pos1,pos2 
			print >> ofp, pos1.chrom+'\t'+str(pos1.start)+'\t'+str(pos1.end)+'\t'+pos2.chrom+'\t'+str(pos2.start)+'\t'+str(pos2.end)+'\t'+record.ID+'\tTYPE:'+sv_dict[record.var_subtype]
