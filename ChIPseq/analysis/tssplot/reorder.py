#!/usr/bin/python
#-*- coding: utf-8 -*-
#===============================================================================
#
#         FILE: reorder.py
#
#        USAGE: ./reorder.py bamfile bedfile genome.size halfwinwidth filetype 
#
#  DESCRIPTION: sort the region bed file based on read coverage for TSS or peak score for peak summit bed file 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: TSS bed format: chr start end name strand; peak summit from macs2 format: chr start end name score.
#       AUTHOR: Rendong Yang (cauyrd@gmail.com), 
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: May 13 2014
#     REVISION: ---
#===============================================================================

import HTSeq
import pybedtools
import numpy
import sys 
import getopt
from operator import itemgetter
from pybedtools.featurefuncs import extend_fields

def usage():
	"""showing help information"""
	print 'Usage:'
	print 'python reorder.py -i input.bam -b enrichment_region.bed -g genome.size -f filetype (peak, tss) [opts]'
	print 'Opts:'
	print ' -w <int>	:half window size extension from center  (default:3000)'
	print ' -h      	:produce this menu'

def main():
	# parameters parsing.
	halfwinwidth = 3000
	bamfile = None
	bedfile = None
	genomefile = None
	filetype = None
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'i:f:b:g:w:h')
	except getopt.GetoptError as err:
		print str(err)
		usage()
		sys.exit(2)
	for o,a in opts:
		if o == '-w':
			halfwinwidth = int(a)
		elif o == '-i':
			bamfile = a
		elif o == '-f':
			filetype = a
		elif o == '-b':
			bedfile = a
		elif o == '-g':
			genomefile = a
		elif o == '-h':
			usage()
			sys.exit()
		else:
			assert False, "unhandled option"
	if not genomefile or not bedfile or not filetype:
		usage()
		sys.exit(2)

	# read genome size file
	genomesize = {}
	ifp = open(genomefile)
	for line in ifp:
		item = line.rstrip().split()
		genomesize[item[0]] = int(item[1])
	ifp.close()

	regions = pybedtools.BedTool(bedfile)
	ofp = open(bedfile+'.sorted.bed','w')

	# sort peak summit file
	if filetype == 'peak':
		peaklist = []
		for peak in regions:
			if peak.start - halfwinwidth < 0 or peak.end + halfwinwidth > genomesize[peak.chrom]:
				continue
			item = peak
			item = extend_fields(item,6)
			item.strand = '+'
			peaklist.append((float(item.score), item))
		# sort peak based on socre large --> small
		sortpeaklist = sorted(peaklist,key=itemgetter(0),reverse=True)
		for each in sortpeaklist:
			print >> ofp, '\t'.join(each[1])

	# sort tss file
	elif filetype == 'tss':
		# read bam file for mapped reads
		tsslist = []
		coverage = HTSeq.GenomicArray('auto',stranded=False,typecode='i')
		bamfile = HTSeq.BAM_Reader(bamfile)
		for read in bamfile:
			if read.aligned:
				coverage[ read.iv ] += 1
		for tss in regions:
			if tss.start - halfwinwidth < 0 or tss.end + halfwinwidth > genomesize[tss.chrom]:
				continue
			item = tss 
			item = extend_fields(item,6)
			item.strand = item.score
			window = HTSeq.GenomicInterval(tss.chrom, tss.start-halfwinwidth, tss.start+halfwinwidth, '.')
			count = sum( numpy.fromiter( coverage[window], dtype='i', count=2*halfwinwidth ))
			item.score = str(count)
			tsslist.append((float(item.score), item))
		# sort tss based on read coverage strong --> weak
		sorttsslist = sorted(tsslist,key=itemgetter(0),reverse=True)
		for each in sorttsslist:
			print >> ofp, '\t'.join(each[1])
	ofp.close()

if __name__ == '__main__':
	main()
