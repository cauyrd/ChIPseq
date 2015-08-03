#!/usr/bin/python
#-*- coding: utf-8 -*-
#===============================================================================
#
#         FILE: heatmap.py
#
#        USAGE: ./heatmap.py bamfile bedfile halfwinwidth binwidth 
#
#  DESCRIPTION: generate the CDT file for Java TreeView 
#
#      OPTIONS: ---
# REQUIREMENTS: sort the bedfile based on read coverage accordingly before running
#         BUGS: ---
#        NOTES: bamfile is the reads mapping file, bedfile is the enrichment genomics positions (e.g. peak summit or TSS)
#       AUTHOR: Rendong Yang (cauyrd@gmail.com), 
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: Mon May 12 14:04:49 CDT 2014
#     REVISION: ---
#===============================================================================

import HTSeq
import pysam
import pybedtools
import numpy
import sys 
import getopt

def usage():
	"""showing help information"""
	print 'Usage:'
	print 'python heatmap.py -i input.bam -b enrichment_region.bed -o <output_file_prefix> [opts]'
	print 'Opts:'
	print ' -w <int>	:half window size extension from center  (default:3000)'
	print ' -l <int>	:bin length (default:30)'
	print ' -h      	:produce this menu'

def get_hm_bin(coverage, window, totalbins):
	hm_list = []
	"""this is to cut promoter into totalbins bins and count the coverage in each bin"""
	bins = numpy.linspace(window.start, window.end, totalbins+1)
	for i in range(totalbins):
		bin_range = HTSeq.GenomicInterval(window.chrom, int(bins[i]), int(bins[i+1]), '.')
		count = sum(numpy.fromiter(coverage[bin_range], dtype='i', count=(int(bins[i+1])-int(bins[i]))))
		hm_list.append(count)
	return hm_list

def main():
	# parameters parsing.
	halfwinwidth = 3000
	binwidth = 30
	bamfile = None
	bedfile = None
	cdtfile = None
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'i:b:o:w:l:h')
	except getopt.GetoptError as err:
		print str(err)
		usage()
		sys.exit(2)
	for o,a in opts:
		if o == '-w':
			halfwinwidth = int(a)
		elif o == '-l':
			binwidth = int(a)
		elif o == '-i':
			bamfile = a
		elif o == '-b':
			bedfile = a
		elif o == '-o':
			cdtfile = a
		elif o == '-h':
			usage()
			sys.exit()
		else:
			assert False, "unhandled option"
	if not bamfile or not bedfile:
		usage()
		sys.exit(2)

	# creat genomic coverage profile
	coverage = HTSeq.GenomicArray('auto',stranded=False,typecode='i')
	#almnt = HTSeq.BAM_Reader(bamfile)
	almnt = pysam.Samfile(bamfile,'rb')
	for read in almnt.fetch():
		if not read.is_unmapped:
			iv = HTSeq.GenomicInterval(almnt.getrname(read.tid),read.pos, read.aend,'.')
			coverage[ iv ] += 1

	# generate the cdt file
	bin_num = halfwinwidth * 2 / binwidth
	label= ['bin' + str(i) for i in range(bin_num)]
	colname = 'UNIQID' + '\t' + 'NAME' + '\t' + '\t'.join(label)
	ofp = open(cdtfile+'.cdt', 'w')
	print >> ofp, colname
	regions = pybedtools.BedTool(bedfile) 
	for pos in regions:
		window = HTSeq.GenomicInterval(pos.chrom, pos.start-halfwinwidth, pos.start+halfwinwidth, '.')
		wincvg = get_hm_bin(coverage,window,bin_num)
		if pos.strand == '+':
			print >> ofp, '_'.join(pos) + '\t' + bamfile + '\t' + '\t'.join(map(str, wincvg))
		else:
			print >> ofp, '_'.join(pos) + '\t' + bamfile + '\t' + '\t'.join(map(str, wincvg[::-1]))

if __name__ == '__main__':
	main()
