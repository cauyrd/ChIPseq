#!/usr/bin/python
#-*- coding: utf-8 -*-
#===============================================================================
#
#         FILE: tssplot.py
#
#        USAGE: ./tssplot.py sample.txt bedfile halfwinwidth binwidth 
#
#  DESCRIPTION: generate the TSS plot for different samples
#
#      OPTIONS: ---
# REQUIREMENTS: run reorder.py before running this program
#         BUGS: ---
#        NOTES: bamfile is the reads mapping file, bedfile is the enrichment genomics positions (e.g. peak summit or TSS)
#       AUTHOR: Rendong Yang (cauyrd@gmail.com), 
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: May 14 2014
#     REVISION: ---
#===============================================================================

import HTSeq
import pybedtools
import numpy
import sys 
import getopt
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

def usage():
	"""showing help information"""
	print 'Usage:'
	print 'python heatmap.py -i sample.txt -b enrichment_region.bed -o <output_file_prefix> [opts]'
	print 'Opts:'
	print ' -w <int>	:half window size extension from center  (default:3000)'
	print ' -l <int>	:bin length (default:30)'
	print ' -h      	:produce this menu'
	print ' Note: run reorder.py before run this program!'

def get_hm_bin(coverage, window, totalbins, totalreads):
	hm_list = []
	"""this is to cut promoter into totalbins bins and count the coverage in each bin"""
	bins = numpy.linspace(window.start, window.end, totalbins+1)
	for i in range(totalbins):
		bin_range = HTSeq.GenomicInterval(window.chrom, int(bins[i]), int(bins[i+1]), '.')
		count = numpy.log1p(sum(numpy.fromiter(coverage[bin_range], dtype='i', count=(int(bins[i+1])-int(bins[i]))))/float(totalreads)*1e6)
		hm_list.append(count)
	return hm_list

def main():
	# parameters parsing.
	halfwinwidth = 3000
	binwidth = 30
	samplefile = None
	bedfile = None
	output = None
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
			samplefile = a
		elif o == '-b':
			bedfile = a
		elif o == '-o':
			output = a
		elif o == '-h':
			usage()
			sys.exit()
		else:
			assert False, "unhandled option"
	if not samplefile or not bedfile:
		usage()
		sys.exit(2)
	
	bin_num = halfwinwidth * 2 / binwidth
	x = numpy.linspace(-halfwinwidth, halfwinwidth, bin_num)
	regions = pybedtools.BedTool(bedfile) 
	
	# read sample file which lists the data in BAM format
	ifp = open(samplefile)
	for line in ifp:
		mapped_read_num = 0
		items = line.rstrip().split()

		# creat genomic coverage profile
		coverage = HTSeq.GenomicArray('auto',stranded=False,typecode='i')
		almnt = HTSeq.BAM_Reader(items[0])
		for read in almnt:
			if read.aligned:
				mapped_read_num += 1
				coverage[ read.iv ] += 1

		# create scaning window for visualization
		profile = numpy.zeros(bin_num, dtype='i')
		for pos in regions:
			window = HTSeq.GenomicInterval(pos.chrom, pos.start-halfwinwidth, pos.start+halfwinwidth, '.')
			wincvg = numpy.fromiter( get_hm_bin(coverage,window,bin_num,mapped_read_num), dtype='i' )
			if pos.strand == '+':
				profile += wincvg
			else:
				profile += wincvg[::-1]
		y = profile/float(len(regions))
		pyplot.plot(x, y, items[1]+'-', lw=2, label=items[2])
	pyplot.legend(prop={'size':8}, loc='best')
	pyplot.axvline(0, linestyle=':', color='k')
	pyplot.xlabel('Distance form center (bp)')
	pyplot.ylabel('log normalized average read coverage\n(per million mapped reads)')
	pyplot.savefig(output+'.png')


if __name__ == '__main__':
	main()
