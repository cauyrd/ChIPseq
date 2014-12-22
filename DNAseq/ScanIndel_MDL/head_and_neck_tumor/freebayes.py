#!/usr/bin/python
#-*- coding: utf-8 -*-
#===============================================================================
#
#         FILE: ScanIndel.py
#
#        USAGE: ./ScanIndel.py -i sample.txt -c config.txt [opts]
#
#  DESCRIPTION: Indel detection for targeted NGS data
#
#      OPTIONS: ---
# REQUIREMENTS: See README.md file
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Rendong Yang (yang4414@umn.edu), 
# ORGANIZATION: 
#      VERSION: 1.0
#===============================================================================
import sys
import os
import re
from time import time, strftime
import getopt

def read_config_file(filename):
	path = {}
	ifp = open(filename)
	for line in ifp:
		if line[0] == '#' or line == '\n':
			continue
		item = line.rstrip().split('=')
		path[item[0]]=item[1]
	ifp.close()
	return path

def read_sample_file(filename):
	input = {}
	ifp = open(filename)
	for line in ifp:
		if line[0] == '#' or line == '\n':
			continue
		item = line.rstrip().split()
		input[item[0]] = ' '.join(item[1:])
	ifp.close()
	return input

def usage():
	"""helping information"""
	print 'Usage:'
	print ' python ScanIndel.py -p config.txt -i sample.txt [opts]'
	print 'Opts:'
	print ' -F  :setting min-alternate-fraction for FreeBayes (default 0.2)'
	print ' -C  :setting min-alternate-count for FreeBayes (default 2)'
	print ' -d  :setting min-coverage for Freebayes (default 0)'
	print ' -t  :setting --target for Freebayes to provide a BED file for analysis'
	print ' -s  :softclipping percentage triggering BLAT re-alignment (default 0.2)'
	print ' --bam  :the input file is BAM format'
	print ' --rmdup  :exccute duplicate removal step before realignment'
	print ' -h --help :produce this menu'
	print 'author: Rendong Yang <yang4414@umn.edu>, MSI, University of Minnesota, 2014'
	print 'version: 1.0'

def main():

	# parameters parsing
	sample_file = 'sample.txt'
	config_file = 'config.txt'
	freebayes_F = 0.2
	freebayes_C = 2
	softclip_ratio = 0.2
	depth = 0
	bam = False
	rmdup = False
	bedfile = ''
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'i:p:F:C:s:d:t:h', ['bam','rmdup','help'])
	except getopt.GetoptError as err:
		print str(err)
		usage()
		sys.exit(2)
	
	for o,a in opts:
		if o == '-i':
			sample_file = a
		elif o == '-p':
			config_file = a
		elif o == '-F':
			freebayes_F = float(a)
		elif o == '-C':
			freebayes_C = int(a)
		elif o == '-s':
			softclip_ratio = float(a)
		elif o == '-d':
			depth = int(a)
		elif o == '-t':
			bedfile = a
		elif o == '-v':
			snp_cutoff = float(a)
		elif o == '--bam':
			bam = True
		elif o == '--rmdup':
			rmdup = True
		elif o in ('-h', '--help'):
			usage()
			sys.exit(0)
		else:
			assert False, "unhandled option"

	# timing progrom start
	sample = read_sample_file(sample_file)
	reference = read_config_file(config_file)	
	for each in sample:
		print "Freebayes start variant calling..."
		os.system('freebayes -F '+str(freebayes_F)+' -C '+str(freebayes_C)+' --min-coverage '+str(depth)+' -f '+reference['freebayes']+' '+each+'.softclip_realigned.sorted.bam|vcfallelicprimitives -k -g|vt normalize -r '+reference['freebayes']+' - > '+each+'.scanindel.vcf')

		if bedfile:
			os.system('intersectBed -a '+each+'.scanindel.vcf -b '+bedfile+' -wa -header -u >'+each+'.scanindel.target.vcf')

if __name__ == '__main__':
	main()
