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
import vcf

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
	start = time()
	print 'ScanIndel starts running: '+strftime("%Y-%m-%d %H:%M:%S")
	sample = read_sample_file(sample_file)
	reference = read_config_file(config_file)	

	path = os.path.dirname(os.path.realpath(__file__))
	
	print "Start up BLAT server"
	cwd = os.getcwd()
	os.chdir(reference['blat'])
	os.system('gfServer start localhost 50000 *.2bit &')
	os.chdir(cwd)

	for each in sample:
		print "Analyzing sample:", each
		if not bam:
			print "Aapter_cutting"
			os.system("cutadapt -a file:/soft/trimmomatic/0.32/adapters/TruSeq3-PE-2.fa -o "+each+".r1.trimmed.fastq "+sample[each].split()[0])
			os.system("cutadapt -a file:/soft/trimmomatic/0.32/adapters/TruSeq3-PE-2.fa -o "+each+".r2.trimmed.fastq "+sample[each].split()[1])
			print "BWA-MEM start aligning raw reads..." 
			os.system("bwa mem -M -t8 "+reference['bwa']+" "+each+".r1.trimmed.fastq "+each+".r2.trimmed.fastq >"+each+".sam")
			os.system("samtools view -bS "+each+".sam >"+each+".bam")
			os.system("samtools sort "+each+".bam "+each+".sorted")
			os.system("samtools index "+each+".sorted.bam")
			os.remove(each+".sam")
			os.remove(each+".bam")
			blat_input = each+'.sorted.bam'
			if rmdup:
				os.system("samtools rmdup "+each+".sorted.bam "+each+".rmdup.bam")
				os.system("samtools index "+each+".rmdup.bam")
				os.remove(each+".sorted.bam")
				os.remove(each+".sorted.bam.bai")
				blat_input = each+'.rmdup.bam'
		else:
			blat_input = sample[each]
			os.system('samtools index '+blat_input)

		print "BLAT start re-aligning soft-clipped reads..."
		blat_start = time()
		os.system('python '+path+'/bwa2blat.py '+reference['blat']+' '+str(softclip_ratio)+' '+blat_input+' '+each+'.softclip_realigned.bam')
		os.system('samtools sort '+each+'.softclip_realigned.bam '+each+'.softclip_realigned.sorted')
		os.system('samtools index '+each+'.softclip_realigned.sorted.bam')
		os.remove(each+'.softclip_realigned.bam')

		# profiling BLAT running time
		blat_end = time()
		print 'Running BLAT using: '+str(blat_end-blat_start)+' seconds.'

		print "Freebayes start variant calling..."
		
		os.system('freebayes -F '+str(freebayes_F)+' -C '+str(freebayes_C)+' --min-coverage '+str(depth)+' -f '+reference['freebayes']+' '+each+'.softclip_realigned.sorted.bam|vcfallelicprimitives -k -g|vt normalize -r '+reference['freebayes']+' - > '+each+'.scanindel.vcf')

		if bedfile:
			os.system('intersectBed -a '+each+'.scanindel.vcf -b '+bedfile+' -wa -header -u >'+each+'.scanindel.target.vcf')

		# profiling FreeBayes running time
		freebayes_end = time()
		print 'Running FreeBayes using: '+str(freebayes_end-blat_end)+' seconds.'

	print "Cleanup sever"
	os.system('gfServer stop localhost 50000')
	print "ScanIndel running done: "+strftime("%Y-%m-%d %H:%M:%S")
	end = time()
	print 'Analyzing your data takes: '+str(end-start)+' seconds.'

if __name__ == '__main__':
	main()
