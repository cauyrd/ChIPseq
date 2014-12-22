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

def subsample(filename,cutoff):
	os.system('zcat '+filename+'|wc > wc.txt')
	ifp = open('wc.txt')
	count = int(ifp.readline().split()[0])/4
	if count > cutoff:
		return True
	else:
		return False

def group_variant(raw_vcffile,indel_vcffile,snp_vcffile,cutoff):
	"""splitting raw vcf to indel and snp vcfs"""
	rawvcf = vcf.Reader(open(raw_vcffile))
	indel_vcf = vcf.Writer(open(indel_vcffile,'w'), rawvcf)
	snp_vcf = vcf.Writer(open(snp_vcffile,'w'), rawvcf)
	for record in rawvcf:
		if record.is_indel:
			indel_vcf.write_record(record)
		elif record.is_snp and record.INFO['AO'][0]/float(record.INFO['DP']) >= cutoff:
			snp_vcf.write_record(record)

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
	print ' -F  :setting min-alternate-fraction for FreeBayes (default 0.005)'
	print ' -C  :setting min-alternate-count for FreeBayes (default 2)'
	print ' -d  :setting min-coverage for Freebayes (default 100)'
	print ' -t  :setting --target for Freebayes to provide a BED file for analysis'
	print ' -s  :softclipping percentage triggering BLAT re-alignment (default 0.2)'
	print ' -n  :the number of read pairs cutoff triggering downsampling (default 1e6)'
	print ' -v  :SNP vaf cutoff (default 0.05)'
	print ' --cutadapt  :run adapter trimming for V3 data'
	print ' -h --help :produce this menu'
	print 'author: Rendong Yang <yang4414@umn.edu>, MSI, University of Minnesota, 2014'
	print 'version: 1.0'

def main():

	# parameters parsing
	config_file = 'config.txt'
	sample_file = 'sample.txt'
	freebayes_F = 0.005
	freebayes_C = 2
	softclip_ratio = 0.2
	depth = 100 
	cutadapt = False
	bedfile = ''
	max_sample_size = 2e6
	snp_cutoff = 0.05
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'v:p:i:F:C:s:d:t:n:h', ['cutadapt','help'])
	except getopt.GetoptError as err:
		print str(err)
		usage()
		sys.exit(2)
	
	for o,a in opts:
		if o == '-p':
			config_file = a
		elif o == '-i':
			sample_file = a
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
		elif o == '-n':
			max_sample_size = float(a)
		elif o == '-v':
			snp_cutoff = float(a)
		elif o == '--cutadapt':
			cutadapt = True
		elif o in ('-h', '--help'):
			usage()
			sys.exit(0)
		else:
			assert False, "unhandled option"

	reference = read_config_file(config_file)	
	sample = read_sample_file(sample_file)
	path = os.path.dirname(os.path.realpath(__file__))
	cwd = os.getcwd()
	os.chdir(reference['blat'])
	os.system('gfServer start localhost 50000 *.2bit &')
	os.chdir(cwd)

	# subsampling
	for each in sample:
		read1,read2 = sample[each].split()
		if subsample(read1,max_sample_size):
			os.system('seqtk sample -s100 '+read1+' '+str(max_sample_size/2)+' >sub_r1.fq')
			os.system('seqtk sample -s100 '+read2+' '+str(max_sample_size/2)+' >sub_r2.fq')
			read1 = 'sub_r1.fq'
			read2 = 'sub_r2.fq'
		if cutadapt:
			os.system("cutadapt -a AGACCAAGTCTCTGCTACCGTA -o sub_r1_trimmed.fq "+read1)
			os.system("cutadapt -a TGTAGAACCATGTCGTCAGTGT -o sub_r2_trimmed.fq "+read2)
			read1 = 'sub_r1_trimmed.fq'
			read2 = 'sub_r2_trimmed.fq'
		os.system("bwa mem -M -t8 "+reference['bwa']+" "+read1+" "+read2+" >bwa.sam")
		os.system("samtools view -bS bwa.sam > bwa.bam")
		os.system("samtools sort bwa.bam bwa.sorted")
		os.system("samtools index bwa.sorted.bam")
		blat_input = 'bwa.sorted.bam'

		os.system('python '+path+'/script/bwa2blat.py '+reference['blat']+' '+str(softclip_ratio)+' bwa.sorted.bam blat.bam')
		os.system('samtools sort blat.bam blat.sorted')
		os.system('mv blat.sorted.bam '+each+'.bam')
		os.system('samtools index '+each+'.bam')

		if bedfile:
			os.system('freebayes -F '+str(freebayes_F)+' -C '+str(freebayes_C)+' --min-coverage '+str(depth)+' -f '+reference['freebayes']+' -t '+bedfile+' '+each+'.bam|vcfallelicprimitives -k -g|vcfbreakmulti|vt normalize -r '+reference['freebayes']+' - > freebayes.raw.vcf')
		else:
			os.system('freebayes -F '+str(freebayes_F)+' -C '+str(freebayes_C)+' --min-coverage '+str(depth)+' -f '+reference['freebayes']+' '+each+'.bam|vcfallelicprimitives -k -g|vcfbreakmulti|vt normalize -r '+reference['freebayes']+' - > freebayes.raw.vcf')
		group_variant('freebayes.raw.vcf',each+'.indel.vcf',each+'.snp.vcf',snp_cutoff)

	os.system('gfServer stop localhost 50000')

if __name__ == '__main__':
	main()
