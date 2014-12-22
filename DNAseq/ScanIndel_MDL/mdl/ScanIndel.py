#!/usr/bin/python
#-*- coding: utf-8 -*-
#===============================================================================
#
#         FILE: ScanIndel.py
#
#        USAGE: ./ScanIndel.py -h
#
#  DESCRIPTION: Indel and SNP detection for targeted NGS data
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
import getopt
import vcf

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

def subsample(filename,cutoff):
	os.system('zcat '+filename+'|wc > wc.txt')
	ifp = open('wc.txt')
	count = int(ifp.readline().split()[0])/4
	if count > cutoff:
		return True
	else:
		return False

def usage():
	"""helping information"""
	print 'Usage:'
	print ' python ScanIndel.py -p config.txt -t targeted.bed -l left_read.fastq.gz -r right_read.fastq.gz -b bam_output -x indel_vcf_output -y snp_vcf_output  [opts]'
	print 'Opts:'
	print ' -F  :setting min-alternate-fraction for FreeBayes (default 0.005)'
	print ' -C  :setting min-alternate-count for FreeBayes (default 2)'
	print ' -d  :setting min-coverage for Freebayes (default 100)'
	print ' -s  :softclipping percentage triggering BLAT re-alignment (default 0.2)'
	print ' -n  :the number of read pairs cutoff triggering downsampling to half (default 2e6)'
	print ' -v  :SNP vaf cutoff (default 0.05)'
	print ' -h --help :produce this menu'
	print 'author: Rendong Yang <yang4414@umn.edu>, MSI, University of Minnesota, 2014'
	print 'version: 1.0'

def main():

	# parameters parsing
	config_file = ''
	freebayes_F = 0.005
	freebayes_C = 2
	softclip_ratio = 0.2
	depth = 100 
	bedfile = ''
	bamfile = ''
	indel_vcffile = ''
	snp_vcffile = ''
	r1_file = ''
	r2_file = ''
	max_sample_size = 2e6
	snp_cutoff = 0.05
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'p:l:r:b:x:y:F:C:s:d:t:n:h', ['help'])
	except getopt.GetoptError as err:
		print str(err)
		usage()
		sys.exit(2)
	
	for o,a in opts:
		if o == '-p':
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
		elif o == '-b':
			bamfile = a
		elif o == '-n':
			max_sample_size = float(a)
		elif o == '-v':
			snp_cutoff = float(a)
		elif o == '-l':
			r1_file = a
		elif o == '-r':
			r2_file = a
		elif o == '-x':
			indel_vcffile = a
		elif o == '-y':
			snp_vcffile = a
		elif o in ('-h', '--help'):
			usage()
			sys.exit(0)
		else:
			assert False, "unhandled option"
	if not bamfile or not indel_vcffile or  not snp_vcffile:
		print "Please provide the bam or vcf outfile name."
		usage()
		sys.exit(1)
	if not r1_file or not r2_file:
		print "Please provide the input R1 or R2 fastq file."
		usage()
		sys.exit(1)
	if not bedfile:
		print "Please provide bedfile for regions of interest"
		usage()
		sys.exit(1)
	if not config_file:
		print "Please provide config file showing the PATH of BWA, BLAT and Freebayes index "
		usage()
		sys.exit(1)

	reference = read_config_file(config_file)	
	path = os.path.dirname(os.path.realpath(__file__))
	cwd = os.getcwd()
	os.chdir(reference['blat'])
	code = os.system('gfServer start localhost 50000 *.2bit &')
	if code != 0:
		print 'Execute gfSever start failed!'
		sys.exit(1)
	os.chdir(cwd)

	# subsampling
	if subsample(r1_file,max_sample_size):
		code = os.system('seqtk sample -s100 '+r1_file+' '+str(max_sample_size/2)+' >sub_r1.fq')
		if code != 0:
			print 'Execute seqtk for read1 failed!'
			sys.exit(1)
		code = os.system('seqtk sample -s100 '+r2_file+' '+str(max_sample_size/2)+' >sub_r2.fq')
		if code != 0:
			print 'Execute seqtk for read1 failed!'
			sys.exit(1)
		r1_file = 'sub_r1.fq'
		r2_file = 'sub_r2.fq'

	code = os.system("cutadapt -a AGACCAAGTCTCTGCTACCGTA -o sub_r1_trimmed.fq "+r1_file)
	if code != 0:
		print 'Execute cutadapt for read1 failed!'
		sys.exit(1)
	code = os.system("cutadapt -a TGTAGAACCATGTCGTCAGTGT -o sub_r2_trimmed.fq "+r2_file)
	if code != 0:
		print 'Execute cutadapt for read2 failed!'
		sys.exit(1)
	r1_file = 'sub_r1_trimmed.fq'
	r2_file = 'sub_r2_trimmed.fq'

	code = os.system("bwa mem -M -t8 "+reference['bwa']+" "+r1_file+" "+r2_file+" >bwa.sam")
	if code != 0:
		print 'Execute BWA mem failed!'
		sys.exit(1)
	code = os.system("samtools view -bS bwa.sam > bwa.bam")
	if code != 0:
		print 'Execute samtools view failed!'
		sys.exit(1)
	code = os.system("samtools sort bwa.bam bwa.sorted")
	if code != 0:
		print 'Execute samtools sort failed!'
		sys.exit(1)
	code = os.system("samtools index bwa.sorted.bam")
	if code != 0:
		print 'Execute samtools index failed!'
		sys.exit(1)
	blat_input = 'bwa.sorted.bam'

	code = os.system('python '+path+'/script/bwa2blat.py '+reference['blat']+' '+str(softclip_ratio)+' bwa.sorted.bam blat.bam')
	if code != 0:
		print 'Execute bwa2blat.py failed!'
		sys.exit(1)
	os.system('samtools sort blat.bam blat.sorted')
	os.system('mv blat.sorted.bam '+bamfile)
	os.system('samtools index '+bamfile)
	os.system('gfServer stop localhost 50000')

	code = os.system('freebayes -F '+str(freebayes_F)+' -C '+str(freebayes_C)+' --min-coverage '+str(depth)+' -f '+reference['freebayes']+' -t '+bedfile+' '+bamfile+'|vcfallelicprimitives -k -g|vcfbreakmulti|vt normalize -r '+reference['freebayes']+' - > freebayes.raw.vcf')
	if code != 0:
		print 'Execute Freebayes failed!'
		sys.exit(1)
	group_variant('freebayes.raw.vcf',indel_vcffile,snp_vcffile,snp_cutoff)

if __name__ == '__main__':
	main()
