Introduction
------------
ScanIndel finds indels (insertions and deletions) smaller than the length of a short read by re-align soft clipped reads. The workflow of ScanIndel is BWA-mem -> BLAT remap softclipped read -> FreeBayes for indel calling. In addition, ScanIndel also predict SNPs (single-nucleotide polymorphisms) in a seperate VCF output.

Pre-installation
----------------
Softwares:
* SAMtools/0.1.18 (http://samtools.sourceforge.net/)
* BWA/0.7.10 (http://bio-bwa.sourceforge.net/) 
* BLAT/34 [gfServer and gfClient] (http://genome.ucsc.edu/FAQ/FAQblat.html)
* freebayes/0.9.16 (https://github.com/ekg/freebayes)
* vcflib [vcfallelicprimitives and vcfbreakmulti] (https://github.com/ekg/vcflib) 
* vt/0.5 (https://github.com/atks/vt)
* seqtk/1.0-r63-dirty (https://github.com/lh3/seqtk)
* cutadapt/1.7dev (https://code.google.com/p/cutadapt/)
* pysam/0.7.7 (https://code.google.com/p/pysam/)
* pyvcf/0.6.7 (https://github.com/jamescasbon/PyVCF)
* biopython/1.64 (http://biopython.org/wiki/Main_Page)



All softwares above are assumed to be installed in your searching path. Ask your admistrator for assistance if necessary. 

Running ScanIndel
-----------------
### command-line usage
	python ScanIndel.py -t common_file/targeted.bed -p common_file/config.txt -l example/read1.fastq.gz -r example/read2.fastq.gz -b output.bam -x output.indel.vcf -y output.snp.vcf [options]
#### Options:
	 -F				:min-alternate-fraction for FreeBayes (default 0.2)
	 -C				:setting min-alternate-count for FreeBayes (default 2)
	 -s  			:softclipping percentage triggering BLAT re-alignment (default 0.2)
	 -d  			:minimal sequencing depth to indetify variants (default 100)
	 -n  			:the number of read pairs cutoff triggering downsampling to half (default 2e6)
	 -v				:SNP vaf cutoff (default 0.05)
	 -h --help 		:produce this menu
#### Input:
	config.txt    	:this file contains the path of reference file for each BWA, BLAT and Freebayes (default name is config.txt)
	targeted.bed	:bed format file for regions to scan indel and snp
	*.fastq.gz		:raw read file
#### Output:
The output files include one VCF file for indel, one VCF file for SNP and BAM files for BWA-MEM and BLAT mapping.

Example Data
------------
The folder example contains the examples of input fastq file, VCF and BAM output by ScanIndel and running ScanIndel pbs script.
The folder common_file contains the examples of config.txt and targeted.bed.
