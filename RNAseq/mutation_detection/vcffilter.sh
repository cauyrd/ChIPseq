#!/bin/bash -l
DPVALUE=20
vcffilter -f " DP > $DPVALUE & QUAL > 20 " wt_rep1.vcf > wt_rep1.filtered.vcf
vcffilter -f " DP > $DPVALUE & QUAL > 20 " wt_rep2.vcf > wt_rep2.filtered.vcf
vcffilter -f " DP > $DPVALUE & QUAL > 20 " wt_rep3.vcf > wt_rep3.filtered.vcf
vcffilter -f " DP > $DPVALUE & QUAL > 20 " pe_rep1.vcf > pe_rep1.filtered.vcf
vcffilter -f " DP > $DPVALUE & QUAL > 20 " pe_rep2.vcf > pe_rep2.filtered.vcf
vcffilter -f " DP > $DPVALUE & QUAL > 20 " pe_rep3.vcf > pe_rep3.filtered.vcf
vcffilter -f " DP > $DPVALUE & QUAL > 20 " pe_rep4.vcf > pe_rep4.filtered.vcf
vcffilter -f " DP > $DPVALUE & QUAL > 20 " pe_rep5.vcf > pe_rep5.filtered.vcf

