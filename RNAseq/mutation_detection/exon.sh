DP_VALUE=1000
vcffilter -f " DP > $DP_VALUE " pe_rep1_not_in_wt.vcf.no_dbSNP.vcf >pe_rep1_not_in_wt.vcf.no_dbSNP.filter.vcf
grep -v '#' pe_rep1_not_in_wt.vcf.no_dbSNP.filter.vcf|wc
intersectBed -a pe_rep1_not_in_wt.vcf.no_dbSNP.filter.vcf -b ~/mm10_exon_merged.bed -wa -wb >pe_rep1_not_in_wt.vcf.no_dbSNP.exon.vcf

vcffilter -f " DP > $DP_VALUE " pe_rep2_not_in_wt.vcf.no_dbSNP.vcf >pe_rep2_not_in_wt.vcf.no_dbSNP.filter.vcf
grep -v '#' pe_rep2_not_in_wt.vcf.no_dbSNP.filter.vcf|wc
intersectBed -a pe_rep2_not_in_wt.vcf.no_dbSNP.filter.vcf -b ~/mm10_exon_merged.bed -wa -wb >pe_rep2_not_in_wt.vcf.no_dbSNP.exon.vcf

vcffilter -f " DP > $DP_VALUE " pe_rep3_not_in_wt.vcf.no_dbSNP.vcf >pe_rep3_not_in_wt.vcf.no_dbSNP.filter.vcf
grep -v '#' pe_rep3_not_in_wt.vcf.no_dbSNP.filter.vcf|wc
intersectBed -a pe_rep3_not_in_wt.vcf.no_dbSNP.filter.vcf -b ~/mm10_exon_merged.bed -wa -wb >pe_rep3_not_in_wt.vcf.no_dbSNP.exon.vcf

vcffilter -f " DP > $DP_VALUE " pe_rep4_not_in_wt.vcf.no_dbSNP.vcf >pe_rep4_not_in_wt.vcf.no_dbSNP.filter.vcf
grep -v '#' pe_rep4_not_in_wt.vcf.no_dbSNP.filter.vcf|wc
intersectBed -a pe_rep4_not_in_wt.vcf.no_dbSNP.filter.vcf -b ~/mm10_exon_merged.bed -wa -wb >pe_rep4_not_in_wt.vcf.no_dbSNP.exon.vcf

vcffilter -f " DP > $DP_VALUE " pe_rep5_not_in_wt.vcf.no_dbSNP.vcf >pe_rep5_not_in_wt.vcf.no_dbSNP.filter.vcf
grep -v '#' pe_rep5_not_in_wt.vcf.no_dbSNP.filter.vcf|wc
intersectBed -a pe_rep5_not_in_wt.vcf.no_dbSNP.filter.vcf -b ~/mm10_exon_merged.bed -wa -wb >pe_rep5_not_in_wt.vcf.no_dbSNP.exon.vcf
