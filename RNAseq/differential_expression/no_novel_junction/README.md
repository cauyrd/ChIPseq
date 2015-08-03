Scripts
1. Initial fastqc on all files to check for the quality encoding version as well as quality issues, see /home/orrharry/shared/riss/wk28/Orr_sequence_initial_fastqc.pbs
2. First I double checked that the reads were properly synced (since there was trouble with these datasets when they were originally sequenced) and then I trimmed the reads with Trimmomatic, then did another fastqc to check that everything looks ok, see /home/orrharry/shared/riss/wk28/Orr_sequence_cleanup.pbs
3. Tophat/Cufflinks: I ran each sample separately, so there is a separate script for each, for one example, see /home/orrharry/shared/riss/wk28/tophat_cuffdiff/D30_rep1.pbs
4. Cuffquant & Cuffdiff: /home/orrharry/shared/riss/wk28/tophat_cuffdiff/Cuffquant_Cuffdiff_D30_B05_FVB.pbs 

Then I used cummeRbund in R to take a look at the data and pull out the significant genes, although you don't have to do it that way.
