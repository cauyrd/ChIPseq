# de novo motif search
module load bedtools
slopBed -i stat5.pairend_summits.bed -b 250 -g /panfs/roc/rissdb/genomes/Mus_musculus/mm10/seq/mm10.len >stat5.pairend_summits_motif.bed
slopBed -i stat5.R1_summits.bed -b 250 -g /panfs/roc/rissdb/genomes/Mus_musculus/mm10/seq/mm10.len >stat5.R1_summits_motif.bed
fastaFromBed -fi /panfs/roc/rissdb/genomes/Mus_musculus/mm10/seq/mm10.fa -bed stat5.pairend_summits_motif.bed -fo stat5.pairend.macs2.fa
fastaFromBed -fi /panfs/roc/rissdb/genomes/Mus_musculus/mm10/seq/mm10.fa -bed stat5.R1_summits_motif.bed -fo stat5.R1.macs2.fa
# summary of venndiagram 
# pairend uniq 1079, common 4829, R1 uniq 315, pairend total 5908, R1 total 5172.

