Introduction
============
The rename.py program is used to rename and index bam files in a specific folder.

Input
=====
an input file named 'list.txt' is needed to list the names of sample (folder) to be renamed. An example input file is included in the same folder with rename.py

Output
======
following files at sample_name/data folder will be renamed with sample_name_file_name.

gene_panel_cvrg.tab
run_count
all_gene_cvrg.tab
full_bam.bam (bam file will be indexed first)
gene_panel.vcf
plot.pdf
sample_sheet.tab

Running the program
===================
Two modes:
1. python rename.py  (this program will rename the files at /home/thyagara/shared/ngsdp_new/projects/)

2. python rename.py  yohes (this program will rename the files at /home/yohes/shared/ngsdp_new/)

Questions?
==========
Contact Rendong Yang (yang4414@umn.edu)
