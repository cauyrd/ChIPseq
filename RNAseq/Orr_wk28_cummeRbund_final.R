#If you haven't installed cummeRbund yet:
#source("http://bioconductor.org/biocLite.R")
#biocLite("cummeRbund")

library(cummeRbund)
setwd("/home/orrharry/shared/riss/wk12/tophat_cuffdiff/cuffdiff_output")

#Read in cuffdiff output
cuff<-readCufflinks(rebuild=T, gtf="/panfs/roc/rissdb/igenomes/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-current/Genes/genes.gtf",genome="/panfs/roc/rissdb/genomes/Mus_musculus/mm10_canonical/seq/mm10_canonical.fa")

#Basic quality control plots
#NOTE: all of these are done at the gene level, but can be repeated at other levels
#Plot dispersion (model fit)
#d<-dispersionPlot(genes(cuff))
#d
##Looking for outliers
#pBoxRep<-csBoxplot(genes(cuff),replicates=T)
#pBoxRep
#pDendro<-csDendro(genes(cuff),replicates=T)
#pDendro
#pBox<-csBoxplot(genes(cuff))
#pBox
#
##Density plots to look at FPKM distribution across samples/replicates
#dens<-csDensity(genes(cuff))
#dens+coord_cartesian(xlim = c(-3,6))
#densRep<-csDensity(genes(cuff),replicates=T)
#densRep+coord_cartesian(xlim = c(-3,6)) 
#
##Pairwise scatterplot
#s<-csScatterMatrix(genes(cuff))+coord_cartesian(xlim = c(-3,6), ylim = c(-3,6)) 
#
##Comparing the CV of the FPKM across samples
#gene.scv<-fpkmSCVPlot(genes(cuff))
#isoforms.scv<-fpkmSCVPlot(isoforms(cuff))
#
##Volcano plots
#v<-csVolcanoMatrix(genes(cuff))
#v+coord_cartesian(xlim = c(-6,6)) 
##MA plots
#m1<-MAplot(genes(cuff),"D30","B05")
#m2<-MAplot(genes(cuff),"D30","FVB")
#m3<-MAplot(genes(cuff),"B05","FVB")




#GENE LEVEL
#Differentially expressed genes

#Rough idea of results
mySigMatgenes<-sigMatrix(cuff,level="genes",alpha=0.05)

#Genes differentially expressed overall (between any pair, all pairs considered together)
sigGeneIds<-getSig(cuff,alpha=0.05,level="genes")
 head(sigGeneIds)
 length(sigGeneIds)
 
#Genes DE between D30 and B05 samples
D30vsB05.sigGeneIds<-getSig(cuff,"D30","B05",alpha=0.05,level="genes")
head(D30vsB05.sigGeneIds)
length(D30vsB05.sigGeneIds)

#Genes DE between D30 and FVB samples
D30vsFVB.sigGeneIds<-getSig(cuff,"D30","FVB",alpha=0.05,level="genes")
head(D30vsFVB.sigGeneIds)
length(D30vsFVB.sigGeneIds)

#Genes DE between B05 and FVB samples
B05vsFVB.sigGeneIds<-getSig(cuff,"B05","FVB",alpha=0.05,level="genes")
head(B05vsFVB.sigGeneIds)
length(B05vsFVB.sigGeneIds)

#Write significant gene lists for each comparison
write.table(D30vsB05.sigGeneIds, file="D30vsB05.sigGeneIds.txt",quote=F, col.names=F, row.names=F, sep="\t")
write.table(D30vsFVB.sigGeneIds, file="D30vsFVB.sigGeneIds.txt",quote=F, col.names=F, row.names=F, sep="\t")
write.table(B05vsFVB.sigGeneIds, file="B05vsFVB.sigGeneIds.txt",quote=F, col.names=F, row.names=F, sep="\t")


#Pull the log fold change, p, q-values, etc. from the gene expression data for each list of significantly DE genes
gene.diff<-diffData(genes(cuff))
write.table(gene.diff[gene.diff$gene_id %in% D30vsB05.sigGeneIds & gene.diff$sample_1=="D30" & gene.diff$sample_2=="B05",], "D30_vs_B05_significant_genes.txt",quote=F, sep="\t",row.names=F, col.names=T)
write.table(gene.diff[gene.diff$gene_id %in% D30vsFVB.sigGeneIds & gene.diff$sample_1=="D30" & gene.diff$sample_2=="FVB",], "D30_vs_FVB_significant_genes.txt",quote=F, sep="\t",row.names=F, col.names=T)
write.table(gene.diff[gene.diff$gene_id %in% B05vsFVB.sigGeneIds & gene.diff$sample_1=="B05" & gene.diff$sample_2=="FVB",], "B05_vs_FVB_significant_genes.txt",quote=F, sep="\t",row.names=F, col.names=T)


#ISOFORM LEVEL
#Differentially expressed isoforms

#Rough idea of results
mySigMatiso<-sigMatrix(cuff,level="isoforms",alpha=0.05)

sigGeneIds<-getSig(cuff,alpha=0.05,level="isoforms")
head(sigGeneIds)
length(sigGeneIds)

