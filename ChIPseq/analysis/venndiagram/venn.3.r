args <- commandArgs(TRUE)
#par(mar = c(2.5, 2.5, 1.6, 1.1), mgp = c(1.5, 0.5, 0), mfrow=c(1,1))
library(VennDiagram)
aa = read.table(args[1])[,1]
bb = read.table(args[2])[,1]
cc = read.table(args[3])[,1]
x = list(aa,bb,cc)
names(x) = args[4:6]
venn.diagram(
			 x = x,
             euler.d = TRUE,
             filename = paste(args[7], "tiff", sep="."), 
             cex = 1.3,
             cat.cex = 1.3,
             cat.pos = c(-20,20,180),
			 fill = c("red", "green", "blue"),
			 col = "transparent",
			 margin = 0.2,
			 fontfamily = "serif",
			 cat.fontfamily = "serif",
   			 cat.dist = 0.05,
  			 alpha = 0.70
             )

