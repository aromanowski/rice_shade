library(tidyverse)
library(edgeR)

cg <- read.delim(file = "../data/Master_counts.txt", sep = "\t")
samples <- read_tsv("../data/samples.txt")

# make samplesID consistent with counts names
table(gsub("\\.", "-", colnames(cg)) == samples[,1]) 
rownames(samples) <- colnames(cg) <- gsub("\\.", "-", colnames(cg))

design <- model.matrix(~rep + treatment * variety, data=samples)
rownames(design) <- colnames(cg)
head(design)
colnames(design)

# normalize the counts (TMM) and estimeate variation and fit the model to the data
dge <- DGEList(counts=cg,
               group=paste(samples$variety, samples$treatment, sep="."),
               genes=rownames(cg))

dge <- dge[filterByExpr( y = dge), , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge)

aa <- plotMDS(dge, gene.selection = "common", top=2000)
mdsdata <- aa$cmdscale.out

par(mar=c(4, 5, 3, 9))
plot(aa, cex=1.5,
     col=rep(c("blue", "red", "purple", "dark green", "brown", "orange"), each=8),
     pch=rep(c(16,17), each=4, time=6))
