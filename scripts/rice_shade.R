############################################################################
# rice_shade.R                                                             #
# --------------                                                           #
# This script analyses gene expression of O. sativa plants grown           #
# under high R:FR or a low R:FR treatment. The analysis looks at 6         #
# rice varieties.                                                          #
# It was used to analyse the data described in Huber et al., 2023          #
#                                                                          #
############################################################################


### Notes
# Many lines of codes are commented. Uncomment them as needed for your usage

###############################################
#             Requires                        #
###############################################
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install()}
if(!require(ggpubr)) {
  if(!require(devtools))install.packages("devtools")
  devtools::install_github("kassambara/ggpubr")}
if(!require(edgeR)) BiocManager::install("edgeR")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(ggrepel)) install.packages("ggrepel")
if(!require(gplots)) install.packages("gplots")
if(!require(gprofiler2)) install.packages("gprofiler2")
if(!require(pheatmap)) install.packages("pheatmap")
if(!require(RColorBrewer)) install.packages("RcolorBrewer")
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(UpSetR)) install.packages("UpSetR")

###############################################
#             Includes                        #
###############################################
library(edgeR)
library(gplots)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
# Other libraries are loaded in the required section

###############################################
#             Begin Analysis                  #
###############################################

# set base directory  ### modify accordingly ###
basedir <- getwd()
basedir <- (paste0(basedir,"/results"))
setwd(basedir)

# load auxiliary functions
source("../scripts/common_functions.R")

# do MDS by variety to look at them individually
source(file="../scripts/get_MDS_by_group.R")

##########################################################
# 1. Import Rice_shade (raw reads)                       #
##########################################################

# Change to project results directory
setwd(basedir)

# Load target files
targets <- read.table("../data/samples_shortxLuk.txt", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
targets$treatment[which(targets$treatment == "C")] <- "WL"

# Load raw counts data
cg <- read.delim("../data/Master_counts_swap.txt", sep = "\t")

# Define groups for EdgeR
# Here I combine the variety and treatment factor so that I will later compare just between treatment within each variety.
# I convert them into a single factor
# See https://bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html for more information on this
lev <- levels(factor(str_c(targets$variety, "_", targets$treatment)))
group <- factor(str_c(targets$variety, "_", targets$treatment), levels = lev)

targets <- targets %>% 
            add_column(population = "indica")
targets$population[which(targets$variety == "LukTakhar")] <- "japonica"
targets$population[which(targets$variety == "MBlatec")] <- "japonica"
write_delim(targets, file = "../data/samples_shortxLuk_populations.txt")

# This allows us to fit a means model to the data, using a design matrix coded as
design <- model.matrix(~0+group)
colnames(design) <- levels(group)

# To check for redundancy of model parameters, one can compare between the number of columns in the design matrix with ncol(design)
# to the rank of the matrix with qr(design)$rank
ncol(design)
qr(design)$rank

# Create the contrasts matrix
# Define contrasts for comparisons of interest using the makeContrasts function. The contrasts are coded as:

################################################################################# 
# Example contrasts                                                             #
#                                                                               #
# contrasts <- makeContrasts(                                                   #
#   BVsT=(groupLUNG_B+groupBRAIN_B)/2-(groupLUNG_T+groupBRAIN_T)/2,             #
#   LungVsBrain=(groupLUNG_B+groupLUNG_T)/2-(groupBRAIN_B+groupBRAIN_T)/2,      #
#   BVsT_Lung=groupLUNG_B-groupLUNG_T,                                          #
#   BVsT_Brain=groupBRAIN_B-groupBRAIN_T,                                       #
#   LungVsBrain_B=groupLUNG_B-groupBRAIN_B,                                     #
#   LungVsBrain_T=groupLUNG_T-groupBRAIN_T,                                     #
#   levels=colnames(design))                                                    #
# rownames(contrasts) <- gsub("group", "", rownames(contrasts))                 #
# contrasts                                                                     #
################################################################################# 

contrasts <- makeContrasts(
  shade_IR64 = (IR64_FR-IR64_WL),
  shade_LukTakhar = (LukTakhar_FR-LukTakhar_WL),
  shade_MBlatec = (MBlatec_FR-MBlatec_WL),
  shade_Mudgo = (Mudgo_FR-Mudgo_WL),
  shade_Sabharaj = (Sabharaj_FR-Sabharaj_WL),
  shade_Zhenshan = (Zhenshan_FR-Zhenshan_WL),
  general_shade = (IR64_FR + LukTakhar_FR + MBlatec_FR + Mudgo_FR + Sabharaj_FR + Zhenshan_FR)/6 - 
    (IR64_WL + LukTakhar_WL + MBlatec_WL + Mudgo_WL + Sabharaj_WL + Zhenshan_WL)/6,
  levels = colnames(design))
contrasts

# Create DGEList element with norm counts
# This object is easy to use as it can be manipulated like an ordinary list in R, and it can also be subsetted like a matrix. The main components of
# a DGEList object are a matrix of read counts, sample information in the data.frame format and optional gene annotation. We enter the counts into a
# DGEList object using the function DGEList in edgeR:

# I removed sample Luk.C.R1 according to the target data file
df <- as.matrix(cg[,setdiff(colnames(cg), c("Luk.C.R1"))])
rownames(df) <- rownames(cg)

y <- DGEList(counts= df ,
             group=group)
colnames(y)=rownames(targets) # Set Column names

y$samples$lib.size = colSums(y$counts) # Recalculate library sizes


# Save the library sizes plot
png(file="raw-library_sizes.png",    # create PNG for the library sizes        
     width = 8*600,        #  8 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
par(mar=c(9,4.1,4.1,2.1))
barplot(y$samples$lib.size, names=colnames(y), las=2)
# Add a title to the plot
title("Barplot of library sizes")
dev.off()

logcounts <- cpm(y,log=TRUE)

png(file=paste0("raw-histogram_logCPM.png"),    # create PNG for the histogram      
     width = 8*600,        #  8 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
par(mar=c(9,4.1,4.1,2.1))
hist(logcounts)
dev.off()

png(filename = paste0("raw-BoxPlot.png"),    # create PNG for the BoxPlot        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)
par(mar=c(9,4.1,4.1,2.1))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title(paste0("Boxplots of logCPMs (unnormalised)"))
dev.off()

# Choose different colours for each sample type
col.cell <- c("violet", "purple", "blue", "cyan", "green", "dark gray","orange", "pink", "magenta", "red", "dark red", "brown", "gray", "black")[factor(str_c(targets$variety, "_", targets$treatment))]
data.frame(str_c(targets$variety, "_", targets$treatment),col.cell)

png(filename = paste0("raw-MDS_plot.png"),    # create PNG for the MDS        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size)
plotMDS(y,col=col.cell)
title("MDS by Sample")
dev.off()


# Filter to remove low counts
# Genes that have very low counts across all the libraries should be removed prior to downstream
# analysis. This is justified on both biological and statistical grounds. From biological point of
# view, a gene must be expressed at some minimal level before it is likely to be translated into a
# protein or to be considered biologically important. From a statistical point of view, genes with
# consistently low counts are very unlikely be assessed as significantly DE because low counts do
# not provide enough statistical evidence for a reliable judgement to be made. Such genes can
# therefore be removed from the analysis without any loss of information.

# Find out how many rows have 0 expression among all samples
table(rowSums(y$counts==0)==47)
# FALSE  TRUE 
# 41493   696

# According to the EdgeR manual, as a rule of thumb, genes are dropped if they can’t possibly be expressed in all
# the samples for any of the conditions. So, here we set the cut-off such that cpm has to be > 1 across the size of the
# group with less number of samples. In our case, the smallest group of samples contains 3 libraries.
keep <- rowSums(cpm(y)>1)>=3
table(keep)
# FALSE  TRUE 
# 16264 25925 

y.filtered <- y[keep, , keep.lib.sizes=FALSE] #25925

# Plot the library sizes
# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
par(mar=c(9,4.1,4.1,2.1))
barplot(y.filtered$samples$lib.size, names=colnames(y), las=2)
# Add a title to the plot
title("Barplot of library sizes")

# Save the library sizes plot
png(file="fil-library_sizes.png",    # create PNG for the MDS        
     width = 8*600,        #  8 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
par(mar=c(9,4.1,4.1,2.1))
barplot(y.filtered$samples$lib.size, names=colnames(y), las=2)
# Add a title to the plot
title("Barplot of library sizes")
dev.off()


# Count data is not normally distributed, so if we want to examine the distributions of the raw counts we need to log the counts. 
# We'll use box plots to check the distribution of the read counts on the log2 scale. 
# We can use the cpm function to get log2 counts per million, which are corrected for the different library sizes. 
# The cpm function also adds a small offset to avoid taking log of zero.
logcounts.filtered <- cpm(y.filtered,log=TRUE)

png(file=paste0("fil-histogram_logCPM.png"),    # create PNG for the histogram        
     width = 8*600,        #  8 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
par(mar=c(9,4.1,4.1,2.1))
hist(logcounts.filtered)
dev.off()

png(filename = paste0("fil-BoxPlot.png"),    # create PNG for the BoxPlot        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)
par(mar=c(9,4.1,4.1,2.1))
boxplot(logcounts.filtered, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts.filtered),col="blue")
title(paste0("Boxplots of logCPMs (filtered)"))
dev.off()

# Create a Multi Dimensional Scaling Plot
# The RNA samples can be clustered in two dimensions using multi-dimensional scaling (MDS) plots. This is both an analysis step and a quality
# control step to explore the overall differences between the expression profiles of the different samples.
# In the MDS plot, the distance between each pair of samples can be interpreted as the leading log-fold change between the samples for the genes
# that best distinguish that pair of samples. By default, leading fold-change is defined as the root-mean-square of the largest 500 log2-fold changes
# between that pair of samples.

# We specify the option to let us plot only one plot
par(mfrow=c(1,1))

# Check number and order of the samples
levels(factor(str_c(targets$variety, "_", targets$treatment)))

## Choose different colours for each sample type
col.cell <- c("violet", "purple", "blue", "cyan", "green", "dark gray","orange", "pink", "magenta", "red", "dark red", "brown", "gray", "black")[factor(str_c(targets$variety, "_", targets$treatment))]
data.frame(str_c(targets$variety, "_", targets$treatment),col.cell)

# MDS with sample type colouring
plotMDS(y.filtered,col=col.cell)

png(filename = paste0("fil-MDS_plot.png"),    # create PNG for the MDS        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size)
plotMDS(y.filtered,col=col.cell)
title("MDS by Sample")
dev.off()

# Hierarchical clustering with heatmaps
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts.filtered, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)
# Subset logcounts matrix
highly_variable_lcpm <- logcounts.filtered[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

## Get some nice colours for the variable genes heatmap
mypalette <- brewer.pal(7,"Spectral")
morecols <- colorRampPalette(mypalette)

# Plot the heatmap
par(mar=c(9,4.1,4.1,2.1))
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Variability of genes across samples",
          ColSideColors=col.cell,scale="row",
          cexRow=1,
          cexCol=1,
          margins=c(7,6),
          srtCol=45,
          labRow=FALSE)

# Save the heatmap
png(file="fil-top500_var_genes-heatmap.png",    # create PNG for the MDS        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
par(mar=c(9,4.1,4.1,2.1))
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Variability of genes across samples",
          ColSideColors=col.cell,scale="row",
          cexRow=1,
          cexCol=1,
          margins=c(7,6),
          srtCol=45,
          labRow=FALSE)
dev.off()

##########################################################
# 2. Normalise reads                                     #
##########################################################

# Create design model to define comparison groups
# Linear modeling and differential expression analysis in edgeR requires a design matrix to be specified. The design matrix records which
# treatment conditions were applied to each samples, and it also defines how the experimental effects are parametrized in the linear models.
#
# Similar example from Limma (see user guide), with a timecourse:
# lev <- c("wt.0hr","wt.6hr","wt.24hr","mu.0hr","mu.6hr","mu.24hr")
# f <- factor(targets$Target, levels=lev)
# design <- model.matrix(~0+f)
# colnames(design) <- lev
# fit <- lmFit(eset, design)

lev <- levels(factor(str_c(targets$variety, "_", targets$treatment)))
group <- factor(str_c(targets$variety, "_", targets$treatment), levels = lev)
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
# To check for redundancy of model parameters, one can compare between the number of columns in the design matrix with ncol(design)
# to the rank of the matrix with qr(design)$rank
ncol(design)
qr(design)$rank
save(design, file = "design.RData")

# Calculate normalization factors by library
# Normalization by trimmed mean of M values (TMM) (Robinson and Oshlack 2010) is performed by using the calcNormFactors function, which returns
# the DGEList argument with only the  norm.factors changed. It calculates a set of normalization factors, one for each sample, to eliminate 
# composition biases between libraries. The product of these factors and the library sizes defines the effective library size, which replaces 
# the original library size in all downstream analyses.
y.filtered.norm<-calcNormFactors(y.filtered, method = "TMM")

# Check before and after TMM normalization effect
# The expression profiles of individual samples can be explored more closely with mean-difference (MD) plots. An MD plot visualizes the library
# size-adjusted log-fold change between two libraries (the difference) against the average log-expression across those libraries (the mean).
par(mfrow=c(2,2)) # plot 2 by 2 graphs in the same layout
plotMD(logcounts.filtered,column = 7) # before normalization
abline(h=0,col="blue")
plotMD(logcounts.filtered,column = 8) # before normalization
abline(h=0,col="blue")
plotMD(y.filtered.norm,column = 7) # after TMM normalization
abline(h=0,col="red")
plotMD(y.filtered.norm,column = 8) # after TMM normalization
abline(h=0,col="red")

# Normalization graphs (one per sample side by side)
pdf("rice_samples_before_and_after_TMM_normalization.pdf",
    paper = "a4r")
for (i in seq(from = 0, to = length(rownames(targets)), by= 6)) {
  par(mfrow=c(2,3), mar=c(9,4.1,4.1,2.1))
  
  print(i)
  for (j in seq(c(0:5))) {
    if ((i+j) <= length(rownames(targets))) {  
      print(paste0("j =", j, " and i=", i+j))
      plotMD(y.filtered,column = i+j) # before normalization
      abline(h=0,col="blue")
      plotMD(y.filtered.norm,column = i+j) # after TMM normalization
      abline(h=0,col="red")
    }
  }
  if (i > length(rownames(targets))) { i = length(rownames(targets))}
  
}
dev.off()

# get normalised logCPMs
logcounts.filtered.norm <- cpm(y.filtered.norm, log = TRUE)
logcounts.filtered.norm.df <- as_tibble(logcounts.filtered.norm, rownames = "Gene_ID")

png(file=paste0("norm-histogram_logCPM.png"),    # create PNG for the histogram        
     width = 8*600,        #  8 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
par(mar=c(9,4.1,4.1,2.1))
hist(logcounts.filtered.norm)
dev.off()

png(filename = paste0("norm-BoxPlot.png"),    # create PNG for the BoxPlot        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)
par(mar=c(9,4.1,4.1,2.1))
boxplot(logcounts.filtered.norm, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts.filtered.norm),col="blue")
title(paste0("Boxplots of logCPMs (normalised) "))
dev.off()

png(filename = paste0("norm-MDS_plot.png"),    # create PNG for the MDS        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size)
plotMDS(y.filtered.norm,col=col.cell)
title("MDS by Sample")
dev.off()


############################################################
# 3. Hierarchical clustering (filtered, normalised reads)  #
############################################################

# Hierarchical clustering ---------------
#hierarchical clustering can only work on a data matrix, not a data frame
distance <- dist(t(logcounts.filtered.norm), method = "euclidean") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters <- hclust(distance, method = "complete") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
par(mfrow=c(1,1))
plot(clusters, labels=rownames(targets))

png(filename = paste0("norm-clustering.png"),    # create PNG for the clustering        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size)
plot(clusters, labels=rownames(y.filtered.norm$samples))
dev.off()

# Correlation pheatmap
# use pearson for CPMs hich looks for linearly dependent variables
# use spearman for log CPMs, as When you log transform them, you change the relationships between genes, ceasing to be linear.

rld_cor <- cor(logcounts.filtered.norm, method = "spearman")
rld_cor2 <- cor(cpm(y.filtered.norm, log = FALSE), method = "pearson") 

heat.colors <- brewer.pal(9, "Blues")

ph1 <- pheatmap(rld_cor, color = heat.colors, border_color=NA, fontsize = 10, 
         fontsize_row = 10, height=20)
ph2 <- pheatmap(rld_cor2, color = heat.colors, border_color=NA, fontsize = 10, 
         fontsize_row = 10, height=20)


#################################################
# 4. PCA analysis (filtered, normalised reads)  #
#################################################

# Principal component analysis (PCA) -------------
pca.res <- prcomp(t(logcounts.filtered.norm), scale.=F, retx=T)
#look at the PCA result (pca.res) that you just created
summary(pca.res) # Prints variance summary for all principal components.
pca.res$rotation #$rotation shows you how much each gene influenced each PC (called 'scores')
pca.res$x # 'x' shows you how much each sample influenced each PC (called 'loadings')
#note that these have a magnitude and a direction (this is the basis for making a PCA plot)

stats::screeplot(pca.res) # A screeplot is a standard way to view eigenvalues for each PCA
pc.var<-pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var/sum(pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC

png(filename = paste0("norm-PCA_percentage_of_variance_by_PC.png"),    # create PNG for the clustering        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size)
plot(pc.per)
title("Percentage variance explained by each PC")
dev.off()


# Visualize your PCA result ------------------
#lets first plot any two PCs against each other
#We know how much each sample contributes to each PC (loadings), so let's plot
pca.res.df <- as_tibble(pca.res$x)

sampleLabels <- rownames(y.filtered.norm$samples)

pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = group) +
  geom_point(size = 4) +
  # geom_text_repel(size = 4) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  guides(color=guide_legend(title="Variety/Treatment")) +
  theme_bw()

png(filename = paste0("Supp_fig_5A-norm-PCA-by_", "group", ".png"),    # create PNG for the PCA    
    width = 7.55*600,        # 7.55 x 600 pixels
    height = 6.8*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size)
print(pca.plot)
dev.off()


# plot PCA plots for PC1 vs PC2 by co-variate
for (i in seq_along(colnames(targets))){
   pca.plot <- ggplot(pca.res.df) +
   aes(x=PC1, y=PC2, label=sampleLabels, color = targets[,i]) +
   geom_point(size = 4) +
   # geom_text_repel(size = 4) +
   xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
   ylab(paste0("PC2 (",pc.per[2],"%",")")) +
   labs(title="PCA plot",
        caption=paste0("produced on ", Sys.time())) +
   coord_fixed() +
     guides(color=guide_legend(title=colnames(targets)[i])) +
   theme_bw()

 png(filename = paste0("Supp_fig_5A-norm-PCA-by_", colnames(targets)[i], ".png"),    # create PNG for the PCA    
     width = 7.55*600,        # 7.55 x 600 pixels
     height = 6.8*600,
     res = 600,            # 600 pixels per inch
     pointsize = 8)        # smaller font size)
  print(pca.plot)
  dev.off()
}

# Create a PCA 'small multiples' chart ----
# this is another way to view PCA laodings to understand impact of each sample on each principal component
pca.res.df <- pca.res$x[,1:20] %>% # note that this is the first time you've seen the 'pipe' operator from the magrittr package
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group,
             variety = targets$variety,
             treatment = targets$treatment,
             rep = targets$rep,
             subpop = targets$subpopulation,
             shade.rank = targets$ShadingRank,
             culm.resp = targets$culm_response,
             internode.resp = targets$internode_response,
             leaf.resp = targets$leaf_response)

pca.pivot <- pivot_longer(pca.res.df, # dataframe to be pivoted
                          cols = PC1:PC20, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)
pca.pivot$PC <- factor(pca.pivot$PC, 
                       levels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10",
                                  "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "P20"))

small.mult.pca <- ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill = group) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  xlab(paste0("Sample")) + 
  ylab(paste0("PCA loadings")) +
  theme_bw() +
  theme(axis.text.y= element_text(size = "4")) +
  coord_flip()
png(filename = "norm-PCA_small_multiples-by_group.png",    # create PNG for the PCA 
    width = 7.55*600,        # 7.55 x 600 pixels
    height = 6.8*600,
    res = 600,            # 600 pixels per inch
    pointsize = 4)        # smaller font size)
print(small.mult.pca)
dev.off()

small.mult.pca <- ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill = treatment) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  xlab(paste0("Sample")) + 
  ylab(paste0("PCA loadings")) +
  theme_bw() +
  theme(axis.text.y= element_text(size = "4")) +
  coord_flip()
png(filename = "norm-PCA_small_multiples-by_treatment.png",    # create PNG for the PCA 
    width = 7.55*600,        # 7.55 x 600 pixels
    height = 6.8*600,
    res = 600,            # 600 pixels per inch
    pointsize = 4)        # smaller font size)
print(small.mult.pca)
dev.off()

small.mult.pca <- ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill = variety) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  xlab(paste0("Sample")) + 
  ylab(paste0("PCA loadings")) +
  theme_bw() +
  theme(axis.text.y= element_text(size = "4")) +
  coord_flip()
png(filename = "norm-PCA_small_multiples-by_variety.png",    # create PNG for the PCA 
    width = 7.55*600,        # 7.55 x 600 pixels
    height = 6.8*600,
    res = 600,            # 600 pixels per inch
    pointsize = 4)        # smaller font size)
print(small.mult.pca)
dev.off()

small.mult.pca <- ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill = subpop) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  xlab(paste0("Sample")) + 
  ylab(paste0("PCA loadings")) +
  theme_bw() +
  theme(axis.text.y= element_text(size = "4")) +
  coord_flip()
png(filename = "norm-PCA_small_multiples-by_subpopulation.png",    # create PNG for the PCA 
    width = 7.55*600,        # 7.55 x 600 pixels
    height = 6.8*600,
    res = 600,            # 600 pixels per inch
    pointsize = 4)        # smaller font size)
print(small.mult.pca)
dev.off()


############################################
# 5. Fit GLM (filtered, normalised reads)  #
############################################

# Estimate dispersions:
# EdgeR estimates an empirical Bayes moderated dispersion for each individual gene. It also estimates a common dispersion, which is a
# global dispersion estimate averaged over all genes, and a trended dispersion where the dispersion of a gene is predicted from its abundance
save(y.filtered.norm, file = "y.filtered.norm.RData")
y.disp <- estimateDisp(y.filtered.norm, design = design, robust = TRUE) #25925
save(y.disp, file = "y.disp.RData")

# This returns a DGEList object with additional components (common.dispersion,  trended.dispersion and tagwise.dispersion) added to hold the 
# estimated dispersions. Here robust=TRUE has been used to protect the empirical Bayes estimates against the possibility of outlier genes with 
# exceptionally large or small individual dispersions (Phipson et al. 2016).

# Plot the dispersions:
# The vertical axis of the plotBCV plot shows square-root dispersion, also known as biological coefficient of variation (BCV) 
# (McCarthy, Chen, and Smyth 2012). For RNA-seq studies, the NB dispersions tend to be higher for genes with very low counts. 
# The dispersion trend tends to decrease smoothly with abundance and to asymptotic to a constant value for genes with larger counts.
par(mfrow=c(1,1)) # Plot only one graph
plotBCV(y.disp) # Plot the dispersion of the data

# Save the dispersion plot
png(file="rice_shade_data_dispersions.png",    # create PNG for the MDS        
    width = 9*600,        # 9 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
plotBCV(y.disp) # Plot the dispersion of the data
dev.off()

# Calculate fitted.values counts (pseudo counts) and perform ANOVA like function (Quasi Likelihood F-test)

# The NB model can be extended with quasi-likelihood (QL) methods to account for gene-specific variability from both biological and technical
# sources (Lund et al. 2012; Lun, Chen, and Smyth 2016). Under the QL framework, the NB dispersion trend is used to describe the overall biological 
# variability across all genes, and gene-specific variability above and below the overall level is picked up by the QL dispersion. In the QL approach,
# the individual (tagwise) NB dispersions are not used.
# The estimation of QL dispersions is performed using the glmQLFit function.
# Setting robust=TRUE in glmQLFit is usually recommended (Phipson et al. 2016). This allows gene-specific prior df estimates, with lower values for
# outlier genes and higher values for the main body of genes. This reduces the chance of getting false positives from genes with extremely high or low
# raw dispersions, while at the same time increasing statistical power to detect differential expression for the main body of genes.
fit <- glmQLFit(y.disp, design = design, robust = TRUE)
save(fit, file = "fit.RData")

# This returns a DGEGLM object with the estimated values of the GLM coefficients for each gene. It also contains a number of empirical Bayes (EB)
# statistics including the QL dispersion trend, the squeezed QL dispersion estimates and the prior degrees of freedom (df). The QL dispersions can be
# visualized by plotQLDisp function.
plotQLDisp(fit)

# Save the QL dispersion plot
png(file="QL_dispersions.png",    # create PNG for the MDS        
    width = 9*600,        # 9 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
plotQLDisp(fit) # Plot the QL dispersion of the data
dev.off()

# We use QL F-tests instead of the more usual likelihood ratio tests (LRT) as they give stricter error rate control by accounting for the uncertainty
# in dispersion estimation
qlf <- glmQLFTest(fit, contrast = contrasts)

fdr<-p.adjust(qlf$table$PValue, method="BH") # Calculate FDR
tt<-cbind(qlf$table, fdr) 
final<-list(qlf, tt)
names(final)=c("qlf", "full")

# Genes that pass read density filter and likelihood test
dim(final$qlf)
genes_qlf <- dim(final$qlf)[1] #25925

# Calculate normalized logCPM data
logcounts2 <- cpm(y.disp, normalized.lib.sizes = TRUE, log=TRUE)

# Calculate normalized counts (cpm)
norm_counts <- cpm(y.disp, normalized.lib.sizes = TRUE) 

# Plot some example genes
par(mfrow=c(2,2), mar=c(9,4.1,4.1,2.1))
barplot(norm_counts["LOC_Os09g20940", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("PAR1 - LOC_Os09g20940")
barplot(norm_counts["LOC_Os03g06654", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("YUC3 - LOC_Os03g06654")
barplot(norm_counts["LOC_Os01g38530", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("ELF3 - LOC_Os01g38530")
barplot(norm_counts["LOC_Os08g06110", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("CCA1 - LOC_Os08g06110")

# Plot example shade response genes
png(file="gene_controls_rice.png",    # create PNG for the MDS        
    width = 9*600,        # 9 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
par(mfrow=c(2,2), mar=c(9,4.1,4.1,2.1))
barplot(norm_counts["LOC_Os09g20940", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("PAR1 - LOC_Os09g20940")
barplot(norm_counts["LOC_Os03g06654", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("YUC3 - LOC_Os03g06654")
barplot(norm_counts["LOC_Os01g38530", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("ELF3 - LOC_Os01g38530")
barplot(norm_counts["LOC_Os08g06110", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("CCA1 - LOC_Os08g06110")
dev.off()

# Plot example shade response genes
png(file="shade_controls_rice.png",    # create PNG for the MDS        
    width = 9*600,        # 9 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
par(mfrow=c(2,2), mar=c(9,4.1,4.1,2.1))
barplot(norm_counts["LOC_Os07g47450", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("FLP1 - LOC_Os07g47450")
barplot(norm_counts["LOC_Os05g04740", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("PIF3 - LOC_Os05g04740")
barplot(norm_counts["LOC_Os07g08460", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("IAA3 - LOC_Os07g08460")
barplot(norm_counts["LOC_Os08g03310", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("Zinc finger CCCH domain-containing protein 54 - LOC_Os08g03310")
dev.off()

# Hierarchical clustering with heatmaps
# We estimate the variance for each row in the logcounts matrix (normalized values)
var_genes <- apply(logcounts2, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)
# Subset logcounts2 matrix
highly_variable_lcpm <- logcounts2[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

# Plot the heatmap
library(gplots)
par(mar=c(9,4.1,4.1,2.1))
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Variability of genes across samples (normalized data)",
          ColSideColors=col.cell,scale="row",
          cexRow=1,
          cexCol=1,
          margins=c(7,6),
          srtCol=45,
          labRow=FALSE)

# Save the heatmap
png(file="norm-top500_var_genes-heatmap.png",    # create PNG for the MDS        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
par(mar=c(9,4.1,4.1,2.1))
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Variability of genes across samples (normalized data)",
          ColSideColors=col.cell,scale="row",
          cexRow=1,
          cexCol=1,
          margins=c(7,6),
          srtCol=45,
          labRow=FALSE)
dev.off()

# Save all count data
write.table(qlf$fitted.values, file="rice.genes.fitted.values.tab", sep="\t", row.names = TRUE, col.names = NA, na = "NA")
write.table(cg, file="rice.genes.raw.values.tab", sep="\t", row.names = TRUE, col.names = NA, na = "NA")
write.table(cg[keep,], file="rice.genes.filtered.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
write.table(logcounts.filtered, file="rice.genes.logCPM.unnormalised.values.tab", sep="\t", row.names = TRUE, col.names = NA, na = "NA")
write.table(logcounts2, file="rice.genes.logCPM.normalised.values.tab", sep="\t", row.names = TRUE, col.names = NA, na = "NA")
write.table(norm_counts, file="rice.genes.CPM.normalised.values.tab", sep="\t", row.names = TRUE, col.names = NA, na = "NA")

# Add annotation data to normalised counts and logCPM data
md <- read.csv(file = "../data/Osativa_323_v7.0.annotation_info.txt", na.strings=c("", "NA"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
md <- md[,setdiff(colnames(md), c("X.pacId"))]
colnames(md)[1] <- "Gene_ID"
md <- md[!duplicated(md$Gene_ID),]
write_delim(md, file = "rice_clean_metadata_msu7.txt", na ="NA", delim = "\t")

# get MSU7 to RAP-DB mappings
msu_to_rap <- read_tsv(file = "../bioData/RAP-DB/MSU_to_RAP_2022-03-11_unique_IDs_clean.tab")

# get Arabidopsis orthologs from BioMart
arab_ortho <- read_tsv(file = "../bioData/RAP-DB/mart_export_Arabidopsis_thaliana_orthologs.txt")
# clean up NAs
arab_ortho$`Arabidopsis thaliana gene name`[which(is.na(arab_ortho$`Arabidopsis thaliana gene name`))] <- arab_ortho$`Arabidopsis thaliana gene stable ID`[which(is.na(arab_ortho$`Arabidopsis thaliana gene name`))]

# Merge orthologs and mappings
colnames(arab_ortho)[1] <- colnames(msu_to_rap)[2]
msu_to_rap_arab_ortho <- merge(msu_to_rap, arab_ortho,"RAP-DB", all = TRUE)[,c(2,1,5,3,4,6)]

# Create a full metadata file for rice
full_md <- merge(md, msu_to_rap_arab_ortho, "Gene_ID", all = TRUE)[,c(1,13:14,17,15:16,10:12,2:9)]
# clean up NAs
full_md$arabi.symbol[which(is.na(full_md$arabi.symbol))] <- full_md$Best.hit.arabi.name[which(is.na(full_md$arabi.symbol))]
full_md$`RAP-DB`[which(is.na(full_md$`RAP-DB`))] <- "None"
full_md[is.na(full_md)] <- "-"
md <- full_md
rm(full_md)

# write down new md file
md <- md[!duplicated(md$Gene_ID),]
write_delim(md, file = "rice_clean_metadata_msu7.txt", na ="NA", delim = "\t")

logcounts2.annot <- getAnnnot(logcounts2, md)
norm_counts.annot <- getAnnnot(norm_counts, md)
write_delim(logcounts2.annot, file="annot.rice.genes.logCPM.normalised.values.tab", delim = "\t", na = "NA")
write_delim(norm_counts.annot, file="annot.rice.genes.CPM.normalised.values.tab", delim = "\t", na = "NA")

##############################
# 6. Sample comparisons      #
##############################

# The QL framework provides more accurate type I error rate control than an exact test, as it accounts for the uncertainty of the dispersion estimates.
# In contrast, the exact test assumes that the estimated dispersion is the true value, which can result in some inaccuracy. (The "exact" refers to the 
# fact that the p-value is calculated exactly rather than relying on approximations; however, this only results in exact type I error control when the 
# estimated and true dispersions are the same.) For this reason, I prefer using the QL methods whenever I apply edgeR.
# 
# The QL methods (and GLM methods) are also more flexible with respect to the experimental design. For example, if you got a second batch of samples, 
# all you would need to do in a GLM framework would be to change the design matrix, while the exact test methods can't handle an extra blocking factor
# for the batch.
#
# In summary, while both of the methods will work for your data set, the QL F-test is probably the better choice. There are some situations where the 
# QL F-test doesn't work well - for example, if you don't have replicates, you'd have to supply a fixed dispersion, which defeats the whole point of 
# modelling estimation uncertainty. Another situation is where the dispersions are very large and the counts are very small, whereby some of the 
# approximations in the QL framework seem to fail. In such cases, I usually switch to the LRT rather than using the exact test, for the reasons of 
# experimental flexibility that I mentioned above.
# Source: https://support.bioconductor.org/p/84291/

setwd(basedir)
dir.create("DE")
setwd(paste0(basedir,"/DE", collapse = NULL))
md <- read.csv(file = "../rice_clean_metadata_msu7.txt", na.strings=c("", "NA"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Set cut-off values
p.value <- 0.05
fdr <- 0.05

# Test for DE genes of IR64 low R-FR vs high R-FR
et <- glmQLFTest(fit, contrast=makeContrasts(shade_IR64 = (IR64_FR-IR64_WL),levels=design))
topTags(et)
fdr.gen <- p.adjust(et$table$PValue, method="BH")
et_merge<-cbind(et$table[,setdiff(colnames(et$table),c("logCPM"))], fdr.gen)
# Get O. sativa best hit and its description for each gene
# First column of metadata file must contain Gene IDs
et_merge <- getAnnnot(et_merge, md = md)
et_merge <- et_merge[order(et_merge$logFC, et_merge$PValue),] # Order by logFC and p value
top_DE<-which(et_merge$PValue < p.value & et_merge$fdr.gen < fdr)
length(unique(et_merge$Gene_ID[top_DE])) #29
detags <- rownames(y)[top_DE]
png(file="smear_IR64_FRvsCtrl.png",    # create PNG for the MDS        
    width = 5*600,        # 5 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
plotSmear(et, de.tags=detags)
title("IR64 low R-FR vs high R-FR")
dev.off()
write.table(et_merge[top_DE,], file = "00-IR64_FRvsCtrl-DE.tab", sep = "\t", row.names = FALSE)
write.table(et_merge, file = "00-IR64_FRvsCtrl-DE_full.tab", sep = "\t", row.names = FALSE)

# Test for DE genes of LukTakhar low R-FR vs high R-FR
et <- glmQLFTest(fit, contrast=makeContrasts(shade_LukTakhar = (LukTakhar_FR-LukTakhar_WL),levels=design))
topTags(et)
fdr.gen <- p.adjust(et$table$PValue, method="BH")
et_merge<-cbind(et$table[,setdiff(colnames(et$table),c("logCPM"))], fdr.gen)
# Get O. sativa best hit and its description for each gene
# First column of metadata file must contain Gene IDs
et_merge <- getAnnnot(et_merge, md = md)
et_merge <- et_merge[order(et_merge$logFC, et_merge$PValue),] # Order by logFC and p value
top_DE<-which(et_merge$PValue < p.value & et_merge$fdr.gen < fdr)
length(unique(et_merge$Gene_ID[top_DE])) #2
detags <- rownames(y)[top_DE]
png(file="smear_LukTakhar_FRvsCtrl.png",    # create PNG for the MDS        
    width = 5*600,        # 5 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
plotSmear(et, de.tags=detags)
title("LukTakhar low R-FR vs high R-FR")
dev.off()
write.table(et_merge[top_DE,], file = "01-LukTakhar_FRvsCtrl-DE.tab", sep = "\t", row.names = FALSE)
write.table(et_merge, file = "01-LukTakhar_FRvsCtrl-DE_full.tab", sep = "\t", row.names = FALSE)

# Test for DE genes of MBlatec low R-FR vs high R-FR
et <- glmQLFTest(fit, contrast=makeContrasts(shade_MBlatec = (MBlatec_FR-MBlatec_WL),levels=design))
topTags(et)
fdr.gen <- p.adjust(et$table$PValue, method="BH")
et_merge<-cbind(et$table[,setdiff(colnames(et$table),c("logCPM"))], fdr.gen)
# Get O. sativa best hit and its description for each gene
# First column of metadata file must contain Gene IDs
et_merge <- getAnnnot(et_merge, md = md)
et_merge <- et_merge[order(et_merge$logFC, et_merge$PValue),] # Order by logFC and p value
top_DE<-which(et_merge$PValue < p.value & et_merge$fdr.gen < fdr)
length(unique(et_merge$Gene_ID[top_DE])) #9
detags <- rownames(y)[top_DE]
png(file="smear_MBlatec_FRvsCtrl.png",    # create PNG for the MDS        
    width = 5*600,        # 5 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
plotSmear(et, de.tags=detags)
title("MBlatec low R-FR vs high R-FR")
dev.off()
write.table(et_merge[top_DE,], file = "02-MBlatec_FRvsCtrl-DE.tab", sep = "\t", row.names = FALSE)
write.table(et_merge, file = "02-MBlatec_FRvsCtrl-DE_full.tab", sep = "\t", row.names = FALSE)

# Test for DE genes of Mudgo low R-FR vs high R-FR
et <- glmQLFTest(fit, contrast=makeContrasts(shade_Mudgo = (Mudgo_FR-Mudgo_WL),levels=design))
topTags(et)
fdr.gen <- p.adjust(et$table$PValue, method="BH")
et_merge<-cbind(et$table[,setdiff(colnames(et$table),c("logCPM"))], fdr.gen)
# Get O. sativa best hit and its description for each gene
# First column of metadata file must contain Gene IDs
et_merge <- getAnnnot(et_merge, md = md)
et_merge <- et_merge[order(et_merge$logFC, et_merge$PValue),] # Order by logFC and p value
top_DE<-which(et_merge$PValue < p.value & et_merge$fdr.gen < fdr)
length(unique(et_merge$Gene_ID[top_DE])) #27
detags <- rownames(y)[top_DE]
png(file="smear_Mudgo_FRvsCtrl.png",    # create PNG for the MDS        
    width = 5*600,        # 5 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
plotSmear(et, de.tags=detags)
title("Mudgo low R-FR vs high R-FR")
dev.off()
write.table(et_merge[top_DE,], file = "03-Mudgo_FRvsCtrl-DE.tab", sep = "\t", row.names = FALSE)
write.table(et_merge, file = "03-Mudgo_FRvsCtrl-DE_full.tab", sep = "\t", row.names = FALSE)

# Test for DE genes of Sabharaj low R-FR vs high R-FR
et <- glmQLFTest(fit, contrast=makeContrasts(shade_Sabharaj = (Sabharaj_FR-Sabharaj_WL),levels=design))
topTags(et)
fdr.gen <- p.adjust(et$table$PValue, method="BH")
et_merge<-cbind(et$table[,setdiff(colnames(et$table),c("logCPM"))], fdr.gen)
# Get O. sativa best hit and its description for each gene
# First column of metadata file must contain Gene IDs
et_merge <- getAnnnot(et_merge, md = md)
et_merge <- et_merge[order(et_merge$logFC, et_merge$PValue),] # Order by logFC and p value
top_DE<-which(et_merge$PValue < p.value & et_merge$fdr.gen < fdr)
length(unique(et_merge$Gene_ID[top_DE])) #35
detags <- rownames(y)[top_DE]
png(file="smear_Sabharaj_FRvsCtrl.png",    # create PNG for the MDS        
    width = 5*600,        # 5 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
plotSmear(et, de.tags=detags)
title("Sabharaj low R-FR vs high R-FR")
dev.off()
write.table(et_merge[top_DE,], file = "04-Sabharaj_FRvsCtrl-DE.tab", sep = "\t", row.names = FALSE)
write.table(et_merge, file = "04-Sabharaj_FRvsCtrl-DE_full.tab", sep = "\t", row.names = FALSE)

# Test for DE genes of Zhenshan low R-FR vs high R-FR
et <- glmQLFTest(fit, contrast=makeContrasts(shade_Zhenshan = (Zhenshan_FR-Zhenshan_WL),levels=design))
topTags(et)
fdr.gen <- p.adjust(et$table$PValue, method="BH")
et_merge<-cbind(et$table[,setdiff(colnames(et$table),c("logCPM"))], fdr.gen)
# Get O. sativa best hit and its description for each gene
# First column of metadata file must contain Gene IDs
et_merge <- getAnnnot(et_merge, md = md)
et_merge <- et_merge[order(et_merge$logFC, et_merge$PValue),] # Order by logFC and p value
top_DE<-which(et_merge$PValue < p.value & et_merge$fdr.gen < fdr)
length(unique(et_merge$Gene_ID[top_DE])) #31
detags <- rownames(y)[top_DE]
png(file="smear_Zhenshan_FRvsCtrl.png",    # create PNG for the MDS        
    width = 5*600,        # 5 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
plotSmear(et, de.tags=detags)
title("Zhenshan low R-FR vs high R-FR")
dev.off()
write.table(et_merge[top_DE,], file = "05-Zhenshan_FRvsCtrl-DE.tab", sep = "\t", row.names = FALSE)
write.table(et_merge, file = "05-Zhenshan_FRvsCtrl-DE_full.tab", sep = "\t", row.names = FALSE)

# Test for DE genes of low R-FR vs high R-FR independent of varieties
et <- glmQLFTest(fit, contrast=makeContrasts(general_shade = ((IR64_FR + LukTakhar_FR + MBlatec_FR + Mudgo_FR + Sabharaj_FR + Zhenshan_FR)/6 - 
                                                              (IR64_WL + LukTakhar_WL + MBlatec_WL + Mudgo_WL + Sabharaj_WL + Zhenshan_WL)/6),
                                                              levels=design)
                                                            )
topTags(et)
fdr.gen <- p.adjust(et$table$PValue, method="BH")
et_merge<-cbind(et$table[,setdiff(colnames(et$table),c("logCPM"))], fdr.gen)
# Get O. sativa best hit and its description for each gene
# First column of metadata file must contain Gene IDs
et_merge <- getAnnnot(et_merge, md = md)
et_merge <- et_merge[order(et_merge$logFC, et_merge$PValue),] # Order by logFC and p value
top_DE<-which(et_merge$PValue < p.value & et_merge$fdr.gen < fdr)
length(unique(et_merge$Gene_ID[top_DE])) #379
detags <- rownames(y)[top_DE]
png(file="smear_all_varieties-FRvsCtrl.png",    # create PNG for the MDS        
    width = 5*600,        # 5 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
plotSmear(et, de.tags=detags)
title("All varieties low R-FR vs high R-FR")
dev.off()
write.table(et_merge[top_DE,], file = "06-all_varieties-FRvsCtrl-DE.tab", sep = "\t", row.names = FALSE)
write.table(et_merge, file = "06-all_varieties-FRvsCtrl-DE_full.tab", sep = "\t", row.names = FALSE)

# Test for DE genes of low R-FR vs high R-FR for indica varieties
et <- glmQLFTest(fit, contrast=makeContrasts(indica_shade = ((IR64_FR + Mudgo_FR + Sabharaj_FR + Zhenshan_FR)/4 - 
                                                                (IR64_WL + Mudgo_WL + Sabharaj_WL + Zhenshan_WL)/4),
                                             levels=design)
)
topTags(et)
fdr.gen <- p.adjust(et$table$PValue, method="BH")
et_merge<-cbind(et$table[,setdiff(colnames(et$table),c("logCPM"))], fdr.gen)
# Get O. sativa best hit and its description for each gene
# First column of metadata file must contain Gene IDs
et_merge <- getAnnnot(et_merge, md = md)
et_merge <- et_merge[order(et_merge$logFC, et_merge$PValue),] # Order by logFC and p value
top_DE<-which(et_merge$PValue < p.value & et_merge$fdr.gen < fdr)
length(unique(et_merge$Gene_ID[top_DE])) #319
detags <- rownames(y)[top_DE]
png(file="smear_indica_varieties-FRvsCtrl.png",    # create PNG for the MDS        
    width = 5*600,        # 5 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
plotSmear(et, de.tags=detags)
title("Indica varieties low R-FR vs high R-FR")
dev.off()
write.table(et_merge[top_DE,], file = "07-indica_varieties-FRvsCtrl-DE.tab", sep = "\t", row.names = FALSE)
write.table(et_merge, file = "07-indica_varieties-FRvsCtrl-DE_full.tab", sep = "\t", row.names = FALSE)


# Test for DE genes of low R-FR vs high R-FR for japonica varieties
et <- glmQLFTest(fit, contrast=makeContrasts(japonica_shade = ((LukTakhar_FR + MBlatec_FR)/2 - 
                                                                (LukTakhar_WL + MBlatec_WL)/2),
                                             levels=design)
)
topTags(et)
fdr.gen <- p.adjust(et$table$PValue, method="BH")
et_merge<-cbind(et$table[,setdiff(colnames(et$table),c("logCPM"))], fdr.gen)
# Get O. sativa best hit and its description for each gene
# First column of metadata file must contain Gene IDs
et_merge <- getAnnnot(et_merge, md = md)
et_merge <- et_merge[order(et_merge$logFC, et_merge$PValue),] # Order by logFC and p value
top_DE<-which(et_merge$PValue < p.value & et_merge$fdr.gen < fdr)
length(unique(et_merge$Gene_ID[top_DE])) #21
detags <- rownames(y)[top_DE]
png(file="smear_japonica_varieties-FRvsCtrl.png",    # create PNG for the MDS        
    width = 5*600,        # 5 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
plotSmear(et, de.tags=detags)
title("Japonica varieties low R-FR vs high R-FR")
dev.off()
write.table(et_merge[top_DE,], file = "08-japonica_varieties-FRvsCtrl-DE.tab", sep = "\t", row.names = FALSE)
write.table(et_merge, file = "08-japonica_varieties-FRvsCtrl-DE_full.tab", sep = "\t", row.names = FALSE)


# Merge all into a single spreadsheet and get DEGs values:
setwd(paste0(basedir,"/DE"))
gene.data.IR64 <- read.delim(file="00-IR64_FRvsCtrl-DE_full.tab", sep="\t", stringsAsFactors = FALSE, header = TRUE, dec= ".")
gene.data.LukTakhar <- read.delim(file="01-LukTakhar_FRvsCtrl-DE_full.tab", sep="\t", stringsAsFactors = FALSE, header = TRUE, dec= ".")
gene.data.MBlatec <- read.delim(file="02-MBlatec_FRvsCtrl-DE_full.tab", sep="\t", stringsAsFactors = FALSE, header = TRUE, dec= ".")
gene.data.Mudgo <- read.delim(file="03-Mudgo_FRvsCtrl-DE_full.tab", sep="\t", stringsAsFactors = FALSE, header = TRUE, dec= ".")
gene.data.Sabharaj <- read.delim(file="04-Sabharaj_FRvsCtrl-DE_full.tab", sep="\t", stringsAsFactors = FALSE, header = TRUE, dec= ".")
gene.data.Zhenshan <- read.delim(file="05-Zhenshan_FRvsCtrl-DE_full.tab", sep="\t", stringsAsFactors = FALSE, header = TRUE, dec= ".")
gene.data.general_shade <- read.delim(file="06-all_varieties-FRvsCtrl-DE_full.tab", sep="\t", stringsAsFactors = FALSE, header = TRUE, dec= ".")
gene.data.indica_shade <- read.delim(file="07-indica_varieties-FRvsCtrl-DE_full.tab", sep="\t", stringsAsFactors = FALSE, header = TRUE, dec= ".")
gene.data.japonica_shade <- read.delim(file="08-japonica_varieties-FRvsCtrl-DE_full.tab", sep="\t", stringsAsFactors = FALSE, header = TRUE, dec= ".")

full.degs.data <- merge(gene.data.IR64[!duplicated(gene.data.IR64$Gene_ID),c(1,18,20,21)], gene.data.LukTakhar[!duplicated(gene.data.LukTakhar$Gene_ID),c(1,18,20,21)], by = "Gene_ID", all = TRUE)
colnames(full.degs.data)[2:7] <- c("IR64.logFC","IR64.PValue","IR64.fdr.gen",
                                   "LukTakhar.logFC","LukTakhar.PValue","LukTakhar.fdr.gen")

full.degs.data <- merge(full.degs.data, gene.data.MBlatec[!duplicated(gene.data.MBlatec$Gene_ID),c(1,18,20,21)], by = "Gene_ID", all = TRUE)
colnames(full.degs.data)[8:10] <- c("MBlatec.logFC","MBlatec.PValue","MBlatec.fdr.gen")

full.degs.data <- merge(full.degs.data, gene.data.Mudgo[!duplicated(gene.data.Mudgo$Gene_ID),c(1,18,20,21)], by = "Gene_ID", all = TRUE)
colnames(full.degs.data)[11:13] <- c("Mudgo.logFC","Mudgo.PValue","Mudgo.fdr.gen")

full.degs.data <- merge(full.degs.data, gene.data.Sabharaj[!duplicated(gene.data.Sabharaj$Gene_ID),c(1,18,20,21)], by = "Gene_ID", all = TRUE)
colnames(full.degs.data)[14:16] <- c("Sabharaj.logFC","Sabharaj.PValue","Sabharaj.fdr.gen")

full.degs.data <- merge(full.degs.data, gene.data.Zhenshan[!duplicated(gene.data.Zhenshan$Gene_ID),c(1,18,20,21)], by = "Gene_ID", all = TRUE)
colnames(full.degs.data)[17:19] <- c("Zhenshan.logFC","Zhenshan.PValue","Zhenshan.fdr.gen")

full.degs.data <- merge(full.degs.data, gene.data.indica_shade[!duplicated(gene.data.indica_shade$Gene_ID),c(1,18,20,21)], by = "Gene_ID", all = TRUE)
colnames(full.degs.data)[20:22] <- c("Indica.logFC","Indica.PValue","Indica.fdr.gen")

full.degs.data <- merge(full.degs.data, gene.data.japonica_shade[!duplicated(gene.data.japonica_shade$Gene_ID),c(1,18,20,21)], by = "Gene_ID", all = TRUE)
colnames(full.degs.data)[23:25] <- c("Japonica.logFC","Japonica.PValue","Japonica.fdr.gen")

full.degs.data <- merge(full.degs.data, gene.data.general_shade[!duplicated(gene.data.general_shade$Gene_ID),c(1,18,20,21)], by = "Gene_ID", all = TRUE)
colnames(full.degs.data)[26:28] <- c("Shade.logFC","Shade.PValue","Shade.fdr.gen")

full.degs.data.annot <- getAnnnot(full.degs.data, md)


## Get number and Gene_ID of all degs
# Set cut-off values
p.value <- 0.05
fdr <- 0.05
logFC <- 0

degs.IR64 <- gene.data.IR64[which(!duplicated(gene.data.IR64$Gene_ID) & abs(gene.data.IR64$logFC) > logFC & gene.data.IR64$PValue < p.value & gene.data.IR64$fdr.gen < fdr),c(1,18,20,21)] #29
degs.LukTakhar <- gene.data.LukTakhar[which(!duplicated(gene.data.LukTakhar$Gene_ID) & abs(gene.data.LukTakhar$logFC) > logFC & gene.data.LukTakhar$PValue < p.value & gene.data.LukTakhar$fdr.gen < fdr),c(1,18,20,21)] #2
degs.MBlatec <- gene.data.MBlatec[which(!duplicated(gene.data.MBlatec$Gene_ID) & abs(gene.data.MBlatec$logFC) > logFC & gene.data.MBlatec$PValue < p.value & gene.data.MBlatec$fdr.gen < fdr),c(1,18,20,21)] #9
degs.Mudgo <- gene.data.Mudgo[which(!duplicated(gene.data.Mudgo$Gene_ID) & abs(gene.data.Mudgo$logFC) > logFC & gene.data.Mudgo$PValue < p.value & gene.data.Mudgo$fdr.gen < fdr),c(1,18,20,21)] #27
degs.Sabharaj <- gene.data.Sabharaj[which(!duplicated(gene.data.Sabharaj$Gene_ID) & abs(gene.data.Sabharaj$logFC) > logFC & gene.data.Sabharaj$PValue < p.value & gene.data.Sabharaj$fdr.gen < fdr),c(1,18,20,21)] #35
degs.Zhenshan <- gene.data.Zhenshan[which(!duplicated(gene.data.Zhenshan$Gene_ID) & abs(gene.data.Zhenshan$logFC) > logFC & gene.data.Zhenshan$PValue < p.value & gene.data.Zhenshan$fdr.gen < fdr),c(1,18,20,21)] #31
degs.Indica_Shade <- gene.data.indica_shade[which(!duplicated(gene.data.indica_shade$Gene_ID) & abs(gene.data.indica_shade$logFC) > logFC & gene.data.indica_shade$PValue < p.value & gene.data.indica_shade$fdr.gen < fdr),c(1,18,20,21)] #319
degs.Japonica_Shade <- gene.data.japonica_shade[which(!duplicated(gene.data.japonica_shade$Gene_ID) & abs(gene.data.japonica_shade$logFC) > logFC & gene.data.japonica_shade$PValue < p.value & gene.data.japonica_shade$fdr.gen < fdr),c(1,18,20,21)] #21
degs.general_Shade <- gene.data.general_shade[which(!duplicated(gene.data.general_shade$Gene_ID) & abs(gene.data.general_shade$logFC) > logFC & gene.data.general_shade$PValue < p.value & gene.data.general_shade$fdr.gen < fdr),c(1,18,20,21)] #379

up.down.degs <- data.frame(sample = c("IR64",
                                      "LukTakhar",
                                      "MBlatec",
                                      "Mudgo",
                                      "Sabharaj",
                                      "Zhenshan",
                                      "Indica",
                                      "Japonica",
                                      "General"),
                           
                           up = c(length(gene.data.IR64$Gene_ID[which(!duplicated(gene.data.IR64$Gene_ID) & gene.data.IR64$logFC > logFC & gene.data.IR64$fdr.gen < fdr)]),
                                  length(gene.data.LukTakhar$Gene_ID[which(!duplicated(gene.data.LukTakhar$Gene_ID) & gene.data.LukTakhar$logFC > logFC & gene.data.LukTakhar$fdr.gen < fdr)]),
                                  length(gene.data.MBlatec$Gene_ID[which(!duplicated(gene.data.MBlatec$Gene_ID) & gene.data.MBlatec$logFC > logFC & gene.data.MBlatec$fdr.gen < fdr)]),
                                  length(gene.data.Mudgo$Gene_ID[which(!duplicated(gene.data.Mudgo$Gene_ID) & gene.data.Mudgo$logFC > logFC & gene.data.Mudgo$fdr.gen < fdr)]),
                                  length(gene.data.Sabharaj$Gene_ID[which(!duplicated(gene.data.Sabharaj$Gene_ID) & gene.data.Sabharaj$logFC > logFC & gene.data.Sabharaj$fdr.gen < fdr)]),
                                  length(gene.data.Zhenshan$Gene_ID[which(!duplicated(gene.data.Zhenshan$Gene_ID) & gene.data.Zhenshan$logFC > logFC & gene.data.Zhenshan$fdr.gen < fdr)]),
                                  length(gene.data.indica_shade$Gene_ID[which(!duplicated(gene.data.indica_shade$Gene_ID) & gene.data.indica_shade$logFC > logFC & gene.data.indica_shade$fdr.gen < fdr)]),
                                  length(gene.data.japonica_shade$Gene_ID[which(!duplicated(gene.data.japonica_shade$Gene_ID) & gene.data.japonica_shade$logFC > logFC & gene.data.japonica_shade$fdr.gen < fdr)]),
                                  length(gene.data.general_shade$Gene_ID[which(!duplicated(gene.data.general_shade$Gene_ID) & gene.data.general_shade$logFC > logFC & gene.data.general_shade$fdr.gen < fdr)])),
                           
                           down = c(length(gene.data.IR64$Gene_ID[which(!duplicated(gene.data.IR64$Gene_ID) & gene.data.IR64$logFC < logFC & gene.data.IR64$fdr.gen < fdr)]),
                                    length(gene.data.LukTakhar$Gene_ID[which(!duplicated(gene.data.LukTakhar$Gene_ID) & gene.data.LukTakhar$logFC < logFC & gene.data.LukTakhar$fdr.gen < fdr)]),
                                    length(gene.data.MBlatec$Gene_ID[which(!duplicated(gene.data.MBlatec$Gene_ID) & gene.data.MBlatec$logFC < logFC & gene.data.MBlatec$fdr.gen < fdr)]),
                                    length(gene.data.Mudgo$Gene_ID[which(!duplicated(gene.data.Mudgo$Gene_ID) & gene.data.Mudgo$logFC < logFC & gene.data.Mudgo$fdr.gen < fdr)]),
                                    length(gene.data.Sabharaj$Gene_ID[which(!duplicated(gene.data.Sabharaj$Gene_ID) & gene.data.Sabharaj$logFC < logFC & gene.data.Sabharaj$fdr.gen < fdr)]),
                                    length(gene.data.Zhenshan$Gene_ID[which(!duplicated(gene.data.Zhenshan$Gene_ID) & gene.data.Zhenshan$logFC < logFC & gene.data.Zhenshan$fdr.gen < fdr)]),
                                    length(gene.data.indica_shade$Gene_ID[which(!duplicated(gene.data.indica_shade$Gene_ID) & gene.data.indica_shade$logFC < logFC & gene.data.indica_shade$fdr.gen < fdr)]),
                                    length(gene.data.japonica_shade$Gene_ID[which(!duplicated(gene.data.japonica_shade$Gene_ID) & gene.data.japonica_shade$logFC < logFC & gene.data.japonica_shade$fdr.gen < fdr)]),
                                    length(gene.data.general_shade$Gene_ID[which(!duplicated(gene.data.general_shade$Gene_ID) & gene.data.general_shade$logFC < logFC & gene.data.general_shade$fdr.gen < fdr)]))
                                    
                             )


up.down.degs <- data.frame(up.down.degs, total = rowSums(up.down.degs[,2:3]))

# plot up/down genes
# adapted from a script by Dr. Hans van Veene

png(file="Fig_3D-up_and_down_DEGs.png",    # create PNG for the MDS        
    width = 7.55*600,        # 7.55 x 600 pixels
    height = 6.8*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))

maxy <- 1.10 * max(cbind(up.down.degs$up, up.down.degs$down)) # range of the y-axis
par(tcl=-0.2, mgp=c(2, 0.5, 0)) # parameters that determin tick size and location of axis labels
par(mar=c(0, 5, 3, 1), mfrow=c(2,1)) # set the margins around the plot (mar) and define two plots in this graph (mfrow, multiple figures by row))
bc <- barplot(up.down.degs[,2], col="dark green", ylim=c(0, maxy),
              axes=F)
axis(side=2); title(ylab="Number of upreg DEGs")
title(main = "Misregulated genes")
text(x = bc, y = up.down.degs$up + 12, label = up.down.degs$sample, pos = 3, cex = 0.8, col = "black")
text(x = bc, y = up.down.degs$up , label = up.down.degs$up, pos = 3, cex = 0.8, col = "dark green")
par(mar=c(3, 5, 0, 1)) # change the margins around the plot to suit the downregulated genes
bc <- barplot(up.down.degs[,3], col="dark red", ylim=c(maxy, 0), axes=F)
axis(side=2)
title(ylab="Number of downreg DEGs")
text(x = bc, y = up.down.degs$down , label = up.down.degs$down, pos = 1, cex = 0.8, col = "red")
dev.off()


all.degs <- unique(c(degs.IR64$Gene_ID,
                     degs.LukTakhar$Gene_ID,
                     degs.MBlatec$Gene_ID,
                     degs.Mudgo$Gene_ID,
                     degs.Sabharaj$Gene_ID,
                     degs.Zhenshan$Gene_ID,
                     degs.Indica_Shade$Gene_ID,
                     degs.Japonica_Shade$Gene_ID,
                     degs.general_Shade$Gene_ID)) # 434

unique.degs <- unique(c(degs.IR64$Gene_ID,
                        degs.LukTakhar$Gene_ID,
                        degs.MBlatec$Gene_ID,
                        degs.Mudgo$Gene_ID,
                        degs.Sabharaj$Gene_ID,
                        degs.Zhenshan$Gene_ID,
                        degs.general_Shade$Gene_ID)) # 396

# Save data
write.table(full.degs.data, file="rice.genes.full.degs.tab", sep="\t", row.names = FALSE, na = "NA")
write.table(full.degs.data.annot, file="rice.genes.full.degs.annot.tab", sep="\t", row.names = FALSE, na = "NA")
write.table(all.degs, file="all.unique.degs.list.tab", sep="\t", row.names = FALSE, col.names = "Gene_ID", na = "NA")
write.table(unique.degs, file = "unique.degs.list.tab", sep="\t", row.names = FALSE, col.names = "Gene_ID", na = "NA")
write.table(degs.general_Shade$Gene_ID, file = "general.degs.list.tab", sep="\t", row.names = FALSE, col.names = "Gene_ID", na = "NA")

# Subset CPM data for clustering
norm_counts_all_DEGs <- norm_counts[which(rownames(norm_counts) %in% all.degs),]
write.table(norm_counts_all_DEGs, file="all.unique.degs.CPM.tab", sep="\t", row.names = TRUE, col.names = NA, na = "NA")
logcounts2_all_DEGs <- logcounts2[which(rownames(logcounts2) %in% all.degs),]
write.table(logcounts2_all_DEGs, file="all.unique.degs.logCPM.tab", sep="\t", row.names = TRUE, col.names = NA, na = "NA")

norm_counts_DEGs <- norm_counts[which(rownames(norm_counts) %in% unique.degs),]
write.table(norm_counts_DEGs, file="unique.degs.CPM.tab", sep="\t", row.names = TRUE, col.names = NA, na = "NA")
logcounts2_DEGs <- logcounts2[which(rownames(logcounts2) %in% unique.degs),]
write.table(logcounts2_DEGs, file="unique.degs.logCPM.tab", sep="\t", row.names = TRUE, col.names = NA, na = "NA")

norm_counts_general_DEGs <- norm_counts[which(rownames(norm_counts) %in% degs.general_Shade$Gene_ID),]
write.table(norm_counts_general_DEGs, file="general.degs.CPM.tab", sep="\t", row.names = TRUE, col.names = NA, na = "NA")
logcounts2_general_DEGs <- logcounts2[which(rownames(logcounts2) %in% degs.general_Shade$Gene_ID),]
write.table(logcounts2_general_DEGs, file="general.degs.logCPM.tab", sep="\t", row.names = TRUE, col.names = NA, na = "NA")

# upsetR graph
if (!(require(UpSetR))) {install.packages("UpSetR")}
if (!(require(ggpubr))) {install.packages("ggpubr")}


all.degs.sets <- list(IR64 = c(degs.IR64$Gene_ID),
                        LukTakhar = c(degs.LukTakhar$Gene_ID),
                        MBlatec = c(degs.MBlatec$Gene_ID),
                        Mudgo = c(degs.Mudgo$Gene_ID),
                        Sabharaj = c(degs.Sabharaj$Gene_ID),
                        Zhenshan = c(degs.Zhenshan$Gene_ID),
                        Indica = c(degs.Indica_Shade$Gene_ID),
                        Japonica = c(degs.Japonica_Shade$Gene_ID),
                        General = c(degs.general_Shade$Gene_ID))

upset(fromList(all.degs.sets), sets = c("IR64", "LukTakhar", "MBlatec", "Mudgo", "Sabharaj", "Zhenshan", "Indica", "Japonica", "General"), sets.bar.color = "#56B4E9",
      order.by = "freq",  set_size.show = TRUE, set_size.scale_max = 600)

upset(fromList(all.degs.sets), sets = c("IR64", "LukTakhar", "MBlatec", "Mudgo", "Sabharaj", "Zhenshan"), sets.bar.color = "#56B4E9",
      order.by = "freq",  set_size.show = TRUE, set_size.scale_max = 50)


png(file="upsetR_plot.png",    # create PNG for the MDS        
    width = 5*600,        # 5 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
upset(fromList(all.degs.sets), sets = c("IR64", "LukTakhar", "MBlatec", "Mudgo", "Sabharaj", "Zhenshan", "Indica", "Japonica", "General"), sets.bar.color = "#56B4E9",
      order.by = "freq",  set_size.show = TRUE, set_size.scale_max = 600)
dev.off()

png(file="Fig_3C-upsetR_plot_per_variety-fdr_0.05.png",    # create PNG for the MDS        
    width = 7.55*600,        # 7.55 x 600 pixels
    height = 6.8*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
upset(fromList(all.degs.sets), sets = c("IR64", "LukTakhar", "MBlatec", "Mudgo", "Sabharaj", "Zhenshan"), sets.bar.color = "#56B4E9",
      order.by = "freq",  set_size.show = TRUE, set_size.scale_max = 50)
dev.off()

png(file="upsetR_plot_per_variety_and_general_shade.png",    # create PNG for the MDS        
    width = 5*600,        # 5 x 600 pixels
    height = 5*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size))
upset(fromList(all.degs.sets), sets = c("IR64", "LukTakhar", "MBlatec", "Mudgo", "Sabharaj", "Zhenshan", "General"), sets.bar.color = "#56B4E9",
      order.by = "freq",  set_size.show = TRUE, set_size.scale_max = 450)
dev.off()

##############################################
# 7. PCA of unique degs (DEGs with q < 0.05) #
##############################################

# PCA of unique degs

# Choose the appropriate subset of unique genes for graphs

## all possible unique DEGs
# logcounts.all.degs <- logcounts2[which(rownames(logcounts2) %in% all.degs),] # all unique degs
# label.degs <- "unique_degs_all"

## variety and general shade unique DEGs
# logcounts.all.degs <- logcounts2[which(rownames(logcounts2) %in% unique.degs),] # degs in all varieties + general shade
# label.degs <- "unique_degs_varieties_and_general"

# general shade unique DEGs
logcounts.all.degs <- logcounts2[which(rownames(logcounts2) %in% degs.general_Shade$Gene_ID),] # degs in general shade
label.degs <- "unique_degs_general"

pca.res <- prcomp(t(logcounts.all.degs), scale.=F, retx=T)

#look at the PCA result (pca.res) that you just created
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
pca.res$rotation #$rotation shows you how much each gene influenced each PC (called 'scores')
pca.res$x # 'x' shows you how much each sample influenced each PC (called 'loadings')
#note that these have a magnitude and a direction (this is the basis for making a PCA plot)
stats::screeplot(pca.res) # A screeplot is a standard way to view eigenvalues for each PCA
pc.var<-pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var/sum(pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC

png(filename = paste0(label.degs,"-PCA_percentage_of_variance_by_PC.png"),    # create PNG for the clustering        
    width = 7.55*600,        # 7.55 x 600 pixels
    height = 6.8*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size)
plot(pc.per)
title("Percentage variance explained by each PC")
dev.off()

# Visualize your PCA result ------------------
#lets first plot any two PCs against each other
#We know how much each sample contributes to each PC (loadings), so let's plot
pca.res.df <- as_tibble(pca.res$x)

sampleLabels <- rownames(y.filtered.norm$samples)

pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = group) +
  geom_point(size = 4) +
  geom_text_repel(size = 4) +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  guides(color=guide_legend(title="group")) +
  theme_bw()

png(filename = paste0(label.degs,"-PCA-by_", "group", ".png"),    # create PNG for the PCA    
    width = 20*600,        # 20 x 600 pixels
    height = 17*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size)
print(pca.plot)
dev.off()

pca.plot2 <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label = c(rep("",47)), color = group) +
  geom_point(size = 6) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  guides(color=guide_legend(title="group")) +
  theme_bw() +
  theme(axis.text.x= element_text(size = "12")) +
  theme(axis.text.y= element_text(size = "12")) 

png(filename = paste0(label.degs,"-PCA-by_", "group", ".png"),    # create PNG for the PCA    
    width = 7.55*600,        # 7.55 x 600 pixels
    height = 6.8*600,
    res = 600,            # 600 pixels per inch
    pointsize = 8)        # smaller font size)
print(pca.plot2)
dev.off()


# plot PCA plots for PC1 vs PC2 by co-variate
for (i in seq_along(colnames(targets))){
  
  pca.plot <- ggplot(pca.res.df) +
    aes(x=PC1, y=PC2, label=c(rep("",47)), color = targets[,i]) +
    geom_point(size = 6) +
    xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
    ylab(paste0("PC2 (",pc.per[2],"%",")")) +
    labs(title="PCA plot",
         caption=paste0("produced on ", Sys.time())) +
    coord_fixed() +
    guides(color=guide_legend(title=colnames(targets)[i])) +
    theme_bw() +
    theme(axis.text.x= element_text(size = "12")) +
    theme(axis.text.y= element_text(size = "12")) 
  
  png(filename = paste0(label.degs,"-PCA-by_", colnames(targets)[i], ".png"),    # create PNG for the PCA    
      width = 7.55*600,        # 7.55 x 600 pixels
      height = 6.8*600,
      res = 600,            # 600 pixels per inch
      pointsize = 8)        # smaller font size)
  print(pca.plot)
  dev.off()
}

# Create a PCA 'small multiples' chart ----
# this is another way to view PCA loadings to understand impact of each sample on each principal component
pca.res.df <- pca.res$x[,1:20] %>% # note that this is the first time you've seen the 'pipe' operator from the magrittr package
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group,
             variety = targets$variety,
             treatment = targets$treatment,
             rep = targets$rep,
             subpop = targets$subpopulation,
             shade.rank = targets$ShadingRank,
             culm.resp = targets$culm_response,
             internode.resp = targets$internode_response,
             leaf.resp = targets$leaf_response)

pca.pivot <- pivot_longer(pca.res.df, # dataframe to be pivoted
                          cols = PC1:PC20, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)

pca.pivot$PC <- factor(pca.pivot$PC, 
                       levels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10",
                                  "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "P20"))

small.mult.pca <- ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill = group) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  xlab(paste0("Sample")) + 
  ylab(paste0("PCA loadings")) +
  theme_bw() +
  theme(axis.text.y= element_text(size = "4")) +
  coord_flip()

png(filename = paste0(label.degs,"-PCA_small_multiples-by_group.png"),    # create PNG for the PCA 
    width = 7.55*600,        # 7.55 x 600 pixels
    height = 6.8*600,     
    res = 600,            # 600 pixels per inch
    pointsize = 4)        # smaller font size)
print(small.mult.pca)
dev.off()

small.mult.pca <- ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill = treatment) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  xlab(paste0("Sample")) + 
  ylab(paste0("PCA loadings")) +
  theme_bw() +
  coord_flip()
png(filename = paste0(label.degs,"-PCA_small_multiples-by_treatment.png"),    # create PNG for the PCA 
    width = 7.55*600,        # 7.55 x 600 pixels
    height = 6.8*600, 
    res = 600,            # 600 pixels per inch
    pointsize = 4)        # smaller font size)
print(small.mult.pca)
dev.off()

data <- pca.pivot %>% filter(PC == "PC2") %>% mutate(PC = paste0("PC2 ", pc.per[2],"%"))
group.colours <- c("gray", "dark red")[factor(data$treatment, levels = c("WL", "FR"))]
data.frame(data$treatment,group.colours)

small.mult.pca <- ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill = variety) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()
png(filename = paste0(label.degs,"-PCA_small_multiples-by_variety.png"),    # create PNG for the PCA 
    width = 7.55*600,        # 7.55 x 600 pixels
    height = 6.8*600,
    res = 600,            # 600 pixels per inch
    pointsize = 4)        # smaller font size)
print(small.mult.pca)
dev.off()

small.mult.pca <- ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill = subpop) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()
png(filename = paste0(label.degs,"-PCA_small_multiples-by_subpopulation.png"),    # create PNG for the PCA 
    width = 7.55*600,        # 7.55 x 600 pixels
    height = 6.8*600,
    res = 600,            # 600 pixels per inch
    pointsize = 4)        # smaller font size)
print(small.mult.pca)
dev.off()

# Generate PC1 plot coloured by variety for Fig. 3C
small.mult.pca <- ggplot(pca.pivot %>% filter(PC == "PC1") %>% mutate(PC = paste0("PC1 ", pc.per[1],"%"))) +
  aes(x=sample, y=loadings, fill = variety) + 
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(caption=paste0("produced on ", Sys.time())) +
  xlab(paste0("Sample")) + 
  ylab(paste0("PCA loadings")) +
  theme_bw(base_size = 12) +
  theme(axis.text.x= element_text(size = "10")) +
  theme(axis.text.y= element_text(size = "10")) +
  coord_flip()
png(filename = paste0("Fig_3A-", label.degs,"-PCA_small_multiples-by_variety-fdr_0.05.png"),    # create PNG for the PCA 
    width = 7.55*600,        # 7.55 x 600 pixels
    height = 6.8*600, 
    res = 600,            # 600 pixels per inch
    pointsize = 4)        # smaller font size)
print(small.mult.pca)
dev.off()

# Generate PC2 plot coloured by treatment for Fig. 3C
small.mult.pca <- ggplot(pca.pivot %>% filter(PC == "PC2") %>% mutate(PC = paste0("PC2 ", pc.per[2],"%"))) +
  aes(x=sample, y=loadings, fill = treatment) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("dark red", "gray")) +
  facet_wrap(~PC) +
  labs(caption=paste0("produced on ", Sys.time())) +
  xlab(paste0("Sample")) + 
  ylab(paste0("PCA loadings")) +
  theme_bw(base_size = 12) +
  theme(axis.text.x= element_text(size = "10")) +
  theme(axis.text.y= element_text(size = "10")) +
  coord_flip()
png(filename = paste0("Fig_3A-", label.degs,"-PCA_small_multiples-by_treatment-fdr_0.05.png"),    # create PNG for the PCA 
    width = 7.55*600,        # 7.55 x 600 pixels
    height = 6.8*600, 
    res = 600,            # 600 pixels per inch
    pointsize = 4)        # smaller font size)
print(small.mult.pca)
dev.off()

################################################
# 8. Plot heatmaps of selected lists of genes  #
################################################
# This list of genes was compiled based on all DEGs found in the dataset

if (!(require(tidyverse))) {install.packages("tidyverse")}
if (!(require(gplots))) {install.packages("gplots")}
if (!(require(readr))) {install.packages("readr")}

#Set wd
setwd(basedir)
setwd("../")

# Load DE data, normalised CPMS and logCPMs
full.degs.data <- read_tsv(file="results/DE/rice.genes.full.degs.tab")
norm_counts <- read_tsv(file="results/rice.genes.CPM.normalised.values.tab")
logcounts2 <- read_tsv(file="results/rice.genes.logCPM.normalised.values.tab")
colnames(full.degs.data)[1] <- "Gene_ID"
colnames(norm_counts)[1] <- "Gene_ID"
colnames(logcounts2)[1] <- "Gene_ID"

# Set cut-off values (modify as desired)
p.value <- 0.05
fdr <- 0.05
logFC <- 0

# Get DEGs IDs
degs.IR64 <- full.degs.data$Gene_ID[which(abs(full.degs.data$IR64.logFC) >= logFC &full.degs.data$IR64.PValue < p.value & full.degs.data$IR64.fdr.gen < fdr)] #29
degs.LukTakhar <- full.degs.data$Gene_ID[which(abs(full.degs.data$LukTakhar.logFC) >= logFC & full.degs.data$LukTakhar.PValue < p.value & full.degs.data$LukTakhar.fdr.gen < fdr)] #2
degs.MBlatec <- full.degs.data$Gene_ID[which(abs(full.degs.data$MBlatec.logFC) >= logFC & full.degs.data$MBlatec.PValue < p.value & full.degs.data$MBlatec.fdr.gen < fdr)] #9
degs.Mudgo <- full.degs.data$Gene_ID[which(abs(full.degs.data$Mudgo.logFC) >= logFC &full.degs.data$Mudgo.PValue < p.value & full.degs.data$Mudgo.fdr.gen < fdr)] #27
degs.Sabharaj <- full.degs.data$Gene_ID[which(abs(full.degs.data$Sabharaj.logFC) >= logFC & full.degs.data$Sabharaj.PValue < p.value & full.degs.data$Sabharaj.fdr.gen < fdr)] #35
degs.Zhenshan <- full.degs.data$Gene_ID[which(abs(full.degs.data$Zhenshan.logFC) >= logFC & full.degs.data$Zhenshan.PValue < p.value & full.degs.data$Zhenshan.fdr.gen < fdr)] #31
degs.Indica <- full.degs.data$Gene_ID[which(abs(full.degs.data$Indica.logFC) >= logFC & full.degs.data$Indica.PValue < p.value & full.degs.data$Indica.fdr.gen < fdr)] #319
degs.Japonica <- full.degs.data$Gene_ID[which(abs(full.degs.data$Japonica.logFC) >= logFC & full.degs.data$Japonica.PValue < p.value & full.degs.data$Japonica.fdr.gen < fdr)] #21
degs.general <- full.degs.data$Gene_ID[which(abs(full.degs.data$Shade.logFC) >= logFC & full.degs.data$Shade.PValue < p.value & full.degs.data$Shade.fdr.gen < fdr)] #379

# Create lists of genes to plot (modify as desired)
filtered.gene.lists <- list("all unique DEGs (varieties, indica, japonica and general shade)" = as.list(unique(c(degs.IR64, degs.LukTakhar, degs.MBlatec, degs.Mudgo, degs.Sabharaj, degs.Zhenshan, degs.Indica, degs.Japonica, degs.general))),
                            "unique DEGs (varieties and general shade)" = as.list(unique(c(degs.IR64, degs.LukTakhar, degs.MBlatec, degs.Mudgo, degs.Sabharaj, degs.Zhenshan, degs.general))),
                            "unique DEGs (varieties)" = as.list(unique(c(degs.IR64, degs.LukTakhar, degs.MBlatec, degs.Mudgo, degs.Sabharaj, degs.Zhenshan))),
                            "general response DEGs" = as.list(degs.general),
                            "indica response DEGs" = as.list(degs.Indica),
                            "japonica response DEGs" = as.list(degs.Japonica),
                            "all unique DEGs (indica, japonica and general shade)"  = as.list(unique(c(degs.Indica, degs.Japonica, degs.general))),
                            "all expressed genes" = as.list(c(norm_counts$Gene_ID))
)
                         

# Get dataframes of the subset of filtered genes of logCPM (df2) and mean logCPM (df4) data

# Scale logCPM normalised dataset to use for the heatmaps
df2 <- logcounts2
# get an index of rows from the full data that match the subset of genes
indexFD <- which(df2$Gene_ID %in% unique(c(unlist(filtered.gene.lists))))  # 952
# scale by row (gene)
df2.filtered <- t(scale(t(data.matrix(df2[indexFD,2:length(colnames(df2))]))))
rownames(df2.filtered) <- df2$Gene_ID[indexFD]
df2.filtered.unscaled <- data.matrix(df2[indexFD,2:length(colnames(df2))])
rownames(df2.filtered.unscaled) <- df2$Gene_ID[indexFD]

# labels for the columns for the heatmap
labels.all.replicates <- c(
  rep("IR-WL",4),   
  rep("IR-WL+FR",4),
  rep("Luk-WL",3),
  rep("Luk-WL+FR",4),
  rep("MB-WL",4),
  rep("MB-WL+FR",4),
  rep("Mudgo-WL",4),   
  rep("Mudgo-WL+FR",4),
  rep("Sab-WL",4),
  rep("Sab-WL+FR",4),
  rep("Zhen-WL",4),
  rep("Zhen-WL+FR",4)
)

# Scale calculated logCPM from rowmeans of CPM data of filtered genes
df3 <- norm_counts
# get an index of rows from the full data that match the subset of genes
indexFD3 <- which(df3$Gene_ID %in% unique(c(unlist(filtered.gene.lists))))  
df3.filtered <- df3[indexFD,2:length(colnames(df3))]
rownames(df3.filtered) <-df3$Gene_ID[indexFD]

# get the average CPM values by genotype and treatment (group)
df3.filtered.means <- data.frame(
  "IR-WL" = rowMeans(cbind(df3.filtered$`IR-C-R1`, df3.filtered$`IR-C-R2`, df3.filtered$`IR-C-R3`, df3.filtered$`IR-C-R4`)),
  "IR-FR" = rowMeans(cbind(df3.filtered$`IR-FR-R1`, df3.filtered$`IR-FR-R2`, df3.filtered$`IR-FR-R3`, df3.filtered$`IR-FR-R4`)),
  
  "Luk-WL" = rowMeans(cbind(df3.filtered$`Luk-C-R2`, df3.filtered$`Luk-C-R3`, df3.filtered$`Luk-C-R4`)),
  "Luk-FR" = rowMeans(cbind(df3.filtered$`Luk-FR-R1`, df3.filtered$`Luk-FR-R2`, df3.filtered$`Luk-FR-R3`, df3.filtered$`Luk-FR-R4`)),
  
  "MB-WL" = rowMeans(cbind(df3.filtered$`MB-C-R1`, df3.filtered$`MB-C-R2`, df3.filtered$`MB-C-R3`, df3.filtered$`MB-C-R4`)),
  "MB-FR"= rowMeans(cbind(df3.filtered$`MB-FR-R1`, df3.filtered$`MB-FR-R2`, df3.filtered$`MB-FR-R3`, df3.filtered$`MB-FR-R4`)),
  
  "Mudgo-WL" = rowMeans(cbind(df3.filtered$`Mudgo-C-R1`, df3.filtered$`Mudgo-C-R2`, df3.filtered$`Mudgo-C-R3`, df3.filtered$`Mudgo-C-R4`)),
  "Mudgo-FR" = rowMeans(cbind(df3.filtered$`Mudgo-FR-R1`, df3.filtered$`Mudgo-FR-R2`, df3.filtered$`Mudgo-FR-R3`, df3.filtered$`Mudgo-FR-R4`)),
  
  "Sab-WL" = rowMeans(cbind(df3.filtered$`Sab-C-R1`, df3.filtered$`Sab-C-R2`, df3.filtered$`Sab-C-R3`, df3.filtered$`Sab-C-R4`)),
  "Sab-FR" = rowMeans(cbind(df3.filtered$`Sab-FR-R1`, df3.filtered$`Sab-FR-R2`, df3.filtered$`Sab-FR-R3`, df3.filtered$`Sab-FR-R4`)),
  
  "Zhen-WL" = rowMeans(cbind(df3.filtered$`Zhen-C-R1`, df3.filtered$`Zhen-C-R2`, df3.filtered$`Zhen-C-R3`, df3.filtered$`Zhen-C-R4`)),
  "Zhen-FR" = rowMeans(cbind(df3.filtered$`Zhen-FR-R1`, df3.filtered$`Zhen-FR-R2`, df3.filtered$`Zhen-FR-R3`, df3.filtered$`Zhen-FR-R4`))
)
rownames(df3.filtered.means) <- rownames(df3.filtered)

# transform mean CPMs to log2 values
df4 <- log2(df3.filtered.means+1/(mean(y.filtered.norm$samples$lib.size) * 1e-6))
rownames(df4) <- rownames(df3.filtered)

# scale by row (gene)
df4.filtered <- t(scale(t(data.matrix(df4))))
df4.filtered.unscaled <- data.matrix(df4)

# labels for the columns for the mean values heatmap
labels.mean.values <- c("IR-WL", "IR-WL+FR", "Luk-WL", "Luk-WL+FR", "MB-WL", "MB-WL+FR",
  "Mudgo-WL", "Mudgo-WL+FR", "Sab-WL", "Sab-WL+FR", "Zhen-WL","Zhen-WL+FR")


# Get logFC values
df5 <- full.degs.data[,c(2,5,8,11,14,17)]
rownames(df5) <- full.degs.data$Gene_ID
# get an index of rows from the full data that match the subset of genes
indexFD5 <- which(rownames(df5) %in% unique(c(unlist(filtered.gene.lists))))  # 952
df5.filtered <- df5[indexFD,]
rownames(df5.filtered) <- rownames(df5)[indexFD]
# # scale by row (gene)
df5.filtered.unscaled <- data.matrix(df5.filtered)
df5.filtered <- t(scale(t(data.matrix(df5.filtered))))


# labels for the columns for the mean values heatmap
labels.logFC.values <- c("IR", "Luk", "MB", "Mudgo", "Sab", "Zhen")

# Graph each of the genes within the vector of gene lists
dir.create("results/heatmaps") 
dir.create("results/heatmaps/selected_genes/") 

# This piece of code also creates the necessary heatmaps for Fig. 3D and Supp. Fig. 4
for (i in seq_along(filtered.gene.lists)) {
  genes_to_plot <- unlist(filtered.gene.lists[i])
  print(genes_to_plot)
  
  selected_dir <- paste0("results/heatmaps/selected_genes/")
  
  gene_subset <- unique(unlist(genes_to_plot)) # obtain all genes belonging to that category
  gene_subset <- gene_subset[which(gene_subset %in% rownames(df2.filtered))] # obtain the subset that has been differentially regulated

  if (length(gene_subset) >= 2){
    # Select subset of genes to be graphed
    plotdf <- df2.filtered[which(rownames(df2.filtered) %in% gene_subset),1:47]
    
    if (names(filtered.gene.lists)[i] != "all expressed genes") {
      plotdf.qvalues <-full.degs.data[which(full.degs.data$Gene_ID %in% gene_subset),c(1,4,7,10,13,16,19,22,25,28)]
      significant_degs <- as.data.frame(plotdf.qvalues) %>% filter_all(any_vars(. < fdr))
      plotdf <- plotdf[which(rownames(plotdf) %in% significant_degs$Gene_ID),]
    }
    
    # Remove rows with zero values
    # Go through each row and determine if a value is zero
    row_sub = apply(plotdf, 1, function(row) all(row = 0 ))
    # Subset rows with non-zero values
    plotdf <- plotdf[!row_sub,]
  
    # Hierarchical clustering of rearranged logCPM data
    d <- dist(plotdf, method = "euclidean") # calculate euclidean distance
    clust <- hclust(d, method = "complete") # hierarchical cluster of the data
    dend <- as.dendrogram(clust) # as dendogram (for clustering rows, Rowv parameter of heatmap.2)
  
    # creates a own color palette from red to green
    my_palette <- colorRampPalette(c("blue", "#ffc000"))(n = 199)
  
    # (optional) defines the color breaks manually for a "skewed" color transition
    col_breaks = c(seq(min(plotdf),median(plotdf),length=100),  # for orange
                 seq(median(plotdf) + 0.1, max(plotdf),length=100))  # for blue
  
    
    # Plot heatmap (without group clustering)
    png(paste0(selected_dir,
               "/",
               i,
               "a-",
               gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
               "_fdr_",
               fdr,
               ".png"),    # create png for the heat map        
        width = 6.5*600,        
        if (length(rownames(plotdf)) > 5) {height = 5*600
        } else {height = 3.5*600},
        res = 600,
        pointsize = 5)        # smaller font size
    heatmap.2(plotdf,
              main = paste0(gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", names(filtered.gene.lists)[i])), # heat map title
              notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(5,1),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier
              breaks=col_breaks,    # enable color transition at specified limits
              dendrogram="none",    # only draw a row dendrogram
              Colv="NA",            # turn off column clustering 
              Rowv = TRUE,
              #scale = "row",
              useRaster=TRUE,
              labRow = NA,
              labCol = labels.all.replicates,
              cexRow = 0.2 + 1/log10(length(rownames(plotdf))),
              cexCol = 0.2 + 1/log10(length(labels.all.replicates)),
              ColSideColors = c(rep("gray",4),
                                rep("#A00000",4),
                                rep("gray",3),
                                rep("#A00000",4),
                                rep(c(rep("gray",4), rep("#A00000",4)),4)
              ),
              key.xlab = "Scaled logCPM",
              keysize = 1,
              lhei = c(1,7))
    dev.off()
    write_delim(cbind("Gene_ID" = rownames(plotdf), as.data.frame(plotdf)),
                file = paste0(selected_dir,
                              "/",
                              i,
                              "a-",
                              gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
                              "_fdr_",
                              fdr,               
                              ".tab"),
                delim = "\t", na = "NA")
    
    # Plot heatmap (with group clustering)
    png(paste0(selected_dir,
               "/",
               i,
               "b-",
               gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
               "_with_group_clustering",
               "_fdr_",
               fdr,               
               ".png"),    # create png for the heat map        
        width = 6.5*600,        
        if (length(rownames(plotdf)) > 5) {height = 5*600
        } else {height = 3.5*600},
        res = 600,
        pointsize = 5)        # smaller font size
    heatmap.2(plotdf,
              main = paste0(gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", names(filtered.gene.lists)[i])), # heat map title
              notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(5,1),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier
              breaks=col_breaks,    # enable color transition at specified limits
              dendrogram="column",    # only draw a row dendrogram
              Colv=TRUE,            # turn off column clustering 
              Rowv = TRUE,
              #scale = "row",
              useRaster=TRUE,
              labRow = NA,
              labCol = labels.all.replicates,
              cexRow = 0.2 + 1/log10(length(rownames(plotdf))),
              cexCol = 0.2 + 1/log10(length(labels.all.replicates)),
              ColSideColors = c(rep("gray",4),
                                rep("#A00000",4),
                                rep("gray",3),
                                rep("#A00000",4),
                                rep(c(rep("gray",4), rep("#A00000",4)),4)
              ),
              key.xlab = "Scaled logCPM",
              keysize = 1,
              lhei = c(1,7))
    dev.off()
    write_delim(cbind("Gene_ID" = rownames(plotdf), as.data.frame(plotdf)),
                file = paste0(selected_dir,
                              "/",
                              i,
                              "b-",
                              gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
                              "_with_group_clustering",
                              "_fdr_",
                              fdr,               
                              ".tab"),
                delim = "\t", na = "NA")
  }
  
  # Redo to graph only mean values by group
  gene_subset <- unique(unlist(genes_to_plot)) # obtain all genes belonging to that category
  gene_subset <- gene_subset[which(gene_subset %in% rownames(df4.filtered))] # obtain the subset that has been differentially regulated
  
  if (length(gene_subset) >= 2){
  # Select subset of genes to be graphed
    plotdf <- df4.filtered[which(rownames(df4.filtered) %in% gene_subset),1:12]
    
    if (names(filtered.gene.lists)[i] != "all expressed genes") {
      plotdf.qvalues <-full.degs.data[which(full.degs.data$Gene_ID %in% gene_subset),c(1,4,7,10,13,16,19,22,25,28)]
      significant_degs <- as.data.frame(plotdf.qvalues) %>% filter_all(any_vars(. < fdr))
      plotdf <- plotdf[which(rownames(plotdf) %in% significant_degs$Gene_ID),]
    }
    # Remove rows with zero values
    # Go through each row and determine if a value is zero
    row_sub = apply(plotdf, 1, function(row) all(row = 0 ))
    # Subset rows with non-zero values
    plotdf <- plotdf[!row_sub,]
    
    # Hierarchical clustering of rearranged logCPM data
    d <- dist(plotdf, method = "euclidean") # calculate euclidean distance
    clust <- hclust(d, method = "complete") # hierarchical cluster of the data
    dend <- as.dendrogram(clust) # as dendogram (for clustering rows, Rowv parameter of heatmap.2)
    
    # creates a own color palette from red to green
    my_palette <- colorRampPalette(c("blue", "#ffc000"))(n = 199)
  
    # (optional) defines the color breaks manually for a "skewed" color transition
    col_breaks = c(seq(min(plotdf),median(plotdf),length=100),  # for orange
                   seq(median(plotdf) + 0.1, max(plotdf),length=100))  # for blue
  
    # Plot means heatmap (without group clustering)
    png(paste0(selected_dir,
               "/",
               i,
               "c-",
               gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
               "_group_means",
               "_fdr_",
               fdr,               
               ".png"),    # create png for the heat map        
        width = 6.5*600,        # 6.5 x 600 pixels
        if (length(rownames(plotdf)) > 5) {height = 5*600
        } else {height = 3.5*600},
        res = 600,            # 600 pixels per inch
        #paper = "a4r",
        pointsize = 7)        # smaller font size
    heatmap.2(plotdf,
              main = paste0(gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", names(filtered.gene.lists)[i])), # heat map title
              notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(8,1),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier
              breaks=col_breaks,    # enable color transition at specified limits
              dendrogram="none",    # only draw a row dendrogram
              Colv = "NA",            # turn off column clustering 
              Rowv = TRUE,
              #scale = "row",
              useRaster=TRUE,
              labRow = NA,
              labCol = labels.mean.values,
              cexRow = 0.2 + 1/log10(length(rownames(plotdf))),
              cexCol = 0.2 + 1/log10(length(colnames(plotdf))),
              ColSideColors = c(rep(c(rep("gray",1), rep("#A00000",1)),6)),
              key.xlab = "Scaled logCPM",
              keysize = 1,
              lhei = c(2,7))
    dev.off()
    write_delim(cbind("Gene_ID" = rownames(plotdf), as.data.frame(plotdf)),
                file = paste0(selected_dir,
               "/",
               i,
               "c-",
               gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
               "_group_means",
               "_fdr_",
               fdr,               
               ".tab"),
               delim = "\t", na = "NA")

    # Plot means heatmap (with group clustering)
    png(paste0(selected_dir,
               "/",
               i,
               "d-",
               gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
               "_group_means_with_group_clustering",
               "_fdr_",
               fdr,               
               ".png"),    # create png for the heat map        
        width = 6.5*600,        # 6.5 x 600 pixels
        if (length(rownames(plotdf)) > 5) {height = 5*600
        } else {height = 3.5*600},
        res = 600,            # 600 pixels per inch
        #paper = "a4r",
        pointsize = 7)        # smaller font size
    heatmap.2(plotdf,
              main = paste0(gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", names(filtered.gene.lists)[i])), # heat map title
              notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(8,1),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier
              breaks=col_breaks,    # enable color transition at specified limits
              dendrogram="column",    # only draw a row dendrogram
              Colv = TRUE,            # turn off column clustering 
              Rowv = TRUE,
              #scale = "row",
              useRaster=TRUE,
              labRow = NA,
              labCol = labels.mean.values,
              cexRow = 0.2 + 1/log10(length(rownames(plotdf))),
              cexCol = 0.2 + 1/log10(length(colnames(plotdf))),
              ColSideColors = c(rep(c(rep("gray",1), rep("#A00000",1)),6)),
              key.xlab = "Scaled logCPM",
              keysize = 1,
              lhei = c(2,7))
    dev.off()
    write_delim(cbind("Gene_ID" = rownames(plotdf), as.data.frame(plotdf)),
                file = paste0(selected_dir,
                              "/",
                              i,
                              "d-",
                              gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
                              "_group_means_with_group_clustering",
                              "_fdr_",
                              fdr,               
                              ".tab"),
                delim = "\t", na = "NA")
  }
  
  # Redo to graph only logFC values by group
  gene_subset <- unique(unlist(genes_to_plot)) # obtain all genes belonging to that category
  gene_subset <- gene_subset[which(gene_subset %in% rownames(df5.filtered))] # obtain the subset that has been differentially regulated
  
  if (length(gene_subset) >= 2){
    # Select subset of genes to be graphed
    plotdf <- df5.filtered[which(rownames(df5.filtered) %in% gene_subset),]
    
    if (names(filtered.gene.lists)[i] != "all expressed genes") {
      plotdf.qvalues <-full.degs.data[which(full.degs.data$Gene_ID %in% gene_subset),c(1,4,7,10,13,16,19,22,25,28)]
      significant_degs <- as.data.frame(plotdf.qvalues) %>% filter_all(any_vars(. < fdr))
      plotdf <- plotdf[which(rownames(plotdf) %in% significant_degs$Gene_ID),]
    }
    # Remove rows with zero values
    # Go through each row and determine if a value is zero
    row_sub = apply(plotdf, 1, function(row) all(row = 0 ))
    # Subset rows with non-zero values
    plotdf <- plotdf[!row_sub,]
    
    # Hierarchical clustering of rearranged logCPM data
    d <- dist(plotdf, method = "euclidean") # calculate euclidean distance
    clust <- hclust(d, method = "complete") # hierarchical cluster of the data
    dend <- as.dendrogram(clust) # as dendogram (for clustering rows, Rowv parameter of heatmap.2)
    
    # creates a own color palette from red to green
    my_palette <- colorRampPalette(c("blue", "#ffc000"))(n = 199)
    
    # (optional) defines the color breaks manually for a "skewed" color transition
    col_breaks = c(seq(min(plotdf),median(plotdf),length=100),  # for orange
                   seq(median(plotdf) + 0.1, max(plotdf),length=100))  # for blue
    
    # Plot means heatmap (without group clustering)
    png(paste0(selected_dir,
               "/",
               i,
               "e-",
               gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
               "_logFC",
               "_fdr_",
               fdr,               
               ".png"),    # create png for the heat map        
        width = 6.5*600,        # 6.5 x 600 pixels
        if (length(rownames(plotdf)) > 5) {height = 5*600
        } else {height = 3.5*600},
        res = 600,            # 600 pixels per inch
        #paper = "a4r",
        pointsize = 7)        # smaller font size
    heatmap.2(plotdf,
              main = paste0(gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", names(filtered.gene.lists)[i])), # heat map title
              notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(8,1),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier
              breaks=col_breaks,    # enable color transition at specified limits
              dendrogram="none",    # only draw a row dendrogram
              Colv = "NA",            # turn off column clustering 
              Rowv = TRUE,
              #scale = "row",
              useRaster=TRUE,
              labRow = NA,
              labCol = labels.logFC.values,
              cexRow = 0.2 + 1/log10(length(rownames(plotdf))),
              cexCol = 0.2 + 1/log10(length(colnames(plotdf))),
              # ColSideColors = NA,
              key.xlab = "Scaled logFC",
              keysize = 1,
              lhei = c(2,7))
    dev.off()
    write_delim(cbind("Gene_ID" = rownames(plotdf), as.data.frame(plotdf)),
                file = paste0(selected_dir,
                              "/",
                              i,
                              "e-",
                              gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
                              "_logFC",
                              "_fdr_",
                              fdr,               
                              ".tab"),
                delim = "\t", na = "NA")
    
    # Plot means heatmap (with group clustering)
    png(paste0(selected_dir,
               "/",
               i,
               "f-",
               gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
               "_logFC_with_group_clustering",
               "_fdr_",
               fdr,               
               ".png"),    # create png for the heat map        
        width = 6.5*600,        # 6.5 x 600 pixels
        if (length(rownames(plotdf)) > 5) {height = 5*600
        } else {height = 3.5*600},
        res = 600,            # 600 pixels per inch
        #paper = "a4r",
        pointsize = 7)        # smaller font size
    heatmap.2(plotdf,
              main = paste0(gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", names(filtered.gene.lists)[i])), # heat map title
              notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(8,1),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier
              breaks=col_breaks,    # enable color transition at specified limits
              dendrogram="column",    # only draw a row dendrogram
              Colv = TRUE,            # turn off column clustering 
              Rowv = TRUE,
              #scale = "row",
              useRaster=TRUE,
              labRow = NA,
              labCol = labels.logFC.values,
              cexRow = 0.2 + 1/log10(length(rownames(plotdf))),
              cexCol = 0.2 + 1/log10(length(colnames(plotdf))),
              # ColSideColors = c(rep(c(rep("gray",1), rep("#A00000",1)),6)),
              key.xlab = "Scaled logFC",
              keysize = 1,
              lhei = c(2,7))
    dev.off()
    write_delim(cbind("Gene_ID" = rownames(plotdf), as.data.frame(plotdf)),
                file = paste0(selected_dir,
                              "/",
                              i,
                              "f-",
                              gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
                              "_logFC _with_group_clustering",
                              "_fdr_",
                              fdr,               
                              ".tab"),
                delim = "\t", na = "NA")
  }
  
}

## repeat but unscaled
dir.create("results/heatmaps/selected_genes/unscaled")

for (i in seq_along(filtered.gene.lists)) {
  genes_to_plot <- unlist(filtered.gene.lists[i])
  print(genes_to_plot)
  
  selected_dir <- paste0("results/heatmaps/selected_genes/unscaled")
  
  gene_subset <- unique(unlist(genes_to_plot)) # obtain all genes belonging to that category
  gene_subset <- gene_subset[which(gene_subset %in% rownames(df2.filtered.unscaled))] # obtain the subset that has been differentially regulated
  
  if (length(gene_subset) >= 2){
    # Select subset of genes to be graphed
    plotdf <- df2.filtered.unscaled[which(rownames(df2.filtered.unscaled) %in% gene_subset),1:47]
    
    if (names(filtered.gene.lists)[i] != "all expressed genes") {
      plotdf.qvalues <-full.degs.data[which(full.degs.data$Gene_ID %in% gene_subset),c(1,4,7,10,13,16,19,22,25,28)]
      significant_degs <- as.data.frame(plotdf.qvalues) %>% filter_all(any_vars(. < fdr))
      plotdf <- plotdf[which(rownames(plotdf) %in% significant_degs$Gene_ID),]
    }
    
    # Remove rows with zero values
    # Go through each row and determine if a value is zero
    row_sub = apply(plotdf, 1, function(row) all(row = 0 ))
    # Subset rows with non-zero values
    plotdf <- plotdf[!row_sub,]
    
    # Hierarchical clustering of rearranged logCPM data
    d <- dist(plotdf, method = "euclidean") # calculate euclidean distance
    clust <- hclust(d, method = "complete") # hierarchical cluster of the data
    dend <- as.dendrogram(clust) # as dendogram (for clustering rows, Rowv parameter of heatmap.2)
    
    # creates a own color palette from red to green
    my_palette <- colorRampPalette(c("blue", "#ffc000"))(n = 199)
    
    # (optional) defines the color breaks manually for a "skewed" color transition
    col_breaks = c(seq(min(plotdf),median(plotdf),length=100),  # for orange
                   seq(median(plotdf) + 0.1, max(plotdf),length=100))  # for blue
    
    
    # Plot heatmap (without group clustering)
    png(paste0(selected_dir,
               "/",
               i,
               "a-",
               gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
               "_fdr_",
               fdr,
               ".png"),    # create png for the heat map        
        width = 6.5*600,        
        if (length(rownames(plotdf)) > 5) {height = 5*600
        } else {height = 3.5*600},
        res = 600,
        pointsize = 5)        # smaller font size
    heatmap.2(plotdf,
              main = paste0(gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", names(filtered.gene.lists)[i])), # heat map title
              notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(5,1),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier
              breaks=col_breaks,    # enable color transition at specified limits
              dendrogram="none",    # only draw a row dendrogram
              Colv="NA",            # turn off column clustering 
              Rowv = TRUE,
              #scale = "row",
              useRaster=TRUE,
              labRow = NA,
              labCol = labels.all.replicates,
              cexRow = 0.2 + 1/log10(length(rownames(plotdf))),
              cexCol = 0.2 + 1/log10(length(labels.all.replicates)),
              ColSideColors = c(rep("gray",4),
                                rep("#A00000",4),
                                rep("gray",3),
                                rep("#A00000",4),
                                rep(c(rep("gray",4), rep("#A00000",4)),4)
              ),
              key.xlab = "Scaled logCPM",
              keysize = 1,
              lhei = c(1,7))
    dev.off()
    write_delim(cbind("Gene_ID" = rownames(plotdf), as.data.frame(plotdf)),
                file = paste0(selected_dir,
                              "/",
                              i,
                              "a-",
                              gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
                              "_fdr_",
                              fdr,               
                              ".tab"),
                delim = "\t", na = "NA")
    
    # Plot heatmap (with group clustering)
    png(paste0(selected_dir,
               "/",
               i,
               "b-",
               gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
               "_with_group_clustering",
               "_fdr_",
               fdr,               
               ".png"),    # create png for the heat map        
        width = 6.5*600,        
        if (length(rownames(plotdf)) > 5) {height = 5*600
        } else {height = 3.5*600},
        res = 600,
        pointsize = 5)        # smaller font size
    heatmap.2(plotdf,
              main = paste0(gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", names(filtered.gene.lists)[i])), # heat map title
              notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(5,1),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier
              breaks=col_breaks,    # enable color transition at specified limits
              dendrogram="column",    # only draw a row dendrogram
              Colv=TRUE,            # turn off column clustering 
              Rowv = TRUE,
              #scale = "row",
              useRaster=TRUE,
              labRow = NA,
              labCol = labels.all.replicates,
              cexRow = 0.2 + 1/log10(length(rownames(plotdf))),
              cexCol = 0.2 + 1/log10(length(labels.all.replicates)),
              ColSideColors = c(rep("gray",4),
                                rep("#A00000",4),
                                rep("gray",3),
                                rep("#A00000",4),
                                rep(c(rep("gray",4), rep("#A00000",4)),4)
              ),
              key.xlab = "Scaled logCPM",
              keysize = 1,
              lhei = c(1,7))
    dev.off()
    write_delim(cbind("Gene_ID" = rownames(plotdf), as.data.frame(plotdf)),
                file = paste0(selected_dir,
                              "/",
                              i,
                              "b-",
                              gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
                              "_with_group_clustering",
                              "_fdr_",
                              fdr,               
                              ".tab"),
                delim = "\t", na = "NA")
  }
  
  # Redo to graph only mean values by group
  gene_subset <- unique(unlist(genes_to_plot)) # obtain all genes belonging to that category
  gene_subset <- gene_subset[which(gene_subset %in% rownames(df4.filtered.unscaled))] # obtain the subset that has been differentially regulated
  
  if (length(gene_subset) >= 2){
    # Select subset of genes to be graphed
    plotdf <- df4.filtered.unscaled[which(rownames(df4.filtered.unscaled) %in% gene_subset),1:12]
    
    if (names(filtered.gene.lists)[i] != "all expressed genes") {
      plotdf.qvalues <-full.degs.data[which(full.degs.data$Gene_ID %in% gene_subset),c(1,4,7,10,13,16,19,22,25,28)]
      significant_degs <- as.data.frame(plotdf.qvalues) %>% filter_all(any_vars(. < fdr))
      plotdf <- plotdf[which(rownames(plotdf) %in% significant_degs$Gene_ID),]
    }
    
    # Remove rows with zero values
    # Go through each row and determine if a value is zero
    row_sub = apply(plotdf, 1, function(row) all(row = 0 ))
    # Subset rows with non-zero values
    plotdf <- plotdf[!row_sub,]
    
    # Hierarchical clustering of rearranged logCPM data
    d <- dist(plotdf, method = "euclidean") # calculate euclidean distance
    clust <- hclust(d, method = "complete") # hierarchical cluster of the data
    dend <- as.dendrogram(clust) # as dendogram (for clustering rows, Rowv parameter of heatmap.2)
    
    # creates a own color palette from red to green
    my_palette <- colorRampPalette(c("blue", "#ffc000"))(n = 199)
    
    # (optional) defines the color breaks manually for a "skewed" color transition
    col_breaks = c(seq(min(plotdf),median(plotdf),length=100),  # for orange
                   seq(median(plotdf) + 0.1, max(plotdf),length=100))  # for blue
    
    # Plot means heatmap (without group clustering)
    png(paste0(selected_dir,
               "/",
               i,
               "c-",
               gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
               "_group_means",
               "_fdr_",
               fdr,               
               ".png"),    # create png for the heat map        
        width = 6.5*600,        # 6.5 x 600 pixels
        if (length(rownames(plotdf)) > 5) {height = 5*600
        } else {height = 3.5*600},
        res = 600,            # 600 pixels per inch
        #paper = "a4r",
        pointsize = 7)        # smaller font size
    heatmap.2(plotdf,
              main = paste0(gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", names(filtered.gene.lists)[i])), # heat map title
              notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(8,1),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier
              breaks=col_breaks,    # enable color transition at specified limits
              dendrogram="none",    # only draw a row dendrogram
              Colv = "NA",            # turn off column clustering 
              Rowv = TRUE,
              #scale = "row",
              useRaster=TRUE,
              labRow = NA,
              labCol = labels.mean.values,
              cexRow = 0.2 + 1/log10(length(rownames(plotdf))),
              cexCol = 0.2 + 1/log10(length(colnames(plotdf))),
              ColSideColors = c(rep(c(rep("gray",1), rep("#A00000",1)),6)),
              key.xlab = "Scaled logCPM",
              keysize = 1,
              lhei = c(2,7))
    dev.off()
    write_delim(cbind("Gene_ID" = rownames(plotdf), as.data.frame(plotdf)),
                file = paste0(selected_dir,
                              "/",
                              i,
                              "c-",
                              gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
                              "_group_means",
                              "_fdr_",
                              fdr,               
                              ".tab"),
                delim = "\t", na = "NA")
    
    # Plot means heatmap (with group clustering)
    png(paste0(selected_dir,
               "/",
               i,
               "d-",
               gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
               "_group_means_with_group_clustering",
               "_fdr_",
               fdr,               
               ".png"),    # create png for the heat map        
        width = 6.5*600,        # 6.5 x 600 pixels
        if (length(rownames(plotdf)) > 5) {height = 5*600
        } else {height = 3.5*600},
        res = 600,            # 600 pixels per inch
        #paper = "a4r",
        pointsize = 7)        # smaller font size
    heatmap.2(plotdf,
              main = paste0(gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", names(filtered.gene.lists)[i])), # heat map title
              notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(8,1),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier
              breaks=col_breaks,    # enable color transition at specified limits
              dendrogram="column",    # only draw a row dendrogram
              Colv = TRUE,            # turn off column clustering 
              Rowv = TRUE,
              #scale = "row",
              useRaster=TRUE,
              labRow = NA,
              labCol = labels.mean.values,
              cexRow = 0.2 + 1/log10(length(rownames(plotdf))),
              cexCol = 0.2 + 1/log10(length(colnames(plotdf))),
              ColSideColors = c(rep(c(rep("gray",1), rep("#A00000",1)),6)),
              key.xlab = "Scaled logCPM",
              keysize = 1,
              lhei = c(2,7))
    dev.off()
    write_delim(cbind("Gene_ID" = rownames(plotdf), as.data.frame(plotdf)),
                file = paste0(selected_dir,
                              "/",
                              i,
                              "d-",
                              gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
                              "_group_means_with_group_clustering",
                              "_fdr_",
                              fdr,               
                              ".tab"),
                delim = "\t", na = "NA")
  }
  
  # Redo to graph only logFC values by group
  gene_subset <- unique(unlist(genes_to_plot)) # obtain all genes belonging to that category
  gene_subset <- gene_subset[which(gene_subset %in% rownames(df5.filtered.unscaled))] # obtain the subset that has been differentially regulated
  
  if (length(gene_subset) >= 2){
    # Select subset of genes to be graphed
    plotdf <- df5.filtered.unscaled[which(rownames(df5.filtered.unscaled) %in% gene_subset),]
    
    if (names(filtered.gene.lists)[i] != "all expressed genes") {
      plotdf.qvalues <-full.degs.data[which(full.degs.data$Gene_ID %in% gene_subset),c(1,4,7,10,13,16,19,22,25,28)]
      significant_degs <- as.data.frame(plotdf.qvalues) %>% filter_all(any_vars(. < fdr))
      plotdf <- plotdf[which(rownames(plotdf) %in% significant_degs$Gene_ID),]
    }    
    
    # Remove rows with zero values
    # Go through each row and determine if a value is zero
    row_sub = apply(plotdf, 1, function(row) all(row = 0 ))
    # Subset rows with non-zero values
    plotdf <- plotdf[!row_sub,]
    
    # Hierarchical clustering of rearranged logCPM data
    d <- dist(plotdf, method = "euclidean") # calculate euclidean distance
    clust <- hclust(d, method = "complete") # hierarchical cluster of the data
    dend <- as.dendrogram(clust) # as dendogram (for clustering rows, Rowv parameter of heatmap.2)
    
    # creates a own color palette from red to green
    my_palette <- colorRampPalette(c("blue", "#ffc000"))(n = 199)
    
    # (optional) defines the color breaks manually for a "skewed" color transition
    col_breaks = c(seq(min(plotdf),median(plotdf),length=100),  # for orange
                   seq(median(plotdf) + 0.1, max(plotdf),length=100))  # for blue
    
    # Plot means heatmap (without group clustering)
    png(paste0(selected_dir,
               "/",
               i,
               "e-",
               gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
               "_logFC",
               "_fdr_",
               fdr,               
               ".png"),    # create png for the heat map        
        width = 6.5*600,        # 6.5 x 600 pixels
        if (length(rownames(plotdf)) > 5) {height = 5*600
        } else {height = 3.5*600},
        res = 600,            # 600 pixels per inch
        #paper = "a4r",
        pointsize = 7)        # smaller font size
    heatmap.2(plotdf,
              main = paste0(gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", names(filtered.gene.lists)[i])), # heat map title
              notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(8,1),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier
              breaks=col_breaks,    # enable color transition at specified limits
              dendrogram="none",    # only draw a row dendrogram
              Colv = "NA",            # turn off column clustering 
              Rowv = TRUE,
              #scale = "row",
              useRaster=TRUE,
              labRow = NA,
              labCol = labels.logFC.values,
              cexRow = 0.2 + 1/log10(length(rownames(plotdf))),
              cexCol = 0.2 + 1/log10(length(colnames(plotdf))),
              # ColSideColors = NA,
              key.xlab = "Scaled logFC",
              keysize = 1,
              lhei = c(2,7))
    dev.off()
    write_delim(cbind("Gene_ID" = rownames(plotdf), as.data.frame(plotdf)),
                file = paste0(selected_dir,
                              "/",
                              i,
                              "e-",
                              gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
                              "_logFC",
                              "_fdr_",
                              fdr,               
                              ".tab"),
                delim = "\t", na = "NA")
    
    # Plot means heatmap (with group clustering)
    png(paste0(selected_dir,
               "/",
               i,
               "f-",
               gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
               "_logFC_with_group_clustering",
               "_fdr_",
               fdr,               
               ".png"),    # create png for the heat map        
        width = 6.5*600,        # 6.5 x 600 pixels
        if (length(rownames(plotdf)) > 5) {height = 5*600
        } else {height = 3.5*600},
        res = 600,            # 600 pixels per inch
        #paper = "a4r",
        pointsize = 7)        # smaller font size
    heatmap.2(plotdf,
              main = paste0(gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", names(filtered.gene.lists)[i])), # heat map title
              notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(8,1),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier
              breaks=col_breaks,    # enable color transition at specified limits
              dendrogram="column",    # only draw a row dendrogram
              Colv = TRUE,            # turn off column clustering 
              Rowv = TRUE,
              #scale = "row",
              useRaster=TRUE,
              labRow = NA,
              labCol = labels.logFC.values,
              cexRow = 0.2 + 1/log10(length(rownames(plotdf))),
              cexCol = 0.2 + 1/log10(length(colnames(plotdf))),
              # ColSideColors = c(rep(c(rep("gray",1), rep("#A00000",1)),6)),
              key.xlab = "Scaled logFC",
              keysize = 1,
              lhei = c(2,7))
    dev.off()
    write_delim(cbind("Gene_ID" = rownames(plotdf), as.data.frame(plotdf)),
                file = paste0(selected_dir,
                              "/",
                              i,
                              "f-",
                              gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", gsub(" ", "_", names(filtered.gene.lists)[i])),
                              "_logFC _with_group_clustering",
                              "_fdr_",
                              fdr,               
                              ".tab"),
                delim = "\t", na = "NA")
  }
  
}


#################################################
# 9. Get G profiler GO terms                    #
#################################################

library(gprofiler2)
for (i in seq_along(all.degs.sets)){
  print(i)
  gostres <- gost(query = c(md$RAP.DB[which(md$Gene_ID %in% all.degs.sets[[i]])]),
                  custom_bg = c(md$RAP.DB[which(md$Gene_ID %in% logcounts2$Gene_ID)]),
                  organism = "osativa",
                  significant = FALSE)
  write_tsv(gostres$result, file = paste0("DE/GO/gprofiler_", names(all.degs.sets[i]),"_",format(Sys.Date(), "%Y-%m-%e"),".tab"))
}


## Save Session Information
setwd(basedir)
sessionInfo()
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

#####################################################################################
#                               End of analysis                                     #
#####################################################################################