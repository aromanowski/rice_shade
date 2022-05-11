if(!require(ggrepel)) install.packages(ggrepel)
  
library(edgeR)
library(ggrepel)

# set base directory  ### modify accordingly ###
setwd(basedir)

targets <- read.table("../data/samples_shortxLuk.txt", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
cg <- read.delim("../data/Master_counts_swap.txt", sep = "\t")
lev_full <- levels(factor(str_c(targets$variety, "_", targets$treatment)))
group_full <- factor(str_c(targets$variety, "_", targets$treatment), levels = lev_full)
design_full <- model.matrix(~0+group_full)
colnames(design_full) <- levels(group_full)
ncol(design_full)
qr(design_full)$rank


df_full <- as.matrix(cg[,setdiff(colnames(cg), c("Luk.C.R1"))])

col.cell <- c("violet", "purple", "blue", "cyan", "green", "dark gray","orange", "pink", "magenta", "red", "dark red", "brown", "gray", "black")[factor(str_c(targets$variety, "_", targets$treatment))]
# data.frame(str_c(targets$variety, "_", targets$treatment),col.cell)

dir.create(paste0(basedir,"/multivariate_analysis_by_variety"))
setwd(paste0(basedir,"/multivariate_analysis_by_variety"))

for (i in c(1,9,16,24,32,40)) {
  
  print(paste0("col start: ", i))
  
  if (i == 9){
    df <- df_full[,i:(i+6)]
    group <-group_full[i:(i+6)]
  }
  else {
  df <- df_full[,i:(i+7)]
  group <-group_full[i:(i+7)]
  }
  
  rownames(df) <- rownames(cg)
  
  variety_name <- targets$variety[i]
  
  # print(colnames(df))
        
  y <- DGEList(counts= df ,
               group=group) # make a DGEList object with the subgroup count data
  if (i == 9){
  colnames(y)=rownames(targets)[i:(i+6)] # Rename the columns
  }
  else {
    colnames(y)=rownames(targets)[i:(i+7)] # Rename the columns
  }
    
  y$samples$lib.size = colSums(y$counts) # Recalculate library sizes
  
  logcounts <- cpm(y,log=TRUE)
  
  png(file=paste0(variety_name,"-raw-histogram_logCPM.png"),    # create PNG for the histogram      
      width = 8*600,        # 10 x 600 pixels
      height = 5*600,
      res = 600,            # 300 pixels per inch
      pointsize = 8)        # smaller font size))
  par(mar=c(9,4.1,4.1,2.1))
  hist(logcounts)
  dev.off()
  
  png(filename = paste0(variety_name,"-raw-BoxPlot.png"),    # create PNG for the BoxPlot        
      width = 10*600,        # 10 x 600 pixels
      height = 5*600,
      res = 600,            # 300 pixels per inch
      pointsize = 8)
  par(mar=c(9,4.1,4.1,2.1))
  boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
  # Let's add a blue horizontal line that corresponds to the median logCPM
  abline(h=median(logcounts),col="blue")
  title(paste0("Boxplots of logCPMs (unnormalised) - ", variety_name))
  dev.off()
  
  png(filename = paste0(variety_name,"-raw-MDS_plot.png"),    # create PNG for the MDS        
      width = 10*600,        # 10 x 600 pixels
      height = 5*600,
      res = 600,            # 300 pixels per inch
      pointsize = 8)        # smaller font size)
  plotMDS(y,col=col.cell)
  title("MDS by Sample")
  dev.off()
  
  ### filter by 10 CPM in 0.7 within one group (at least) 
  keep <- filterByExpr(y, group = group)
  table(keep)
  y.filtered <- y[keep, , keep.lib.sizes=FALSE]
  
  y.filtered$samples$lib.size = colSums(y$counts) # Recalculate library sizes
  
  logcounts.filtered <- cpm(y.filtered,log=TRUE)
  
  png(file=paste0(variety_name,"-fil-histogram_logCPM.png"),    # create PNG for the histogram        
      width = 8*600,        # 10 x 600 pixels
      height = 5*600,
      res = 600,            # 300 pixels per inch
      pointsize = 8)        # smaller font size))
  par(mar=c(9,4.1,4.1,2.1))
  hist(logcounts.filtered)
  dev.off()
  
  png(filename = paste0(variety_name,"-fil-BoxPlot.png"),    # create PNG for the BoxPlot        
      width = 10*600,        # 10 x 600 pixels
      height = 5*600,
      res = 600,            # 300 pixels per inch
      pointsize = 8)
  par(mar=c(9,4.1,4.1,2.1))
  boxplot(logcounts.filtered, xlab="", ylab="Log2 counts per million",las=2)
  # Let's add a blue horizontal line that corresponds to the median logCPM
  abline(h=median(logcounts.filtered),col="blue")
  title(paste0("Boxplots of logCPMs (filtered) - ", variety_name))
  dev.off()
  
  png(filename = paste0(variety_name,"-fil-MDS_plot.png"),    # create PNG for the MDS        
      width = 10*600,        # 10 x 600 pixels
      height = 5*600,
      res = 600,            # 300 pixels per inch
      pointsize = 8)        # smaller font size)
  plotMDS(y.filtered,col=col.cell)
  title("MDS by Sample")
  dev.off()
  
  y.filtered.norm <-  calcNormFactors(y.filtered, method = "TMM")
  logcounts.filtered.norm <- cpm(y.filtered.norm,log=TRUE)
  
  png(file=paste0(variety_name,"-norm-histogram_logCPM.png"),    # create PNG for the histogram        
      width = 8*600,        # 10 x 600 pixels
      height = 5*600,
      res = 600,            # 300 pixels per inch
      pointsize = 8)        # smaller font size))
  par(mar=c(9,4.1,4.1,2.1))
  hist(logcounts.filtered.norm)
  dev.off()
  
  png(filename = paste0(variety_name,"-norm-BoxPlot.png"),    # create PNG for the BoxPlot        
      width = 10*600,        # 10 x 600 pixels
      height = 5*600,
      res = 600,            # 300 pixels per inch
      pointsize = 8)
  par(mar=c(9,4.1,4.1,2.1))
  boxplot(logcounts.filtered.norm, xlab="", ylab="Log2 counts per million",las=2)
  # Let's add a blue horizontal line that corresponds to the median logCPM
  abline(h=median(logcounts.filtered.norm),col="blue")
  title(paste0("Boxplots of logCPMs (normalised) - ", variety_name))
  dev.off()
  
  png(filename = paste0(variety_name,"-norm-MDS_plot.png"),    # create PNG for the MDS        
      width = 10*600,        # 10 x 600 pixels
      height = 5*600,
      res = 600,            # 300 pixels per inch
      pointsize = 8)        # smaller font size)
  plotMDS(y.filtered.norm,col=col.cell)
  title("MDS by Sample")
  dev.off()
  
  # Hierarchical clustering
  distance <- dist(t(cpm(y.filtered.norm, log = TRUE)), method = "euclidean") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
  clusters <- hclust(distance, method = "complete") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
  png(filename = paste0(variety_name,"-norm-clustering.png"),    # create PNG for the clustering        
      width = 10*600,        # 10 x 600 pixels
      height = 5*600,
      res = 600,            # 300 pixels per inch
      pointsize = 8)        # smaller font size)
  plot(clusters, labels=rownames(y.filtered.norm$samples))
  dev.off()
  
  # Principal component analysis (PCA) -------------
  pca.res <- prcomp(t(logcounts.filtered.norm), scale.=F, retx=T)
  pca.res.df <- as_tibble(pca.res$x)
  pc.var<-pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
  pc.per<-round(pc.var/sum(pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
  
  png(filename = paste0(variety_name,"-norm-PCA_percentage_of_variance_by_PC.png"),    # create PNG for the clustering        
      width = 10*600,        # 10 x 600 pixels
      height = 5*600,
      res = 600,            # 300 pixels per inch
      pointsize = 8)        # smaller font size)
  plot(pc.per)
  title(paste0("Percentage variance explained by each PC - ", variety_name))
  dev.off()
  
  sampleLabels <- rownames(y.filtered.norm$samples)
  
  pca.plot <- ggplot(pca.res.df) +
    aes(x=PC1, y=PC2, label=sampleLabels, color = group) +
    geom_point(size = 4) +
    geom_text_repel(size = 4) +
    # stat_ellipse() +
    xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
    ylab(paste0("PC2 (",pc.per[2],"%",")")) +
    labs(title=paste0("PCA plot - ", variety_name),
         caption=paste0("produced on ", Sys.time())) +
    coord_fixed() +
    theme_bw()
  
  png(filename = paste0(variety_name,"-norm-PCA.png"),    # create PNG for the clustering        
      width = 10*600,        # 10 x 600 pixels
      height = 5*600,
      res = 600,            # 300 pixels per inch
      pointsize = 8)        # smaller font size)
  print(pca.plot)
  dev.off()
  
  # Create a PCA 'small multiples' chart ----
  # this is another way to view PCA laodings to understand impact of each sample on each principal component
  pca.res.df <- pca.res$x[,1:6] %>% # note that this is the first time you've seen the 'pipe' operator from the magrittr package
    as_tibble() %>%
    add_column(sample = sampleLabels,
               group = group)
  
  pca.pivot <- pivot_longer(pca.res.df, # dataframe to be pivoted
                            cols = PC1:PC6, # column names to be stored as a SINGLE variable
                            names_to = "PC", # name of that new variable (column)
                            values_to = "loadings") # name of new variable (column) storing all the values (data)
  
  
  pca.small.multiples.plot <- ggplot(pca.pivot) +
    aes(x=sample, y=loadings, fill=group) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
    geom_bar(stat="identity") +
    facet_wrap(~PC) +
    labs(title="PCA 'small multiples' plot",
         caption=paste0("produced on ", Sys.time())) +
    theme_bw() +
    coord_flip()
  
  png(filename = paste0(variety_name,"-norm-PCA_small_multiples.png"),    # create PNG for the clustering        
      width = 15*600,        # 10 x 600 pixels
      height = 10*600,
      res = 600,            # 300 pixels per inch
      pointsize = 8)        # smaller font size)
  print(pca.small.multiples.plot)
  dev.off()
  
}

setwd(basedir)













