# Function to check for NaN values in data frames
is.nan.data.frame <- function(x){
  do.call(cbind, lapply(x, is.nan))}

# Function to check for infinite values in data frames
is.infinite.data.frame <- function(x){
  do.call(cbind, lapply(x, is.infinite))}

# Rounding functions
floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)
round_dec <- function(x, level=1) round(x+ 5*10^(-level-1), level)


# Filter a gene counts dataframe by read density
#
#' \dontrun{
#' filterByRd(df, targets, min, type)
#' df is a dataframe containing counts data
#' targets is the file describing sample, bam and conditions
#' min is the minimum read density that the gene is required to have
#' type can be "any" (at least one sample) or "all" (all samples)
#' }
filterByRd<-function(df0, targets, min, type) 
{
  a<-table(targets$condition) 
  start=ncol(df0)-sum(a)+1
  genes.rd=df0[,start:ncol(df0)]/df0$effective_length
  write.table(genes.rd, file="genes.rd.txt")
  colnames(genes.rd)<-targets$condition
  #independiente del numero de condiciones
  list<-matrix(unlist(
    lapply( unique(colnames(genes.rd)), 
            function(x) rowMeans(genes.rd[,colnames(genes.rd) == x])  > min )),  
    nrow=nrow(genes.rd), 
    byrow=FALSE)
  
  if (type=="any")  { 
    ii=rowSums(list)>0  
    df=df0[ii,]
  }else  { 
    ii=rowSums(list)==ncol(list)
    df=df0[ii,]
  }
  return (df)
  
}


# Function to get annotation data
getAnnnot <- function(df, md){
  #' \dontrun{
  #' getAnnot(df = dataframe where first column is named "Gene ID" and contains genes,
  #'          md = metadata where first column is named "Gene ID" and contains genes and has extra columns with annotations  
  #'          to add to the dataframe)
  #' }
  if (colnames(md)[1] != "Gene_ID") {
    colnames(md)[1] <- "Gene_ID"
  }
  if (colnames(df)[1] != "Gene_ID") {
    df <- data.frame(cbind(rownames(df)),df)
    colnames(df)[1] <- "Gene_ID"
  }
  res <- merge(x = md, y = df, by = "Gene_ID", all.y = TRUE) # Left inner join of df and getAnnnot <- function(df, md){
  res[is.na(res)] <- "-"
  return(res)
}


# Function to plot a heatmap using heatmap.2
plotHM <- function(df,                    # dataframe to plot
                   title = "",            # Graph title
                   oma=c(2,2,2,12),       # (bottom,left,top,right)
                   margins =c(10,4),      # margins for column and row names
                   my_palette = "",       # Color palette to use
                   dendrogram = "none",   # Dendrogram 
                   Colv = NULL,          # Column clustering
                   labCol = "",           # Column labels
                   colSideColors = "",    # Column colours
                   Rowv = NULL,          # Row clustering
                   labRow = NULL)         # Now row labels
{
  if (!require("gplots")) {
    install.packages("gplots", dependencies = TRUE)
    library(dplyr)
  }
  if (!require("RColorBrewer")) {
    install.packages("RColorBrewer", dependencies = TRUE)
    library(RColorBrewer)
  }
  
  if (my_palette == "") { my_palette = colorRampPalette(c("blue", "yellow"))(n = 199)} # Color palette to use
  
  if (is.na(labCol)){ labCol = colnames(df)} # Set column names
  
  par(oma=oma) # Set outer margins
  HMplot <-heatmap.2(df,
                     main = title,         # Set title
                     # density.info="none",  # turns off density plot inside color legend
                     trace="none",         # turns off trace lines inside the heat map
                     margins = margins,    # widens margins around plot
                     col=my_palette,       # use on color palette defined earlier
                     dendrogram = dendrogram,     # draw a row/column dendrogram
                     Colv=Colv,            # turn on/off column clustering 
                     useRaster=TRUE,
                     labCol = labCol,      # Set column names
                     Rowv = Rowv,
                     labRow = labRow,
                     ColSideColors = colSideColors)    # Set column colours
  return(HMplot)
}

# Function to check the operating system
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd+
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

writeGSEA <- function(subset, universe, PValue = 0.05, fdr = 0.05, annotation = "org.At.tair.db", ontology = "BP", description = "my", project = getwd(), sub_project = "", all = FALSE)
{
  
  ## Remove genes that have no entrezGene id
  subset <- as.character(subset)
  names(subset) <- subset
  entrezIds <- mget(subset, envir=org.At.tairENTREZID, ifnotfound=NA)
  haveEntrezId <- names(entrezIds)[sapply(entrezIds, function(x) !is.na(x))]
  numNoEntrezId <- length(subset) - length(haveEntrezId)
  subset <- subset[haveEntrezId]
  
  ## Remove genes from Universe that have no entrezGene id
  universe <- as.character(universe)
  names(universe) <- universe
  entrezIds <- mget(universe, envir=org.At.tairENTREZID, ifnotfound=NA)
  haveEntrezId <- names(entrezIds)[sapply(entrezIds, function(x) !is.na(x))]
  numNoEntrezId <- length(universe) - length(haveEntrezId)
  universe <- universe[haveEntrezId]
  
  ## Remove genes with no GO mapping
  haveGo <- sapply(mget(subset, org.At.tairGO),
                   function(x) {
                     if (length(x) == 1 && is.na(x)) 
                       FALSE 
                     else TRUE
                   })
  numNoGO <- sum(!haveGo)
  subset <- subset[haveGo]
  
  ## Remove universe genes with no GO mapping
  haveGo <- sapply(mget(universe, org.At.tairGO),
                   function(x) {
                     if (length(x) == 1 && is.na(x)) 
                       FALSE 
                     else TRUE
                   })
  numNoGO <- sum(!haveGo)
  universe <- universe[haveGo]
  
  params <- new("GOHyperGParams",
                geneIds=subset,
                universeGeneIds=universe,
                annotation="org.At.tair.db",
                ontology= ontology,
                pvalueCutoff=1,
                conditional=FALSE,
                testDirection="over")
  # paramsCond <- params
  # conditional(paramsCond) <- TRUE
  hgOver <- hyperGTest(params)
  # hgCondOver <- hyperGTest(paramsCond)
  hg <- hgOver
  
  ## Get the p-values of the test
  hg.pv <- pvalues(hg)
  length(hg.pv)
  ## Adjust p-values for multiple test (FDR)
  hg.pv.fdr <- p.adjust(hg.pv, method="BH")
  length(hg.pv.fdr)
  ## select the GO terms with adjusted p-value less than the cut off
  sigGO.ID <- names(hg.pv.fdr[hg.pv.fdr <= fdr])
  length(sigGO.ID)
  # get table from HyperG test result
  df <- summary(hg)
  length(df[,1])
  # keep only significant GO terms in the table
  # GOannot.table <- df[df[,1] %in% sigGO.ID,]
  GOannot.table <- cbind(df[, 1:2],
                         "FDR" = hg.pv.fdr[names(hg.pv.fdr) %in% df[,1]],
                         df[, 3:5],
                         "Subset" = geneMappedCount(hg),
                         "Expected" = df[, 6],
                         "Universe" =  universeMappedCount(hg),
                         "Enrichment" = (df[,5]/geneMappedCount(hg))/(df[,6]/universeMappedCount(hg)),
                         "Term" = df[, 7])
  
  if (all == FALSE) {
    GOannot.table <- GOannot.table[GOannot.table[,1] %in% sigGO.ID,]
  }
  
  # Create text report of the significantly over-represented GO terms
  write.table(GOannot.table,file=paste0(project, "/", sub_project, "/", description, "_GOterms_OverRep_",ontology,".txt"),sep="\t",row.names=F)
  # Create html report of all over-represented GO terms
  htmlReport(hg, file=paste0(project, "/", sub_project, "/", description, "_GOterms_OverRep_", ontology,".html"))
  print(hg)
  print(paste0("Only ", length(sigGO.ID), " with FDR less than ", fdr))
  print(paste0("Saving output to ", project, "/", sub_project, "/", description, "_GOterms_OverRep_", ontology, ".txt"))
}

drawKEGG <- function(genes_KEGGresult, description, destdir, gene.data, pathways){
  old.dir <- getwd()
  dir.create(destdir)
  setwd(destdir)
  
  for (i in c(as.character(genes_KEGGresult$pathway.code))) {
    if (i != "ath01100") {
      print(pathways[paste0("path:",i)])
      map <- i
      pv.out <- pathview(gene.data = gene.data,
                         gene.idtype = "KEGG",
                         pathway.id = map,
                         species = "ath",
                         out.suffix = map,
                         keys.align = "y",
                         kegg.native = T,
                         match.data = T, key.pos = "topright",
                         map.symbol = TRUE,
                         low = "#0070C0", mid = "white", high = "#FFC000")
      plot.name <- paste0(map,".", map, ".", description," - ", gsub("([.|()\\^{}+$*?]|\\[|\\]|/|)", "",(pathways[paste0("path:ath",i)]),2), ".png")
      if (!is.atomic(pv.out)){
        cat(capture.output(print(pv.out$plot.data.gene)), file = paste0(map,".", map, ".", description," - ", gsub("([.|()\\^{}+$*?]|\\[|\\]|/|)", "",(pathways[paste0("path:",i)]),2), ".txt"))
      }
    }
  }
  setwd(old.dir)
}

drawKEGGpathway <- function(map = "", description = "", destdir = "", gene.data = NULL, pathways = ""){
  old.dir <- getwd()
  dir.create(destdir)
  setwd(destdir)
  if (substr(map,1,3) == "ath"){
    map <- substr(map,4,str_length(map))
  }
  
  print(pathways[paste0("path:ath",map)])
  pv.out <- pathview(gene.data = gene.data,
                     gene.idtype = "KEGG",
                     pathway.id = paste0("ath",map),
                     species = "ath",
                     out.suffix = paste0(map, c(description)),
                     keys.align = "y",
                     kegg.native = T,
                     match.data = T, 
                     key.pos = "topright",
                     map.symbol = TRUE,
                     low = "#0070C0", mid = "white", high = "#FFC000")
  plot.name <- paste0(map,".", map, ".", description," - ", gsub("([.|()\\^{}+$*?]|\\[|\\]|/|)", "",(pathways[paste0("path:ath", map)]),2), ".png")
  if (!is.atomic(pv.out)){
    cat(capture.output(print(pv.out$plot.data.gene)), file = paste0(map,".", map, ".", description," - ", gsub("([.|()\\^{}+$*?]|\\[|\\]|/|)", "",(pathways[paste0("path:ath", map)]),2), ".txt"))
  }
  
  setwd(old.dir)
}

pathwayEnrichment <- function(gene.data, genes.by.pathway, pathways.list, universe){
  pVals.by.pathway <- t(sapply(names(genes.by.pathway),
                               function(pathway) {
                                 pathway.genes <- intersect(universe, genes.by.pathway[[pathway]])
                                 list.genes.in.pathway <- intersect(names(gene.data), pathway.genes)
                                 list.genes.not.in.pathway <- setdiff(names(gene.data), list.genes.in.pathway)
                                 scores.in.pathway <- gene.data[list.genes.in.pathway]
                                 scores.not.in.pathway <- gene.data[list.genes.not.in.pathway]
                                 if (length(scores.in.pathway) > 0){
                                   p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
                                 } else{
                                   p.value <- NA
                                 }
                                 return(c(p.value = p.value, Annotated = length(list.genes.in.pathway), Total = length(genes.by.pathway[[pathway]]) ))
                               }
  ))
  # Assemble output table
  outdat <- data.frame(pathway.code = rownames(pVals.by.pathway))
  outdat$pathway.name <- pathways.list[paste0("path:", outdat$pathway.code)]
  outdat$p.value <- pVals.by.pathway[,"p.value"]
  # outdat$fdr <- p.adjust(outdat$p.value, method="BH") # Calculate FDR
  # outdat$p.adj.holm <- p.adjust(outdat$p.value, method="holm") # Calculate FDR
  outdat$annotated <- pVals.by.pathway[,"Annotated"]
  outdat$total <- pVals.by.pathway[,"Total"]
  outdat <- outdat[order(outdat$p.value),]
  
  return(outdat)
}

getGenesInPathway <- function(gene.data, genes.by.pathway, pathways.list){
  res <- data.frame()
  h <- 0
  for (i in c(pathways.list)) {
    pathway <- i 
    pathway.genes <- genes.by.pathway[[pathway]]
    list.genes.in.pathway <- intersect(names(gene.data), pathway.genes)
    h <- h + 1
    for (j in seq_along(list.genes.in.pathway)){
      res[j,h] <- list.genes.in.pathway[j]
    }
    colnames(res)[h] <- i    
  }
  res[is.na(res)] <- ""
  
  return(res)
}


