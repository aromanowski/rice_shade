###############################################
#             Requires                        #
###############################################
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(Cairo)) install.packages("Cairo")
if(!require(tidyverse)) install.packages("tidyverse")

###############################################
#             Includes                        #
###############################################
library(ggplot2)
library(Cairo)
library(tidyverse)
library(plyr)

###############################################
# Begin Bubble plots                          #
###############################################

# set base directory  ### modify accordingly ###
basedir <- getwd()
basedir <- (paste0(basedir,"/results"))
setwd(basedir)

# load GO terms
general_GO <- read_tsv("DE/GO/gprofiler_General_2023-05-19.tab")
IR64_GO<- read_tsv("DE/GO/gprofiler_IR64_2023-05-19.tab")
LukTakhar_GO<- read_tsv("DE/GO/gprofiler_LukTakhar_2023-05-19.tab")
MBlatec_GO<- read_tsv("DE/GO/gprofiler_MBlatec_2023-05-19.tab")
Mudgo_GO<- read_tsv("DE/GO/gprofiler_Mudgo_2023-05-19.tab")
Sabharaj_GO<- read_tsv("DE/GO/gprofiler_Sabharaj_2023-05-19.tab")
Zhenshan_GO<- read_tsv("DE/GO/gprofiler_Zhenshan_2023-05-19.tab")

# Get all processes/pathways that are significant in at least one variety
interesting_list <- Reduce(union, list(general_GO$term_id[which(general_GO$p_value < 0.05)],
                                           IR64_GO$term_id[which(IR64_GO$p_value < 0.05)],
                                           LukTakhar_GO$term_id[which(LukTakhar_GO$p_value < 0.05)],
                                           MBlatec_GO$term_id[which(MBlatec_GO$p_value < 0.05)],
                                           Mudgo_GO$term_id[which(Mudgo_GO$p_value < 0.05)],
                                           Sabharaj_GO$term_id[which(Sabharaj_GO$p_value < 0.05)],
                                           Zhenshan_GO$term_id[which(Zhenshan_GO$p_value < 0.05)]))

# Generate a dataframe with the required format for the bubble plots with the general shade enrichment dataset
general_shade_GO <- data.frame(Type = general_GO$source,
                               GO_Id = general_GO$term_id,
                               Term = general_GO$term_name,
                               x_var = "General shade",
                               y_var = c(1:1933),
                               intersection_size = general_GO$intersection_size[match(general_GO$term_id, general_GO$term_id)],
                               query_size = general_GO$query_size[match(general_GO$term_id, general_GO$term_id)],
                               term_size = general_GO$term_size[match(general_GO$term_id, general_GO$term_id)],
                               universe_size = general_GO$effective_domain_size[match(general_GO$term_id, general_GO$term_id)],
                               RF = ((general_GO$intersection_size/general_GO$query_size) / (general_GO$term_size/general_GO$effective_domain_size)),
                               p.value = general_GO$p_value)

# Get the significant list of BPs and CCs for the entire enrichment dataset
sig_list_general_BP <- unique(general_shade_GO$GO_Id[which(general_shade_GO$p.value < 0.05 & general_shade_GO$Type == "GO:BP")])
sig_list_general_CC <- unique(general_shade_GO$GO_Id[which(general_shade_GO$p.value < 0.05 & general_shade_GO$Type == "GO:CC")])

# Generate a dataframe with the required format for the bubble plots with all enrichment datasets
merged_GO <- rbind.fill(data.frame(Type = general_GO$source,
                                   GO_Id = general_GO$term_id,
                                   Term = general_GO$term_name,
                                   x_var = "General shade",
                                   y_var = c(1:1933),
                                   intersection_size = general_GO$intersection_size[match(general_GO$term_id, general_GO$term_id)],
                                   query_size = general_GO$query_size[match(general_GO$term_id, general_GO$term_id)],
                                   term_size = general_GO$term_size[match(general_GO$term_id, general_GO$term_id)],
                                   universe_size = general_GO$effective_domain_size[match(general_GO$term_id, general_GO$term_id)],
                                   RF = ((general_GO$intersection_size/general_GO$query_size) / (general_GO$term_size/general_GO$effective_domain_size)),
                                   p.value = general_GO$p_value),
                   
                        data.frame(Type = general_GO$source,
                              GO_Id = general_GO$term_id,
                              Term = general_GO$term_name,
                              x_var = "IR64",
                              y_var = c(1:1933),
                              intersection_size = IR64_GO$intersection_size[match(general_GO$term_id, IR64_GO$term_id)],
                              query_size = IR64_GO$query_size[match(general_GO$term_id, IR64_GO$term_id)],
                              term_size = IR64_GO$term_size[match(general_GO$term_id, IR64_GO$term_id)],
                              universe_size = IR64_GO$effective_domain_size[match(general_GO$term_id, IR64_GO$term_id)],
                              RF = ((IR64_GO$intersection_size[match(general_GO$term_id, IR64_GO$term_id)]/IR64_GO$query_size[match(general_GO$term_id, IR64_GO$term_id)]) / (IR64_GO$term_size[match(general_GO$term_id, IR64_GO$term_id)]/general_GO$effective_domain_size)),
                              p.value = IR64_GO$p_value[match(general_GO$term_id, IR64_GO$term_id)]),
                        
                        data.frame(Type = general_GO$source,
                                   GO_Id = general_GO$term_id,
                                   Term = general_GO$term_name,
                                   x_var = "LukTakhar",
                                   y_var = c(1:1933),
                                   intersection_size = LukTakhar_GO$intersection_size[match(general_GO$term_id, LukTakhar_GO$term_id)],
                                   query_size = LukTakhar_GO$query_size[match(general_GO$term_id, LukTakhar_GO$term_id)],
                                   term_size = LukTakhar_GO$term_size[match(general_GO$term_id, LukTakhar_GO$term_id)],
                                   universe_size = LukTakhar_GO$effective_domain_size[match(general_GO$term_id, LukTakhar_GO$term_id)],
                                   RF = ((LukTakhar_GO$intersection_size[match(general_GO$term_id, LukTakhar_GO$term_id)]/LukTakhar_GO$query_size[match(general_GO$term_id, LukTakhar_GO$term_id)]) / (LukTakhar_GO$term_size[match(general_GO$term_id, LukTakhar_GO$term_id)]/general_GO$effective_domain_size)),
                                   p.value = LukTakhar_GO$p_value[match(general_GO$term_id, LukTakhar_GO$term_id)]),
                        
                        data.frame(Type = general_GO$source,
                                   GO_Id = general_GO$term_id,
                                   Term = general_GO$term_name,
                                   x_var = "MBlatec",
                                   y_var = c(1:1933),
                                   intersection_size = MBlatec_GO$intersection_size[match(general_GO$term_id, MBlatec_GO$term_id)],
                                   query_size = MBlatec_GO$query_size[match(general_GO$term_id, MBlatec_GO$term_id)],
                                   term_size = MBlatec_GO$term_size[match(general_GO$term_id, MBlatec_GO$term_id)],
                                   universe_size = MBlatec_GO$effective_domain_size[match(general_GO$term_id, MBlatec_GO$term_id)],
                                   RF = ((MBlatec_GO$intersection_size[match(general_GO$term_id, MBlatec_GO$term_id)]/MBlatec_GO$query_size[match(general_GO$term_id, MBlatec_GO$term_id)]) / (MBlatec_GO$term_size[match(general_GO$term_id, MBlatec_GO$term_id)]/general_GO$effective_domain_size)),
                                   p.value = MBlatec_GO$p_value[match(general_GO$term_id, MBlatec_GO$term_id)]),
                        
                        data.frame(Type = general_GO$source,
                                   GO_Id = general_GO$term_id,
                                   Term = general_GO$term_name,
                                   x_var = "Mudgo",
                                   y_var = c(1:1933),
                                   intersection_size = Mudgo_GO$intersection_size[match(general_GO$term_id, Mudgo_GO$term_id)],
                                   query_size = Mudgo_GO$query_size[match(general_GO$term_id, Mudgo_GO$term_id)],
                                   term_size = Mudgo_GO$term_size[match(general_GO$term_id, Mudgo_GO$term_id)],
                                   universe_size = Mudgo_GO$effective_domain_size[match(general_GO$term_id, Mudgo_GO$term_id)],
                                   RF = ((Mudgo_GO$intersection_size[match(general_GO$term_id, Mudgo_GO$term_id)]/Mudgo_GO$query_size[match(general_GO$term_id, Mudgo_GO$term_id)]) / (Mudgo_GO$term_size[match(general_GO$term_id, Mudgo_GO$term_id)]/general_GO$effective_domain_size)),
                                   p.value = Mudgo_GO$p_value[match(general_GO$term_id, Mudgo_GO$term_id)]),
                        
                        data.frame(Type = general_GO$source,
                                   GO_Id = general_GO$term_id,
                                   Term = general_GO$term_name,
                                   x_var = "Sabharaj",
                                   y_var = c(1:1933),
                                   intersection_size = Sabharaj_GO$intersection_size[match(general_GO$term_id, Sabharaj_GO$term_id)],
                                   query_size = Sabharaj_GO$query_size[match(general_GO$term_id, Sabharaj_GO$term_id)],
                                   term_size = Sabharaj_GO$term_size[match(general_GO$term_id, Sabharaj_GO$term_id)],
                                   universe_size = Sabharaj_GO$effective_domain_size[match(general_GO$term_id, Sabharaj_GO$term_id)],
                                   RF = ((Sabharaj_GO$intersection_size[match(general_GO$term_id, Sabharaj_GO$term_id)]/Sabharaj_GO$query_size[match(general_GO$term_id, Sabharaj_GO$term_id)]) / (Sabharaj_GO$term_size[match(general_GO$term_id, Sabharaj_GO$term_id)]/general_GO$effective_domain_size)),
                                   p.value = Sabharaj_GO$p_value[match(general_GO$term_id, Sabharaj_GO$term_id)]),
                        
                        data.frame(Type = general_GO$source,
                                   GO_Id = general_GO$term_id,
                                   Term = general_GO$term_name,
                                   x_var = "Zhenshan",
                                   y_var = c(1:1933),
                                   intersection_size = Zhenshan_GO$intersection_size[match(general_GO$term_id, Zhenshan_GO$term_id)],
                                   query_size = Zhenshan_GO$query_size[match(general_GO$term_id, Zhenshan_GO$term_id)],
                                   term_size = Zhenshan_GO$term_size[match(general_GO$term_id, Zhenshan_GO$term_id)],
                                   universe_size = Zhenshan_GO$effective_domain_size[match(general_GO$term_id, Zhenshan_GO$term_id)],
                                   RF = ((Zhenshan_GO$intersection_size[match(general_GO$term_id, Zhenshan_GO$term_id)]/Zhenshan_GO$query_size[match(general_GO$term_id, Zhenshan_GO$term_id)]) / (Zhenshan_GO$term_size[match(general_GO$term_id, Zhenshan_GO$term_id)]/general_GO$effective_domain_size)),
                                   p.value = Zhenshan_GO$p_value[match(general_GO$term_id, Zhenshan_GO$term_id)])
                   )

# Get the significant list of BPs and CCs for the entire enrichment dataset
sig_list_BP <- unique(merged_GO$GO_Id[which(merged_GO$p.value < 0.05 & merged_GO$Type == "GO:BP")])
sig_list_CC <- unique(merged_GO$GO_Id[which(merged_GO$p.value < 0.05 & merged_GO$Type == "GO:CC")])

# Graph general shade BP
df.bp<-general_shade_GO %>% select(everything()) %>% filter(GO_Id %in% c(sig_list_general_BP))
df.bp$y_var <- c(1:(length(df.bp$GO_Id)))

# assign a representation factor of 0.01 to those terms not present in the data
df.bp[which(is.na(df.bp$p.value)),"RF"] <- 0.01 
# assign NA to non-significant terms
df.bp[which(df.bp$p.value > 0.05),"p.value"] <- NA

# obtain the p.values and assign a score
score <- df.bp$p.value
score[is.na(df.bp$p.value)]<-NA
score[df.bp$p.value < 0.05]<-1
score[df.bp$p.value <= 0.01]<-2
score[df.bp$p.value <= 0.001]<-3
score[df.bp$p.value <= 0.0001]<-4

go.bp <-ggplot(df.bp, aes(x = x_var, y = y_var, labels = Term, size =RF, col = score)) +
  geom_point(alpha = 0.7) +
  scale_colour_gradient(limits= c(1, 4), low = "#ffc000",  high = "red", na.value = "light gray") +     # pink "#f6e7e4"
  scale_size_continuous(range = c(2,20)) +
  scale_y_continuous("GO Terms (Biological Process)", breaks = df.bp$y_var, labels = df.bp$Term) +
  scale_x_discrete("") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey"), text=element_text(size=12))

go.bp

# export the plot to a png file
png(filename = "DE/GO/Supp_fig_5A-bubble_plot_selected_BP.png",
    width = 600 * 7.55,        
    height = 600 * 6.8,
    units = "px",
    res = 600,
    pointsize = 8)  
plot(go.bp)
dev.off()

# Graph general shade CC
df.cc<-general_shade_GO %>% select(everything()) %>% filter(GO_Id %in% c(sig_list_general_CC))
df.cc$y_var <- c(1:(length(df.cc$GO_Id)))

# assign a representation factor of 0.01 to those terms not present in the data
df.cc[which(is.na(df.cc$p.value)),"rf"] <- 0.01 
# assign NA to non-significant terms
df.cc[which(df.cc$p.value > 0.05),"p.value"] <- NA

# obtain the p.values and assign a score
score <- df.cc$p.value
score[is.na(df.cc$p.value)]<-NA
score[df.cc$p.value < 0.05]<-1
score[df.cc$p.value <= 0.01]<-2
score[df.cc$p.value <= 0.001]<-3
score[df.cc$p.value <= 0.0001]<-4

go.cc <-ggplot(df.cc, aes(x = x_var, y = y_var, labels = Term, size =RF, col = score)) +
  geom_point(alpha = 0.7) +
  scale_colour_gradient(limits= c(1, 4), low = "#ffc000",  high = "red", na.value = "light gray") +     # pink "#f6e7e4"
  scale_size_continuous(range = c(2,20)) +
  scale_y_continuous("GO Terms (Cellular Component)", breaks = df.cc$y_var, labels = df.cc$Term) +
  scale_x_discrete("") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey"), text=element_text(size=12))

go.cc

# export the plot to a png file
png(filename = "DE/GO/Supp_fig_5B-bubble_plot_selected_CC.png",
    width = 600 * 5.55,        
    height = 600 * 6.8,
    units = "px",
    res = 600,
    pointsize = 8)    
plot(go.cc)
dev.off()

#######################
# End of bubble plots #
#######################