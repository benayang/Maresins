#####################################################################################
# Set up environment and Var2s
#####################################################################################
library(readxl)
library(ggplot2)
library(MetaboDiff)
library(RColorBrewer)
library(reshape2)
library(ggpubr)
library(gplots)
library(ggrepel)

projdir <- "/nas/homes/benyang/Maresins"

qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
Day0 = "grey60"
Day3 = col_vector[1]
Day7 = col_vector[2]
Day14 = col_vector[3]

#####################################################################################
# Load and pre-process data
#####################################################################################
sample_data <- read_excel(file.path(projdir,"20210818JCAMeicosanomics - pg per mg.xlsx"), 
                          sheet = 'pg per mg', skip = 1, col_names = TRUE) %>% as.data.frame()
metabolites <- sample_data[,1]
#colnames(sample_data) <- sample_data[1,]
sample_data <- sample_data[,-1]
assay <- data.frame(lapply(sample_data,as.numeric))
rownames(assay) <- metabolites
assay[1:5,1:5]

meta_data <- read_excel(file.path(projdir, "TA VML injury_LC-MS_08-13-21 extract JAL.xlsx"), sheet = 1, skip = 3) %>% as.data.frame()
#colnames(meta_data) <- meta_data[3,]
#meta_data <- meta_data[-c(1:3),]
meta_data$`Code/mouse#` <- as.numeric(meta_data$`Code/mouse#`)
rownames(meta_data) <- meta_data$`LC-MS/MS code`
meta_data$Condition <- NA
for (i in 1:length(meta_data$Condition)) {
  meta_data$Condition[i] <- sprintf("%s_%s",meta_data$Time[i],meta_data$Injury[i])
}
head(meta_data)

rowData <- rownames(assay) %>% as.data.frame()
colnames(rowData) <- "BIOCHEMICAL"
rownames(rowData) <- rowData$BIOCHEMICAL
rowData$SUPER_PATHWAY <- rep("Lipid",length(rowData))
met <- create_mae(assay,rowData,meta_data)

# Impute missing data
na_heatmap(met, group_factor="Condition", label_colors=col_vector[1:7])
met <- knn_impute(met,cutoff=0.4) # use knn to impute missing values as long as metabolite is expressed in >60% of samples
met <- normalize_met(met) # variance stabilizing normalization
met

met <- met %>%
    diss_matrix %>%
    identify_modules(min_module_size=5) %>%
    name_modules(pathway_annotation="SUB_PATHWAY") %>%
    calculate_MS(group_factors=c("Time","Injury"))

MS_plot(met,
        group_factor="Condition",
        p_value_cutoff=0.05,
        p_adjust=FALSE)

# saveRDS(met, file.path(projdir,"lipids_met_obj.RDS"))

# QC Plot
png(file.path(projdir,"Plots","Quality_Plot.png"), res = 300, height = 6, width = 8, units = "in")
quality_plot(met, group_factor="Condition",label_colors=col_vector[1:8])
dev.off()

# Remove samples 25 and 19
png(file.path(projdir,"Plots","Outlier_Heatmap.png"), res = 300, height = 6, width = 8, units = "in")
outlier_heatmap(met, group_factor="Condition",
                label_colors=col_vector[1:7],
                k=7)
dev.off()


#####################################################################################
# Create Data Frame for statistical analysis and plotting
#####################################################################################
data <- assays(met)
data.norm <- data$norm_imputed
#data.norm <- data.norm[,-c(19,25)]
data.raw <- data$imputed
#data.raw <- data.raw[,-c(19,25)]
data.df <- melt(data.raw)
colnames(data.df) <- c('Lipid', 'Sample', 'Concentration')
#data.df$Time <- c(rep('Day 0',790), rep('Day 3', 79*9), rep('Day 7', 79*9), rep('Day 14', 790)) %>% as.factor() #if removing samples 19 and 25
data.df$Time <- c(rep('Day 0',790), rep('Day 3', 790), rep('Day 7', 790), rep('Day 14', 790)) %>% as.factor()
#data.df$Injury <- c(rep('Uninjured',790), rep('1mm', 79*5), rep('2mm',79*4), rep('1mm', 79*4), rep('2mm',79*5), rep('1mm', 79*5), rep('2mm',79*5)) %>% as.factor() #if removing samples 19, 25
data.df$Injury <- c(rep('Uninjured',790), rep('1mm', 79*5), rep('2mm',79*5), rep('1mm', 79*5), rep('2mm',79*5), rep('1mm', 79*5), rep('2mm',79*5)) %>% as.factor()
data.df$Condition <- NA
for (i in 1:length(data.df$Condition)) {
  data.df$Condition[i] <- sprintf("%s_%s",data.df$Time[i],data.df$Injury[i])
}
data.df$Condition <- as.factor(data.df$Condition)

table(data.df$Injury, data.df$Time)
data.df %>%
    mutate(log10conc = log10(Concentration)) %>%
    ggboxplot(x = "Time", y = "log10conc", color = "Injury", 
          palette = c("#00AFBB", "#E7B800", "grey50"))
ggsave(file.path(projdir,"Plots","log10concentration.png"), dpi = 300, height = 6, width = 8)

source("http://peterhaschke.com/Code/multiplot.R")
pca_plot(met, group_factor="Condition", label_colors=brewer.pal(name="Paired", n=7))
ggsave(file.path(projdir, "Plots", "condition_pca.png"), dpi=300, width=6, height=5)
# multiplot(
# pca_plot(met,
#          group_factor="Condition"),
# tsne_plot(met,
#           group_factor="Condition",
#              label_colors=c("darkseagreen","dodgerblue")),
# cols=2)

#####################################################################################
# Heatmap of top 25 DE analytes
#####################################################################################
data.stats <- data.frame()
for (lipid in levels(data.df$Lipid)) {
  res.aov2 <- aov(Concentration ~ Injury + Time, data = data.df[data.df$Lipid == lipid,])
  summary(res.aov2)
  de <- data.frame("Lipid" = lipid, "pval.time" = summary(res.aov2)[[1]][["Pr(>F)"]][2], "pval.injury" = summary(res.aov2)[[1]][["Pr(>F)"]][1])
  data.stats <- rbind(data.stats, de)
}

top.25 <- data.stats %>% top_n(n = -25, wt = pval.injury)
to.plot <- data.raw[top.25$Lipid,]
#to.plot <- to.plot[,c(1:15,20:23,29:33,16:19,24:28,34:38)] # reorder columns, if removing samples
to.plot <- to.plot[,c(1:15,21:25,31:35,16:20,26:30,36:40)] # reorder columns


# Time Series Heatmap of top 25-50 DEGs (z-score)
# annotation.df <- data.frame(Time = factor(c(rep('Day 0',10), 
#                                             rep('Day 3', 5), 
#                                             rep('Day 7', 4), 
#                                             rep('Day 14', 5), 
#                                             rep('Day 3', 4), 
#                                             rep('Day 7', 5), 
#                                             rep('Day 14', 5))),
#                             Injury = factor(c(rep('Uninjured',10), 
#                                               rep('1mm', 14), 
#                                               rep('2mm', 14))))
annotation.df <- data.frame(Time = factor(c(rep('Day 0',10), 
                                            rep('Day 3', 5), 
                                            rep('Day 7', 5), 
                                            rep('Day 14', 5), 
                                            rep('Day 3', 5), 
                                            rep('Day 7', 5), 
                                            rep('Day 14', 5))),
                            Injury = factor(c(rep('Uninjured',10), 
                                              rep('1mm', 15), 
                                              rep('2mm', 15))))
rownames(annotation.df) <- colnames(to.plot)
mat_breaks <- seq(-2, 2, length.out = 100)

# Specify colors
ann_colors = list(
  Time = c("Day 0" = "grey60", 
           "Day 3" = col_vector[1], 
           "Day 7" = col_vector[2], 
           "Day 14" = col_vector[3]),
  Injury = c("Uninjured" = col_vector[5], 
             "1mm" = col_vector[6], 
             "2mm" = col_vector[7])
)
png(file.path(projdir,"Plots","Heatmap_Top25_DELipids.png"), height = 7.5, width = 7, res = 300, units = 'in')
pheatmap(to.plot,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  annotation_col = annotation.df,
  annotation_colors = ann_colors,
  annotation_legend = TRUE,
  border_color = NA,
  #breaks = mat_breaks,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  legend = TRUE,
  show_rownames = T,
  show_colnames = T,
  main = "Top 25 DE Lipids",
  fontsize = 10,
  angle_col = "90",
  gaps_col = c(10,25)
)
dev.off()

#####################################################################################
# Heatmap for all DE lipids
#####################################################################################

#to.plot <- to.plot[,c(1:15,20:23,29:33,16:19,24:28,34:38)] # reorder columns, if removing samples
all.plot <- data.raw[,c(1:15,21:25,31:35,16:20,26:30,36:40)] # reorder columns
png(file.path(projdir,"Plots","Heatmap_all_DELipids.png"), height = 12, width = 10, res = 300, units = 'in')
pheatmap(all.plot,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(250),
  annotation_col = annotation.df,
  annotation_colors = ann_colors,
  annotation_legend = TRUE,
  border_color = NA,
  #breaks = mat_breaks,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  legend = TRUE,
  show_rownames = T,
  show_colnames = T,
  main = "All DE Lipids",
  fontsize = 10,
  angle_col = "90",
  gaps_col = c(10,25)
)
dev.off()

#####################################################################################
# Export DE and all data for DPGP
#####################################################################################

export_data_for_dpgp <- function(data, meta_data, prefix) {
  # average by condition
  data_avg <- data.frame(d0_left = rowMeans(data[,meta_data %>% filter(Time == "Day 0" & Leg == "Left") %>% rownames()]),
                            d0_right = rowMeans(data[,meta_data %>% filter(Time == "Day 0" & Leg == "Right") %>% rownames()]),
                            d3_1mm = rowMeans(data[,meta_data %>% filter(Condition == "Day 3_1 mm") %>% rownames()]),
                            d3_2mm = rowMeans(data[,meta_data %>% filter(Condition == "Day 3_2 mm") %>% rownames()]),
                            d7_1mm = rowMeans(data[,meta_data %>% filter(Condition == "Day 7_1 mm") %>% rownames()]),
                            d7_2mm = rowMeans(data[,meta_data %>% filter(Condition == "Day 7_2 mm") %>% rownames()]),
                            d14_1mm = rowMeans(data[,meta_data %>% filter(Condition == "Day 14_1 mm") %>% rownames()]),
                            d14_2mm = rowMeans(data[,meta_data %>% filter(Condition == "Day 14_2 mm") %>% rownames()]))
  # calculate fold changes
  data_fc <- data.frame(d0_fc = data_avg$d0_right / data_avg$d0_left,
                      d3_fc = data_avg$d3_2mm / data_avg$d3_1mm,
                      d7_fc = data_avg$d7_2mm / data_avg$d7_1mm,
                      d14_fc = data_avg$d14_2mm / data_avg$d14_1mm,
                      row.names = rownames(data_avg))
  # find z-scores of fold changes
  data_fc_zscore <- t(scale(t(data_fc), center=T, scale=T))
  # write tables
  write.table(data_avg, file.path(projdir, "Tables", sprintf("%s_avg.tsv",prefix)), sep='\t', row.names=T, col.names=T, quote=F)
  write.table(data_fc, file.path(projdir, "Tables", sprintf("%s_fc.tsv",prefix)), sep='\t', row.names=T, col.names=T, quote=F)
  write.table(data_fc_zscore, file.path(projdir, "Tables", sprintf("%s_fc_zscore.tsv",prefix)), sep='\t', row.names=T, col.names=T, quote=F)
  # return lists
  return(list(data_avg = data_avg, data_fc = data_fc, data_fc_zscore = data_fc_zscore))
}

sig_metabolites <- data.stats %>% filter(pval.injury < 0.05) %>% pull(Lipid)
sig_data <- data.norm[sig_metabolites, ]
#sample_subset <- meta_data %>% filter(Injury != "Uninjured") %>% rownames()
#sig_data <- data.norm[sig_metabolites, sample_subset]
sig_data <- export_data_for_dpgp(data = data.norm[sig_metabolites, ],
                                meta_data = as.data.frame(colData(met)),
                                prefix = "sig_data")

all_data <- export_data_for_dpgp(data = data.norm,
                                meta_data = as.data.frame(colData(met)),
                                prefix = "all_data")

saveRDS(sig_data, file.path(projdir,"sig_data_for_dpgp.RDS"))
saveRDS(all_data, file.path(projdir,"all_data_for_dpgp.RDS"))

#####################################################################################
# Mean and standard deviation of DE lipids
#####################################################################################

make_DPGP_plot <- function(data, dpgp_clusters) {
  dpgp_clusters_plot_df <- data[dpgp_clusters$gene, ] %>% 
    tibble::rownames_to_column("gene") %>% 
    pivot_longer(cols = !gene, names_to="day", values_to="zscore") %>%
    left_join(dpgp_clusters, by="gene") %>%
    group_by(cluster, day) %>%
    summarise(mean = mean(zscore), sd = sd(zscore)) %>%
    mutate(ymin = mean - 2*sd, ymax = mean + 2*sd)
  dpgp_clusters_plot_df$day <- factor(dpgp_clusters_plot_df$day, levels=c("d0_fc","d3_fc","d7_fc","d14_fc"))
  levels(dpgp_clusters_plot_df$day) <- c("d0_fc"="d0", "d3_fc"="d3", "d7_fc"="d7", "d14_fc"="d14")
  dpgp_clusters_plot_df$cluster <- factor(dpgp_clusters_plot_df$cluster)

  p = ggplot(data=dpgp_clusters_plot_df, aes(x=day)) + 
    geom_hline(yintercept=0, lty=2) +
    geom_ribbon(aes(ymin=ymin, ymax=ymax, group=cluster), fill="gray", alpha=0.6) +
    geom_line(aes(y=mean, group=cluster), size=1) +
    geom_point(aes(y=mean, color=cluster), size=5, stroke=0) + 
    scale_color_brewer(palette="Paired") +
    theme_bw() +
    labs(x = "Day", y = "Fold-change z-score")
  
  panel_labs <- paste0("Cluster ",unique(dpgp_clusters$cluster), ", N=", table(dpgp_clusters$cluster))
  names(panel_labs) <- unique(dpgp_clusters$cluster)

  p <- p + facet_wrap(~cluster, labeller=as_labeller(panel_labs)) +
      theme(legend.position = "none",
          text = element_text(family="Arial"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.background = element_rect(color="black", fill="gray95"),
          strip.text = element_text(size = 12, face="bold"),
          axis.text = element_text(face="bold", color="black", size=12),
          axis.title = element_text(face="bold", color="black", size=12))
  return(p)
}

sig_dpgp_clusters <- read.table(file.path(projdir, "DPGP_within_timepoint", "sig_within_timepoint_optimal_clustering.txt"), sep='\t', header=T)
sig_DPGP_plt <- make_DPGP_plot(as.data.frame(sig_data$data_fc_zscore), sig_dpgp_clusters)
ggsave(plot=sig_DPGP_plt, filename = file.path(projdir,"Plots","DE_cluster_trajectory.png"), dpi=300, width=6, height=6)

all_dpgp_clusters <- read.table(file.path(projdir, "DPGP_within_timepoint", "all_within_timepoint_optimal_clustering.txt"), sep='\t', header=T)
all_DPGP_plt <- make_DPGP_plot(as.data.frame(all_data$data_fc_zscore), all_dpgp_clusters)
ggsave(plot=all_DPGP_plt, filename = file.path(projdir,"Plots","all_cluster_trajectory.png"), dpi=300, width=6, height=6)
