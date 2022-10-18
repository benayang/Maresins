# Metabolomics High-Level Data Analysis and Plotting
# Input: excel file containing LC-MS/MS analyte concentrations (pg per mg) and csv containing sample metadata.
# Output: Heatmap of top 25 DE, volcano plots comparing injuries at each timepoint, pca plots showing all timepoints for each injury
# Uses MetaboDiff for data imputation and QC, ggplot2, R stats package, and pheatmap for differential analysis and plotting
# Jacqueline Larouche
# August 29, 2021

#####################################################################################
# Set up environment and variables
#####################################################################################
library(readxl)
library(ggplot2)
library(MetaboDiff)
library(RColorBrewer)
library(reshape2)
library(ggpubr)
library(gplots)
library(ggrepel)

setwd("~/Documents/UMichigan/Aguilar_Lab/ResearchProjects/TA_VML/AICAR_Metabolomics")

qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
Day0 = "grey60"
Day3 = col_vector[1]
Day7 = col_vector[2]
Day14 = col_vector[3]

#####################################################################################
# Load and pre-process data
#####################################################################################
sample_data <- read_excel("20210818JCAMeicosanomics - pg per mg.xlsx", 
                          sheet = 'pg per mg', col_names = TRUE) %>% as.data.frame()
metabolites <- sample_data[,1]
colnames(sample_data) <- sample_data[1,]
sample_data <- sample_data[-1,-1]
assay <- data.frame(lapply(sample_data,as.numeric))
rownames(assay) <- metabolites[-1]
assay[1:5,1:5]

meta_data <- read_excel(file.path(getwd(), "TA VML injury_LC-MS_08-13-21 extract JAL.xlsx"), sheet = 1) %>% as.data.frame()
colnames(meta_data) <- meta_data[3,]
meta_data <- meta_data[-c(1:3),]
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

# QC Plot
png("Plots/Quality_Plot.png", res = 300, height = 6, width = 8, units = "in")
quality_plot(met, group_factor="Condition",label_colors=col_vector[1:8])
dev.off()

# Remove samples 25 and 19
png("Plots/Outlier_Heatmap.png", res = 300, height = 6, width = 8, units = "in")
outlier_heatmap(met, group_factor="Condition",
                label_colors=col_vector[1:7],
                k=7)
dev.off()

#####################################################################################
# Create Data Frame for statistical analysis and plotting
#####################################################################################
data <- assays(met)
data.norm <- data$norm_imputed
data.norm <- data.norm[,-c(19,25)]
data.raw <- data$imputed
data.raw <- data.raw[,-c(19,25)]
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
ggboxplot(data.df, x = "Time", y = "Concentration", color = "Injury",
          palette = c("#00AFBB", "#E7B800", "grey50"))


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
png(file = "Plots/Heatmap_Top25_DELipids.png", height = 7.5, width = 7, res = 300, units = 'in')
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

#############################################################################
# Volcano Plots comparing 2mm vs 1mm at each time point
##############################################################################
# Function for running t-test on each analyte and outputting as data frame
ttestRat <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}

## Day 3
data.3dpi <- data.raw[,11:19]
colnames(data.3dpi)
rawpvalue <- apply(data.3dpi, 1, ttestRat, grp1 = c(1:5), grp2 = c(6:9))
hist(rawpvalue)
data.3dpi <- log2(data.3dpi) #transform into log2 base.
control <- apply(data.3dpi[,1:5], 1, mean) #control group (1mm)
test <- apply(data.3dpi[, 6:9], 1, mean) #test group (2mm)
foldchange <- control - test
hist(foldchange, xlab = "log2 Fold Change (Control vs Test)")

results <- cbind(foldchange, rawpvalue)
results <- as.data.frame(results)
results$Lipid <- rownames(results)
results$diffexpressed <- "p>0.05"
results$diffexpressed[results$foldchange > 0.6 & results$rawpvalue < 0.05] <- "p<0.05 & >0.6 fold change"
results$diffexpressed[results$foldchange < -0.6 & results$rawpvalue < 0.05] <- "p<0.05 & <0.6 fold change"
results$delabel <- NA
results$delabel[results$diffexpressed != "p>0.05"] <- results$Lipid[results$diffexpressed != "p>0.05"]

# Transform the p-value (-1*log(p-value)) and create a volcano plot using ggplot2.
volcano <- ggplot(data = results, aes(y = -1*log10(rawpvalue), x = foldchange, col=diffexpressed, label = delabel))
png("Plots/Volcano_3days.png", res = 300, width = 10, height = 6, units = "in")
volcano + 
  geom_point(size = 3) +
  geom_text_repel() +
  theme_classic() + 
  geom_vline(xintercept=c(-0.6, 0.6), col="grey50", linetype='dotted') +
  geom_hline(yintercept=-log10(0.05), col="grey50", linetype='dotted') +
  scale_color_manual(values=c("blue", "red", "grey50"), name = "2mm vs 1mm") +
  labs(title = "Day 3 Post Injury", 
       x = "-Log10 raw p-value", 
       y = "Log2 Fold Change") +
  theme(plot.title = element_text(hjust = 0.5, size = 25)) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.title.x = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 15, colour = "black")) +
  theme(legend.title = element_text(size = 15, colour = "black"), 
        legend.text = element_text(size = 12, colour = "black"))
dev.off()
  
## Day 7
data.7dpi <- data.raw[,20:28]
colnames(data.7dpi)
rawpvalue <- apply(data.7dpi, 1, ttestRat, grp1 = c(1:4), grp2 = c(5:9))
hist(rawpvalue)
data.7dpi <- log2(data.7dpi)
control <- apply(data.7dpi[,1:5], 1, mean) #control group (1mm)
test <- apply(data.7dpi[, 6:10], 1, mean) #test group (2mm)
foldchange <- control - test
hist(foldchange, xlab = "log2 Fold Change (Control vs Test)")

results <- cbind(foldchange, rawpvalue)
results <- as.data.frame(results)
results$Lipid <- rownames(results)
results$diffexpressed <- "p>0.05"
results$diffexpressed[results$foldchange > 0.6 & results$rawpvalue < 0.05] <- "p<0.05 & >0.6 fold change"
results$diffexpressed[results$foldchange < -0.6 & results$rawpvalue < 0.05] <- "p<0.05 & <0.6 fold change"
results$delabel <- NA
results$delabel[results$diffexpressed != "p>0.05"] <- results$Lipid[results$diffexpressed != "p>0.05"]

# Transform the p-value (-1*log(p-value)) and create a volcano plot using ggplot2.
volcano <- ggplot(data = results, aes(y = -1*log10(rawpvalue), x = foldchange, col=diffexpressed, label = delabel))
png("Plots/Volcano_7days.png", res = 300, units = "in", width = 10, height = 6)
volcano + 
  geom_point(size = 3) +
  geom_text_repel() +
  theme_classic() + 
  xlim(-1.7,1.7) +
  ylim(-0.01,1.75) +
  geom_vline(xintercept=c(-0.6, 0.6), col="grey50", linetype='dotted') +
  geom_hline(yintercept=-log10(0.05), col="grey50", linetype='dotted') +
  scale_color_manual(values=c("blue", "grey50"), name = "2mm vs 1mm") +
  labs(title = "Day 7 Post Injury", 
       x = "-Log10 raw p-value", 
       y = "Log2 Fold Change") +
  theme(plot.title = element_text(hjust = 0.5, size = 25)) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.title.x = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 15, colour = "black")) +
  theme(legend.title = element_text(size = 15, colour = "black"), 
        legend.text = element_text(size = 12, colour = "black"))
dev.off()

## Day 14
data.14dpi <- data.raw[,29:38]
colnames(data.14dpi)
rawpvalue <- apply(data.14dpi, 1, ttestRat, grp1 = c(1:5), grp2 = c(6:10))
hist(rawpvalue)
data.14dpi <- log2(data.14dpi)
control <- apply(data.14dpi[,1:5], 1, mean) #control group (1mm)
test <- apply(data.14dpi[, 6:10], 1, mean) #test group (2mm)
foldchange <- control - test
hist(foldchange, xlab = "log2 Fold Change (Control vs Test)")

results <- cbind(foldchange, rawpvalue)
results <- as.data.frame(results)
results$Lipid <- rownames(results)
results$diffexpressed <- "p>0.05"
results$diffexpressed[results$foldchange > 0.6 & results$rawpvalue < 0.05] <- "p<0.05 & >0.6 fold change"
results$diffexpressed[results$foldchange < -0.6 & results$rawpvalue < 0.05] <- "p<0.05 & <0.6 fold change"
results$delabel <- NA
results$delabel[results$diffexpressed != "p>0.05"] <- results$Lipid[results$diffexpressed != "p>0.05"]

# Transform the p-value (-1*log(p-value)) and create a volcano plot using ggplot2.
volcano <- ggplot(data = results, aes(x = -1*log10(rawpvalue), y = foldchange, col=diffexpressed, label = delabel))
png("Plots/Volcano_14days.png", width = 900, height = 600)
volcano + 
  geom_point(size = 3) +
  geom_text_repel() +
  theme_classic() + 
  xlim(-2,2) +
  ylim(-0.01,3) +
  geom_vline(xintercept=c(-0.6, 0.6), col="grey50", linetype='dotted') +
  geom_hline(yintercept=-log10(0.05), col="grey50", linetype='dotted') +
  scale_color_manual(values=c("blue", "grey50"), name = "2mm vs 1mm") +
  labs(title = "Day 14 Post Injury", 
       x = "-Log10 raw p-value", 
       y = "Log2 Fold Change") +
  theme(plot.title = element_text(hjust = 0.5, size = 25)) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.title.x = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 15, colour = "black")) +
  theme(legend.title = element_text(size = 15, colour = "black"), 
        legend.text = element_text(size = 12, colour = "black"))
dev.off()

######################################################################################
# PCA plots (one for 2mm and one for 1mm, color by time-point)
######################################################################################
#pca_2mm <- prcomp(t(data.raw[,c(1:10, 16:19, 24:28, 34:38)]), center = FALSE, scale = FALSE)
pca_2mm <- prcomp(t(data.raw[,c(1:10, 16:20, 26:30, 36:40)]), center = FALSE, scale = FALSE)
pcaresults <- summary(pca_2mm)
scree.data <- as.data.frame(pcaresults$importance)
score.data <- as.data.frame(pcaresults$x)
loadings.data <- as.data.frame(pcaresults$rotation)
write.csv(score.data, file = "pca_score_2mm.csv")
pca_data <- read.csv("pca_score_2mm.csv", header = TRUE)
pca_data <- pca_data[,c(1:3)]
time <- meta_data$Time[(meta_data$Injury[-c(19,25)] %in% c("Uninjured", "2 mm"))]
time <- meta_data$Time[(meta_data$Injury %in% c("Uninjured", "2 mm"))]
pca_data$Group <- time[1:24]
pca_data$Group <- time[1:25]

png("Plots/PCA_score_2mm.png", units = "in", res = 300, width = 8, height = 6)
ggplot(pca_data, aes(PC1, PC2)) +
  geom_point(aes(color = Group, size = 5)) +
  stat_ellipse(aes(color = Group)) +
  scale_color_manual(values = c(Day0, Day14, Day3, Day7), name = "Days Post Injury") +
  theme_classic() +
  labs(title = "2mm PCA Score Plot", 
       x = "Principle Component 1", 
       y = "Principle Component 2") +
  theme(plot.title = element_text(hjust = 0.5, size = 25)) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.title.x = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 15, colour = "black"))
dev.off()


pca_1mm <- prcomp(t(data.raw[,c(1:15, 20:23, 29:33)]), center = FALSE, scale = FALSE)
pcaresults <- summary(pca_1mm)
scree.data <- as.data.frame(pcaresults$importance)
score.data <- as.data.frame(pcaresults$x)
loadings.data <- as.data.frame(pcaresults$rotation)
write.csv(score.data, file = "pca_score_1mm.csv")
pca_data <- read.csv("pca_score_1mm.csv", header = TRUE)
pca_data <- pca_data[,c(1:3)]
time <- meta_data$Time[(meta_data$Injury[-c(19,25)] %in% c("Uninjured", "1 mm"))]
pca_data$Group <- time[-c(16:17)]

png("Plots/PCA_score_1mm.png", units = "in", res = 300, width = 8, height = 6)
ggplot(pca_data, aes(PC1, PC2)) +
  geom_point(aes(color = Group, size = 5)) +
  stat_ellipse(aes(color = Group)) +
  scale_color_manual(values = c(Day0, Day14, Day3, Day7), name = "Days Post Injury") +
  theme_classic() +
  labs(title = "1mm PCA Score Plot", 
       x = "Principle Component 1", 
       y = "Principle Component 2") +
  theme(plot.title = element_text(hjust = 0.5, size = 25)) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.title.x = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 15, colour = "black"))
dev.off()
