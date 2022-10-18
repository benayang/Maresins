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
library(rstatix)
library(ComplexHeatmap)

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
sample_data <- sample_data[,-1]
assay <- data.frame(lapply(sample_data,as.numeric))
rownames(assay) <- metabolites
assay[1:5,1:5]

meta_data <- read_excel(file.path(projdir, "TA VML injury_LC-MS_08-13-21 extract JAL.xlsx"), sheet = 1, skip = 3) %>% as.data.frame()
meta_data$`Code/mouse#` <- as.numeric(meta_data$`Code/mouse#`)
rownames(meta_data) <- meta_data$`LC-MS/MS code`
meta_data$Injury <- factor(meta_data$Injury, levels = c("Uninjured","1 mm","2 mm"))
meta_data$Time <- factor(meta_data$Time, levels=stringr::str_sort(unique(meta_data$Time), numeric=T))
meta_data$Leg <- factor(recode(meta_data$Leg, right = "Right", rght = "Right"))
meta_data <- mutate(meta_data, Condition = factor(paste(Time, Injury, sep="_"), ))

rowData <- rownames(assay) %>% as.data.frame()
colnames(rowData) <- "BIOCHEMICAL"
rownames(rowData) <- rowData$BIOCHEMICAL
rowData$SUPER_PATHWAY <- rep("Lipid",length(rowData))
met <- create_mae(assay,rowData,meta_data)

# Impute missing data
png(file.path(projdir,"Plots","na_heatmap.png"), res=300, units='in', width=6, height=6)
na_heatmap(met, group_factor="Condition", label_colors=col_vector[1:7])
dev.off()

met <- knn_impute(met,cutoff=0.4) # use knn to impute missing values as long as metabolite is expressed in >60% of samples
# [1] raw: SummarizedExperiment with 143 rows and 40 columns
# [2] imputed: SummarizedExperiment with 79 rows and 40 columns
# 64 lipids removed since present in <40% of samples

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
data.imputed <- data$imputed
data.df <- data.imputed %>% 
              as.data.frame() %>%
              tibble::rownames_to_column("Lipid") %>%
              pivot_longer(cols=!Lipid, names_to="Sample", values_to="Concentration") %>%
              left_join(meta_data %>% tibble::rownames_to_column("Sample"), by="Sample") %>%
              mutate_if(is.character, as.factor)

table(data.df$Injury, data.df$Time)

# source("http://peterhaschke.com/Code/multiplot.R")
# pca_plot(met, group_factor="Condition", label_colors=brewer.pal(name="Paired", n=12))
# ggsave(file.path(projdir, "Plots", "condition_pca.png"), dpi=300, width=6, height=5)

# png(file.path(projdir, "Plots", "condition_and_sample_pca.png"), res=300, units="in", width=6, height=10)
# multiplot(pca_plot(met, group_factor="Condition", label_colors=brewer.pal(name="Paired", n=12)),
#           pca_plot(met, group_factor="LC-MS/MS code", 
#                   label_colors=col_vector[seq_along(as.data.frame(colData(met))$Condition)]))
# dev.off()

#####################################################################################
# PCA on all data without removing outliers
#####################################################################################

all_pca <- prcomp(t(data.imputed), scale=T)
summary(all_pca)
var_explained <- all_pca$sdev^2/sum(all_pca$sdev^2)
png(file.path(projdir,"Plots","raw_imputed_PCA_var_prop.png"), res=300, units='in', width=6, height=4)
barplot(var_explained[1:15], names.arg=1:15, xlab="PCs", ylab="Proportion of Variance")
dev.off()

pca_plt_data <- all_pca$x %>% 
  as.data.frame %>%
  rownames_to_column("Sample") %>%
  mutate(i = stringr::str_order(Sample, numeric = T)) %>%
  arrange(i) %>%
  mutate(Sample = factor(Sample, levels=Sample)) %>%
  left_join(meta_data %>% select(Condition, Leg, `LC-MS/MS code`), by=c("Sample"="LC-MS/MS code")) %>%
  group_by(Condition, Leg) %>%
  mutate(Rep = factor(seq_len(n()))) %>%
  ungroup() %>%
  mutate(Condition = factor(Condition, levels=unique(Condition)))

png(file.path(projdir,"Plots","raw_imputed_PCA_PC12.png"), res=300, units='in', width=6, height=6)
ggplot(pca_plt_data, aes(x=PC1,y=PC2)) + 
  geom_point(aes(color=Condition, shape=Rep), size=4) +
  #scale_color_manual(values = col_vector) +
  ggsci::scale_color_jco() +
  theme_bw() + coord_equal() +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(text = element_text(size=12, color="black"))
dev.off()

png(file.path(projdir,"Plots","raw_imputed_PCA_PC23.png"), res=300, units='in', width=6, height=6)
ggplot(pca_plt_data, aes(x=PC2,y=PC3)) + 
  geom_point(aes(color=Condition, shape=Rep), size=4) +
  #scale_color_manual(values = col_vector) +
  ggsci::scale_color_jco() +
  theme_bw() + coord_equal() +
  labs(x=paste0("PC2: ",round(var_explained[2]*100,1),"%"),
       y=paste0("PC3: ",round(var_explained[3]*100,1),"%")) +
  theme(text = element_text(size=12, color="black"))
dev.off()


#####################################################################################
# Heatmap for all lipids and samples
#####################################################################################

#to.plot <- to.plot[,c(1:15,20:23,29:33,16:19,24:28,34:38)] # reorder columns, if removing samples
hmp_sample_order <- colData(met) %>%
  as.data.frame() %>%
  mutate(Injury = factor(Injury, levels=c("Uninjured","1 mm", "2 mm")), 
        Time = factor(Time, levels=stringr::str_sort(unique(Time), numeric=T))) %>% 
  arrange(Injury, Time) %>%
  rownames()

ann_colors <- list(
  Time = c("grey60", col_vector[1:3]),
  Injury = col_vector[5:7],
  Leg = col_vector[8:9]
)
names(ann_colors$Time) <- unique(as.data.frame(colData(met))$Time)
names(ann_colors$Injury) <- unique(as.data.frame(colData(met))$Injury)
names(ann_colors$Leg) <- unique(as.data.frame(colData(met))$Leg)

mat_breaks <- seq(-2, 2, length.out = 100)

#all.plot <- data.imputed[,c(1:15,21:25,31:35,16:20,26:30,36:40)] # reorder columns
all.plot <- data.imputed[,hmp_sample_order]
png(file.path(projdir,"Plots","Heatmap_all_Lipids.png"), height = 12, width = 10, res = 300, units = 'in')
pheatmap(all.plot,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(250),
  annotation_col = as.data.frame(colData(met))[,c("Time","Injury","Leg")],
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
  main = "All Lipids",
  fontsize = 10,
  angle_col = "90",
  gaps_col = c(10,25)
)
dev.off()

#####################################################################################
# Distance matrix
#####################################################################################

sampleDists <- dist(t(assays(met)$imputed))
sampleDistMatrix <- as.matrix(sampleDists)
png(file.path(projdir,"Plots","Heatmap_all_Lipids_distance.png"), height = 12, width = 12, res = 300, units = 'in')
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(250),
          annotation_col = as.data.frame(colData(met))[,c("Time","Injury","Leg")],
          annotation_colors = ann_colors,
          annotation_legend = TRUE,
          border_color = NA,
          show_rownames = T,
          show_colnames = T,
          main = "All Lipids Distance Matrix",
          fontsize = 10,
          angle_col = "90")
dev.off()

#####################################################################################
# Drop outliers
#####################################################################################

outlier_samples <- c("JACM.1","JACM.19","JACM.25","JACM.28","JACM.38")

filt.met <- create_mae(assay = assay[,setdiff(colnames(assay), outlier_samples)],
                      rowData = rowData,
                      colData = meta_data[setdiff(rownames(meta_data), outlier_samples),])
filt.met <- knn_impute(filt.met,cutoff=0.4) # use knn to impute missing values as long as metabolite is expressed in >60% of samples

filt.data <- assays(filt.met)
filt.data.imputed <- filt.data$imputed
filt.data.df <- filt.data.imputed %>% 
              as.data.frame() %>%
              tibble::rownames_to_column("Lipid") %>%
              pivot_longer(cols=!Lipid, names_to="Sample", values_to="Concentration") %>%
              left_join(meta_data %>% tibble::rownames_to_column("Sample"), by="Sample") %>%
              mutate_if(is.character, as.factor)

#####################################################################################
# Check PCA after cleaning data
#####################################################################################

filt_pca <- prcomp(t(filt.data.imputed), scale=T)
summary(filt_pca)
var_explained <- filt_pca$sdev^2/sum(filt_pca$sdev^2)
png(file.path(projdir,"Plots","filt_imputed_PCA_var_prop.png"), res=300, units='in', width=6, height=4)
barplot(var_explained[1:15], names.arg=1:15, xlab="PCs", ylab="Proportion of Variance")
dev.off()

filt_pca_plt_data <- filt_pca$x %>% 
  as.data.frame %>%
  rownames_to_column("Sample") %>%
  mutate(i = stringr::str_order(Sample, numeric = T)) %>%
  arrange(i) %>%
  mutate(Sample = factor(Sample, levels=Sample)) %>%
  left_join(meta_data %>% select(Condition, Leg, `LC-MS/MS code`), by=c("Sample"="LC-MS/MS code")) %>%
  group_by(Condition, Leg) %>%
  mutate(Rep = factor(seq_len(n()))) %>%
  ungroup() %>%
  mutate(Condition = factor(Condition, levels=unique(Condition)))

png(file.path(projdir,"Plots","filt_imputed_PCA_PC12.png"), res=300, units='in', width=6, height=6)
ggplot(filt_pca_plt_data, aes(x=PC1,y=PC2)) + 
  geom_point(aes(color=Condition, shape=Rep), size=4) +
  #scale_color_manual(values = col_vector) +
  ggsci::scale_color_jco() +
  theme_bw() + coord_equal() +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(text = element_text(size=12, color="black"))
dev.off()

png(file.path(projdir,"Plots","filt_imputed_PCA_PC23.png"), res=300, units='in', width=6, height=6)
ggplot(filt_pca_plt_data, aes(x=PC2,y=PC3)) + 
  geom_point(aes(color=Condition, shape=Rep), size=4) +
  #scale_color_manual(values = col_vector) +
  ggsci::scale_color_jco() +
  theme_bw() + coord_equal() +
  labs(x=paste0("PC2: ",round(var_explained[2]*100,1),"%"),
       y=paste0("PC3: ",round(var_explained[3]*100,1),"%")) +
  theme(text = element_text(size=12, color="black"))
dev.off()

#####################################################################################
# Heatmap for all lipids and samples after removing outliers
#####################################################################################

gaps_col_idx <- colData(filt.met) %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  group_by(Injury) %>% 
  summarise(first_value = first(Sample),
            last_value = last(Sample)) %>% 
  mutate(last_idx = na.omit(match(last_value, rownames(as.data.frame(colData(filt.met))))),
        first_idx = na.omit(match(first_value, rownames(as.data.frame(colData(filt.met)))))) 

hmp_sample_order <- colData(filt.met) %>%
  as.data.frame() %>%
  mutate(Injury = factor(Injury, levels=c("Uninjured","1 mm", "2 mm")), 
        Time = factor(Time, levels=stringr::str_sort(unique(Time), numeric=T))) %>% 
  arrange(Injury, Time) %>%
  rownames()

filt.plot <- filt.data.imputed[,hmp_sample_order]
png(file.path(projdir,"Plots","Heatmap_all_Lipids_filtered.png"), height = 12, width = 10, res = 300, units = 'in')
pheatmap(filt.plot,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(250),
  annotation_col = as.data.frame(colData(filt.met))[,c("Time","Injury","Leg")],
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
  main = "All Lipids (outliers removed)",
  fontsize = 10,
  angle_col = "90",
  gaps_col = c(9,23)
)
dev.off()

#####################################################################################
# Normalize filtered object
#####################################################################################

filt.met <- normalize_met(filt.met) # variance stabilizing normalization
filt.met

# saveRDS(filt.met, file.path(projdir,"lipids_met_filt_obj.RDS"))

filt.data.df <- assays(filt.met)$norm_imputed %>% 
              as.data.frame() %>%
              tibble::rownames_to_column("Lipid") %>%
              pivot_longer(cols=!Lipid, names_to="Sample", values_to="Concentration") %>%
              left_join(meta_data %>% tibble::rownames_to_column("Sample"), by="Sample") %>%
              mutate_if(is.character, as.factor)

# QC Plot
png(file.path(projdir,"Plots","filtered_Quality_Plot.png"), res = 300, height = 6, width = 8, units = "in")
quality_plot(filt.met, group_factor="Condition",label_colors=col_vector[1:8])
dev.off()

# filt.data.df.test <- filt.data.df %>% mutate(log10conc = log10(Concentration)) %>% group_by(Time) %>% t_test(log10conc ~ Leg)
# filt.data.df.test <- add_xy_position(filt.data.df.test, x="Time", fun="mean_se")

# filt_conc_plt <- filt.data.df %>%
#     ggplot(aes(x=Time, y=log10(Concentration))) +
#     stat_summary(aes(color = Leg), fun.data="mean_se", geom="pointrange") +
#     stat_summary(aes(color = Leg), fun.data="mean_se", geom="point") + 
#     stat_summary(aes(color = Leg, group = Leg), fun="mean", geom="line") +
#     scale_y_continuous(expand=expansion(mult=c(0.1,0.1))) +
#     stat_pvalue_manual(filt.data.df.test, remove.bracket = T, label.size=3) +
#     labs(y = "log10(Imputed Raw Abundance)") +
#     theme_bw() + theme(axis.text = element_text(family = "Arial", color = "black"),
#                       axis.title = element_text(family = "Arial", color = "black", face="bold"))
# ggsave(plot=filt_conc_plt, filename=file.path(projdir,"Plots","filt_log10concentration.png"), dpi = 300, height = 3, width = 5)

#####################################################################################
# Heatmap of imputed and normalized data after outlier removal
#####################################################################################

norm.filt.plot <- assays(filt.met)$norm_imputed[,hmp_sample_order]
png(file.path(projdir,"Plots","Heatmap_all_Lipids_filtered_norm.png"), height = 12, width = 10, res = 300, units = 'in')
pheatmap(norm.filt.plot,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(250),
  annotation_col = as.data.frame(colData(filt.met))[hmp_sample_order,c("Time","Injury","Leg")],
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
  main = "All Lipids (outliers removed)",
  fontsize = 10,
  angle_col = "90",
  gaps_col = c(9,23)
)
dev.off()

#####################################################################################
# PCA of normalized and filtered object
#####################################################################################

norm_filt_pca <- prcomp(t(assays(filt.met)$norm_imputed), scale=T)
summary(norm_filt_pca)
var_explained <- norm_filt_pca$sdev^2/sum(norm_filt_pca$sdev^2)
png(file.path(projdir,"Plots","filt_norm_imputed_PCA_var_prop.png"), res=300, units='in', width=6, height=4)
barplot(var_explained[1:15], names.arg=1:15, xlab="PCs", ylab="Proportion of Variance")
dev.off()

filt_norm_pca_plt_data <- norm_filt_pca$x %>% 
  as.data.frame %>%
  rownames_to_column("Sample") %>%
  mutate(i = stringr::str_order(Sample, numeric = T)) %>%
  arrange(i) %>%
  mutate(Sample = factor(Sample, levels=Sample)) %>%
  left_join(meta_data %>% select(Condition, Leg, `LC-MS/MS code`), by=c("Sample"="LC-MS/MS code")) %>%
  group_by(Condition, Leg) %>%
  mutate(Rep = factor(seq_len(n()))) %>%
  ungroup() %>%
  mutate(Condition = factor(Condition, levels=unique(Condition)))

png(file.path(projdir,"Plots","filt_norm_imputed_PCA_PC12.png"), res=300, units='in', width=6, height=6)
ggplot(filt_norm_pca_plt_data, aes(x=PC1,y=PC2)) + 
  geom_point(aes(color=Condition, shape=Rep), size=4) +
  #scale_color_manual(values = col_vector) +
  ggsci::scale_color_jco() +
  theme_bw() + coord_equal() +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(text = element_text(size=12, color="black"))
dev.off()

#####################################################################################
# Find DE lipids
#####################################################################################

# https://towardsdatascience.com/anovas-three-types-of-estimating-sums-of-squares-don-t-make-the-wrong-choice-91107c77a27a
# aov and lm use Type 1 Sums of Squares, so the order of independent variables matters! Here, we test for effects of Injury followed by those of Time (Injury is the more important variable). We also don't consider interactions here.
data.stats <- data.frame()
for (lipid in levels(filt.data.df$Lipid)) {
  # MODIFIED TO EXCLUDE UNINJURED FROM DE LIPID CALCULATIONS
  no_uninjured_data <- filt.data.df %>% filter(Time != "Day 0" & Injury != "Uninjured") %>% droplevels()
  res.aov2 <- aov(Concentration ~ 0 + Condition, data = no_uninjured_data[no_uninjured_data$Lipid == lipid,])
  de <- data.frame("Lipid" = lipid, 
                  "pval" = summary(res.aov2)[[1]]["Condition", "Pr(>F)"])
  de <- cbind(de, t(res.aov2$coefficients))

  # res.aov2 <- aov(Concentration ~ Condition, data = no_uninjured_data[no_uninjured_data$Lipid == lipid,])
  # de <- data.frame("Lipid" = lipid, 
  #                 "pval.time" = summary(res.aov2)[[1]]["Time", "Pr(>F)"], 
  #                 "pval.injury" = summary(res.aov2)[[1]]["Injury", "Pr(>F)"])
  # de <- cbind(de, t(res.aov2$coefficients))

  # full_model <- lm(Concentration ~ Injury + Time + Injury:Time, data = filt.data.df[filt.data.df$Lipid == lipid,])
  # reduced_model <- lm(Concentration ~ Injury + Time, data = filt.data.df[filt.data.df$Lipid == lipid,])
  # lrt_test <- lmtest::lrtest(full_model, reduced_model)
  # de <- data.frame(Lipid = lipid, pval = na.omit(lrt_test$`Pr(>Chisq)`))

  data.stats <- rbind(data.stats, de)
}
data.stats$padj <- p.adjust(data.stats$pval, "BH")
# data.stats$padj.time <- p.adjust(data.stats$pval.time, "BH")
# data.stats$padj.injury <- p.adjust(data.stats$pval.injury, "BH")

# filt.met <- diff_test(filt.met, group_factors = c("Injury","Time"))
#  
# volcano_plot(filt.met,
#              group_factor="Time",
#              label_colors = c("darkseagreen","orange"),
#              main="Time comparisons")

#####################################################################################
# Find DE lipids using limma-voom
#####################################################################################

design_data <- as.data.frame(colData(filt.met)) %>% 
  filter(Time != "Day 0" & Injury != "Uninjured") %>% 
  mutate(Time = factor(gsub(" ","",Time)),
    Injury = factor(gsub(" ","",Injury)),
    Condition = factor(paste(Time, Injury, sep="_"))) %>%
  droplevels() 
design_matrix <- model.matrix(~ 0 + Condition, data = design_data)
no_uninjured_data <- assays(filt.met)$imputed
no_uninjured_data <- no_uninjured_data[,colData(filt.met)$Time != "Day 0"]

fit <- lmFit(no_uninjured_data, design = design_matrix)
colnames(fit$coefficients)
# [1] "ConditionDay14_1mm" "ConditionDay14_2mm" "ConditionDay3_1mm"
# [4] "ConditionDay3_2mm"  "ConditionDay7_1mm"  "ConditionDay7_2mm"
cont.matrix=makeContrasts(
    d14_2mm_vs_1mm="ConditionDay14_2mm - ConditionDay14_1mm",
    d7_2mm_vs_1mm="ConditionDay7_2mm - ConditionDay7_1mm",
    d3_2mm_vs_1mm="ConditionDay3_2mm - ConditionDay3_1mm",
    levels=design_matrix)
cont.fit <- contrasts.fit(fit, cont.matrix)
y <- voom(no_uninjured_data, design_matrix, plot = T)
fit_ebayes <- eBayes(cont.fit)
top.table <- topTable(fit_ebayes, sort.by = "F", p.value = 0.05, n = Inf)


#####################################################################################
# Heatmap for all DE lipids
#####################################################################################

sig_metabolites <- data.stats %>% filter(padj < 0.05) %>% pull(Lipid)
#sig_metabolites <- data.stats %>% filter(padj.injury < 0.05) %>% pull(Lipid)

#to.plot <- to.plot[,c(1:15,20:23,29:33,16:19,24:28,34:38)] # reorder columns, if removing samples
#all.de.plot <- data.imputed[sig_metabolites,c(1:15,21:25,31:35,16:20,26:30,36:40)] # reorder columns
all.de.plot <- assays(filt.met)$norm_imputed[sig_metabolites, hmp_sample_order]
png(file.path(projdir,"Plots","Heatmap_all_DELipids_filtered.png"), height = 12, width = 10, res = 300, units = 'in')
pheatmap(all.de.plot,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(250),
  annotation_col = as.data.frame(colData(filt.met))[hmp_sample_order, c("Time","Injury")],
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
  gaps_col = c(9,23)
)
dev.off()


#####################################################################################

# Export DE and all data for DPGP
#####################################################################################

export_data_for_dpgp <- function(data, meta_data, fc_type, prefix) {
  # average by condition
  data_avg <- data.frame(d0 = rowMeans(data[,meta_data %>% filter(Condition == "Day 0_Uninjured") %>% rownames()]),
                          d3_1mm = rowMeans(data[,meta_data %>% filter(Condition == "Day 3_1 mm") %>% rownames()]),
                          d3_2mm = rowMeans(data[,meta_data %>% filter(Condition == "Day 3_2 mm") %>% rownames()]),
                          d7_1mm = rowMeans(data[,meta_data %>% filter(Condition == "Day 7_1 mm") %>% rownames()]),
                          d7_2mm = rowMeans(data[,meta_data %>% filter(Condition == "Day 7_2 mm") %>% rownames()]),
                          d14_1mm = rowMeans(data[,meta_data %>% filter(Condition == "Day 14_1 mm") %>% rownames()]),
                          d14_2mm = rowMeans(data[,meta_data %>% filter(Condition == "Day 14_2 mm") %>% rownames()]))
  # calculate fold changes
  data_fc <- data.frame(d3_fc = data_avg$d3_2mm / data_avg$d3_1mm,
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

sig_data <- assays(filt.met)$norm_imputed[sig_metabolites, ]
#sample_subset <- meta_data %>% filter(Injury != "Uninjured") %>% rownames()
#sig_data <- data.norm[sig_metabolites, sample_subset]
sig_data <- export_data_for_dpgp(data = assays(filt.met)$norm_imputed[sig_metabolites, ],
                                meta_data = as.data.frame(colData(filt.met)),
                                prefix = "sig_data")

all_data <- export_data_for_dpgp(data = assays(filt.met)$norm_imputed,
                                meta_data = as.data.frame(colData(filt.met)),
                                prefix = "all_data")

# saveRDS(sig_data, file.path(projdir,"sig_data_for_dpgp.RDS"))
# saveRDS(all_data, file.path(projdir,"all_data_for_dpgp.RDS"))

#####################################################################################
# Make DPGP trajectory plots for all and DE lipid data
#####################################################################################

make_DPGP_plot <- function(data, dpgp_clusters, onecol) {
  dpgp_clusters_plot_df <- data[dpgp_clusters$gene, ] %>% 
    tibble::rownames_to_column("gene") %>% 
    pivot_longer(cols = !gene, names_to="day", values_to="zscore") %>%
    left_join(dpgp_clusters, by="gene") %>%
    group_by(cluster, day) %>%
    summarise(mean = mean(zscore), sd = sd(zscore)) %>%
    mutate(ymin = mean - 2*sd, ymax = mean + 2*sd)
  dpgp_clusters_plot_df$day <- factor(dpgp_clusters_plot_df$day, levels=c("d3_fc","d7_fc","d14_fc"))
  levels(dpgp_clusters_plot_df$day) <- c("d3_fc"="d3", "d7_fc"="d7", "d14_fc"="d14")
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

  if(onecol) {
    p <- p + facet_wrap(~cluster, labeller=as_labeller(panel_labs), ncol=1) +
        theme(legend.position = "none",
            text = element_text(family="Arial"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            strip.background = element_rect(color="black", fill="gray95"),
            strip.text = element_text(size = 12, face="bold"),
            axis.text = element_text(face="bold", color="black", size=12),
            axis.title = element_text(face="bold", color="black", size=12))
  } else {
    p <- p + facet_wrap(~cluster, labeller=as_labeller(panel_labs)) +
        theme(legend.position = "none",
            text = element_text(family="Arial"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            strip.background = element_rect(color="black", fill="gray95"),
            strip.text = element_text(size = 12, face="bold"),
            axis.text = element_text(face="bold", color="black", size=12),
            axis.title = element_text(face="bold", color="black", size=12))
  }
  return(p)
}

sig_dpgp_clusters <- read.table(file.path(projdir, "DPGP_within_timepoint", "sig_within_timepoint_optimal_clustering.txt"), sep='\t', header=T)
sig_DPGP_plt <- make_DPGP_plot(as.data.frame(sig_data$data_fc_zscore), sig_dpgp_clusters, F)
ggsave(plot=sig_DPGP_plt, filename = file.path(projdir,"Plots","DE_cluster_trajectory.png"), dpi=300, width=6, height=4)

all_dpgp_clusters <- read.table(file.path(projdir, "DPGP_within_timepoint", "all_within_timepoint_optimal_clustering.txt"), sep='\t', header=T)
all_DPGP_plt <- make_DPGP_plot(as.data.frame(all_data$data_fc_zscore), all_dpgp_clusters, F)
ggsave(plot=all_DPGP_plt, filename = file.path(projdir,"Plots","all_cluster_trajectory.png"), dpi=300, width=6, height=4)

sig_DPGP_onecol_plt <- make_DPGP_plot(as.data.frame(sig_data$data_fc_zscore), sig_dpgp_clusters, T)
ggsave(plot=sig_DPGP_onecol_plt, filename = file.path(projdir,"Plots","sig_cluster_trajectory_onecol.png"), dpi=300, width=2, height=12)

all_DPGP_onecol_plt <- make_DPGP_plot(as.data.frame(all_data$data_fc_zscore), all_dpgp_clusters, T)
ggsave(plot=all_DPGP_onecol_plt, filename = file.path(projdir,"Plots","all_cluster_trajectory_onecol.png"), dpi=300, width=2, height=12)

#####################################################################################
# Add DPGP clusters to heatmap
#####################################################################################

dpgp_hmp_sample_order <- colData(filt.met) %>%
  as.data.frame() %>%
  mutate(Injury = factor(Injury, levels=c("Uninjured","1 mm", "2 mm")), 
        Time = factor(Time, levels=stringr::str_sort(unique(Time), numeric=T))) %>% 
  arrange(Injury, Time) %>%
  rownames()

dpgp_gaps <- colData(filt.met) %>%
  as.data.frame() %>%
  mutate(Injury = factor(Injury, levels=c("Uninjured","1 mm", "2 mm")), 
        Time = factor(Time, levels=stringr::str_sort(unique(Time), numeric=T))) %>% 
  arrange(Injury, Time) %>%
  group_by(Injury) %>%
  tally() %>%
  mutate(gaps = cumsum(n))

dpgp_clust_ann <- data.frame(Cluster = factor(all_dpgp_clusters$cluster))
rownames(dpgp_clust_ann) <- all_dpgp_clusters$gene

dpgp_ann_colors <- list(
  Time = c("grey60", col_vector[1:3]),
  Injury = col_vector[5:7],
  Cluster = brewer.pal(11,"Paired")[seq_along(unique(all_dpgp_clusters$cluster))]
)
names(dpgp_ann_colors$Time) <- unique(as.data.frame(colData(filt.met))$Time)
names(dpgp_ann_colors$Injury) <- unique(as.data.frame(colData(filt.met))$Injury)
names(dpgp_ann_colors$Cluster) <- unique(all_dpgp_clusters$cluster)

png(file.path(projdir, "Plots", "Heatmap_all_Lipids_filtered_norm_DPGP.png"), height = 12, width = 10, res = 300, units = 'in')
pheatmap(assays(filt.met)$norm_imputed[all_dpgp_clusters$gene, dpgp_hmp_sample_order],
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(250),
  annotation_col = as.data.frame(colData(filt.met))[dpgp_hmp_sample_order, c("Time","Injury")],
  annotation_colors = dpgp_ann_colors,
  annotation_row = dpgp_clust_ann, 
  annotation_legend = TRUE,
  border_color = NA,
  #breaks = mat_breaks,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  legend = TRUE,
  show_rownames = T,
  show_colnames = T,
  main = "All Lipids",
  fontsize = 10,
  angle_col = "90",
  gaps_col = dpgp_gaps$gaps
)
dev.off()

#####################################################################################
# Add DPGP clusters to DE lipid heatmap
#####################################################################################

dpgp_clust_ann <- data.frame(Cluster = factor(sig_dpgp_clusters$cluster))
rownames(dpgp_clust_ann) <- sig_dpgp_clusters$gene

dpgp_ann_colors <- list(
  Time = c("grey60", col_vector[1:3]),
  Injury = col_vector[5:7],
  Cluster = brewer.pal(11,"Paired")[seq_along(unique(sig_dpgp_clusters$cluster))]
)
names(dpgp_ann_colors$Time) <- unique(as.data.frame(colData(filt.met))$Time)
names(dpgp_ann_colors$Injury) <- unique(as.data.frame(colData(filt.met))$Injury)
names(dpgp_ann_colors$Cluster) <- unique(sig_dpgp_clusters$cluster)

png(file.path(projdir, "Plots", "Heatmap_all_DELipids_filtered_norm_DPGP.png"), height = 12, width = 10, res = 300, units = 'in')
pheatmap(assays(filt.met)$norm_imputed[sig_dpgp_clusters$gene, dpgp_hmp_sample_order],
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(250),
  annotation_col = as.data.frame(colData(filt.met))[dpgp_hmp_sample_order, c("Time","Injury")],
  annotation_colors = dpgp_ann_colors,
  annotation_row = dpgp_clust_ann, 
  annotation_legend = TRUE,
  border_color = NA,
  #breaks = mat_breaks,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  legend = TRUE,
  show_rownames = T,
  show_colnames = T,
  main = "All DE Lipids",
  fontsize = 10,
  angle_col = "90",
  gaps_col = dpgp_gaps$gaps
)
dev.off()
