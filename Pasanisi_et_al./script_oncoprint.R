# Oncoprint script
# R version 4.1.2 (2021-11-01)

require(tidyverse)
require(ggplot2)
library(ComplexHeatmap)

################################################################################
### Add clinical data
countsdatapath <-"path"     #Define path on the cluster
clinical.data <- read.table(file.path(countsdatapath,"clinical_data.csv"),sep=",", header = TRUE)
sample_order <- read.table(file.path(countsdatapath,"sample_order.txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
column_order <- sample_order$ID.SCANDARE

### Mutational signatures : major
sig <- read.table(file.path(countsdatapath,"SBS_major_per_sample.tsv"), header=TRUE)

clinical.data <- clinical.data %>%
  left_join(sig, by=c("ID"="PAIR_ID")) %>%
  mutate(RCB = sub('I/II/III', 'I, II, III', RCB),
         Tstage = sub('/', ', ', Tstage)) %>%
  mutate(Age = factor(Age, levels = c(">= 50 years old","< 50 years old")),
         Subcohort = factor(Subcohort, levels = c("Neoadjuvant therapy","Upfront surgery")),
         Tstage = factor(replace_na(as.character(Tstage), "Unknown"), levels = c("T1, T2", "T3, T4","Unknown")),
         Nstage = factor(replace_na(as.character(Nstage), "Unknown"), levels = c("N0","N+","Unknown")),
         RCB = factor(replace_na(as.character(RCB), "Unknown"), levels = c("0","I, II, III", "Unknown")),
         DFS = factor(replace_na(as.character(DFS), "Unknown"), levels = c("Yes","No","Unknown")),
         HRD = factor(replace_na(as.character(HRD), "Unknown"), levels = c("HRD","HRP","Unknown")),
         TMB = factor(TMB, levels=c("High","Low")),
         MSI = factor(MSI, levels=c("MSI","MSS")),
         Max_SBS = factor(Max_SBS)) %>%
  mutate(Statut.methylation.constit.BRCA1 = factor(case_when(Statut.methylation.constit.BRCA1 == "Oui" ~ "Yes",
                                                             Statut.methylation.constit.BRCA1 == "Non" ~ "No",
                                                             is.na(Statut.methylation.constit.BRCA1) ~ "Unknown"), levels =c("Yes","No","Unknown")),
         Statut.methylation.tumeur.BRCA1 = factor(case_when(Statut.methylation.tumeur.BRCA1 == "Oui" ~ "Yes",
                                                            Statut.methylation.tumeur.BRCA1 == "Non" ~ "No",
                                                            is.na(Statut.methylation.tumeur.BRCA1) ~ "Unknown"), levels =c("Yes","No","Unknown")),
         LOH = factor(case_when(LOH == "Oui" ~ "Yes",
                                LOH == "Non" ~ "No",
                                is.na(LOH) ~ "Unknown"), levels =c("Yes","No","Unknown"))) %>%
  column_to_rownames(var="ID")

################################################################################
### Load matrix
input.mat <- read.csv(file.path(countsdatapath,"HRR_gene_matrix.csv"),",",header=TRUE,row.names = 1)
sub_mat <- t(input.mat)

### Add missing patient
missing_patient <- setdiff(rownames(clinical.data),colnames(sub_mat))
empty_mat <- matrix(NA, nrow=nrow(sub_mat), ncol = length(missing_patient))
colnames(empty_mat) <- missing_patient

mat <- cbind(sub_mat, empty_mat)
mat <- as.matrix(mat)
mat[is.na(mat)] <- "" ## replace NAs by ""

################################################################################
### Colors

col = c("background"="#CCCCCC",
        "frameshift_variant" ="#7DB3E2FF",
        "splice_donor_variant"="#DB3D06FF",
        "missense_variant"="#8EB035FF",
        "stop_gained"="gold",
        "constit" = "#3D3E52FF")


################################################################################
### Figure parameters

heatmap_legend_param = list(title = "Alterations",
                            at = c("stop_gained", "frameshift_variant", "missense_variant", "splice_donor_variant", "constit"),
                            labels = c("Stop gain", "Frameshift indel", "Missense", "Splicing", "Constitutional"),
                            ncol = 1,
                            legend_gp = gpar(fontsize = 8))
alter_fun_list = list(
  background = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA)),
  stop_gained = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["stop_gained"], col = NA)),
  frameshift_variant = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.6, gp = gpar(fill = col["frameshift_variant"], col = NA)),
  missense_variant = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["missense_variant"], col = NA)),
  splice_donor_variant = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4,gp = gpar(fill = col["splice_donor_variant"], col = NA)),
  constit = function(x, y, w, h) grid.points(x, y, pch = 16, size = unit(2, "mm"), gp = gpar(col = col["constit"]))
)

################################################################################
### Oncoprint

heatmap.annot <- HeatmapAnnotation(Age = clinical.data[colnames(mat),"Age"],
                                   `cT` = clinical.data[colnames(mat),"Tstage"],
                                   `cN` = clinical.data[colnames(mat),"Nstage"],
                                   `RCB status` = clinical.data[colnames(mat),"RCB"],
                                   HRD = clinical.data[colnames(mat),"HRD"],
                                   `BRCA1 LOH` = clinical.data[colnames(mat),"LOH"],
                                   `Meth-cBRCA1` = clinical.data[colnames(mat),"Statut.methylation.constit.BRCA1"],
                                   `Meth-tBRCA1` = clinical.data[colnames(mat),"Statut.methylation.tumeur.BRCA1"],
                                   annotation_legend_param = list(
                                     by_nrow = FALSE,
                                     # legend_padding = unit(20, 'mm'),
                                     Age = list(title = "Age", at = c(">= 50 years old","< 50 years old"), labels = expression("" >= 50 ~ years ~ old, "" < 50 ~ years ~ old))),
                                   cbar = anno_oncoprint_barplot(height = unit (1.5, "cm")),
                                   annotation_name_gp = gpar(fontsize = 9),
                                   # simple_anno_size = unit(4, "cm"),
                                   col = list(
                                     TMB = c("High" = "#FFB800", "Low" = "#FFF9B1"),
                                     MSI = c("MSI" = "#D97A32", "MSS" = "#FEE8C8"),
                                     HRD=c("HRD"="red3" ,"HRP"="#FFF5EE","Unknown" = "grey75"),
                                     Age = c("< 50 years old" = "#EB949D", ">= 50 years old" = "#FBD1D6"),
                                     Subcohort = c("Neoadjuvant therapy" = "#C3E6DC", "Upfront surgery" = "#62AA9D"),
                                     `cN` = c("N0" = "thistle", "N+" = "orchid4","Unknown" = "grey75"),
                                     Event = c("Yes" = "#7CB5FF", "No" = "#F0F8FF","Unknown" = "grey75"),
                                     `Mutational signatures` = c(
                                       "Alkyl"="#8DD3C7",
                                       "APOBEC"="#FFED6F",
                                       "Aristol_ac"="#BEBADA",
                                       "Artifact"="#FB8072",
                                       "Aza"="#80B1D3",
                                       "Chemotherapy/Treatment"="#FDB462",
                                       "Colibactin"="#B3DE69",
                                       "HRD"="#FCCDE5",
                                       "Ig_hypermutation"="darkgrey",
                                       "MMR"="#BC80BD",
                                       "MUTYH"="#CCEBC5",
                                       "Platin"="#FFFFB3",
                                       "Unknown"="grey90"),
                                     `RCB status` = c("0" = "#cbe1e1", "Unknown" = "grey75", "I, II, III"= "#008b8b"),
                                     `cT` = c("T1, T2" = "lightsteelblue2", "T3, T4"= "dodgerblue3","Unknown" = "grey75"),
                                     `BRCA1 LOH` = c("Yes" = "red3", "No"= "#FFF5EE"),
                                     `Meth-cBRCA1` = c("Yes" = "forestgreen", "No"= "#C7E9C0"),
                                     `Meth-tBRCA1` = c("Yes" = "forestgreen", "No"= "#C7E9C0")),
                                   show_legend = TRUE)

col2 <- col
col2["constit"] <- NA

op1 <- oncoPrint(
  mat,
  alter_fun = alter_fun_list,
  col = col2,   
  heatmap_legend_param = heatmap_legend_param,
  column_order = column_order,
  top_annotation = heatmap.annot,
  remove_empty_columns = FALSE,
  show_column_names = FALSE,
  alter_fun_is_vectorized = TRUE
  )


pdf(file.path(countsdatapath,"oncoprint.pdf"), width=16, height=10)
draw(op1, heatmap_legend_side = "left", annotation_legend_side = "bottom")
dev.off()