# Script to generate dendrogram on 1000 most variable coding genes
# Cohort: SCANDARE TNBC

### load R packages
require(dplyr)
require(stringr)
require(DESeq2)
library(ComplexHeatmap)

################################################################################
### clinical data
################################################################################
splan.table <- read.table("clinical_data.tsv", sep='\t', header = TRUE, row.names = "NGS_ID")
splan <- splan.table %>%
  mutate(run = factor(str_replace(run, "-","\\."))) %>% 
  mutate(age.cat = factor(age.cat, levels=c(">=50","<50"), labels=c(">= 50 yo", "< 50 yo"))) %>% 
  mutate(INC_PROT_I = factor(INC_PROT_I)) %>% 
  mutate(INC_BR_K1T.cat = factor(INC_BR_K1T.cat)) %>% 
  mutate(INC_BR_K1N.cat = factor(INC_BR_K1N.cat)) %>% 
  mutate(DFS = factor(DFS, levels = c(1, 0), labels=c("Yes","No"))) %>% 
  mutate(RCB.status.cat  = factor(ifelse(is.na(RCB.status), NA, 
                                         ifelse(RCB.status == "0", "0", "I, II, III"))))

################################################################################
### Counts
################################################################################
d.coding <- readRDS("d.coding.RDS")
d.tpm.coding <- readRDS("d.tpm.coding.RDS")

d.coding <- d.coding[,rownames(splan)]
d.tpm.coding <- d.tpm.coding[,rownames(splan)]

################################################################################
### Normalization
################################################################################

dds <- DESeqDataSetFromMatrix(countData=d.coding, DataFrame(condition=splan$run), ~ condition)
dds <- estimateSizeFactors(dds)         ## Estimate size factors
dds <- dds[rowSums(counts(dds)) > 0, ]  ## Remove lines with (only zeros) 
rld <- vst(dds, blind=TRUE)             ## Run the vst normalization
d.rlog <- assay(rld)                    ## extract matrix
gvar <- apply(d.rlog, 1, var) 

################################################################################
### Figure: dendrogram
################################################################################
mostvargenes <- order(gvar, decreasing=TRUE)[1:1000]    ## Select 1000 most variable genes

heatmap.annot <- rowAnnotation(Age = splan$age.cat,
                               Subcohort = splan$INC_PROT_I,
                               `T stage` = splan$INC_BR_K1T.cat,
                               `N stage` = splan$INC_BR_K1N.cat,
                               `RCB status` = splan$RCB.status.cat,
                               Event = splan$DFS,
                               
                               annotation_legend_param = list(Age = list(title = "Age", 
                                                                         at = c("≥ 50 yo","< 50 yo"), 
                                                                         labels = expression("" >= 50 ~ years ~ old, "" < 50 ~ years ~ old)),
                                                              title_gp = gpar(fontsize = 10, fontface = "bold"),   # Title font size
                                                              labels_gp = gpar(fontsize = 9)),   # Label font size),
                               
                               annotation_name_gp = gpar(fontsize = 10),
                               annotation_name_side = "top",
                               annotation_name_rot = 45,
                               
                               col = list(Age = c("< 50 yo" = "#EB949D", "≥ 50 yo" = "#FBD1D6"),
                                          Subcohort = c("Neoadjuvant therapy" = "#C3E6DC", "Upfront surgery" = "#62AA9D"),
                                          `N stage` = c("N0" = "thistle", "N+" = "orchid4","Unknown" = "grey75"),
                                          Event = c("Yes" = "#7CB5FF", "No" = "#F0F8FF","Unknown" = "grey75"),
                                          `RCB status` = c("0" = "#cbe1e1", "Unknown" = "grey75", "I, II, III"= "#008b8b"),                                   
                                          `T stage` = c("T1, T2" = "lightsteelblue2", "T3, T4"= "dodgerblue3","Unknown" = "grey75")),
                               
                               show_legend = TRUE)

scaled_data <- t(scale(t(d.rlog[mostvargenes,])))  # row scaling as in pheatmap(scale = "row")

hc <- ComplexHeatmap::Heatmap(t(scaled_data),
                              cluster_columns = FALSE,
                              clustering_distance_rows = "pearson",
                              clustering_method_rows = "ward.D2", 
                              col = circlize::colorRamp2(c(0, 1), c("white", "white")), #color_fun,
                              show_row_names = FALSE,
                              show_column_names = FALSE,
                              column_names_gp = gpar(fontsize = 7),
                              column_names_rot = 45,
                              left_annotation = heatmap.annot,
                              width = unit(1, "cm"),
                              row_dend_width = unit(3, "cm"),
                              show_heatmap_legend = FALSE)


# pdf("SCAN-TNBC_dendrogram.pdf", height=8)
draw(hc, annotation_legend_side = "bottom")
# dev.off()