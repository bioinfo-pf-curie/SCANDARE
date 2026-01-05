require(tidyverse)
require(ggplot2)
library(ComplexHeatmap)

################################################################################
### Load SNV
################################################################################

in.dir <- "table/" ## to be defined
files <- dir(path = Sys.glob(in.dir), pattern = "*_snv_oncodriver.tsv", recursive = TRUE, full.names=TRUE)
snv.data <- files %>%
  map_dfr(~ {
    sample <- gsub("_Mutect2.*","", basename(.))
    read.table(., header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE) %>% 
      mutate(across(.fns = as.character)) %>% 
      mutate(ID = sample) %>%
      rename("gene" = "ANN[0].GENE", "variant"="ANN[0].EFFECT") %>%
      mutate(pos = as.numeric(POS)) %>%
      mutate(end = pos + 1) %>%
      mutate(variant = sub("&.*","", variant)) %>%  ## keep only the first annotation in case of multiple
      rename(chrom = CHROM, start  = pos) %>%  relocate(chrom, start, end) %>% select(-c("POS"))
  })

snv.data <- snv.data %>% 
  mutate(variant = case_when(variant == "missense_variant" ~ "missense",
                             variant == "conservative_inframe_deletion" ~ "inframe indel",
                             variant == "conservative_inframe_insertion" ~ "inframe indel",
                             variant == "disruptive_inframe_deletion" ~ "inframe indel",
                             variant == "disruptive_inframe_insertion" ~ "inframe indel",
                             variant == "frameshift_variant" ~ "frameshift indel",
                             variant == "splice_acceptor_variant" ~ "splicing",
                             variant == "splice_donor_variant" ~ "splicing",
                             variant == "stop_gained" ~ "stop gain",
                             variant == "upstream_gene_variant" ~ "promoter"
  ))


################################################################################
### Load CNV
################################################################################

in.dir <- "table/" ## to be defined
files <- dir(path = Sys.glob(in.dir), pattern = "*_cnv_oncodriver.csv", recursive=TRUE, full.names=TRUE)

cnv.data <- files %>%
  map_dfr(~ {
    df <- read.table(., header = FALSE, skip = 1, sep = ",", 
                     check.names = FALSE, stringsAsFactors = FALSE, 
                     col.names = c("chrom","loc.start","loc.end","ID","CNt","Geno","logratio","ploidy","call","LOH","gene","driver_status","gene_type"))
    if (nrow(df) == 0) return(NULL)  # Skip empty files
    df %>% 
      rename("start" = "loc.start", "end" = "loc.end", "variant"="call") %>% 
      select(-c("CNt", "Geno", "logratio", "ploidy", "LOH")) %>%
      mutate(
        chrom = as.character(chrom),
        start = as.numeric(start),
        end = as.numeric(end),
        length = end - start
      )
  })

cnv.data <- cnv.data %>% mutate(isOncogene = ifelse(gene_type == "oncogene", "Yes","No"),
                                isTSG = ifelse(gene_type == "tsg", "Yes","No"),
                                isBoth = ifelse(gene_type == "both", "Yes","No"),
                                isUnknown = ifelse(gene_type == "unknown", "Yes","No"))


#########################################################################################
### Convert to matrix 
#########################################################################################

mat <- snv.data %>% 
  full_join(., y=cnv.data, by = c("gene", "ID", "chrom", "start", "end", "variant")) %>%
  select(gene, ID, variant) %>%
  group_by(gene, ID) %>%
  distinct() %>% 
  mutate(variant = paste0(variant, collapse = ";")) %>%
  distinct() %>%
  ungroup() %>% 
  pivot_wider(.,names_from = ID, values_from = variant) %>% 
  column_to_rownames(var = "gene") %>% 
  as.data.frame()

mat <- as.matrix(mat)
mat[is.na(mat)] <- ""


################################################################################
### Colors 
################################################################################

col = c("background"="#CCCCCC",
        "frameshift indel" ="orange", 
        "inframe indel"="yellowgreen", 
        "splicing"="maroon3", #"pink", 
        "unknown_splice" = "darkorchid4",
        "missense"="forestgreen", 
        "synonymous_variant"="mediumorchid2",
        "stop gain"="gold", 
        "start_lost" = "yellow2", 
        "stop_lost" = "red",
        "unknown"="darkgrey",
        "DEL" = 'deepskyblue',
        "AMP" = 'red3',
        "promoter" = "dodgerblue4", 
        "fusion" = "black") 

################################################################################
### Labelling with oncoKB
################################################################################
oncogene <- read.csv("oncoKB/cancerGeneList.tsv", header = TRUE, sep = "\t")

oncogene <- mutate(oncogene, class = case_when(
  Is.Oncogene == "Yes" & Is.Tumor.Suppressor.Gene == "Yes"  ~ "Both",
  Is.Oncogene == "Yes" & Is.Tumor.Suppressor.Gene == "No"  ~ "Oncogene",
  Is.Oncogene == "No" & Is.Tumor.Suppressor.Gene == "Yes"  ~ "TSG",
  Is.Oncogene == "No" & Is.Tumor.Suppressor.Gene == "No"  ~ "Unknown")) %>%
  filter(class != "Unknown") %>%     # Remove unknown genes
  mutate(class = factor(class, levels = c("Oncogene", "TSG", "Both")))


################################################################################
### Processing for oncoprint 
################################################################################
### Alteration matrix

fqseuil <- 0.04
fq <- apply(mat, 1, function (x) sum((x != ""), na.rm = TRUE) / ncol(mat))

mat.filtered <- mat.filtered[rownames(mat.filtered) %in% oncogene$Hugo.Symbol,] ## remove unknown genes and genes not included in oncoKb

fq.filtered <- apply(mat.filtered, 1, function (x) sum((x != ""), na.rm = TRUE) / ncol(mat.filtered))
fqOrdered <- fq.filtered[order(fq.filtered, decreasing = TRUE)]
oncogene_order <- names(fqOrdered[which(names(fqOrdered) %in% rownames(mat.filtered))])


################################################################################
### Figure parameters 
################################################################################
heatmap_legend_param = list(title = "Alterations", 
                            at = c("AMP", "DEL", "stop gain", "frameshift indel","inframe indel", "missense", "splicing", "promoter","fusion","constit"), 
                            labels = c("AMP", "DEL", "Stop gain", "Frameshift indel", "Inframe indel", "Missense", "Splicing", "Promoter", "Fusion", "Constitutional"),
                            ncol = 1,
                            legend_gp = gpar(fontsize = 8))

alter_fun_list = list(
  background = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA)),
  AMP = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["AMP"], col = NA)), 
  DEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["DEL"], col = NA)),
  `stop gain` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,  gp = gpar(fill = col["stop gain"], col = NA)),
  `frameshift indel` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.6, gp = gpar(fill = col["frameshift indel"], col = NA)),
  `inframe indel` = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.6, gp = gpar(fill = col["inframe indel"], col = NA)),
  missense = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["missense"], col = NA)),
  splicing = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4,gp = gpar(fill = col["splicing"], col = NA)),
  stop_lost = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.2, gp = gpar(fill = col["stop_lost"], col = NA)),
  start_lost = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.2, gp = gpar(fill = col["start_lost"], col = NA)),
  promoter = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.2, gp = gpar(fill = col["promoter"], col = NA)),
  fusion = function(x, y, w, h) {
    grid.segments(x - w*0.3, y - h*0.3, x + w*0.3, y + h*0.3, gp = gpar(lwd = 1, col = col["fusion"]))
    grid.segments(x + w*0.3, y - h*0.3, x - w*0.3, y + h*0.3, gp = gpar(lwd = 1, col = col["fusion"]))},
  constit = function(x, y, w, h) 
    grid.points(x, y, pch = 16, size = unit(2, "mm")) 
)


################################################################################
### Add clinical data
################################################################################

clinical.data <- read.table("clinical_data.csv",sep=",", header = TRUE)

clinical.data <- clinical.data %>% 
  mutate(RCB = sub('I/II/III', 'I, II, III', RCB),
         Tstage = sub('/', ', ', Tstage)) %>% 
  mutate(Age = factor(Age, levels = c(">= 50 years old","< 50 years old")),
         Subcohort = factor(Subcohort, levels = c("Neoadjuvant therapy","Upfront surgery")), 
         Tstage = factor(replace_na(as.character(Tstage), "Unknown"), levels = c("T1, T2", "T3, T4","Unknown")),
         Nstage = factor(replace_na(as.character(Nstage), "Unknown"), levels = c("N0","N+","Unknown")),
         RCB = factor(replace_na(as.character(RCB), "Unknown"), levels = c("0","I, II, III", "Unknown")),
         DFS =  factor(replace_na(as.character(DFS), "Unknown"), levels = c("Yes","No","Unknown")),
         HRD = factor(replace_na(as.character(HRD), "Unknown"), levels = c("HRD","HRP","Unknown")),
         TMB = factor(TMB, levels=c("High","Low")),
         MSI = factor(MSI, levels=c("MSI","MSS"))) %>%
  mutate(RNA.analysis = factor(ifelse(ID %in% splan.rna, "Yes","No"),levels=c("Yes","No"))) %>% 
  column_to_rownames(var="ID")



################################################################################
### Oncoprint
################################################################################

heatmap.annot <- HeatmapAnnotation(Age = clinical.data[colnames(mat.filtered),"Age"],
                                   Subcohort = clinical.data[colnames(mat.filtered),"Subcohort"],
                                   `T stage` = clinical.data[colnames(mat.filtered),"Tstage"],
                                   `N stage` = clinical.data[colnames(mat.filtered),"Nstage"],
                                   Event = clinical.data[colnames(mat.filtered),"DFS"],
                                   `RCB status` = clinical.data[colnames(mat.filtered),"RCB"],
                                   TMB = clinical.data[colnames(mat.filtered),"TMB"],
                                   MSI = clinical.data[colnames(mat.filtered),"MSI"],
                                   HRD = clinical.data[colnames(mat.filtered),"HRD"],
                                   `Tumor with RNAseq analysis` = clinical.data[colnames(mat.filtered),"RNA.analysis"],
                                   
                                   annotation_legend_param = list(
                                     Age = list(title = "Age", at = c(">= 50 years old","< 50 years old"), labels = expression("" >= 50 ~ years ~ old, "" < 50 ~ years ~ old))),
                                   
                                   cbar = anno_oncoprint_barplot(height = unit (2, "cm")),
                                   
                                   annotation_name_gp = gpar(fontsize = 9),

                                   col = list(
                                     TMB = c("High" = "#FFB800", "Low" = "#FFF9B1"),
                                     MSI = c("MSI" = "#D97A32", "MSS" = "#FEE8C8"),
                                     HRD=c("HRD"="red3" ,"HRP"="#FFF5EE","Unknown" = "grey75"),
                                     
                                     
                                     Age = c("< 50 years old" = "#EB949D", ">= 50 years old" = "#FBD1D6"),
                                     Subcohort = c("Neoadjuvant therapy" = "#C3E6DC", "Upfront surgery" = "#62AA9D"),
                                     `N stage` = c("N0" = "thistle", "N+" = "orchid4","Unknown" = "grey75"),
                                     Event = c("Yes" = "#7CB5FF", "No" = "#F0F8FF","Unknown" = "grey75"),
                                     
                                     `RCB status` = c("0" = "#cbe1e1", "Unknown" = "grey75", "I, II, III"= "#008b8b"),
                                     `T stage` = c("T1, T2" = "lightsteelblue2", "T3, T4"= "dodgerblue3","Unknown" = "grey75"),
                                     `Tumor with RNAseq analysis` = c("Yes"="#A8D5BA", "No"="darkgreen")
                                   ),
                                   
                                   show_legend = TRUE)



op1 <- oncoPrint(mat.filtered,
                 alter_fun = alter_fun_list,
                 col = col,
                 heatmap_legend_param = heatmap_legend_param,
                 column_names_rot = 45,
                 column_gap = unit(20,"mm"),
                 row_order = oncogene_order,
                 row_gap = unit(2, "mm"),
                 row_names_gp = gpar(fontsize = 9), pct_gp = gpar(fontsize = 8),
                 row_split = oncogene$class[match(rownames(mat.filtered), oncogene$Hugo.Symbol)],
                 gap = unit(c(10),"mm"),
                 top_annotation = heatmap.annot,
                 remove_empty_columns = FALSE, show_column_names = FALSE, alter_fun_is_vectorized = TRUE)


# pdf("SCAN-TNBC_oncoprint.pdf", width=16, height=10)
draw(op1, heatmap_legend_side = "left", annotation_legend_side = "right")
# dev.off()




