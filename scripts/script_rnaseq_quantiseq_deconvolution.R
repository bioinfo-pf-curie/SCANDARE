### Script to generate barplot of cellular fraction per tumor
# Cohort : SCANDARE ovarian

### load R packages
library(dplyr)
library(immunedeconv)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(grid)

splan.table <- read.table("clinical.data.csv", header = TRUE, sep = ",")
splan <- splan.table %>% 
  mutate(ID.cBioportal = case_when(Timepoint == "Baseline" ~ paste0(ID.cBioportal, '_B'),
                                   Timepoint == "Post-NAC" ~ paste0(ID.cBioportal, '_PN'),
                                   Timepoint == "Rechute" ~ paste0(ID.cBioportal, '_R'))) %>% 
  mutate(PlatFS = factor(PlatFS, levels=c("Platin resistant","Platin intermediate","Platin sensitive"))) %>%
  mutate(INC_OV_HISTYP_GRADE = factor(INC_OV_HISTYP_GRADE, levels=c("Serous High","Serous Low","Other"))) %>% 
  mutate(FIGO.cat = factor(FIGO.cat, levels=c("I.II","III.IV"), labels=c("I, II", "III, IV"))) %>% 
  mutate(AGE.cat = factor(AGE.cat, levels=c("<65", ">=65"), labels=c("< 65 years old",">= 65 years old"))) %>% 
  mutate(ID.cBioportal = sub('Patient_OV_','P',ID.cBioportal)) %>% 
  tibble::column_to_rownames("ID.cBioportal") 

################################################################################
## Table ENSEMBL gene id to Hugo
gene.annot <- read.csv("tableannot.csv", header = TRUE, row.names = 1)

## TPM count table
tpm.counts <- read.csv("tablecounts_tpm.csv", row.names = 1, check.names = FALSE) 
d.tpm <- tpm.counts[,splan$SCANDARE.ID] %>% 
  tibble::column_to_rownames("gene") 

## Prepare matrix for deconvolution
exp.mat <- d.tpm %>% tibble::rownames_to_column("gene_id") %>% 
  left_join(gene.annot, by=c("gene_id")) %>% 
  select(-gene_id) %>% 
  distinct(gene_name, .keep_all=TRUE) %>% 
  tibble::column_to_rownames("gene_name") %>% 
  as.matrix()

## Deconvolution with quantiseq using immunedeconv R package
# res_quantiseq = deconvolute(exp.mat, "quantiseq", tumor = TRUE)
# saveRDS(res_quantiseq, "SCAN-ovarian_res_quantiseq.rds")
res_quantiseq <- readRDS(file = "SCAN-ovarian_res_quantiseq.rds.rds")
colnames(res_quantiseq) <- sub("Patient_OV_","P",colnames(res_quantiseq))

res_quantiseq$cell_type <- factor(res_quantiseq$cell_type, levels = res_quantiseq$cell_type)
color_order = c(brewer.pal(n = 11, name = 'PuOr')[-7], brewer.pal(n = 11, name = 'PuOr')[7])

################################################################################
### Baseline
################################################################################
bssample <- splan %>% filter(Timepoint == "Baseline") %>% rownames()

res_quantiseq_baseline <- res_quantiseq[,c("cell_type",bssample)]

res_quantiseq_baseline %>% 
  tidyr::gather(sample, fraction, -cell_type) %>%
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values = color_order) +
  scale_x_discrete(limits = rev(levels(res_quantiseq))) +
  xlab("") + 
  ylab("Cell fractions") +
  labs(fill="Cell type") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust= 1))


res_quantiseq_baseline_tumor <- res_quantiseq_baseline %>% filter(cell_type == "uncharacterized cell") %>% tibble::column_to_rownames("cell_type") 
res_quantiseq_baseline_tumor_t <- res_quantiseq_baseline_tumor %>% as.data.frame() %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("sample")
sample_order <- res_quantiseq_baseline_tumor_t[order(res_quantiseq_baseline_tumor_t$`uncharacterized cell`),"sample"]

## Cellular fraction barplot ordered according to "uncharacterized cell" cell type
res_quantiseq_baseline %>%
  tidyr::gather(sample, fraction, -cell_type) %>% 
  mutate(sample = factor(sample, levels = sample_order)) %>% 
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values = color_order) +
  scale_x_discrete(limits = rev(levels(res_quantiseq))) +
  ylab("Cell fractions") +
  labs(fill="Cell type") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust= 1))


################################################################################
### Longitudinal 
################################################################################

dup <- splan %>% filter(duplicated(SCANDARE.ID)) %>% select(SCANDARE.ID) %>% unique() %>% pull(SCANDARE.ID)          ## Select duplicated patients
doublesample <- splan %>% filter(SCANDARE.ID %in% dup) %>% rownames()           ## Extract duplicated samples

res_quantiseq_longitudinal <- res_quantiseq[,c("cell_type",doublesample)]

res_quantiseq_longitudinal %>% 
  tidyr::gather(sample, fraction, -cell_type) %>%
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values = color_order) +
  scale_x_discrete(limits = rev(levels(res_quantiseq))) +
  xlab("") + 
  ylab("Cell fractions") +
  labs(fill="Cell type") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust= 1))

res_quantiseq_longitudinal_tumor <- res_quantiseq_longitudinal %>% filter(cell_type == "uncharacterized cell") %>% tibble::column_to_rownames("cell_type") 
res_quantiseq_longitudinal_tumor_t <- res_quantiseq_longitudinal_tumor %>% as.data.frame() %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("sample")
sample_order <- res_quantiseq_longitudinal_tumor_t[order(res_quantiseq_longitudinal_tumor_t$`uncharacterized cell`),"sample"]

## Cellular fraction barplot ordered according to "uncharacterized cell" cell type
res_quantiseq_longitudinal %>%
  tidyr::gather(sample, fraction, -cell_type) %>% 
  mutate(sample = factor(sample, levels = sample_order)) %>% 
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values = color_order) +
  scale_x_discrete(limits = rev(levels(res_quantiseq))) +
  ylab("Cell fractions") +
  labs(fill="Cell type") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust= 1))


################################################################################
### STROMA only
################################################################################
res_quantiseq_stroma <- res_quantiseq %>% tibble::column_to_rownames("cell_type") %>% t() %>% 
  as.data.frame() %>% 
  mutate(total_stroma = 1-`uncharacterized cell`) %>% 
  select(-`uncharacterized cell`) %>% 
  mutate(`B cell` = `B cell` / total_stroma,
         `Macrophage M1` = `Macrophage M1` / total_stroma,
         `Macrophage M2` = `Macrophage M2` / total_stroma,
         Monocyte = Monocyte / total_stroma,
         Neutrophil = Neutrophil / total_stroma,
         `NK cell` = `NK cell` / total_stroma,
         `T cell CD4+ (non-regulatory)` = `T cell CD4+ (non-regulatory)` / total_stroma,
         `T cell CD8+` = `T cell CD8+` / total_stroma,
         `T cell regulatory (Tregs)` = `T cell regulatory (Tregs)` / total_stroma,
         `Myeloid dendritic cell` = `Myeloid dendritic cell` / total_stroma) %>% 
  select(-total_stroma) %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("cell_type") 

res_quantiseq_stroma$cell_type <- factor(res_quantiseq_stroma$cell_type) 
color_order = c(brewer.pal(n = 11, name = 'PuOr')[c(-7)])

################################################################################
### STROMA only - Baseline samples
################################################################################
res_quantiseq_stroma_baseline <- res_quantiseq_stroma[,c("cell_type",bssample)]

res_quantiseq_stroma_baseline_Treg <- res_quantiseq_stroma_baseline %>% filter(cell_type == "T cell regulatory (Tregs)") %>% tibble::column_to_rownames("cell_type") 
res_quantiseq_stroma_baseline_Treg_t <- res_quantiseq_stroma_baseline_Treg %>% as.data.frame() %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("sample")
sample_order <- res_quantiseq_stroma_baseline_Treg_t[order(res_quantiseq_stroma_baseline_Treg_t$`T cell regulatory (Tregs)`),"sample"]

# pdf("SCAN-ovarian_deconvolution_quantiseq_baseline_stroma_ordered_by_Tregs.pdf", width = 16, height = 10)
res_quantiseq_stroma_baseline %>%
  tidyr::gather(sample, fraction, -cell_type) %>% 
  mutate(sample = factor(sample, levels = sample_order)) %>% 
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values = color_order) +
  xlab("") + 
  ylab("Cell fractions") +
  labs(fill="Cell type") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust= 1),
        plot.margin = unit(c(0.1, 0.1, 0.1, 1), "cm"))
# dev.off()

################################################################################
### STROMA only - Longitudinal samples
################################################################################
res_quantiseq_stroma_longitudinal <- res_quantiseq_stroma[,c("cell_type",doublesample)]

res_quantiseq_stroma_longitudinal_tumor <- res_quantiseq_stroma_longitudinal %>% filter(cell_type == "T cell regulatory (Tregs)") %>% tibble::column_to_rownames("cell_type") 
res_quantiseq_stroma_longitudinal_tumor_t <- res_quantiseq_stroma_longitudinal_tumor %>% as.data.frame() %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("sample")
sample_order <- res_quantiseq_stroma_longitudinal_tumor_t[order(res_quantiseq_stroma_longitudinal_tumor_t$`T cell regulatory (Tregs)`),"sample"]

# pdf("SCAN-ovarian_deconvolution_quantiseq_longitudinal_stroma_ordered_by_Tregs.pdf", width = 16, height = 10)
res_quantiseq_stroma_longitudinal %>%
  tidyr::gather(sample, fraction, -cell_type) %>% 
  mutate(sample = factor(sample, levels = sample_order)) %>% 
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values = color_order) +
  xlab("") + 
  ylab("Cell fractions") +
  labs(fill="Cell type") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust= 1),
        plot.margin = unit(c(0.1, 0.1, 0.1, 1), "cm"))
# dev.off()


