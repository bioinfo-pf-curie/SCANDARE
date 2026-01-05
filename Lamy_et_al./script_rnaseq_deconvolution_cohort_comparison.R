### Script to generate boxplots of the cellular fraction per cohort
# Cohort : SCANDARE ovarian, SCANDARE TNBC, SCANDARE HNSCC

### load R packages
library(dplyr)
library(immunedeconv)
library(RColorBrewer)
library(ggplot2)


################################################################################
### load quantiseq deconvolution result (see script_rnaseq_quantiseq_deconvolution.R)
###  Cohort: ovarian
res_quantiseq_ovarian <- readRDS(file = "SCAN-ovarian_res_quantiseq.rds")

res_quantiseq_ovarian$cell_type <- factor(res_quantiseq_ovarian$cell_type, levels = res_quantiseq_ovarian$cell_type)
color_order = c(brewer.pal(n = 11, name = 'PuOr')[-7], brewer.pal(n = 11, name = 'PuOr')[7])

### Focus on Baseline tumor
splan.ovarian.table <- read.table("clinical_data_SCANDARE_ovarian.csv", header = TRUE, sep = ",")
bssample <- splan.ovarian.table %>% filter(Timepoint == "Baseline") %>% rownames()
res_quantiseq_baseline <- res_quantiseq_ovarian[,c("cell_type",bssample)]

################################################################################
### load quantiseq deconvolution result (see script_rnaseq_quantiseq_deconvolution.R)
### Cohort: TNBC
res_quantiseq_TNBC <- readRDS(file = "SCAN-TNBC_res_quantiseq.rds")

################################################################################
### load quantiseq deconvolution result (see script_rnaseq_quantiseq_deconvolution.R)
### Cohort: HNSCC
res_quantiseq_HNSCC <- readRDS(file = "SCAN-HNSCC_res_quantiseq.rds")

################################################################################
df.ovarian <- res_quantiseq_baseline %>% 
  tibble::column_to_rownames("cell_type") %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("sampleID") %>% 
  mutate(cohort = "Ovarian")

df.tnbc <- res_quantiseq_TNBC %>% 
  tibble::column_to_rownames("cell_type") %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("sampleID") %>% 
  mutate(cohort = "TNBC")

df.HNSCC <- res_quantiseq_HNSCC %>% 
  tibble::column_to_rownames("cell_type") %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("sampleID") %>% 
  mutate(cohort = "HNSCC")

df.full <- full_join(df.ovarian, df.tnbc) %>% full_join(df.HNSCC) %>% 
  mutate(cohort = factor(cohort)) %>% 
  tidyr::pivot_longer(cols = -c(sampleID, cohort), 
                      names_to = "CellType", 
                      values_to = "Value") %>% 
  mutate(CellType = factor(CellType)) %>% 
  filter(CellType != "uncharacterized cell")  %>% 
  as.data.frame()

################################################################################
### Boxplot figure

# pdf("Deconvolution_boxplot_per_cohort_stroma_padj.pdf", width = 12)
ggplot(df.full, aes(x=cohort, y=Value, fill=cohort)) +
  geom_boxplot(width=0.9, outlier.size = 0.5 ) +
  ggpubr::stat_compare_means(aes(group = cohort),
                             comparisons = list(c("HNSCC","Ovarian"), c("Ovarian","TNBC"),c("TNBC","HNSCC")),
                             method = "wilcox.test", label = "p.adj", paired=FALSE) +
  facet_wrap(~ CellType, nrow = 2) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +  # add more space on top
  labs(x="",y="Cell fractions") +
  theme_classic() + 
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))  +
  scale_fill_manual(values= c("#239D76","#3B6CAF","#EF077E"))
# dev.off()

