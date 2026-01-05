# Survival curve script
# Author: Sophie Vacher
# R version 4.1.2 (2021-11-01)

library(survival)
library(survminer)

# Define path on the cluster
countsdatapath <-"path"

# Load data
data <- read.table(file.path(countsdatapath,"PourFig.csv"),sep=",", header = TRUE)
head(data)
dim(data)
unique(data$Statut_HRD)

# Create survfit object
surv_obj <- Surv(time = data$DEL_OS, event = data$VITAL_STATUS)
surv_obj.2 <- Surv(time = data$DEL_DFS, event = data$DFS)

# Adjust the survival curve by risk group
fit <- survfit(surv_obj ~ Statut_HRD, data = data)
fit.2 <- survfit(surv_obj.2 ~ Statut_HRD, data = data)

# Survival curve with risk profile
# OS
ggsurvplot(
  fit, 
  data = data,
  risk.table = TRUE,       # Affiche le profil de risque en dessous
  pval = TRUE,             # Affiche la valeur p du log-rank test
  conf.int = FALSE,         # Intervalle de confiance
  palette = c("red","blue"), # Couleurs par groupe "#E7B800","#2E9FDF"
  legend.labs = c("HRD", "HRP"),
  legend.title = NULL,
  risk.table.title = NULL,
  xlab = "Time (months)", 
  ylab = "OS (%)",
  risk.table.height = 0.25 # Hauteur de la table de risque
)

# DFS
ggsurvplot(
  fit.2, 
  data = data,
  risk.table = TRUE,       # Affiche le profil de risque en dessous
  pval = TRUE,             # Affiche la valeur p du log-rank test
  conf.int = FALSE,         # Intervalle de confiance
  palette = c("red","blue"), # Couleurs par groupe "#E7B800","#2E9FDF"
  legend.labs = c("HRD", "HRP"),
  xlab = "Time (months)", 
  ylab = "DFS (%)",
  risk.table.height = 0.25 # Hauteur de la table de risque
)

