# Load packages
library(STAAR)

# Parse arguments
args <- commandArgs(T)
# ID names (header)
sample_id <- args[1]
group_id <- args[2]
# File inputs
pheno_file <- args[3]
kinship_file <- args[6]
# Variable names
outcome <- args[4]
covariates <- args[5]

# Read in files: phenotypes, kinship
pheno <- as.data.frame(read.table(pheno_file, sep="", stringsAsFactors=FALSE, head=T))
kinship <- readRDS(kinship_file)

# Null model formula
covar_splt <- strsplit(covariates, split=" ")[[1]]
null_model_char <- paste0(outcome, "~", paste(covar_splt, collapse="+"))

# Fit, save null model
null_model <- STAAR::fit_null_glmmkin(as.formula(null_model_char), id = sample_id, use_sparse = TRUE,
                               groups = group_id, data = pheno, kins = kinship, family = gaussian(link = "identity"))
saveRDS(null_model, file="null_model.Rds")
