# Description: Generate a null model using the STAAR package using a wrapper of the GMMAT package.
# Inputs:
# pheno_file: delimited text file that has been fully prepared for analysis, includes (transformed) outcome and covariates (.txt or .csv or .tsv)
# kinship_file: file containing sparse relatedness matrix for all samples (.Rdata)
# sample_id: column name for sample ids (string)
# group_id: column name for group for heteroscedastic errors (string)
# outcome: column name for outcome variable (string)
# outcome_type: type of outcome: 'continuous' or 'binary'
# covariates: column names for covariates, separated by commas (string)

# Load packages
suppressMessages(library(STAAR))

# Parse arguments
args <- commandArgs(T)
# File inputs
pheno_file <- args[1]
kinship_file <- args[2]
# ID names (header)
sample_id <- args[3]
group_id <- args[4]
# Variable names
outcome <- args[5]
outcome_type <- args[6]
covariates <- args[7]

# Read in files: phenotypes, kinship
pheno <- as.data.frame(read.table(pheno_file, stringsAsFactors=FALSE, header=T))
id_include <- pheno[,which(names(pheno)==sample_id)]
kinship <- load(kinship_file)
sub_kinship <- skm[which(row.names(skm) %in% id_include),which(colnames(skm) %in% id_include)]
sub_kinship <- sub_kinship[match(id_include,row.names(sub_kinship)),match(id_include,colnames(sub_kinship))]

# Null model formula
if (covariates=='NA'){
  covar_splt <- 1
} else {
  covar_splt <- strsplit(covariates, split=",")[[1]]
}
null_model_char <- paste0(outcome, "~", paste(covar_splt, collapse="+"))

# Fit, save null model
if (tolower(outcome_type)=='continuous'){
  if (group_id=='NA'){
    null_model <- STAAR::fit_null_glmmkin(as.formula(null_model_char), id = sample_id, use_sparse = TRUE,
                                          data = pheno, kins = sub_kinship, family = gaussian(link = "identity"))
  } else {
    null_model <- STAAR::fit_null_glmmkin(as.formula(null_model_char), id = sample_id, use_sparse = TRUE,
                                          groups = group_id, data = pheno, kins = sub_kinship, family = gaussian(link = "identity"))
  }
} 
if (tolower(outcome_type)=='binary'){
   null_model <- STAAR::fit_null_glmmkin(as.formula(null_model_char), id = sample_id, use_sparse = TRUE,
                                        data = pheno, kins = sub_kinship, family = binomial(link = "logit"))
  
}
saveRDS(null_model, file="null_model.Rds")
