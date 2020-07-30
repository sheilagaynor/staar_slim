# Description: Generate a null model using the STAAR package using a wrapper of the GMMAT package.
# Inputs:
# null_file : file containing output from null model fitting via STAAR (.Rds)
# geno_file : annotated GDS file containing the given annotation channels (.gds)
# agg_file : file variant masks, with a column for gene name and a column for unique variant identifier (.Rds)
# num_cores : number of cores for multi-threading when looping tests (numeric)

# Load packages
library(gdsfmt); library(SeqArray); library(STAAR)
library(SeqVarTools); library(dplyr); library(doMC)

# Parse arguments
args <- commandArgs(T)
# File inputs
null_file <- args[1]
geno_file <- args[2]
agg_file <- args[3]
# Compute inputs
num_cores <- as.numeric(args[4])

# Prepare annotations
prepare_annot <- function(geno){
  CADD.PHRED <- seqGetData(geno, "annotation/info/TOPMedAnnotation/CADD.FULL/PHRED")
  CADD.PHRED[is.na(CADD.PHRED)] <- 0
  LINSIGHT.PHRED <- seqGetData(geno, "annotation/info/TOPMedAnnotation/LINSIGHT.PHRED.rounded")
  FATHMM.XF.PHRED <- seqGetData(geno, "annotation/info/TOPMedAnnotation/FATHMM.XF.PHRED.rounded")
  aPC.EpigeneticActive.PHRED <- seqGetData(geno, "annotation/info/TOPMedAnnotation/APC.PHRED.rounded/aPC.EpigeneticActive")
  aPC.EpigeneticRepressed.PHRED <- seqGetData(geno, "annotation/info/TOPMedAnnotation/APC.PHRED.rounded/aPC.EpigeneticRepressed")
  aPC.EpigeneticTranscription.PHRED <- seqGetData(geno, "annotation/info/TOPMedAnnotation/APC.PHRED.rounded/aPC.EpigeneticTranscription")
  aPC.Conservation.v2.PHRED <- seqGetData(geno, "annotation/info/TOPMedAnnotation/APC.PHRED.rounded/aPC.Conservation.v2")
  aPC.Protein.PHRED <- seqGetData(geno, "annotation/info/TOPMedAnnotation/APC.PHRED.rounded/aPC.Protein")
  aPC.LocalDiversity.v2.PHRED <- seqGetData(geno, "annotation/info/TOPMedAnnotation/APC.PHRED.rounded/aPC.LocalDiversity.v2")
  aPC.RegulatoryDistance.PHRED <- seqGetData(geno, "annotation/info/TOPMedAnnotation/APC.PHRED.rounded/aPC.RegulatoryDistance")
  aPC.TF.PHRED <- seqGetData(geno, "annotation/info/TOPMedAnnotation/APC.PHRED.rounded/aPC.TF")
  aPC.Dist2TSSTES.PHRED <- seqGetData(geno, "annotation/info/TOPMedAnnotation/APC.PHRED.rounded/aPC.Dist2TSSTES")
  aPC.MicroRNA.PHRED <- seqGetData(geno, "annotation/info/TOPMedAnnotation/APC.PHRED.rounded/aPC.MicroRNA")
  phred_df <- data.frame(CADD.PHRED,LINSIGHT.PHRED,FATHMM.XF.PHRED,
                         aPC.EpigeneticActive.PHRED,aPC.EpigeneticRepressed.PHRED,aPC.EpigeneticTranscription.PHRED,
                         aPC.Conservation.v2.PHRED,aPC.Protein.PHRED,aPC.LocalDiversity.v2.PHRED,
                         aPC.RegulatoryDistance.PHRED,aPC.TF.PHRED,aPC.Dist2TSSTES.PHRED,aPC.MicroRNA.PHRED)
  return(phred_df) }


# Read in files: null model, genotypes
null_model <- readRDS(null_file)
geno <- seqOpen(geno_file)
agg_units <- readRDS(agg_file)

# Define sample and variants of interest
pheno_id <- as.character(null_model$id_include)
variant_id <- seqGetData(geno, "variant.id")

# Set up potential multi-core analysis
n_cores <- min(c(num_cores, parallel::detectCores(logical = TRUE)))
# Break into larger chunks to read in annotations, faster than doing per window
gene_names <- unique(agg_units$gene)
n_groups <- length(gene_names)
group_chunks <- split(1:n_groups, ceiling(seq_along(1:n_groups)/50))
n_group_chunks <- length(group_chunks)

if(n_cores > 1) {
  doMC::registerDoMC(cores = n_cores)
  mc_options <- list(preschedule=FALSE, set.seed=FALSE)
  out <- foreach(i=1:n_group_chunks, .combine=rbind, .inorder=FALSE, .options.multicore = mc_options) %dopar% {
      results <- c()
      # Subset to this chunk of the chromosome
      agg_units_chunk <- agg_units[agg_units$gene %in% gene_names[group_chunks[[i]]],]
      is.in <- (variant_id %in% as.character(unique(agg_units_chunk$variant)))
      seqSetFilter(geno,variant.id=variant_id[is.in],sample.id=pheno_id)
      # Extract data from gds, in order aligning with phenotype
      match_genopheno <- match(pheno_id,seqGetData(geno,"sample.id"))
      geno_chunk <- seqGetData(geno, "$dosage")
      geno_chunk <- geno_chunk[match_genopheno,]
      annot_df <- prepare_annot(geno)
      variant_chunk <- seqGetData(geno, "variant.id")
      for (j in 1:length(unique(agg_units_chunk$gene))){
        set_id <- variant_chunk %in% agg_units_chunk$variant[agg_units_chunk$gene == unique(agg_units_chunk$gene)[j]]
        geno_set <- geno_chunk[,set_id]
        results_temp <- c(chr, unique(agg_units_chunk$gene)[j])
        phred_set <- annot_df[set_id,]
        pvalues <- 0
        try(pvalues <- STAAR(geno_set, null_model, phred_set))
        if(class(pvalues)=="list"){
          results_temp <- c(results_temp,pvalues$num_variant,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                            pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
          results <- rbind(results,results_temp)
        }
      results
    }
  }
} else {
  out <- c()
  for (i in 1:n_group_chunks){
    # Subset to this chunk of the chromosome
    agg_units_chunk <- agg_units[agg_units$gene %in% gene_names[group_chunks[[i]]],]
    is.in <- (variant_id %in% as.character(unique(agg_units_chunk$variant)))
    seqSetFilter(geno,variant.id=variant_id[is.in],sample.id=pheno_id)
      # Extract data from gds, in order aligning with phenotype
      match_genopheno <- match(pheno_id,seqGetData(geno,"sample.id"))
      geno_chunk <- seqGetData(geno, "$dosage")
      geno_chunk <- geno_chunk[match_genopheno,]
      annot_df <- prepare_annot(geno)
      variant_chunk <- seqGetData(geno, "variant.id")
      for (j in 1:length(unique(agg_units_chunk$gene))){
        set_id <- variant_chunk %in% agg_units_chunk$variant[agg_units_chunk$gene == unique(agg_units_chunk$gene)[j]]
        geno_set <- geno_chunk[,set_id]
        out_temp <- c(chrom, unique(agg_units_chunk$gene)[j])
        phred_set <- annot_df[set_id,]
        pvalues <- 0
        try(pvalues <- STAAR(geno_set, null_model, phred_set))
        if(class(pvalues)=="list"){
          out_temp <- c(out_temp,pvalues$num_variant,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
          out <- rbind(out,out_temp)
        }
      }
    
    seqResetFilter(geno)
  }
}

# Clean up results
if(!is.null(out))
{
  colnames(out) <- colnames(out, do.NULL = FALSE, prefix = "col")
  colnames(out)[(length(colnames(out))-1):length(colnames(out))] <- c("ACAT-O","STAAR-O")
  colnames(out)[1:3] <- c('chr', 'gene', 'n_snv')
}
seqClose(geno)

# Save output
saveRDS(cbind(out,substr(agg_file,1,nchar(agg_file)-14)), 'gene_centric.Rds')
