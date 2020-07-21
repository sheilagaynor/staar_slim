# Load packages
library(gdsfmt); library(SeqArray); library(STAAR)
library(SeqVarTools); library(dplyr); library(doMC)

# Parse arguments
args <- commandArgs(T)
# Compute inputs
num_cores <- as.numeric(args[1])
# File inputs
null_file <- args[2]
geno_file <- args[3]
# Window inputs
step_length <- as.numeric(args[4])
window_length <- as.numeric(args[5])

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

# Define sample and variants of interest
pheno_id <- as.character(null_model$id_include)
filter <- seqGetData(geno, "annotation/filter")
AVGDP <- seqGetData(geno, "annotation/info/AVGDP")
SNVlist <- filter == "PASS" & AVGDP > 10 & isSNV(geno)
variant_id <- seqGetData(geno, "variant.id")
rm(filter); rm(AVGDP)

# Define windows by position, parameters
position <- as.numeric(seqGetData(geno, "position"))
chrom <- unique(as.numeric(seqGetData(geno, "chromosome")))
variant_range <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr=chrom, start=min(position), end=max(position)))
windows <- GenomicRanges::slidingWindows(variant_range, window_length, step_length)
n_windows <- length(windows@unlistData@ranges)

# Set up potential multi-core analysis
n_cores <- min(c(num_cores, parallel::detectCores(logical = TRUE)))
# Break into larger chunks to read in annotations, faster than doing per window
window_chunks <- split(1:n_windows, ceiling(seq_along(1:n_windows)/75))
n_window_chunks <- length(window_chunks)

if(n_cores > 1) {
  doMC::registerDoMC(cores = n_cores)
  mc_options <- list(preschedule=FALSE, set.seed=FALSE)
  out <- foreach(i=1:n_window_chunks, .combine=rbind, .inorder=FALSE, .options.multicore = mc_options) %dopar% {
    results <- c()
    # Subset to this chunk of the chromosome
    range_chunk <- windows[[1]][window_chunks[[i]]]
    is.in <- (SNVlist) & (position>=min(range_chunk@ranges@start)) & (position<=max(range_chunk@ranges@start))
    if (sum(is.in) >= 2){
      seqSetFilter(geno,variant.id=variant_id[is.in],sample.id=pheno_id)
      # Extract data from gds, in order aligning with phenotype
      match_genopheno <- match(pheno_id,seqGetData(geno,"sample.id"))
      geno_chunk <- seqGetData(geno, "$dosage")
      geno_chunk <- geno_chunk[match_genopheno,]
      annot_df <- prepare_annot(geno)
      for (j in 1:length(range_chunk)){
        chunk_range <- GenomicRanges::granges(geno)
        set_id <- GenomicRanges::findOverlaps(chunk_range, range_chunk[j])@from
        geno_set <- geno_chunk[,set_id]
        results_temp <- c(chrom, range_chunk[j]@ranges@start, range_chunk[j]@ranges@start+range_chunk[j]@ranges@width-1)
        phred_set <- annot_df[set_id,]
        pvalues <- 0
        try(pvalues <- STAAR(geno_set, null_model, phred_set))
        if(class(pvalues)=="list"){
          results_temp <- c(results_temp,pvalues$num_variant,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                            pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
          results <- rbind(results,results_temp)
        }
      }
      results
     }
   }
 } else {
   out <- c()
   for (i in 1:n_window_chunks){
     # Subset to this chunk of the chromosome
     range_chunk <- windows[[1]][window_chunks[[i]]]
     is.in <- (SNVlist) & (position>=min(range_chunk@ranges@start)) & (position<=max(range_chunk@ranges@start))
     if (sum(is.in) >= 2){
       seqSetFilter(geno,variant.id=variant_id[is.in],sample.id=pheno_id)
       # Extract data from gds, in order aligning with phenotype
       match_genopheno <- match(pheno_id,seqGetData(geno,"sample.id"))
       geno_chunk <- seqGetData(geno, "$dosage")
       geno_chunk <- geno_chunk[match_genopheno,]
       annot_df <- prepare_annot(geno)
       for (j in 1:length(range_chunk)){
         chunk_range <- GenomicRanges::granges(geno)
         set_id <- GenomicRanges::findOverlaps(chunk_range, range_chunk[j])@from
         geno_set <- geno_chunk[,set_id]
         out_temp <- c(chrom, range_chunk[j]@ranges@start, range_chunk[j]@ranges@start+range_chunk[j]@ranges@width-1)
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
   }
  seqResetFilter(geno)
 }
}

# Clean up results
if(!is.null(out))
{
  colnames(out) <- colnames(out, do.NULL = FALSE, prefix = "col")
  colnames(out)[(length(colnames(out))-1):length(colnames(out))] <- c("ACAT-O","STAAR-O")
  colnames(out)[1:4] <- c('chr', 'start', 'end', 'n_snv')
}
seqClose(geno)

# Save output
saveRDS(out, 'genetic_region.Rds')
