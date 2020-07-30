version 1.0

workflow STAAR_TOPMed_genome_wide {
    input {
        File pheno_file
        File kinship_file
        String sample_id
        String group_id = "none"
        String outcome
        String outcome_type = "continuous"
        String covariates = "none"
        Array[File] geno_files
        File? agg_file
        Int num_cores = 1
        Int step_length = 1000
        Int window_length = 2000
    }

    call run_null_model {
        input:
            sample_id = sample_id,
            group_id = group_id,
            pheno_file = pheno_file,
            outcome = outcome,
            outcome_type = outcome_type,
            covariates = covariates,
            kinship_file = kinship_file
    }

    scatter (geno_file in geno_files) {
        call run_genetic_region {
            input:
                null_file = run_null_model.null_model,
                geno_file = geno_file,
                num_cores = num_cores,
                step_length = step_length,
                window_length = window_length
        }
    }
    call run_summarize {
        input:
            results = run_genetic_region.results
    }

    output {
        File null_model = run_null_model.null_model
        File summarized_result = run_summarize.summarized_result
    }

    parameter_meta {
        null_file_precompute: "Optional precomputed null model."
        pheno_file: "Phenotypic data for sample set."
        kinship_file: "Optional relatedness matrix for sample set."
        sample_id: "Column header name of samples."
        group_id: "Optional column header name of group for groups in heteroscedastic linear model."
        outcome: "Column header name of outcome."
        outcome_type: "Type of outcome: continuous or binary."
        covariates: "Optional comma-separated column header names of covariates."
        null_memory: "Requested memory for null model (GB)."
        null_disk: "Requested disk space for null model (GB)."
        geno_files: "Array of genotypes in annotated GDS files."
        agg_file: "Optional file of aggregation units for gene-centric tests otherwise genetic region performed."
        num_cores: "Number of cores for multi-threading during test task."
        test_memory: "Requested memory for aggregate tests (GB)."
        test_disk: "Requested disk space for aggregate tests (GB)."
        step_length: "Step size for sliding windows (bp)."
        window_length: "Window length for sliding windows (bp)."
    }

    meta {
        author: "Sheila Gaynor"
        email: "sheilagaynor@hsph.harvard.edu"
        description: "Run and summarize genome-wide STAAR rare variant tests on annotated TOPMed Freeze 8 data."
    }
}

task run_null_model {
    input {
        File pheno_file
        File kinship_file
        String sample_id
        String group_id
        String outcome
        String outcome_type
        String covariates
    }
    command {
        Rscript /Users/sheilagaynor/Desktop/staar_topmed/STAAR_null_model.R  ${pheno_file} ${kinship_file} ${sample_id} ${group_id} ${outcome} ${outcome_type} ${covariates}
    }
    output {
        File null_model = "null_model.Rds"
    }
}

task run_genetic_region {
    input {
        File null_file
        File geno_file
        Int num_cores
        Int step_length
        Int window_length
    }
    command {
        Rscript /Users/sheilagaynor/Desktop/staar_topmed/STAAR_genetic_region.R ${null_file} ${geno_file} ${num_cores} ${step_length} ${window_length}
    }
    output {
        File results = "genetic_region.Rds"
    }
}

task run_gene_centric {
    input {
        File null_file
        File geno_file
        File agg_file
        Int num_cores
    }
    command {
        Rscript /Users/sheilagaynor/Desktop/staar_topmed/STAAR_gene_centric.R ${null_file} ${geno_file} ${agg_file} ${num_cores}
    }
    output {
        File results = "gene_centric.Rds"
    }
}

task run_summarize {
    input {
        Array[File] results
    }
    command {
        Rscript /Users/sheilagaynor/Desktop/staar_topmed/STAAR_summarize.R ${sep = ',' results}
    }
    output {
        File summarized_result = "result_out.Rds"
    }
}
