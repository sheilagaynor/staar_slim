version 1.0

workflow STAAR_TOPMed_genome_wide {
    File? null_file_precompute
    File pheno_file
    File kinship_file
    String sample_id
    String group_id = "none"
    String outcome
    String outcome_type = "continuous"
    String covariates = "none"
    Int null_memory = 25
    Int null_disk = 50

    # run_genetic_region or run_gene_centric inputs
    Array[File] geno_files
    File? agg_file
    Int num_cores = 1
    Int test_memory = 25
    Int test_disk = 50
    Int step_length = 1000
    Int window_length = 2000

    if (!defined(null_file_precompute)) {
        call run_null_model {
            input:
                sample_id = sample_id,
                group_id = group_id,
                pheno_file = pheno_file,
                outcome = outcome,
                covariates = covariates,
                kinship_file = kinship_file
        }
    }

    File? null_file = if (defined(null_file_precompute)) then null_file_precompute else run_null_model.null_model

    if (!defined(agg_file)) {
        scatter (geno_file in geno_files) {
            call run_genetic_region {
                input:
                    null_file = null_file,
                    geno_file = geno_file,
                    num_cores = num_cores,
                    test_memory = test_memory,
                    test_disk = test_disk,
                    step_length = step_length,
                    window_length = window_length
            }
        }
        call run_summarize {
            input:
                results = run_genetic_region.results
        }
    }

    if (defined(agg_file)) {
        scatter (geno_file in geno_files) {
            call run_gene_centric {
                input:
                    null_file = null_file,
                    geno_file = geno_file,
                    agg_file = agg_file,
                    num_cores = num_cores,
                    test_memory = test_memory,
                    test_disk = test_disk
            }
        }
        call run_summarize {
            input:
                results = run_gene_centric.results
        }
    }

    output {
        File null_model = null_file
        File result_out = run_summarize.result_out
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
        Int null_memory
        Int null_disk
    }
    command {
        Rscript /STAAR_null_model.R  ${pheno_file} ${kinship_file} ${sample_id} ${group_id} ${outcome} ${outcome_type} ${covariates}
    }
    runtime {
        docker: "quay.io/sheilagaynor/staar_slim"
        memory: "${null_memory} GB"
        disks: "local-disk ${null_disk} HDD"
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
        Int test_memory
        Int test_disk
        Int step_length
        Int window_length
    }
    command {
        Rscript /STAAR_genetic_region.R ${null_file} ${geno_file} ${num_cores} ${step_length} ${window_length}
    }
    runtime {
        docker: "quay.io/sheilagaynor/staar_slim"
        memory: "${test_memory} GB"
        disks: "local-disk ${test_disk} HDD"
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
        Int test_memory
        Int test_disk
    }
    command {
        Rscript /STAAR_gene_centric.R ${null_file} ${geno_file} ${agg_file} ${num_cores}
    }
    runtime {
        docker: "quay.io/sheilagaynor/staar_slim"
        memory: "${test_memory} GB"
        disks: "local-disk ${test_disk} HDD"
    }
    output {
        File results = "gene_centric.Rds"
    }
}

task run_summarize {
    input {
        Array[File] results
        Int summary_memory
        Int summary_disk
    }
    command {
        Rscript /Summarize.R ${sep = ',' assoc}
    }
    runtime {
        docker: "quay.io/sheilagaynor/staar_slim"
        memory: "${summary_memory} GB"
        disks: "local-disk ${summary_disk} HDD"
    }
    output {
        File results = "result_out.Rds"
    }
}
