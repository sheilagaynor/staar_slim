workflow STAAR_SLIM {

    # run_null_model inputs
    File? null_file_precompute
    File? pheno_file
    String? sample_id
    String? group_id = "NA"
    String? outcome
    String? outcome_type = "continuous"
    String? covariates = "NA"
    File? kinship_file
    Int? null_memory = 25
    Int? null_disk = 50

    # run_genetic_region or run_gene_centric inputs
    File? agg_file
    Int num_cores
    Array[File] geno_files
    Int step_length = 1000
    Int window_length = 2000
    Int test_memory = 25
    Int test_disk = 50

    # run_summarize inputs
    Int summary_memory = 25
    Int summary_disk = 50

    if (!defined(null_file_precompute)) {
        call run_null_model {
            input:
                #File inputs
                pheno_file = pheno_file,
                kinship_file = kinship_file,
                #ID names
                sample_id = sample_id,
                group_id = group_id,
                #Variable names
                outcome = outcome,
                outcome_type = outcome_type,
                covariates = covariates
        }
    }

    File null_file = select_first([null_file_precompute, run_null_model.null_model])

    scatter (geno_file in geno_files) {
        call run_genetic_region {
            input:
                #File inputs
                null_file=null_file,
                geno_file=geno_file,
                #Compute inputs
                num_cores=num_cores,
                #Window inputs
                step_length=step_length,
                window_length=window_length
        }
    }

    call run_summarize {
        input:
            results = run_genetic_region.results
    }

    output {
        File null_model = null_file
        #File summarized_results = run_summarize.summarized_result
        #Array[File] summarized_results = run_genetic_region.results
    }
}

task run_null_model {

    File pheno_file
    File kinship_file
    String sample_id
    String group_id
    String outcome
    String outcome_type
    String covariates

    command {
        Rscript /STAAR_null_model.R ${pheno_file} ${kinship_file} ${sample_id} ${group_id} ${outcome} ${outcome_type} ${covariates}
    }
    runtime {
        docker: "quay.io/sheilagaynor/staar_slim"
        memory: "20 GB"
        disks: "local-disk 1 HDD"
    }
    output {
        File null_model = "null_model.Rds"
    }
}

task run_genetic_region {

    File null_file
    File geno_file
    Int num_cores
    Int step_length
    Int window_length

    command {
        Rscript /STAAR_genetic_region.R ${null_file} ${geno_file} ${num_cores} ${step_length} ${window_length}
    }
    runtime {
        docker: "quay.io/sheilagaynor/staar_slim"
        memory: "20 GB"
        disks: "local-disk 4 HDD"
    }
    output {
        File results = "genetic_region.Rds"
    }
}

task run_summarize {

    Array[File] results

    command {
        Rscript /STAAR_summarize.R ${sep = ',' results}
    }
    runtime {
        docker: "quay.io/sheilagaynor/staar_slim"
        memory: "8 GB"
        disks: "local-disk 1 HDD"
    }
    output {
        File summarized_result = "staar_out.Rds"
    }
}
