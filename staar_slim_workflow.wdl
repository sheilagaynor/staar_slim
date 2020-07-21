version 1.0

workflow STAAR_SLIM {
    input {
        String sample_id
        String group_id
        File pheno_file
        String outcome
        String covariates
        File kinship_file
        Int num_cores
        File geno_file
        Int step_length
        Int window_length
    }
    call run_null_model {
        input:
            sample_id = sample_id,
            group_id = group_id,
            pheno_file = pheno_file,
            outcome = outcome,
            covariates = covariates,
            kinship_file = kinship_file
    }
    call run_genetic_region {
        input:
            num_cores=num_cores,
            null_file=run_null_model.null_model,
            geno_file=geno_file,
            step_length=step_length,
            window_length=window_length
    }
    output {
        File null_model = run_null_model.null_model
        File results = run_genetic_region.results
    }
}

task run_null_model {
    input {
        String sample_id
        String group_id
        File pheno_file
        String outcome
        String covariates
        File kinship_file
    }
    command {
        Rscript /Users/sheilagaynor/Desktop/staar-workflow/STAAR_null_model.R ${sample_id} ${group_id} ${pheno_file} ${outcome} ${covariates} ${kinship_file}
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
    input {
        Int num_cores
        File null_file
        File geno_file
        Int step_length
        Int window_length
    }
    command {
        Rscript /Users/sheilagaynor/Desktop/staar-workflow/STAAR_genetic_region.R ${num_cores} ${null_file} ${geno_file} ${step_length} ${window_length}
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
