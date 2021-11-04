# Process haplotype blocks data and
rule create_gwas_input:
    input:
        os.path.join(config["data_store_dir"], "recombination_blocks/{date_of_processing}/F2_{bin_length}.txt"),
    output:
        os.path.join(config["data_store_dir"], "association_testing/{date_of_assoc_test}_true/{bin_length}.rds"),
    log:
        os.path.join(config["working_dir"], "logs/create_gwas_input/{bin_length}.log")
    params:
        bin_length = lambda wildcards: wildcards.bin_length    
    resources:
        mem_mb = 50
    container:
        config["R"]
    script:
        "../scripts/create_gwas_input.R"    