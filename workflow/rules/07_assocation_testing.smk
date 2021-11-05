# Process haplotype blocks data and create GWAS input
rule create_gwas_input:
    input:
        genotypes = os.path.join(config["data_store_dir"], "recombination_blocks/20211027/F2_{bin_length}.txt")
    output:
        os.path.join(config["data_store_dir"], "association_testing/{date_of_assoc_test}_true/inputs/{bin_length}.rds"),
    log:
        os.path.join(config["working_dir"], "logs/create_gwas_input/{date_of_assoc_test}/{bin_length}.log")
    params:
        bin_length = lambda wildcards: wildcards.bin_length,
    resources:
        mem_mb = 50000
    container:
        config["R"]
    script:
        "../scripts/create_gwas_input.R"

rule test_gwls:

rule run_gwls:
    input:
        gt_pos_list = os.path.join(config["data_store_dir"], "association_testing/{date_of_assoc_test}_true/inputs/{bin_length}.rds"),
        phenotypes_file = config["phenotypes_file"],
        source_file = "workflow/scripts/run_gwas_source.R"
    output:
        os.path.join(config["data_store_dir"], "association_testing/{date_of_assoc_test}_true/results/{target_phenotype}/{bin_length}.rds"),
    log:
        os.path.join(config["working_dir"], "logs/run_gwas/{date_of_assoc_test}/{target_phenotype}/{bin_length}.log")
    params:
        bin_length = lambda wildcards: wildcards.bin_length,
        target_phenotype = "{target_phenotype}"
    resources:
        mem_mb = 20000
    container:
        config["R"]
    script:
        "../scripts/run_gwas.R"    
