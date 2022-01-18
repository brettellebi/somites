# Process haplotype blocks data and create GWAS input

# Note: this step filters out the samples with low coverage (config["low_cov_samples"])
rule create_gwas_input:
    input:
        genotypes = os.path.join(config["data_store_dir"], "recombination_blocks/F2/{site_filter}/{bin_length}.txt"),
    output:
        os.path.join(config["data_store_dir"], "association_testing/{date_of_assoc_test}/{site_filter}/inputs/{bin_length}.rds"),
    log:
        os.path.join(config["working_dir"], "logs/create_gwas_input/{date_of_assoc_test}/{site_filter}/{bin_length}.log")
    params:
        bin_length = "{bin_length}",
        low_cov_samples = config["low_cov_samples"]
    resources:
        mem_mb = 50000
    container:
        config["R"]
    script:
        "../scripts/create_gwas_input.R"

rule simulate_phenotypes:
    input:
        os.path.join(config["data_store_dir"], "association_testing/{date_of_assoc_test}_true/{site_filter}/inputs/{bin_length}.rds"),
    output:
        sample_genos = os.path.join(config["data_store_dir"], "association_testing/{date_of_assoc_test}/{site_filter}/sample_genos/{bin_length}.csv"),
        sim_phenos = os.path.join(config["data_store_dir"], "association_testing/{date_of_assoc_test}/{site_filter}/sim_phenos/{bin_length}.xlsx"),
    log:
        os.path.join(config["working_dir"], "logs/simulate_phenotypes/{date_of_assoc_test}/{site_filter}/{bin_length}.log"),
    params:
        n_sample_gts = config["n_sample_gts"]
    container:
        config["PhenotypeSimulator"]
    script:
        "../scripts/simulate_phenotypes.R"    

rule test_gwls:
    input:
        gt_pos_list = os.path.join(config["data_store_dir"], "association_testing/{date_of_assoc_test}/{site_filter}/inputs/{bin_length}.rds"),
        phenotypes_file = os.path.join(config["data_store_dir"], "association_testing/{date_of_assoc_test}/{site_filter}/sim_phenos/{bin_length}.xlsx"),
        source_file = "workflow/scripts/run_gwls_source.R"
    output:
        os.path.join(config["data_store_dir"], "association_testing/{date_of_assoc_test}/{site_filter}/test_results/{bin_length}.rds"),
    log:
        os.path.join(config["working_dir"], "logs/test_gwls/{date_of_assoc_test}/{site_filter}/{bin_length}.log")
    params:
        bin_length = lambda wildcards: wildcards.bin_length,
        target_phenotype = "Y"
    resources:
        mem_mb = 20000
    container:
        config["R"]
    script:
        "../scripts/run_gwls.R"   

rule run_gwls:
    input:
        gt_pos_list = os.path.join(config["data_store_dir"], "association_testing/{date_of_assoc_test}/{site_filter}/inputs/{bin_length}.rds"),
        phenotypes_file = config["phenotypes_file"],
        source_file = "workflow/scripts/run_gwls_source.R"
    output:
        os.path.join(config["data_store_dir"], "association_testing/{date_of_assoc_test}/{site_filter}/true_results/{target_phenotype}/{bin_length}.rds"),
    log:
        os.path.join(config["working_dir"], "logs/run_gwls/{date_of_assoc_test}/{site_filter}/{target_phenotype}/{bin_length}.log")
    params:
        bin_length = "{bin_length}",
        target_phenotype = "{target_phenotype}"
    resources:
        mem_mb = 20000
    container:
        config["R"]
    script:
        "../scripts/run_gwls.R"

rule run_permutations:
    input:
        gt_pos_list = os.path.join(config["data_store_dir"], "association_testing/{date_of_assoc_test}/{site_filter}/inputs/{bin_length}.rds"),
        phenotypes_file = config["phenotypes_file"],
        source_file = "workflow/scripts/run_gwls_source.R"
    output:
        os.path.join(config["data_store_dir"], "association_testing/{date_of_assoc_test}/{site_filter}/permutations/{target_phenotype}/{bin_length}/{permutation_seed}.rds"),
    log:
        os.path.join(config["working_dir"], "logs/run_permutations/{date_of_assoc_test}/{site_filter}/{target_phenotype}/{bin_length}/{permutation_seed}.log")
    params:
        bin_length = "{bin_length}",
        target_phenotype = "{target_phenotype}",
        permutation_seed = "{permutation_seed}"
    resources:
        mem_mb = 10000
    container:
        config["R"]
    script:
        "../scripts/run_gwls_permutation.R"

