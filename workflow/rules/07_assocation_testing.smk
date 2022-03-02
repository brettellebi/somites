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
    resources:
        mem_mb = 2000
    container:
        config["PhenotypeSimulator"]
    script:
        "../scripts/simulate_phenotypes.R"    

rule test_gwls:
    input:
        gt_pos_list = rules.create_gwas_input.output,
        phenotypes_file = rules.simulate_phenotypes.output.sim_phenos,
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

def get_mem_mb(wildcards, attempt):
    return attempt * 50000

rule run_gwls:
    input:
        gt_pos_list = rules.create_gwas_input.output,
        phenotypes_file = config["phenotypes_file"],
        source_file = "workflow/scripts/run_gwls_source.R"
    output:
        os.path.join(config["working_dir"], "association_testing/{date_of_assoc_test}/{site_filter}/true_results/{target_phenotype}/{covariates}/{inverse_norm}/{bin_length}.rds"),
    log:
        os.path.join(config["working_dir"], "logs/run_gwls/{date_of_assoc_test}/{site_filter}/{target_phenotype}/{covariates}/{inverse_norm}/{bin_length}.log")
    params:
        target_phenotype = "{target_phenotype}",
        covariates = "{covariates}",
        inverse_norm = "{inverse_norm}",
        bin_length = "{bin_length}"
    resources:
        mem_mb = get_mem_mb
    container:
        config["R"]
    script:
        "../scripts/run_gwls.R"

rule create_permuted_phenotypes:
    input:
        config["phenotypes_file"]
    output:
        os.path.join(config["working_dir"], "permuted_phenos/{date_of_assoc_test}/{target_phenotype}/{permutation_seed}.xlsx")
    log:
        os.path.join(config["working_dir"], "logs/create_permuted_phenotypes/{date_of_assoc_test}/{target_phenotype}/{permutation_seed}.log")
    params:
        target_phenotype = "{target_phenotype}",
        permutation_seed = "{permutation_seed}"
    resources:
        mem_mb = 5000
    container:
        config["R"]
    script:
        "../scripts/create_permuted_phenotypes.R"

rule run_permutations:
    input:
        gt_pos_list = rules.create_gwas_input.output,
        phenotypes_file = rules.create_permuted_phenotypes.output,
        source_file = "workflow/scripts/run_gwls_source.R"
    output:
        os.path.join(config["working_dir"], "association_testing/{date_of_assoc_test}/{site_filter}/permutations/{target_phenotype}/{covariates}/{inverse_norm}/{bin_length}/{permutation_seed}.rds"),
    log:
        os.path.join(config["working_dir"], "logs/run_permutations/{date_of_assoc_test}/{site_filter}/{target_phenotype}/{covariates}/{inverse_norm}/{bin_length}/{permutation_seed}.log")
    params:
        target_phenotype = "{target_phenotype}",
        covariates = "{covariates}",
        inverse_norm = "{inverse_norm}",
        bin_length = "{bin_length}"
    resources:
        mem_mb = 10000
    container:
        config["R"]
    script:
        "../scripts/run_gwls.R"

rule get_permutations_min_p:
    input:
        expand(os.path.join(
                config["working_dir"],
                "association_testing/{date_of_assoc_test}/{site_filter}/permutations/{target_phenotype}/{bin_length}/{permutation_seed}.rds"
                ),
            date_of_assoc_test = config["date_of_assoc_test"],
            site_filter = config["site_filter"],
            target_phenotype = config["target_phenotypes"],
            bin_length = config["bin_length"],
            permutation_seed = PERM_SEEDS,
        ),
    output:
        csv = "data/{date_of_assoc_test}_permutation_mins.csv",
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_permutations_min_p/{date_of_assoc_test}/get_permutations_min_p.log"
        ),
    resources: 
        mem_mb = 10000
    container:
        config["R"]
    script:
        "../get_permutations_min_p.R"

rule get_manhattan:
    input:
        gwas_results = rules.run_gwls.output,
        sig_levels = rules.get_permutations_min_p.output.csv,
        source_file = "workflow/scripts/get_manhattan_source.R"
    output:
        fig = "book/plots/manhattans/{date_of_assoc_test}/{site_filter}/{target_phenotype}/{covariates}/{inverse_norm}/{bin_length}.png"
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_manhattan/{date_of_assoc_test}/{site_filter}/{target_phenotype}/{covariates}/{inverse_norm}/{bin_length}.log"
        ),
    params:
        site_filter = "{site_filter}",
        target_phenotype = "{target_phenotype}",
        covariates = "{covariates}",
        inverse_norm = "{inverse_norm}",
        bin_length = "{bin_length}",
    container:
        config["R"]
    resources:
        mem_mb = 5000
    script:
        "../workflow/scripts/get_manhattan.R"
                
