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

# Create functions to assign different memory and queue requirements for different
# types of GWAS (no covariates require little memory; two covariates require a lot)
def get_mem_mb(wildcards, attempt):
    if wildcards.covariates == "None":
        multiplier = 4000
    elif wildcards.covariates == "Microscope":
        multiplier = 6000
    elif wildcards.covariates == "Microscope-Date_of_imaging":
        multiplier = 80000
    return attempt * multiplier

def get_queue(wildcards, attempt):
    if wildcards.covariates == "None":
        queue = "short"
    elif wildcards.covariates == "Microscope":
        queue = "short"
    elif wildcards.covariates == "Microscope-Date_of_imaging":
        queue = "short"
    return queue

rule run_gwls:
    input:
        gt_pos_list = rules.create_gwas_input.output,
        phenotypes_file = config["phenotypes_file"],
    output:
        os.path.join(config["working_dir"], "association_testing/{date_of_assoc_test}/{site_filter}/true_results/{target_phenotype}/{covariates}/{inverse_norm}/{bin_length}.rds"),
    log:
        os.path.join(config["working_dir"], "logs/run_gwls/{date_of_assoc_test}/{site_filter}/{target_phenotype}/{covariates}/{inverse_norm}/{bin_length}.log")
    params:
        target_phenotype = "{target_phenotype}",
        covariates = "{covariates}",
        inverse_norm = "{inverse_norm}",
        bin_length = "{bin_length}",
        source_file = "workflow/scripts/run_gwls_source.R"
    resources:
        mem_mb = get_mem_mb,
        queue = get_queue
    container:
        config["R"]
    script:
        "../scripts/run_gwls.R"

rule create_permuted_phenotypes:
    input:
        config["phenotypes_file"]
    output:
        os.path.join(config["working_dir"], "permuted_phenos/{date_of_assoc_test}/{permutation_seed}.xlsx")
    log:
        os.path.join(config["working_dir"], "logs/create_permuted_phenotypes/{date_of_assoc_test}/{permutation_seed}.log")
    params:
        permutation_seed = "{permutation_seed}",
    resources:
        mem_mb = 100
    container:
        config["R"]
    script:
        "../scripts/create_permuted_phenotypes.R"

rule run_permutations:
    input:
        gt_pos_list = rules.create_gwas_input.output,
        phenotypes_file = rules.create_permuted_phenotypes.output,
    output:
        os.path.join(
            config["working_dir"],
            "association_testing/{date_of_assoc_test}/{site_filter}/permutations/{target_phenotype}/{covariates}/{inverse_norm}/{bin_length}/{permutation_seed}.rds"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/run_permutations/{date_of_assoc_test}/{site_filter}/{target_phenotype}/{covariates}/{inverse_norm}/{bin_length}/{permutation_seed}.log"
        ),
    params:
        target_phenotype = "{target_phenotype}",
        covariates = "{covariates}",
        inverse_norm = "{inverse_norm}",
        bin_length = "{bin_length}",
        source_file = "workflow/scripts/run_gwls_source.R"
    resources:
        mem_mb = get_mem_mb,
        queue = get_queue
    container:
        config["R"]
    script:
        "../scripts/run_gwls.R"

rule get_manhattan:
    input:
        gwas_results = rules.run_gwls.output,
        perm_results = expand(os.path.join(
            config["working_dir"],
            "association_testing/{{date_of_assoc_test}}/{{site_filter}}/permutations/{{target_phenotype}}/{{covariates}}/{{inverse_norm}}/{{bin_length}}/{permutation_seed}.rds"),
                permutation_seed = PERM_SEEDS
        ),
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
        mem_mb = 2500
    script:
        "../scripts/get_manhattan.R"
                
