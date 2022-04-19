# Convert processed recominbation blocks to Plink format
# NOTE: this step imputes missing genotypes
rule create_ped:
    input:
        genos = rules.true_hmmlearn.output.csv
    output:
        ped = os.path.join(
            config["working_dir"],
            "peds/F2/{ref}/{max_reads}/{bin_length}/{cov}.ped"
        ),
        map = os.path.join(
            config["working_dir"],
            "peds/F2/{ref}/{max_reads}/{bin_length}/{cov}.map"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/create_ped/{ref}/{max_reads}/{bin_length}/{cov}.log"
        ),
    resources:
        mem_mb = 30000,
    container:
        # requires tidyr >= v1.2
        config["tidyverse_4.1.3"]
    script:
        "../scripts/create_ped.R"

# Convert .ped to .bed
rule create_bed:
    input:
        rules.create_ped.output.ped
    output:
        bed = os.path.join(
            config["working_dir"],
            "beds/F2/{ref}/{max_reads}/{bin_length}/{cov}.bed"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/create_bed/{ref}/{max_reads}/{bin_length}/{cov}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input[0].replace(".ped", ""),
        out_pref = lambda wildcards, output: output.bed.replace(".bed", ""),
    resources:
        mem_mb = 100
    container:
        config["plink1.9"]
    shell:
        """
        plink1.9 \
            --make-bed \
            --no-fid \
            --no-parents \
            --no-sex \
            --no-pheno \
            --file {params.in_pref} \
            --out {params.out_pref} \
                2> {log}
        """

# Create .phen files
# as specified here: https://gcta.freeforums.net/thread/247/greml-estimating-variance-explained-snps
rule create_phen:
    input:
        ped = os.path.join(
            config["working_dir"],
            "peds/F2/hdrr/None/5000/1/intercept.ped"
        ),
        phenos = config["phenotypes_file"]
    output:
        os.path.join(
            config["working_dir"],
            "phens/true/{phenotype}.phen"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/create_phen/{phenotype}.log"
        ),
    params:
        phenotype = "{phenotype}"
    resources:
        mem_mb = 2000,
    container:
        # requires tidyr >= v1.2
        config["tidyverse_4.1.3"]
    script:
        "../scripts/create_phen.R"

rule permute_phen:
    input:
        rules.create_phen.output,
    output:
        os.path.join(
            config["working_dir"],
            "phens/permuted/{phenotype}/{seed}.phen"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/permute_phen/{phenotype}/{seed}.log"
        ),
    params:
        phenotype = "{phenotype}",
        seed = "{seed}"
    resources:
        mem_mb = 500,
    container:
        # requires tidyr >= v1.2
        config["tidyverse_4.1.3"]
    script:
        "../scripts/permute_phen.R"

# Set rule order for creating covariate files and running mlma-loco
# Because the covariate files aren't in the mlma rules' inputs
ruleorder: create_covar > permute_covar > run_mlma_loco > run_mlma_loco_permuted

# Create covariate files
# as specified here: https://gcta.freeforums.net/thread/247/greml-estimating-variance-explained-snps
rule create_covar:
    input:
        #genos = rules.process_rc_blocks.output,
        ped = os.path.join(
            config["working_dir"],
            "peds/F2/hdrr/None/5000/1/intercept.ped"
        ),
        phenos = config["phenotypes_file"]
    output:
        os.path.join(
            config["working_dir"],
            "covars/true/{covars}.covar"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/create_covar/{covars}.log"
        ),
    params:
        covars = "{covars}"
    resources:
        mem_mb = 1000,
    container:
        # requires tidyr >= v1.2
        config["tidyverse_4.1.3"]
    script:
        "../scripts/create_covar.R"

rule permute_covar:
    input:
        rules.create_covar.output,
    output:
        os.path.join(
            config["working_dir"],
            "covars/permuted/{covars}/{seed}.covar"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/permute_covar/{covars}/{seed}.log"
        ),
    params:
        covars = "{covars}",
        seed = "{seed}"
    resources:
        mem_mb = 500,
    container:
        # requires tidyr >= v1.2
        config["tidyverse_4.1.3"]
    script:
        "../scripts/permute_covars.R"

rule create_excluded_samples_list:
    output:
        "config/low_cov_samples.list"
    log:
        os.path.join(
            config["working_dir"],
            "logs/create_excluded_samples_list/all.log"
        ),
    resources:
        mem_mb = 100 
    run:
        file = open(output[0], 'w')
        for item in config["low_cov_samples"]:
            file.writelines(str(item)+'\t'+str(item)+'\n')
        file.close()

def set_covars(wildcards):
    if wildcards.covars == "None":
        out = ""
    else:
        covars_file = os.path.join(config["working_dir"], "covars/true/" + wildcards.covars + ".covar")
        out = '--covar ' + covars_file
    return(out)

rule run_mlma_loco:
    input:
        bed = rules.create_bed.output.bed,
        phen = rules.create_phen.output,
        excl_samples = rules.create_excluded_samples_list.output,
    output:
        os.path.join(
            config["working_dir"],
            "gcta/mlma_loco/true/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{covars}.loco.mlma"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/run_mlma_loco/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{covars}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input.bed.replace(".bed", ""),
        out_pref = lambda wildcards, output: output[0].replace(".loco.mlma", ""),
        covars = set_covars,
    resources:
        mem_mb = 500,
        threads = 1
    container:
        config["GCTA"]
    shell:
        """
        gcta64 \
            --mlma-loco \
            --bfile {params.in_pref} \
            --pheno {input.phen} \
            {params.covars} \
            --out {params.out_pref} \
            --remove {input.excl_samples} \
            --autosome-num 24 \
            --thread-num {resources.threads} \
                2> {log}
        """

def set_covars_permuted(wildcards):
    if wildcards.covars == "None":
        out = ""
    else:
        covars_file = os.path.join(config["working_dir"], "covars/permuted/" + wildcards.covars + "/" + wildcards.seed + ".covar")
        out = '--covar ' + covars_file
    return(out)

rule run_mlma_loco_permuted:
    input:
        bed = rules.create_bed.output.bed,
        phen = rules.permute_phen.output,
        excl_samples = rules.create_excluded_samples_list.output,
    output:
        os.path.join(
            config["working_dir"],
            "gcta/mlma_loco/permuted/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{covars}/{seed}.loco.mlma"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/run_mlma_loco_permuted/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{covars}/{seed}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input.bed.replace(".bed", ""),
        out_pref = lambda wildcards, output: output[0].replace(".loco.mlma", ""),
        covars = set_covars_permuted,
    resources:
        mem_mb = 500,
        threads = 1
    container:
        config["GCTA"]
    shell:
        """
        gcta64 \
            --mlma-loco \
            --bfile {params.in_pref} \
            --pheno {input.phen} \
            {params.covars} \
            --out {params.out_pref} \
            --remove {input.excl_samples} \
            --autosome-num 24 \
            --thread-num {resources.threads} \
                2> {log}
        """

rule get_min_p_perms:
    input:
        expand(os.path.join(
            config["working_dir"],
            "gcta/mlma_loco/permuted/{{ref}}/{{max_reads}}/{{bin_length}}/{{cov}}/{{phenotype}}/{{covars}}/{seed}.loco.mlma"
            ),
                seed = PERM_SEEDS         
        ),
    output:
        os.path.join(
            config["working_dir"],
            "gcta/mlma_loco/min_p/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{covars}.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_min_p_perms/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{covars}.log"
        ),
    resources:
        mem_mb = 200
    container:
        config["tidyverse_4.1.3"]
    script:
        "../scripts/get_min_p_perms.R"

rule get_manhattan_gcta:
    input:
        res = rules.run_mlma_loco.output,
        min_p = rules.get_min_p_perms.output,
    output:
        man = "book/plots/gcta/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{covars}.png"
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_manhattan_gcta/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{covars}.log"
        ),    
    params:
        max_reads = "{max_reads}",
        bin_length = "{bin_length}",
        cov = "{cov}",
        phenotype = "{phenotype}",
        covars = "{covars}"
    resources:
        mem_mb = 1000
    container:
        config["R_4.1.3"]
    script:
        "../scripts/get_manhattan_gcta.R"

    