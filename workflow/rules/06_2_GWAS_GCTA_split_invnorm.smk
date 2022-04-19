rule split_inv_norm_pheno:
    input:
        phen = rules.create_phen.output,
        covar = expand(rules.create_covar.output,
            covars = "Microscope"
        ),
    output:
        os.path.join(
            config["working_dir"],
            "phens/inv_norm/{phenotype}.phen"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/split_inv_norm_pheno/{phenotype}.log"
        ),
    resources:
        mem_mb = 300,
    container:
        config["tidyverse_4.1.3"]
    script:
        "../scripts/split_inv_norm_pheno.R"

rule permute_invnorm_phen:
    input:
        rules.split_inv_norm_pheno.output,
    output:
        os.path.join(
            config["working_dir"],
            "phens/inv_norm_permuted/{phenotype}/{seed}.phen"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/permute_invnorm_phen/{phenotype}/{seed}.log"
        ),
    params:
        phenotype = "{phenotype}",
        seed = "{seed}"
    resources:
        mem_mb = 500,
    container:
        config["tidyverse_4.1.3"]
    script:
        "../scripts/permute_phen.R"

rule run_mlma_loco_invnorm:
    input:
        bed = rules.create_bed.output.bed,
        phen = rules.split_inv_norm_pheno.output,
        excl_samples = rules.create_excluded_samples_list.output,
    output:
        os.path.join(
            config["working_dir"],
            "gcta/mlma_loco_invnorm/true/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}.loco.mlma"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/run_mlma_loco_invnorm/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input.bed.replace(".bed", ""),
        out_pref = lambda wildcards, output: output[0].replace(".loco.mlma", "")
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
            --out {params.out_pref} \
            --remove {input.excl_samples} \
            --autosome-num 24 \
            --thread-num {resources.threads} \
                2> {log}
        """

rule run_mlma_loco_invnorm_permuted:
    input:
        bed = rules.create_bed.output.bed,
        phen = rules.permute_invnorm_phen.output,
        excl_samples = rules.create_excluded_samples_list.output,
    output:
        os.path.join(
            config["working_dir"],
            "gcta/mlma_loco_invnorm/permuted/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{seed}.loco.mlma"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/run_mlma_loco_invnorm_permuted/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{seed}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input.bed.replace(".bed", ""),
        out_pref = lambda wildcards, output: output[0].replace(".loco.mlma", ""),
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
            --out {params.out_pref} \
            --remove {input.excl_samples} \
            --autosome-num 24 \
            --thread-num {resources.threads} \
                2> {log}
        """

rule get_min_p_perms_invnorm:
    input:
        expand(os.path.join(
            config["working_dir"],
            "gcta/mlma_loco_invnorm/permuted/{{ref}}/{{max_reads}}/{{bin_length}}/{{cov}}/{{phenotype}}/{seed}.loco.mlma"
            ),
                seed = PERM_SEEDS         
        ),
    output:
        os.path.join(
            config["working_dir"],
            "gcta/mlma_loco_invnorm/min_p/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_min_p_perms/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}.log"
        ),
    resources:
        mem_mb = 200
    container:
        config["tidyverse_4.1.3"]
    script:
        "../scripts/get_min_p_perms.R"

rule get_manhattan_gcta_invnorm:
    input:
        res = rules.run_mlma_loco_invnorm.output,
        min_p = rules.get_min_p_perms_invnorm.output,
    output:
        man = "book/plots/gcta_invnorm/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}.png"
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_manhattan_gcta_invnorm/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}.log"
        ),    
    params:
        max_reads = "{max_reads}",
        bin_length = "{bin_length}",
        cov = "{cov}",
        phenotype = "{phenotype}",
    resources:
        mem_mb = 1000
    container:
        config["R_4.1.3"]
    script:
        "../scripts/get_manhattan_gcta_invnorm.R"
