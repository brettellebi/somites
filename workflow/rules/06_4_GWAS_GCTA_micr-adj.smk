rule create_micr_adj_phen:
    input:
        ped = os.path.join(
            config["working_dir"],
            "peds/F2/hdrr/None/5000/1/intercept.ped"
        ),
        phenos = "config/DF_648_periodNORMALISED.csv"
    output:
        os.path.join(
            config["working_dir"],
            "phens/micr_adj/{phenotype}.phen"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/create_micr_adj_phen/{phenotype}.log"
        ),
    params:
        phenotype = "{phenotype}"
    resources:
        mem_mb = 2000,
    container:
        # requires tidyr >= v1.2
        config["tidyverse_4.1.3"]
    script:
        "../scripts/create_micr_adj_phen.R"

rule permute_micr_adj_phen:
    input:
        rules.create_micr_adj_phen.output,
    output:
        os.path.join(
            config["working_dir"],
            "phens/micr_adj_permuted/{phenotype}/{seed}.phen"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/permute_micr_adj_phen/{phenotype}/{seed}.log"
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

rule run_mlma_loco_micr_adj:
    input:
        bed = rules.create_bed.output.bed,
        phen = rules.create_micr_adj_phen.output,
        excl_samples = rules.create_excluded_samples_list.output,
    output:
        os.path.join(
            config["working_dir"],
            "gcta/mlma_loco_micr_adj/true/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}.loco.mlma"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/run_mlma_loco_micr_adj/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}.log"
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

rule run_mlma_loco_micr_adj_permuted:
    input:
        bed = rules.create_bed.output.bed,
        phen = rules.permute_micr_adj_phen.output,
        excl_samples = rules.create_excluded_samples_list.output,
    output:
        os.path.join(
            config["working_dir"],
            "gcta/mlma_loco_micr_adj/permuted/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{seed}.loco.mlma"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/run_mlma_loco_micr_adj_permuted/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{seed}.log"
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

rule get_min_p_perms_micr_adj:
    input:
        expand(os.path.join(
            config["working_dir"],
            "gcta/mlma_loco_micr_adj/permuted/{{ref}}/{{max_reads}}/{{bin_length}}/{{cov}}/{{phenotype}}/{seed}.loco.mlma"
            ),
                seed = PERM_SEEDS         
        ),
    output:
        os.path.join(
            config["working_dir"],
            "gcta/mlma_loco_micr_adj/min_p/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_min_p_perms_micr_adj/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}.log"
        ),
    resources:
        mem_mb = 500
    container:
        config["tidyverse_4.1.3"]
    script:
        "../scripts/get_min_p_perms.R"

rule get_manhattan_gcta_micr_adj:
    input:
        res = rules.run_mlma_loco_micr_adj.output,
        min_p = rules.get_min_p_perms_micr_adj.output,
    output:
        man = "book/plots/gcta_micr_adj/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}.png"
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_manhattan_gcta_micr_adj/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}.log"
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
        "../scripts/get_manhattan_gcta.R"
