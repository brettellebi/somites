rule split_pheno_by_microscope:
    input:
        phen = rules.create_phen.output,
        covar = expand(rules.create_covar.output,
            covars = "Microscope"
        ),
    output:
        phen = os.path.join(
            config["working_dir"],
            "phens/split_microscope/{phenotype}/{microscope}.phen"
        ),
        ids = os.path.join(
            config["working_dir"],
            "phens/split_microscope/{phenotype}/{microscope}.list"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/split_pheno_by_microscope/{phenotype}/{microscope}.log"
        ),
    params:
        microscope = "{microscope}"
    resources:
        mem_mb = 300,
    container:
        config["tidyverse_4.1.3"]
    script:
        "../scripts/split_pheno_by_microscope.R"

rule split_bed:
    input:
        ped = rules.create_ped.output.ped,
        ids = rules.split_pheno_by_microscope.output.ids,
    output:
        bed = os.path.join(
            config["working_dir"],
            "beds/split_microscope/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{microscope}.bed"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/split_bed/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{microscope}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input.ped.replace(".ped", ""),
        out_pref = lambda wildcards, output: output.bed.replace(".bed", ""),
    resources:
        mem_mb = 200
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
            --chr-set 24 no-xy no-mt \
            --keep {input.ids} \
            --file {params.in_pref} \
            --out {params.out_pref} \
                2> {log}
        """

rule permute_split_microscope_phen:
    input:
        rules.split_pheno_by_microscope.output.phen,
    output:
        os.path.join(
            config["working_dir"],
            "phens/split_microscope_permuted/{phenotype}/{microscope}/{seed}.phen"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/permute_split_microscope_phen/{phenotype}/{microscope}/{seed}.log"
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

rule run_mlma_loco_split_microscope:
    input:
        bed = rules.split_bed.output.bed,
        phen = rules.split_pheno_by_microscope.output.phen,
        excl_samples = rules.create_excluded_samples_list.output,
    output:
        os.path.join(
            config["working_dir"],
            "gcta/mlma_loco_split_microscope/true/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{microscope}.loco.mlma"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/run_mlma_loco_split_microscope/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{microscope}.log"
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

rule run_mlma_loco_split_microscope_permuted:
    input:
        bed = rules.split_bed.output.bed,
        phen = rules.split_pheno_by_microscope.output.phen,
        excl_samples = rules.create_excluded_samples_list.output,
    output:
        os.path.join(
            config["working_dir"],
            "gcta/mlma_loco_split_microscope/permuted/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{microscope}/{seed}.loco.mlma"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/run_mlma_loco_split_microscope_permuted/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{microscope}/{seed}.log"
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

rule get_min_p_perms_split_microscope:
    input:
        expand(os.path.join(
            config["working_dir"],
            "gcta/mlma_loco_split_microscope/permuted/{{ref}}/{{max_reads}}/{{bin_length}}/{{cov}}/{{phenotype}}/{{microscope}}/{seed}.loco.mlma"
            ),
                seed = PERM_SEEDS         
        ),
    output:
        os.path.join(
            config["working_dir"],
            "gcta/mlma_loco_split_microscope/min_p/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{microscope}.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_min_p_perms_split_microscope/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{microscope}.log"
        ),
    resources:
        mem_mb = 200
    container:
        config["tidyverse_4.1.3"]
    script:
        "../scripts/get_min_p_perms.R"

rule get_manhattan_gcta_split_microscope:
    input:
        res = rules.run_mlma_loco_split_microscope.output,
        min_p = rules.get_min_p_perms_split_microscope.output,
    output:
        man = "book/plots/gcta_split_microscope/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{microscope}.png"
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_manhattan_gcta_split_microscope/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{microscope}.log"
        ),    
    params:
        max_reads = "{max_reads}",
        bin_length = "{bin_length}",
        cov = "{cov}",
        phenotype = "{phenotype}",
        microscope = "{microscope}"
    resources:
        mem_mb = 1000
    container:
        config["R_4.1.3"]
    script:
        "../scripts/get_manhattan_gcta_split_microscope.R"
