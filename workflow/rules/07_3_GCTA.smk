# Convert processed recominbation blocks to Plink format
rule create_ped:
    input:
        genos = rules.process_rc_blocks.output,
        phenos = config["phenotypes_file"]
    output:
        ped = os.path.join(
            config["working_dir"],
            "peds/F2/{site_filter}/{bin_length}/{phenotype}.ped"
        ),
        map = os.path.join(
            config["working_dir"],
            "peds/F2/{site_filter}/{bin_length}/{phenotype}.map"
        ),
        phen = os.path.join(
            config["working_dir"],
            "peds/F2/{site_filter}/{bin_length}/{phenotype}.phen"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/create_ped/{site_filter}/{bin_length}/{phenotype}.log"
        ),
    params:
        phenotype = "{phenotype}"
    resources:
        mem_mb = 30000,
    container:
        # requires tidyr >= v1.2
        config["tidyverse_4.1.3"]
    script:
        "../scripts/create_ped.R"

rule create_bed:
    input:
        rules.create_ped.output.ped
    output:
        bed = os.path.join(
            config["working_dir"],
            "beds/F2/{site_filter}/{bin_length}/{phenotype}.bed"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/create_bed/{site_filter}/{bin_length}/{phenotype}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input[0].replace(".ped", ""),
        out_pref = lambda wildcards, output: output.bed.replace(".bed", ""),
    resources:
        mem_mb = 10000
    container:
        config["plink1.9"]
    shell:
        """
        plink1.9 \
            --make-bed \
            --file {params.in_pref} \
            --out {params.out_pref} \
                2> {log}
        """


rule run_mlma_loco:
    input:
        bed = rules.create_bed.output.bed,
        phen = rules.create_ped.output.phen,
    output:
        os.path.join(
            config["working_dir"],
            "gcta/mlma_loco/{site_filter}/{bin_length}/{phenotype}.loco.mlma"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/run_mlma_loco/{site_filter}/{bin_length}/{phenotype}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input.bed.replace(".bed", ""),
        out_pref = lambda wildcards, output: output[0].replace(".loco.mlma", ""),
    resources:
        mem_mb = 50000
    container:
        config["GCTA"]
    shell:
        """
        gcta64 \
            --mlma-loco \
            --bfile {params.in_pref} \
            --pheno {input.phen} \
            --out {params.out_pref} \
            --thread-num 10 \
                > {log}
        """
