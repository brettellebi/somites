rule get_annotations_invnorm:
    input:
        res = rules.run_mlma_loco_invnorm.output,
        min_p = rules.get_min_p_perms_invnorm.output,
        sites = rules.extract_homo_div_snps.output.sites,
    output:
        bins_all = "results/annotations_invnorm/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/bin_overlaps_all.csv",
        bins_unq = "results/annotations_invnorm/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/bin_overlaps_unique.csv",
        snps = "results/annotations_invnorm/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/snps.csv"
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_annotations_invnorm/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}.log"
        ),
    params:
        bin_length = "{bin_length}"
    resources:
        mem_mb = 2000
    container:
        config["R_4.1.3"]
    script:
        "../scripts/get_annotations_gcta.R"

rule get_annotations_psm:
    input:
        res = rules.run_mlma_loco.output,
        min_p = rules.get_min_p_perms.output,
        sites = rules.extract_homo_div_snps.output.sites,
    output:
        bins_all = "results/annotations_psm/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{covars}/bin_overlaps_all.csv",
        bins_unq = "results/annotations_psm/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{covars}/bin_overlaps_unique.csv",
        snps = "results/annotations_psm/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{covars}/snps.csv",
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_annotations_psm/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}/{covars}.log"
        ),
    params:
        bin_length = "{bin_length}"
    resources:
        mem_mb = 2000
    container:
        config["R_4.1.3"]
    script:
        "../scripts/get_annotations_gcta.R"
