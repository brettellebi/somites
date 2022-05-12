rule top_and_bottom_karyos:
    input:
        genos = rules.true_hmmlearn.output.csv,
        phenos = rules.split_inv_norm_pheno.output,
    output:
        prop_sites_total_top = "book/plots/{ref}/F1_het_min_DP/top_and_bottom_karyos/{max_reads}/{bin_length}/{cov}/{phenotype}/prop_sites_total_top.png",
        prop_sites_total_bottom = "book/plots/{ref}/F1_het_min_DP/top_and_bottom_karyos/{max_reads}/{bin_length}/{cov}/{phenotype}/prop_sites_total_bottom.png",
        karyoplot_no_missing_top = "book/plots/{ref}/F1_het_min_DP/top_and_bottom_karyos/{max_reads}/{bin_length}/{cov}/{phenotype}/karyoplot_no_missing_top.png",
        karyoplot_no_missing_bottom = "book/plots/{ref}/F1_het_min_DP/top_and_bottom_karyos/{max_reads}/{bin_length}/{cov}/{phenotype}/karyoplot_no_missing_bottom.png",
        karyoplot_with_missing_top = "book/plots/{ref}/F1_het_min_DP/top_and_bottom_karyos/{max_reads}/{bin_length}/{cov}/{phenotype}/karyoplot_wi_missing_top.png",
        karyoplot_with_missing_bottom = "book/plots/{ref}/F1_het_min_DP/top_and_bottom_karyos/{max_reads}/{bin_length}/{cov}/{phenotype}/karyoplot_wi_missing_bottom.png",
    log:
        os.path.join(
            config["working_dir"],
            "logs/top_and_bottom_karyos/F1_het_min_DP/{ref}/{max_reads}/{bin_length}/{cov}/{phenotype}.log"
        ),
    params:
        cov = "{cov}",
        bin_length = "{bin_length}",
        max_reads = "{max_reads}"
    resources:
        mem_mb = 50000
    container:
        config["R_4.1.3"]
    script:
        "../scripts/top_and_bottom_karyos.R"

rule kaga_F0_chr3:
    input:
        counts = rules.trio_gt_counts_in_bins.output,
        res = rules.run_mlma_loco_invnorm.output,
        min_p = rules.get_min_p_perms_invnorm.output,
    output:
        fig = "book/plots/{ref}/sig_region_zoom/{max_reads}/{cov}/{phenotype}/{sample}_{bin_length}_chr3.png"
    log:
        os.path.join(
            config["working_dir"],
            "logs/kaga_F0_chr3/{ref}_{max_reads}_{cov}_{phenotype}_{sample}_{bin_length}.log"
        ),
    params:
        bin_length = "{bin_length}"
    resources:
        mem_mb = 20000
    container:
        config["R_4.1.3"]
    script:
        "../scripts/kaga_F0_chr3.R"

rule F2_chr3_pheno_by_gt:
    input:
        genos = rules.true_hmmlearn.output.csv,
        phenos = rules.split_inv_norm_pheno.output,
        res = rules.run_mlma_loco_invnorm.output,    
        min_p = rules.get_min_p_perms_invnorm.output,
    output:
        fig = "book/plots/{ref}/pheno_by_geno/{max_reads}/{cov}/{phenotype}/{bin_length}_chr3.png",
    log:
        os.path.join(
            config["working_dir"],
            "logs/F2_chr3_pheno_by_gt/{ref}_{max_reads}_{cov}_{phenotype}_{bin_length}.log"
        ),
    params:
        bin_length = "{bin_length}"
    resources:
        mem_mb = 20000
    container:
        config["R_4.1.3"]
    script:
        "../scripts/F2_chr3_pheno_by_gt.R"        