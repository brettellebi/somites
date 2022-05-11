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
        genos = rules.extract_trio_genos.output,
        res = rules.run_mlma_loco_invnorm.output,
        min_p = rules.get_min_p_perms_invnorm.output,
    output:
        fig = "book/plots/{ref}/F1_het_min_DP/Kaga_F0_chr3.png"
    log:
        os.path.join(
            config["working_dir"],
            "logs/kaga_F0_chr3/F1_het_min_DP/{ref}.log"
        ),
    resources:
        mem_mb = 50000
    container:
        config["R_4.1.3"]
    script:
        "../scripts/kaga_F0_chr3.R"       
