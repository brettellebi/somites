# Send plots to be used for the paper to the Google Drive
rule send_plots_to_google_drive:
    input:
        circos = expand(
            rules.circos_homozygosity.output.plot,
                ref = "hdrr",
                bin_length = 5000,
                sample = ["Cab", "Kaga", "F1"]
        ),
        snp_counts_chr = expand(
            rules.plot_SNP_counts_per_chr.output,
                ref = "hdrr"
        ),
        int_split_invnorm = expand(
            rules.get_manhattan_gcta_invnorm.output,
                ref = "hdrr",
                max_reads = "None",
                bin_length = "5000",
                cov = config["hmm_covariance"],
                phenotype = "intercept"
        ),
    output:
        touch(
            os.path.join(
                config["working_dir"],
                "logs/send_plots_to_google_drive/all.done"
            )
        )
    log:
        os.path.join(
            config["working_dir"],
            "logs/send_plots_to_google_drive/all.log"
        ),
    params:
        drive_dir = config["google_drive_dir"]
    conda:
        "../envs/rclone.yaml"
    resources:
        mem_mb = 1000
    shell:
        """
        for i in $(echo {input}); do \
            rclone copy $i {params.drive_dir}/ ; \
        done
        """
