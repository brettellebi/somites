rule run_LB_impute:
    input:
        os.path.join(config["data_store_dir"], "vcfs/F0_and_F2/final/all.vcf.gz"),
    output:
        os.path.join(config["working_dir"], "vcfs/LB-Impute/out.vcf.gz"),
    