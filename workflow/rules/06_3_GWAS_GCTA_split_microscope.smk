rule split_pheno_by_microscope:
    input:
        phen = rules.create_phen.output,
        covar = expand(rules.create_covar.output,
            covars = "Microscope"
        ),
    output:
        os.path.join(
            config["working_dir"],
            "phens/split_microscope/{phenotype}/{microscope}.phen"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/split_pheno_by_microscope/{phenotype}/{microscope}.log"
        ),
    params:
        microsope = "{microscope}"
    resources:
        mem_mb = 300,
    container:
        config["tidyverse_4.1.3"]
    script:
        "../scripts/split_pheno_by_microscope.R"
