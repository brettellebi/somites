# Scripts from https://github.com/betharowan/TIGER_Scripts-for-distribution

# Create input files and get set up
rule call_bases:
    input:
        os.path.join(config["working_dir"], "dpABs/F2/no_repeat_reads/{F2_sample}.txt"),
    output:
        os.path.join(config["working_dir"], "TIGER/01_base_called/{F2_sample}.txt"),
    log:
        os.path.join(config["working_dir"], "logs/call_bases/{F2_sample}.log")
    params:
        basecaller_script = "workflow/scripts/TIGER/base_caller.jar",
    container:
        config["java"]
    shell:
        """
        java -jar {params.basecaller_script} \
            -r {input} \
            -o {output} \
            -n bi \
                2> {log}
        """

rule estimate_allele_freqs:
    input:
        # uses the same INPUT as `call_bases`
        rules.call_bases.input
    output:
        os.path.join(config["working_dir"], "TIGER/02_af_estimated/{F2_sample}.txt"),
    log:
        os.path.join(config["working_dir"], "logs/estimate_allele_freqs/{F2_sample}.log")
    params:
        af_estimator_script = "workflow/scripts/TIGER/allele_freq_estimator.jar",
        window_size = 1000
    container:
        config["java"]
    shell:
        """
        java -jar {params.af_estimator_script} \
            -r {input} \
            -o {output} \
            -n bi \
            -w {params.window_size} \
                2> {log}
        """    

rule apply_beta_mixture_model:
    input:
        rules.estimate_allele_freqs.output
    output:
        os.path.join(config["working_dir"], "TIGER/03_beta_modelled/{F2_sample}.txt"),
    log:
        os.path.join(config["working_dir"], "logs/apply_beta_mixture_model/{F2_sample}.log")
    params:
        beta_mixture_script = "workflow/scripts/TIGER/beta_mixture_model.R",
    container:
        config["R"]
    shell:
        """
        Rscript --vanilla {params.beta_mixture_script} \
            {input} \
            {output} \
                2> {log}
        """

rule prepare_for_HMM:
    input:
        # uses the same INPUT as `call_bases`
        counts = rules.call_bases.input,
        base_call_output = rules.call_bases.output,
        chr_sizes = config["ref_chr_lengths"]
    output:
        os.path.join(config["working_dir"], "TIGER/04_hmm_prepared/{F2_sample}.txt"),
    log:
        os.path.join(config["working_dir"], "logs/prepare_for_HMM/{F2_sample}.log"),
    params:
        prep_prob_script = "workflow/scripts/TIGER/prep_prob.pl",
        sample = "{F2_sample}"
    container:
        config["perl"]
    shell:
        """
        perl {params.prep_prob_script} \
            -s {params.sample} \
            -m {input.counts} \
            -b {input.base_call_output} \
            -c {input.chr_sizes} \
            -o {output} \
                2> {log}
        """
    
rule calculate_HMM_probs:
    input:
        allele_freqs = rules.estimate_allele_freqs.output,
        file_for_prob = rules.prepare_for_HMM.output,
        bmm_output = rules.apply_beta_mixture_model.output,
        chr_sizes = config["ref_chr_lengths"]
    output:
        os.path.join(config["working_dir"], "TIGER/05_hmm_probs/{F2_sample}.txt"),
    log:
        os.path.join(config["working_dir"], "logs/calculate_HMM_probs/{F2_sample}.log"),
    params:
        hmm_prob_script = "workflow/scripts/TIGER/hmm_prob.pl",
        sample = "{F2_sample}",
        output_prefix = os.path.join(config["working_dir"], "TIGER/05_hmm_probs/{F2_sample}")
    container:
        config["perl"]
    shell:
        """
        perl {params.hmm_prob_script} \
            -s {input.allele_freqs} \
            -p {input.file_for_prob} \
            -o {params.output_prefix} \
            -a {input.bmm_output} \
            -c {input.chr_sizes}
        """