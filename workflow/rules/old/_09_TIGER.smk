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
        model = os.path.join(config["working_dir"], "TIGER/05_hmm_probs/{F2_sample}_hmm_model"),
        breaks = os.path.join(config["working_dir"], "TIGER/05_hmm_probs/{F2_sample}_sliding_window.breaks.txt"),
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
            -c {input.chr_sizes} \
                2> {log}
        """

rule run_HMM:
    input:
        called_bases = rules.call_bases.output,
        hmm_probs = rules.calculate_HMM_probs.output.model
    output:
        os.path.join(config["working_dir"], "TIGER/06_hmm_out/{F2_sample}.txt"),
    log:
        os.path.join(config["working_dir"], "logs/run_HMM/{F2_sample}.log"),
    params:
        hmm_play_script = "workflow/scripts/TIGER/hmm_play.jar",
    container:
        config["java"]
    shell:
        """
        java -jar {params.hmm_play_script} \
            -r {input.called_bases} \
            -o {output} \
            -t bi \
            -z {input.hmm_probs} \
                2> {log}
        """

rule estimate_breakpoints:
    input:
        # uses the same INPUT as `call_bases`
        counts = rules.call_bases.input,
        hmm_out = rules.run_HMM.output,
        chr_sizes = config["ref_chr_lengths"],
    output:
        main = os.path.join(config["working_dir"], "TIGER/07_est_breakpoints/{F2_sample}.txt"),
        breaks = os.path.join(config["working_dir"], "TIGER/07_est_breakpoints/{F2_sample}.breaks.txt"),
    log:
        os.path.join(config["working_dir"], "logs/estimate_breakpoints/{F2_sample}.log"),
    params:
        prepare_break_script = "workflow/scripts/TIGER/prepare_break.pl",
        sample = "{F2_sample}"
    container:
        config["perl"]
    shell:
        """
        perl {params.prepare_break_script} \
            -s {params.sample} \
            -m {input.counts} \
            -b {input.hmm_out} \
            -c {input.chr_sizes} \
            -o {output} \
                2> {log}
        """

rule refine_breakpoints:
    input:
        # note: this should be the "complete" set of markers, but we're using the filtered markers 
        counts = rules.call_bases.input,
        rough_breaks = rules.estimate_breakpoints.output.breaks,
    output:
        recomb_file = os.path.join(config["working_dir"], "TIGER/07_est_breakpoints/{F2_sample}.recomb.txt"),
        refined_breaks = os.path.join(config["working_dir"], "TIGER/07_est_breakpoints/{F2_sample}.refined.breaks.txt"),
        refined_recomb = os.path.join(config["working_dir"], "TIGER/07_est_breakpoints/{F2_sample}.refined.recomb.txt"),
    log:
        os.path.join(config["working_dir"], "logs/refine_breakpoints/{F2_sample}.log"),
    params:
        refine_recombination_break_script = "workflow/scripts/TIGER/refine_recombination_break.pl",
    container:
        config["perl"]
    shell:
        """
        perl {params.refine_recombination_break_script} \
            {input.counts} \
            {input.rough_breaks} \
                2> {log}
        """

rule smooth_breaks:
    input:
        refined_breaks = rules.refine_breakpoints.output.refined_breaks,
    output:
        smoothed_breaks = os.path.join(config["working_dir"], "TIGER/08_final/{F2_sample}.txt"),
    log:
        os.path.join(config["working_dir"], "logs/smooth_breaks/{F2_sample}.log"),
    params:
        breaks_smoother_script = "workflow/scripts/TIGER/breaks_smoother.pl",
    container:
        config["perl"]
    shell:
        """
        perl {params.breaks_smoother_script} \
            -b {input.refined_breaks} \
            -o {output.smoothed_breaks} \
                2> {log}
        """

#rule visualise_tiger_out:
#    input:
#        rough_breaks = rules.estimate_breakpoints.output.breaks,
#        refined_breaks = rules.refine_breakpoints.output.refined_breaks,
#        smoothed_breaks = rules.smooth_breaks.output.smoothed_breaks,
#        allele_freqs = rules.estimate_allele_freqs.output,
#        sliding_windows = rules.calculate_HMM_probs.output.breaks,
#    output:
#        os.path.join(config["working_dir"], "TIGER/09_plots/{F2_sample}.pdf"),
#    log:
#        os.path.join(config["working_dir"], "logs/visualise_tiger_out/{F2_sample}.log"),
#    params:
#        sample = "{F2_sample}",
#        plot_genotyping_script = "workflow/scripts/TIGER/plot_genotyping.R",
#    container:
#        config["R"]
#    shell:
#        """
#        Rscript --vanilla {params.plot_genotyping_script} \
#            {params.sample} \
#            {output} \
#            {input.rough_breaks} \
#            {input.refined_breaks} \
#            {input.smoothed_breaks} \
#            {input.allele_freqs} \
#            {input.sliding_windows} \
#                2> {log}
#        """