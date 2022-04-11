# Get read counts supporting each F0-homozygous-divergent allele
rule bam_readcount_F2:
    input:
        bam = rules.mark_duplicates_F2.output.bam,
        index = rules.samtools_index_F2.output,
        sites_file = rules.extract_homo_div_snps.output.sites,
        ref = rules.get_genome.output,
    output:
        os.path.join(
            config["working_dir"],
            "dp4s/F2/{ref}/F1_het_min_DP/{F2_sample}.dp4.txt"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/bam_readcount_F2_F1_het_min_DP/{ref}/{F2_sample}.log"
        ),
    resources:
        mem_mb = 200
    container:
        config["bam-readcount"]
    shell:
        """
        bam-readcount \
            -l {input.sites_file} \
            -f {input.ref} \
            {input.bam} | \
            cut -f 1,15,28,41,54,67 -d ":" | sed 's/=//g' | sed 's/\\t:/\\t/g' | sed 's/:/\\t/g' \
                > {output} 2> {log}
        """

# Make dpAB files
rule make_dp_AB_F2:
    input:
        dp4 = rules.bam_readcount_F2.output,
        sites_file = rules.extract_homo_div_snps.output.sites,
    output:
        os.path.join(
            config["working_dir"],
            "dpABs/F2/{ref}/F1_het_min_DP/{F2_sample}.txt"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/make_dp_AB_F2/F1_het_min_DP/{ref}/{F2_sample}.log"
        ),
    resources:
        mem_mb = 2000
    script:
        "../scripts/make_dp_AB.py"

# Create HMM inputs
rule make_hmm_input:
    input:
        expand(rules.make_dp_AB_F2.output,
            ref = "hdrr",
            F2_sample = F2_samples['SAMPLE']
        ),
    output:
        os.path.join(
            config["working_dir"],
            "hmm_in/F2/{ref}/F1_het_min_DP/{max_reads}/{bin_length}.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/make_hmm_input/F1_het_min_DP/{ref}/{max_reads}/{bin_length}.log"
        ),
    params:
        bin_length = "{bin_length}",
        max_reads = "{max_reads}"
    resources:
        mem_mb = 10000
    container:
        config["tidyverse_4.1.3"]
    script:
        "../scripts/make_hmm_input.R"


# Test HMM with hmmlearn
rule test_hmmlearn:
    input:
        rules.make_hmm_input.output,
    output:
        os.path.join(
            config["working_dir"],
            "hmm_out/F2/{ref}/hmmlearn/{max_reads}/{bin_length}/{mod}.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/test_hmmlearn/F1_het_min_DP/{ref}/{max_reads}/{bin_length}/{mod}.log"
        ),
    params:
        mod = "{mod}",
        low_cov_samples = lambda wildcards: config["low_cov_samples"]
    resources:
        mem_mb = 5000
    container:
        config["hmmlearn"]
    script:
        "../scripts/test_hmmlearn.py"