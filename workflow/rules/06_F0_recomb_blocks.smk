# Filter reads that overlap repeat regions
rule filter_repeats_from_bams_F0:
    input:
        bam = os.path.join(config["working_dir"], "bams/F0/merged/{F0_sample}.bam"),
        repeats_bed = os.path.join(config["data_store_dir"], "repeats/hdrr_repeats.bed"),
    output:
        keep = os.path.join(config["working_dir"], "bams/F0/no_repeat_keep/{F0_sample}.bam"),
        throw = os.path.join(config["working_dir"], "bams/F0/no_repeat_throw/{F0_sample}.bam"),
    log:
        os.path.join(config["working_dir"], "logs/filter_repeats_from_bams_F0/{F0_sample}.log"),
    container:
        config["samtools"]
    shell:
        """
        samtools view \
            --bam \
            --targets-file {input.repeats_bed} \
            --output-unselected {output.keep} \
            {input.bam} \
            > {output.throw} \
                2> {log}
        """

rule samtools_index_F0_filtered_bams:
    input:
        os.path.join(config["working_dir"], "bams/F0/no_repeat_keep/{F0_sample}.bam")
    output:
        os.path.join(config["working_dir"], "bams/F0/no_repeat_keep/{F0_sample}.bam.bai")
    log:
        os.path.join(config["working_dir"], "logs/samtools_index_F0_filtered_bams/{F0_sample}.log") 
    container:
        config["samtools"]
    shell:
        """
        samtools index \
            {input[0]} \
            {output[0]} \
                2> {log}
        """

# Bam readcounts
    ## all sites
rule bam_readcount_F0_all_sites:
    input:
        bam = os.path.join(config["working_dir"], "bams/F0/merged/{F0_sample}.bam"),
        index = os.path.join(config["working_dir"], "bams/F0/merged/{F0_sample}.bam.bai"),
        sites_file = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/all.txt"),
        ref = config["ref_prefix"] + ".fasta"
    output:
        os.path.join(config["working_dir"], "dp4s/F0/all_sites/{F0_sample}.dp4.txt"),
    log:
        os.path.join(config["working_dir"], "logs/bam_readcount_F0_all_sites/{F0_sample}.log")
    container:
        config["bam-readcount"]
    shell:
        """
        bam-readcount \
            -l {input.sites_file} \
            -f {input.ref} \
            {input.bam} | \
            cut -f 1,15,28,41,54,67 -d ":" | sed 's/=//g' | sed 's/\\t:/\\t/g' | sed 's/:/\\t/g' \
                > {output} \
                2> {log}
        """

    ## Filter out SITES that overlap repeat regions
rule bam_readcount_F0_excl_repeat_sites:
    input:
        bam = os.path.join(config["working_dir"], "bams/F0/merged/{F0_sample}.bam"),
        index = os.path.join(config["working_dir"], "bams/F0/merged/{F0_sample}.bam.bai"),
        # `sites_file` differs from above
        sites_file = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/no_repeats.txt"),
        ref = config["ref_prefix"] + ".fasta",
    output:
        os.path.join(config["working_dir"], "dp4s/F0/no_repeat_sites/{F0_sample}.dp4.txt")
    log:
        os.path.join(config["working_dir"], "logs/bam_readcount_F2_excl_repeat_sites/{F0_sample}.log")
    resources:
        mem_mb = 10000
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

    ## Filter out READS that overlap repeat regions 
rule bam_readcount_F0_excl_repeat_reads:
    input:
        #bam = os.path.join(config["working_dir"], "bams/F0/merged/{F0_sample}.bam"),
        #index = os.path.join(config["working_dir"], "bams/F0/merged/{F0_sample}.bam.bai"),
        bam = os.path.join(config["working_dir"], "bams/F0/no_repeat_keep/{F0_sample}.bam"),
        index = os.path.join(config["working_dir"], "bams/F0/no_repeat_keep/{F0_sample}.bam.bai"),
        sites_file = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/no_repeats.txt"),
        ref = config["ref_prefix"] + ".fasta"
    output:
        os.path.join(config["working_dir"], "dp4s/F0/no_repeat_reads/{F0_sample}.dp4.txt"),
    log:
        os.path.join(config["working_dir"], "logs/bam_readcount_F0_excl_repeat_reads/{F0_sample}.log")
    container:
        config["bam-readcount"]
    shell:
        """
        bam-readcount \
            -l {input.sites_file} \
            -f {input.ref} \
            {input.bam} | \
            cut -f 1,15,28,41,54,67 -d ":" | sed 's/=//g' | sed 's/\\t:/\\t/g' | sed 's/:/\\t/g' \
                > {output} \
                2> {log}
        """

rule make_dp_AB_F0:
    input:
        dp4 = os.path.join(config["working_dir"], "dp4s/F0/{site_filter}/{F0_sample}.dp4.txt"),
        sites_file = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/all.txt"),
    output:
        os.path.join(config["working_dir"], "dpABs/F0/{site_filter}/{F0_sample}.txt"),
    log:
        os.path.join(config["working_dir"], "logs/make_dp_AB_F0/{site_filter}/{F0_sample}.log")
    resources:
        mem_mb = 10000
    script:
        "../scripts/make_dp_AB.py"

rule run_rc_block_F0:
    input:
        dp_files = expand(os.path.join(config["working_dir"], "dpABs/F0/{{site_filter}}/{F0_sample}.txt"),
            F0_sample = config["F0_lines"]
        ),
        source_code = "workflow/scripts/rc_block_hmm.R"
    output:
        os.path.join(config["data_store_dir"], "recombination_blocks/F0/{site_filter}/{bin_length}.txt"),
    log:
        os.path.join(config["working_dir"], "logs/run_rc_block_F0/{site_filter}/{bin_length}.log")
    params:
        bin_length = lambda wildcards: wildcards.bin_length
    resources:
        mem_mb = 50000
    container:
        config["R"]
    script:
        "../scripts/run_rc_block.R"