# Extract all homozygous sites per F0 sample
rule get_homozygous_sites:
    input:
        os.path.join(config["data_store_dir"], "vcfs/F0/final/all.vcf.gz")
    output:
        os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_all/{F0_sample}.csv")
    log:
        os.path.join(config["working_dir"], "logs/get_homozygous_sites/{F0_sample}.log")
    params:
        sample = "{F0_sample}"
    resources:
        mem_mb = 5000
    container:
        config["bcftools"]
    shell:
        """
        bcftools view \
            --samples {params.sample} \
            --types snps \
            --max-alleles 2 \
            --output-type u \
            {input} |\
        bcftools query \
            --include 'GT="hom"' \
            --format '%CHROM,%POS,%POS,%REF,%ALT,[%GT]\\n' \
            --output {output} \
                2> {log}
        """

rule get_divergent_sites:
    input:
    # Note we're using the order of the F0 lines provided in config["F0_lines"].
    # This order is critical for the next steps, so ensure it is correct.
        expand(os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_all/{F0_sample}.csv"),
                F0_sample = config["F0_lines"]
        ),
    output:
        os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/all.txt")
    log:
        os.path.join(config["working_dir"], "logs/get_divergent_sites/all.log")
    resources:
        mem_mb = 5000
    container:
        config["pandas"]
    script:
        "../scripts/get_divergent_sites.py"

# Exclude target sites that overlap repeat regions
rule exclude_repeat_sites_F2:
    input:
        target_sites = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/all.txt"),
        repeats_file = config["hdrr_repeats_gff"]
    output:
        os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/no_repeats.txt")
    log:
        os.path.join(config["working_dir"], "logs/exclude_repeat_sites_F2/all.log")
    resources:
        mem_mb = 5000
    container:
        config["R"]
    script:
        "../scripts/exclude_repeat_sites_F2.R"

# Exclude "black list" sites found to have persistent heterozygosity in the MIKK panel
rule filter_black_list:
    input:
        sites_excl_repeats = rules.exclude_repeat_sites_F2.output,
        black_list = config["het_black_list"],
    output:
        os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/no_repeats_no_persistent_hets.txt"),
    log:
        os.path.join(config["working_dir"], "logs/filter_black_list/all.log")
    resources:
        mem_mb = 5000
    container:
        config["R"]
    script:
        "../scripts/filter_black_list.R"        

# Create a BED file of repeat regions
rule create_repeats_bed:
    input:
        config["hdrr_repeats_gff"]
    output:
        os.path.join(config["data_store_dir"], "repeats/hdrr_repeats.bed")
    log:
        os.path.join(config["working_dir"], "logs/create_repeats_bed/all.log")
    resources:
        mem_mb = 5000
    container:
        config["R"]
    script:
        "../scripts/create_repeats_bed.R"

# Filter reads that overlap repeat regions
rule filter_repeats_from_bams_F2:
    input:
        bam = os.path.join(config["working_dir"], "bams/F2/bwamem2/marked/{F2_sample}.bam"),
        repeats_bed = os.path.join(config["data_store_dir"], "repeats/hdrr_repeats.bed"),
    output:
        keep = os.path.join(config["working_dir"], "bams/F2/bwamem2/no_repeat_keep/{F2_sample}.bam"),
        throw = os.path.join(config["working_dir"], "bams/F2/bwamem2/no_repeat_throw/{F2_sample}.bam"),
    log:
        os.path.join(config["working_dir"], "logs/filter_repeats_from_bams_F2/{F2_sample}.log"),
    resources:
        mem_mb = 5000
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

rule samtools_index_f2_filtered_bams:
    input:
        os.path.join(config["working_dir"], "bams/F2/bwamem2/no_repeat_keep/{F2_sample}.bam")
    output:
        os.path.join(config["working_dir"], "bams/F2/bwamem2/no_repeat_keep/{F2_sample}.bam.bai")
    log:
        os.path.join(config["working_dir"], "logs/samtools_index_f2_filtered_bams/{F2_sample}.log") 
    resources:
        mem_mb = 5000
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
    ## All sites
rule bam_readcount_F2_all:
    input:
        bam = os.path.join(config["working_dir"], "bams/F2/bwamem2/marked/{F2_sample}.bam"),
        index = os.path.join(config["working_dir"], "bams/F2/bwamem2/marked/{F2_sample}.bam.bai"),
        sites_file = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/all.txt"),
        ref = config["ref_prefix"] + ".fasta",
    output:
        os.path.join(config["working_dir"], "dp4s/F2/all_sites/{F2_sample}.dp4.txt")
    log:
        os.path.join(config["working_dir"], "logs/bam_readcount_F2_all/{F2_sample}.log")
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

    ## Excluding SITES overlapping repeat regions
rule bam_readcount_F2_excl_repeat_sites:
    input:
        bam = os.path.join(config["working_dir"], "bams/F2/bwamem2/marked/{F2_sample}.bam"),
        index = os.path.join(config["working_dir"], "bams/F2/bwamem2/marked/{F2_sample}.bam.bai"),
        # `sites_file` differs from above
        sites_file = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/no_repeats.txt"),
        ref = config["ref_prefix"] + ".fasta",
    output:
        os.path.join(config["working_dir"], "dp4s/F2/no_repeat_sites/{F2_sample}.dp4.txt")
    log:
        os.path.join(config["working_dir"], "logs/bam_readcount_F2_excl_repeat_sites/{F2_sample}.log")
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

    ## Excluding READS overlapping repeat regions
rule bam_readcount_F2_excl_repeat_reads:
    input:
        # `bam` and `index` differ from above
        bam = os.path.join(config["working_dir"], "bams/F2/bwamem2/no_repeat_keep/{F2_sample}.bam"),
        index = os.path.join(config["working_dir"], "bams/F2/bwamem2/no_repeat_keep/{F2_sample}.bam.bai"),
        sites_file = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/no_repeats.txt"),
        ref = config["ref_prefix"] + ".fasta",
    output:
        os.path.join(config["working_dir"], "dp4s/F2/no_repeat_reads/{F2_sample}.dp4.txt")
    log:
        os.path.join(config["working_dir"], "logs/bam_readcount_F2_excl_repeat_reads/{F2_sample}.log")
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

    ## Excluding READS overlapping repeat regions and persistent het black list
rule bam_readcount_F2_excl_repeat_reads_and_black_list:
    input:
        # `bam` and `index` differ from above
        bam = os.path.join(config["working_dir"], "bams/F2/bwamem2/no_repeat_keep/{F2_sample}.bam"),
        index = os.path.join(config["working_dir"], "bams/F2/bwamem2/no_repeat_keep/{F2_sample}.bam.bai"),
        sites_file = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/no_repeats_no_persistent_hets.txt"),
        ref = config["ref_prefix"] + ".fasta",
    output:
        os.path.join(config["working_dir"], "dp4s/F2/no_repeat_reads_or_pers_hets/{F2_sample}.dp4.txt")
    log:
        os.path.join(config["working_dir"], "logs/bam_readcount_F2_excl_repeat_reads_and_black_list/{F2_sample}.log")
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
    ## Excluding READS overlapping repeat regions and persistent het black list, and excluding sites with extreme F2-wide read counts and proportion of Cab
rule bam_readcount_F2_excl_repeat_reads_and_black_list_and_extreme_read_count_and_cab_prop:
    input:
        # `bam` and `index` differ from above
        bam = os.path.join(config["working_dir"], "bams/F2/bwamem2/no_repeat_keep/{F2_sample}.bam"),
        index = os.path.join(config["working_dir"], "bams/F2/bwamem2/no_repeat_keep/{F2_sample}.bam.bai"),
        sites_file = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/no_repeats_no_persistent_hets_filtered_for_read_count_and_cab_prop.txt"),
        ref = config["ref_prefix"] + ".fasta",
    output:
        os.path.join(config["working_dir"], "dp4s/F2/no_repeat_reads_or_pers_hets_filtered_for_read_count_and_cab_prop/{F2_sample}.dp4.txt")
    log:
        os.path.join(config["working_dir"], "logs/bam_readcount_F2_excl_repeat_reads_and_black_list/{F2_sample}.log")
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

# Make dpAB files
rule make_dp_AB_F2:
    input:
        dp4 = os.path.join(config["working_dir"], "dp4s/F2/{site_filter}/{F2_sample}.dp4.txt"),
        sites_file = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/all.txt"),
    output:
        os.path.join(config["working_dir"], "dpABs/F2/{site_filter}/{F2_sample}.txt"),
    log:
        os.path.join(config["working_dir"], "logs/make_dp_AB_F2/{site_filter}/{F2_sample}.log")
    resources:
        mem_mb = 10000
    script:
        "../scripts/make_dp_AB.py"

# Run HMM recombination blocks
rule run_rc_block_F2:
    input:
        dp_files = expand(os.path.join(config["working_dir"], "dpABs/F2/{{site_filter}}/{F2_sample}.txt"),
            F2_sample = F2_samples['SAMPLE']
        ),
        source_code = "workflow/scripts/rc_block_hmm.R"
    output:
        os.path.join(config["data_store_dir"], "recombination_blocks/F2/{site_filter}/{bin_length}.txt"),
    log:
        os.path.join(config["working_dir"], "logs/run_rc_block_F2/{site_filter}/{bin_length}.log")
    params:
        bin_length = lambda wildcards: wildcards.bin_length
    resources:
        mem_mb = 50000
    container:
        config["R"]
    script:
        "../scripts/run_rc_block.R"

rule consolidate_dbABs_Ewan:
    input:
        dp_files = expand(os.path.join(config["working_dir"], "dpABs/F2/{{site_filter}}/{F2_sample}.txt"),
            F2_sample = F2_samples['SAMPLE']
        ),
    output:
        os.path.join(config["data_store_dir"], "dpABs/F2_consolidated/{site_filter}.txt"),
    log:
        os.path.join(config["working_dir"], "logs/consolidate_dbABs_Ewan/{site_filter}.log")
    resources:
        mem_mb = 350000
    container:
        config["tidyverse_4"]
    script:
        "../scripts/consolidate_dbABs_Ewan.R"

rule filter_sites_for_read_count_and_cab_prop:
    input:
        dp_files = expand(os.path.join(config["working_dir"], "dpABs/F2/no_repeat_reads_or_pers_hets/{F2_sample}.txt"),
            F2_sample = F2_samples['SAMPLE']
        ),
        sites_file = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/no_repeats_no_persistent_hets.txt"),
    output:
        counts_pre_filter  = os.path.join(config["data_store_dir"], "dpABs/F2_total_read_counts/pre_filter.txt.gz"),
        counts_post_filter = os.path.join(config["data_store_dir"], "dpABs/F2_total_read_counts/post_filter.txt.gz"),
        filtered_sites = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/no_repeats_no_persistent_hets_filtered_for_read_count_and_cab_prop.txt"),
    log:
        os.path.join(config["working_dir"], "logs/filter_sites_for_read_count_and_cab_prop/all.log")
    params:
        min_reads = config["min_reads"],
        max_reads = config["max_reads"],
        min_prop_cab = config["min_prop_cab"],
        max_prop_cab = config["max_prop_cab"]
    resources:
        mem_mb = 80000
    container:
        config["R"]
    script:
        "../scripts/filter_sites_for_read_count_and_cab_prop.R"

rule plot_recombination_blocks:
    input:
        os.path.join(config["data_store_dir"], "recombination_blocks/F2/{site_filter}/{bin_length}.txt"),
    output:
        base_cov_total = "book/plots/snakemake/{site_filter}/{bin_length}/base_cov_total.png",
        base_cov_by_chrom = "book/plots/snakemake/{site_filter}/{bin_length}/base_cov_by_chrom.png",
        prop_sites_total = "book/plots/snakemake/{site_filter}/{bin_length}/prop_sites_total.png",
        prop_sites_by_chrom = "book/plots/snakemake/{site_filter}/{bin_length}/prop_sites_by_chrom.png",
        karyoplot_no_missing = "book/plots/snakemake/{site_filter}/{bin_length}/karyoplot_no_missing.png",
        karyoplot_with_missing = "book/plots/snakemake/{site_filter}/{bin_length}/karyoplot_with_missing.png",
    log:
        os.path.join(config["working_dir"], "logs/plot_recombination_blocks/{site_filter}/{bin_length}.log")
    params:
        site_filter = "{site_filter}",
        bin_length = "{bin_length}"
    resources:
        mem_mb = 40000
    container:
        config["R"]
    script:
        "../scripts/plot_recombination_blocks.R"