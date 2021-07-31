# Extract all homozygous sites per F0 sample
rule get_homozygous_sites:
    input:
        os.path.join(config["data_store_dir"], "vcfs/F0/final/all.vcf.gz")
    output:
        os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_all/{F0_sample}.csv")
    params:
        sample = lambda wildcards: wildcards.sample
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
            --output {output}
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
    container:
        config["pandas"]
    script:
        "../scripts/get_divergent_sites.py"

rule bam_readcount:
    input:
        bam = os.path.join(config["working_dir"], "bams/F2/bwamem2/marked/{F2_sample}.bam"),
        index = os.path.join(config["working_dir"], "bams/F2/bwamem2/marked/{F2_sample}.bam.bai"),
        sites_file = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/all.txt"),
        ref = config["ref_prefix"] + ".fasta",
    output:
        os.path.join(config["working_dir"], "dp4s/batch_01/bwamem2/{F2_sample}.dp4.txt")
    container:
        config["bam-readcount"]
    shell:
        """
        bam-readcount \
            -l {input.sites_file} \
            -f {input.ref} \
            {input.bam} | \
            cut -f 1,15,28,41,54,67 -d ":" | sed 's/=//g' | sed 's/\\t:/\\t/g' | sed 's/:/\\t/g' \
                > {output}
        """

rule make_dp_AB:
    input:
        dp4 = os.path.join(config["working_dir"], "dp4s/batch_01/bwamem2/{F2_sample}.dp4.txt"),
        sites_file = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/final/all.txt"),
    output:
        os.path.join(config["working_dir"], "dpABs/batch_01/bwamem2/{F2_sample}.txt"),
    script:
        "../scripts/make_dp_AB.py"

rule run_rc_block:
    input:
        dp_files = expand(os.path.join(config["working_dir"], "dpABs/batch_01/bwamem2/{F2_sample}.txt"),
            F2_sample = F2_samples['SAMPLE']
        ),
        source_code = "workflow/scripts/rc_block_hmm.R"
    output:
        os.path.join(config["data_store_dir"], "recombination_blocks/20210716_hmm_output_batch_01_binlen_{bin_length}.txt"),
    params:
        bin_length = lambda wildcards: wildcards.bin_length
    container:
        config["R"]
    script:
        "../scripts/run_rc_block.R"
