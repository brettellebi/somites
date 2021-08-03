rule get_ref_alt_sites:
    input:
        os.path.join(config["data_store_dir"], "vcfs/F0/final/all.vcf.gz")
    output:
        os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/ref_alt/all.txt")
    log:
        os.path.join(config["working_dir"], "logs/get_ref_alt_sites/all.log")
    container:
        config["bcftools"]
    shell:
        """
        bcftools view \
            --types snps \
            --max-alleles 2 \
            --output-type u \
            {input[0]} |\
        bcftools query \
            --format '%CHROM\\t%POS\\t%POS\\t%REF\\t%ALT\\t0/0\\t1/1\\n' \
            --output {output[0]} \
                2> {log}
        """

rule bam_readcount_F0:
    input:
        bam = os.path.join(config["working_dir"], "bams/F0/merged/{F0_sample}.bam"),
        index = os.path.join(config["working_dir"], "bams/F0/merged/{F0_sample}.bam.bai"),
        sites_file = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/ref_alt/all.txt"),
        ref = config["ref_prefix"] + ".fasta"
    output:
        os.path.join(config["working_dir"], "dp4s/F0/{F0_sample}.dp4.txt"),
    log:
        os.path.join(config["working_dir"], "logs/bam_readcount_F0/{F0_sample}.log")
    resources:
        mem_mb = 30000
    threads:
        4
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

rule make_dp_AB_F0:
    input:
        dp4 = os.path.join(config["working_dir"], "dp4s/F0/{F0_sample}.dp4.txt"),
        sites_file = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/ref_alt/all.txt"),
    output:
        os.path.join(config["working_dir"], "dpABs/F0/{F0_sample}.txt"),
    log:
        os.path.join(config["working_dir"], "logs/make_dp_AB_F0/{F0_sample}.log")
    resources:
        mem_mb = 50000
    script:
        "../scripts/make_dp_AB.py"

rule run_rc_block_F0:
    input:
        dp_files = expand(os.path.join(config["working_dir"], "dpABs/F0/{F0_sample}.txt"),
            F0_sample = config["F0_lines"]
        ),
        source_code = "workflow/scripts/rc_block_hmm.R"
    output:
        os.path.join(config["data_store_dir"], "recombination_blocks/20210802_hmm_output_F0_binlen_{bin_length_F0}.txt"),
    log:
        os.path.join(config["working_dir"], "logs/run_rc_block_F0/{bin_length_F0}.log")
    params:
        bin_length = lambda wildcards: wildcards.bin_length
    resources:
        mem_mb = 50000
    container:
        config["R"]
    script:
        "../scripts/run_rc_block.R"