rule get_ref_alt_sites:
    input:
        os.path.join(config["data_store_dir"], "vcfs/F0/final/all.vcf.gz")
    output:
        os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/ref_alt/all.txt")
    container:
        config["bcftools"]
    shell:
        """
        bcftools query \
            --format '%CHROM\\t%POS\\t%POS\\t%REF\\t%ALT\\t0/0\\t1/1\\n' \
            --output {output} \
            {input}
        """

rule bam_readcount_F0:
    input:
        bam = os.path.join(config["working_dir"], "bams/F0/merged/{F0_sample}.bam"),
        index = os.path.join(config["working_dir"], "bams/F0/merged/{F0_sample}.bam.bai"),
        sites_file = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/ref_alt/all.txt"),
        ref = config["ref_prefix"] + ".fasta",
    output:
        os.path.join(config["working_dir"], "dp4s/F0/{F0_sample}.dp4.txt")
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