# Create file containing sites that are homozygous for the two F0 parental strains
#rule get_homozygous_sites:
#    input:
#        os.path.join(config["working_dir"], "vcfs/F0/final/all.vcf.gz")
#    output:
#        os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/all_sites/{contig}.csv")
#    params:
#        contig = lambda wildcards: wildcards.contig
#    container:
#        "docker://biocontainers/bcftools:v1.9-1-deb_cv1"
#    shell:
#        """
#        bcftools view \
#            --regions {params.contig} \
#            --types snps \
#            --max-alleles 2 \
#            --output-type u \
#            {input} |\
#        bcftools query \
#            --format '%CHROM,%POS,%POS,%REF,%ALT[,%GT]\\n' \
#            --include 'GT="hom"' \
#            --output {output}
#        """

# Extract all homozygous sites per F0 sample
rule get_homozygous_sites:
    input:
        os.path.join(config["working_dir"], "vcfs/F0/final/all.vcf.gz")
    output:
        os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_all/{sample}.csv")
    params:
        sample = lambda wildcards: wildcards.sample
    container:
        "docker://biocontainers/bcftools:v1.9-1-deb_cv1"
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
        expand(os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_all/{sample}.csv"),
                sample = config["F0_lines"]
        ),
    output:
        os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/all.txt")
    container:
        "docker://quay.io/biocontainers/pandas:1.1.5"
    script:
        "../scripts/get_divergent_sites.py"

#rule get_divergent_sites:
#    input:
#        os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_all/{contig}.csv")
#    output:
#        os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/{contig}.txt")
#    run:
#        # set column names
#        col_names = ["CHROM", "POS_1", "POS_2", "REF", "ALT", "Cab", "Kaga"]
#        # read in DF
#        df = pd.read_csv(input[0], names = col_names)
#        # drop rows with missing data
#        df = df.dropna().reset_index(drop = True)
#        # take last two columns and remove / and |
#        df_new = pd.DataFrame()
#        target_cols = ["Cab", "Kaga"]
#        df_new[target_cols] = df.loc[:, target_cols].replace('/|\|', '', regex = True)
#        # find which rows are divergent
#        out = df.loc[df_new['Cab'] != df_new['Kaga']].copy()
#        # replace | with / in Cab and Kaga columns
#        out.loc[:, target_cols] = out.loc[:, target_cols].replace('\|', '/', regex = True)
#        # write to file
#        out.to_csv(output[0], sep = '\t', header = False, index = False)
#
#
#rule bind_divergent_sites:
#    input:
#        expand(os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/divergent/{contig}.txt"),
#                contig = get_contigs())
#    output:
#        os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/final/all.txt")
#    run:
#        # make empty list
#        dfs = list()
#        # read in each file as data frame and add to list
#        for f in input:
#            df = pd.read_csv(f, sep = "\t", header = None)
#            dfs.append(df)
#        # combine and write to file
#        pd.concat(dfs).to_csv(output[0], sep = '\t', header = False, index = False)


rule bam_readcount:
    input:
        bam = os.path.join(config["working_dir"], "bams/F2/bwamem2/marked/{sample}.bam"),
        index = os.path.join(config["working_dir"], "bams/F2/bwamem2/marked/{sample}.bam.bai"),
        sites_file = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/homo_divergent/all.txt"),
        ref = config["ref_prefix"] + ".fasta",
    output:
        os.path.join(config["working_dir"], "dp4s/batch_01/bwamem2/{sample}.dp4.txt")
    container:
        "docker://quay.io/biocontainers/bam-readcount:0.8--py36pl526h94a8ba4_4"
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
        dp4 = os.path.join(config["working_dir"], "dp4s/batch_01/bwamem2/{sample}.dp4.txt"),
        sites_file = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/final/all.txt"),
    output:
        os.path.join(config["working_dir"], "dpABs/batch_01/bwamem2/{sample}.txt"),
    script:
        "../scripts/make_dp_AB.py"


rule run_rc_block:
    input:
        dp_files = expand(os.path.join(config["working_dir"], "dpABs/batch_01/bwamem2/{sample}.txt"),
            sample = F2_samples['SAMPLE']
        ),
        source_code = "workflow/scripts/rc_block_hmm.R"
    output:
        os.path.join(config["data_store_dir"], "hmm_output_batch_01_corrected.txt"),
    container:
        "docker://brettellebi/somite_f2:latest"
    script:
        "../scripts/run_rc_block.R"
