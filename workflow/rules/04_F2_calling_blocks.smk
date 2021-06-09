# Create file containing sites that are homozygous for the two F0 parental strains
rule get_homozygous_sites:
    input:
        os.path.join(config["working_dir"], "vcfs/F0/final/all.vcf.gz")
    output:
        os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/all_sites/{contig}.csv")
    params:
        contig = lambda wildcards: wildcards.contig
    container:
        "docker://biocontainers/bcftools:v1.9-1-deb_cv1"
    shell:
        """
        bcftools view \
            --regions {params.contig} \
            --types snps \
            --max-alleles 2 \
            --output-type u \
            {input} |\
        bcftools query \
            --format '%CHROM,%POS,%POS,%REF,%ALT[,%GT]\\n' \
            --include 'GT="hom"' \
            --output {output}
        """

rule get_divergent_sites:
    input:
        os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/all_sites/{contig}.csv")
    output:
        os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/divergent/{contig}.txt")
    run:
        # set column names
        col_names = ["CHROM", "POS_1", "POS_2", "REF", "ALT", "Cab", "Kaga"]
        # read in DF
        df = pd.read_csv(input[0], names = col_names)
        # drop rows with missing data
        df = df.dropna().reset_index(drop = True)
        # take last two columns and remove / and |
        df_new = pd.DataFrame()
        df_new['Cab'] = df.Cab.str.replace('/|\|', '', regex = True)
        df_new['Kaga'] = df.Kaga.str.replace('/|\|', '', regex = True)
        # find which rows are the same
        out = df.loc[df_new['Cab'] != df_new['Kaga']]
        # write to file
        out.to_csv(output[0], sep = '\t', header = False, index = False)

rule bind_divergent_sites:
    input:
        expand(os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/divergent/{contig}.txt"),
                contig = get_contigs())
    output:
        os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/final/all.txt")
    run:
        # make empty list
        dfs = list()
        # read in each file as data frame and add to list
        for f in input:
            df = pd.read_csv(f, sep = "\t", header = None)
            dfs.append(df)
        # combine and write to file
        pd.concat(dfs).to_csv(output[0], sep = '\t', header = False, index = False)

rule extract_dp4:
    input:
        