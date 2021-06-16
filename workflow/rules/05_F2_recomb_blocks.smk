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
        target_cols = ["Cab", "Kaga"]
        df_new[target_cols] = df.loc[:, target_cols].replace('/|\|', '', regex = True)
        # find which rows are divergent
        out = df.loc[df_new['Cab'] != df_new['Kaga']].copy()
        # replace | with / in Cab and Kaga columns
        out.loc[:, target_cols] = out.loc[:, target_cols].replace('\|', '/', regex = True)
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


rule bam_readcount:
    input:
        bam = os.path.join(config["working_dir"], "bams/F2/bwamem2/marked/{sample}.bam"),
        sites_file = os.path.join(config["working_dir"], "data/sites_files/F0_Cab_Kaga/final/all.txt"),
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
    run:
        # Read in dp4
        dp4_cols = ["CHR", "POS", "REF", "TOTAL", "A", "C", "G", "T", "N"]
#        dp4 = pd.read_csv("/hps/nobackup/birney/users/ian/somites/dp4s/batch_01/bwamem2/11.dp4.txt", sep = "\t", names = dp4_cols, index_col = ["CHR", "POS"])
        dp4 = pd.read_csv(input.dp4, sep = "\t", names = dp4_cols, index_col = ["CHR", "POS"])
        # Read in sites
        sites_cols = ["CHR", "POS", "POS_2", "REF", "ALT", "F0_1", "F0_2"]
#        sites = pd.read_csv("/hps/nobackup/birney/users/ian/somites/data/sites_files/F0_Cab_Kaga/final/all.txt", sep = "\t", names = sites_cols, index_col = ["CHR", "POS"])
        sites = pd.read_csv(input.sites_file, sep = "\t", names = sites_cols, index_col = ["CHR", "POS"])
        # Add columns with F0 and F1 alleles
        sites = sites.assign(F0_1_ALLELE =lambda x: np.where(x["F0_1"] == "0/0", x["REF"], x["ALT"]))
        sites = sites.assign(F0_2_ALLELE =lambda x: np.where(x["F0_2"] == "0/0", x["REF"], x["ALT"]))
        # Filter dp4 for necessary columns
        dp4 = dp4[["A", "C", "G", "T"]]
        # Join dp4 to sites
        joined = dp4.join(sites)
        # Create new columns with counts for F0_1 and F0_2 alleles
        def rules(row, column):
            if row[column] == "A":
                return row["A"]
            elif row[column] == "C":
                return row["C"]
            elif row[column] == "G":
                return row["G"]            
            elif row[column] == "T":
                return row["T"]               
            else:
                return None
        joined['F0_1_COUNT'] = joined.apply(lambda x: rules(x, "F0_1_ALLELE"), 1)
        joined['F0_2_COUNT'] = joined.apply(lambda x: rules(x, "F0_2_ALLELE"), 1)
        # Write to file
#        joined[["F0_1_ALLELE", "F0_1_COUNT", "F0_2_ALLELE", "F0_2_COUNT"]].to_csv("tmp.txt", sep = "\t", header = False)
        joined[["F0_1_ALLELE", "F0_1_COUNT", "F0_2_ALLELE", "F0_2_COUNT"]].to_csv(output[0], sep = "\t", header = False)

rule run_rc_block:
    input:
        os.path.join(config["working_dir"], "dpABs/batch_01/bwamem2/{sample}.txt"),
    output:
        
