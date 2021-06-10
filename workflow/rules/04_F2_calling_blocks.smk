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

rule bwa_mem2_mem:
    input:
        reads=get_fastq_F2,
        idx=rules.bwa_index.output,
    output:
        os.path.join(config["working_dir"], "bams/F2/bwamem/sorted/{sample}.bam"),
    params:
        index=lambda wildcards: config["ref_prefix"] + ".fasta",
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="samtools",
        sort_order="coordinate",
        sort_extra=""
    threads: 8
    wrapper:
        "v0.75.0/bio/bwa-mem2/mem"

rule map_F2:
    input:
        target = config["ref_prefix"] + ".fasta",
        query = get_fastq_F2,
    output:
        os.path.join(config["working_dir"], "sams/F2/mapped/{sample}.sam")
    params:
        extra="-ax sr"
    threads: 4
    wrapper:
        "0.74.0/bio/minimap2/aligner"

## Note
# sample 121 fails with the following error:
# from /hps/nobackup/birney/users/ian/somites/logs/map_F2/2432668_sample=121.err
#minimap2: sketch.c:84: mm_sketch: Assertion `len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)' failed.
#/usr/local/bin/bash: line 1: 2362261 Aborted                 (core dumped) ( minimap2 -t 4 -ax sr -o /hps/nobackup/birney/users/ian/somites/sams/F2/mapped/121.sam /hps/nobackup/birney/users/ian/somites/refs/hdrr.fasta /hps/nobackup/birney/projects/indigene/raw_data/Ali/Kaga-Cab_F2_First200WGS/HKNVLBBXY_Pool_1_21s001372-1-1_Seleit_lane121_1_sequence.txt.gz /hps/nobackup/birney/projects/indigene/raw_data/Ali/Kaga-Cab_F2_First200WGS/HKNVLBBXY_Pool_1_21s001372-1-1_Seleit_lane121_2_sequence.txt.gz )

rule replace_rg_f2:
    input:
        os.path.join(config["working_dir"], "sams/F2/mapped/{sample}.sam")
    output:
        os.path.join(config["working_dir"], "sams/F2/grouped/{sample}.sam")
    container:
        "docker://quay.io/biocontainers/picard:2.9.2--2"
    params:
        "RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM={sample}"
    shell:
        """
        picard AddOrReplaceReadGroups \
            -Xmx1024M \
            {params} \
            I={input} \
            O={output}
        """

# Many failed jobs with this rule, e.g.
# 365 fails the above stage:
# Exception in thread "main" htsjdk.samtools.SAMFormatException: Error parsing text SAM file. Invalid character in read bases; File /hps/nobackup/birney/users/ian/somites/sams/F2/mapped/365.sam; Line 7975375
# Line: ST-K00119:220:HKNVLBBXY:3:1212:252AJJJJJJ5TTGJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJFJFJJJJJJJJJJJJJJJFFJJJJJJJFJGGCA+TJJJJJJJFJGGCA+TJJJJJJJFJGCTGT84444JJJJCATCTATTACCATCTGACTTCCACCTCCTCC 4       *       0       0       **0       0       GGJFJ<J7JJAAJ<AJF-F-FJAJJ77FJJFJAJ-GGTGGTCTCTG5TTGJGGCTCGGTGAACATGJJJJ<CTTCCTCAAJJJJJJJJJJJJJ:3:1212:2FJJJJJJ<FJJFTJTFJJFJJFJY:3::::::TCTTACGCAGTCGGGTCAJ *       rl:i:0

# 5119:
# Exception in thread "main" htsjdk.samtools.SAMFormatException: Error parsing text SAM file. Empty field at position 9 (zero-based); File /hps/nobackup/birney/users/ian/somites/sams/F2/mapped/5119.sam; Line 8872706
#Line: ST-K00119:220:HKNVLBBXY:5:1JJJJJJJJJJJJJJJA1:1CGTAJJJJJJJJJJJA1:1GCCAATCACAGACATGCGGGAACCGGCCTTCGAACCGJJJJJJCCAGACATGCJ01   4       *       0       0       *       *       0       0               *       rl:i:0

# So try using bwa-mem instead

rule sort_sam_f2:
    input:
        os.path.join(config["working_dir"], "sams/F2/grouped/{sample}.sam")
    output:
        os.path.join(config["working_dir"], "bams/F2/sorted/{sample}.bam")
    params:
        sort_order="coordinate",
        extra=lambda wildcards: "VALIDATION_STRINGENCY=LENIENT TMP_DIR=" + config["tmp_dir"]
    resources:
        mem_mb=1024
    wrapper:
        "0.74.0/bio/picard/sortsam"

rule mark_duplicates_f2:
    input:
        os.path.join(config["working_dir"], "bams/F2/sorted/{sample}.bam")
    output:
        bam=os.path.join(config["working_dir"], "bams/F2/marked/{sample}.bam"),
        metrics=os.path.join(config["working_dir"], "bams/F2/marked/{sample}.metrics.txt")
    params:
        lambda wildcards: "REMOVE_DUPLICATES=true TMP_DIR=" + config["tmp_dir"]
    resources:
        mem_mb=1024
    wrapper:
        "0.74.0/bio/picard/markduplicates"

rule samtools_index_f2:
    input:
        os.path.join(config["working_dir"], "bams/F2/marked/{sample}.bam")
    output:
        os.path.join(config["working_dir"], "bams/F2/marked/{sample}.bam.bai")
    wrapper:
        "0.74.0/bio/samtools/index"
