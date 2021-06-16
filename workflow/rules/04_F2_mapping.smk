#rule map_F2:
#    input:
#        target = config["ref_prefix"] + ".fasta",
#        query = get_fastq_F2,
#    output:
#        os.path.join(config["working_dir"], "sams/F2/mapped/{sample}.sam")
#    params:
#        extra="-ax sr"
#    threads: 4
#    wrapper:
#        "0.74.0/bio/minimap2/aligner"

## Note
# sample 121 fails with the following error:
# from /hps/nobackup/birney/users/ian/somites/logs/map_F2/2432668_sample=121.err
#minimap2: sketch.c:84: mm_sketch: Assertion `len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)' failed.
#/usr/local/bin/bash: line 1: 2362261 Aborted                 (core dumped) ( minimap2 -t 4 -ax sr -o /hps/nobackup/birney/users/ian/somites/sams/F2/mapped/121.sam /hps/nobackup/birney/users/ian/somites/refs/hdrr.fasta /hps/nobackup/birney/projects/indigene/raw_data/Ali/Kaga-Cab_F2_First200WGS/HKNVLBBXY_Pool_1_21s001372-1-1_Seleit_lane121_1_sequence.txt.gz /hps/nobackup/birney/projects/indigene/raw_data/Ali/Kaga-Cab_F2_First200WGS/HKNVLBBXY_Pool_1_21s001372-1-1_Seleit_lane121_2_sequence.txt.gz )

#rule replace_rg_f2:
#    input:
#        os.path.join(config["working_dir"], "sams/F2/mapped/{sample}.sam")
#    output:
#        os.path.join(config["working_dir"], "sams/F2/grouped/{sample}.sam")
#    container:
#        "docker://quay.io/biocontainers/picard:2.9.2--2"
#    params:
#        "RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM={sample}"
#    shell:
#        """
#        picard AddOrReplaceReadGroups \
#            -Xmx1024M \
#            {params} \
#            I={input} \
#            O={output}
#        """

# Many failed jobs with this rule, e.g.
# 365 fails the above stage:
# Exception in thread "main" htsjdk.samtools.SAMFormatException: Error parsing text SAM file. Invalid character in read bases; File /hps/nobackup/birney/users/ian/somites/sams/F2/mapped/365.sam; Line 7975375
# Line: ST-K00119:220:HKNVLBBXY:3:1212:252AJJJJJJ5TTGJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJFJFJJJJJJJJJJJJJJJFFJJJJJJJFJGGCA+TJJJJJJJFJGGCA+TJJJJJJJFJGCTGT84444JJJJCATCTATTACCATCTGACTTCCACCTCCTCC 4       *       0       0       **0       0       GGJFJ<J7JJAAJ<AJF-F-FJAJJ77FJJFJAJ-GGTGGTCTCTG5TTGJGGCTCGGTGAACATGJJJJ<CTTCCTCAAJJJJJJJJJJJJJ:3:1212:2FJJJJJJ<FJJFTJTFJJFJJFJY:3::::::TCTTACGCAGTCGGGTCAJ *       rl:i:0

# 5119:
# Exception in thread "main" htsjdk.samtools.SAMFormatException: Error parsing text SAM file. Empty field at position 9 (zero-based); File /hps/nobackup/birney/users/ian/somites/sams/F2/mapped/5119.sam; Line 8872706
#Line: ST-K00119:220:HKNVLBBXY:5:1JJJJJJJJJJJJJJJA1:1CGTAJJJJJJJJJJJA1:1GCCAATCACAGACATGCGGGAACCGGCCTTCGAACCGJJJJJJCCAGACATGCJ01   4       *       0       0       *       *       0       0               *       rl:i:0

# So try using bwa-mem instead

#rule bwa_mem2_mem:
#    input:
#        reads=get_fastq_F2,
#        idx=rules.bwa_mem2_index.output,
#    output:
#        os.path.join(config["working_dir"], "bams/F2/bwamem2/sorted/{sample}.bam"),
#    params:
#        index=lambda wildcards: config["ref_prefix"] + ".fasta",
#        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
#        sort="samtools",
#        sort_order="coordinate",
#        sort_extra=""
#    threads: 8
#    wrapper:
#        "v0.75.0/bio/bwa-mem2/mem"

# Wrapper doesn't work with many files, so use container instead

rule bwa_mem2_mem:
    input:
        reads=get_fastq_F2,
        idx=rules.bwa_mem2_index.output,
    output:
        os.path.join(config["working_dir"], "sams/F2/bwamem2/mapped/{sample}.sam"),
    params:
        index=lambda wildcards: config["ref_prefix"] + ".fasta",
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
    container:
        "docker://quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0"
    threads: 8
    shell:
        """
        bwa-mem2 mem \
            -t {threads} \
            {params.extra} \
            {params.index} \
            {input.reads} \
                > {output}
        """

rule sort_sam_f2:
    input:
        os.path.join(config["working_dir"], "sams/F2/bwamem2/mapped/{sample}.sam"),
    output:
        os.path.join(config["working_dir"], "bams/F2/bwamem2/sorted/{sample}.bam"),
    params:
        sort_order="coordinate",
        extra=lambda wildcards: "VALIDATION_STRINGENCY=LENIENT TMP_DIR=" + config["tmp_dir"],
    container:
        "docker://quay.io/biocontainers/picard:2.9.2--2"
    resources:
        mem_mb=1024
    shell:
        """
        picard SortSam \
            -Xmx{resources.mem_mb}M \
            {params.extra} \
            INPUT={input[0]} \
            OUTPUT={output[0]} \
            SORT_ORDER={params.sort_order}
        """

rule mark_duplicates_f2:
    input:
        os.path.join(config["working_dir"], "bams/F2/bwamem2/sorted/{sample}.bam")
    output:
        bam=os.path.join(config["working_dir"], "bams/F2/bwamem2/marked/{sample}.bam"),
        metrics=os.path.join(config["working_dir"], "bams/F2/bwamem2/marked/{sample}.metrics.txt")
    params:
        extra = lambda wildcards: "REMOVE_DUPLICATES=true TMP_DIR=" + config["tmp_dir"]
    container:
        "docker://quay.io/biocontainers/picard:2.9.2--2"
    resources:
        mem_mb=1024
    shell:
        """
        picard MarkDuplicates \
            -Xmx{resources.mem_mb}M \
            {params.extra} \
            INPUT={input[0]} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metrics}
        """

rule samtools_index_f2:
    input:
        os.path.join(config["working_dir"], "bams/F2/bwamem2/marked/{sample}.bam")
    output:
        os.path.join(config["working_dir"], "bams/F2/bwamem2/marked/{sample}.bam.bai")
    container:
        "docker://quay.io/biocontainers/samtools:0.1.19--h94a8ba4_5"
    shell:
        """
        samtools index \
            {input[0]} \
            {output[0]}
        """
