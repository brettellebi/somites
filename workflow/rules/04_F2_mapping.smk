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

rule copy_f2_seq_data:
    input:
        get_fastq_F2,
    output:
        expand(os.path.join(config["working_dir"], "fastqs/F2/{{sample}}_{pair}.txt.gz"),
            pair = PAIRS
        ),
    shell:
        """
        cp {input[0]} {output[0]} ;
        cp {input[1]} {output[1]}
        """

rule bwa_mem2_mem:
    input:
        reads=expand(os.path.join(config["working_dir"], "fastqs/F2/{{sample}}_{pair}.txt.gz"),
            pair = PAIRS
        ),
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
# One error for sample (LANE) 7171
# log file here: /hps/nobackup/birney/users/ian/somites/logs/bwa_mem2_mem/2604307_sample\=7171.err
# [mem_sam_pe_batch_post] paired reads have different names: "�0�����O&�����x�7�Z�HĻ�@�%Q�ü��", "ST-K00119:220:HKNVLBBXY:7:1101:1661:1138"
# That sample is hashed out of `config/F2_samples.txt` until the issue is resolved.

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
