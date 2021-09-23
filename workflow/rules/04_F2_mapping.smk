rule copy_f2_seq_data:
    input:
        lambda wildcards: F2_samples.loc[wildcards.F2_sample, ["fq1", "fq2"]].dropna().tolist(),
    output:
        expand(os.path.join(config["working_dir"], "fastqs/F2/{{F2_sample}}_{pair}.txt.gz"),
            F2_sample = F2_samples['SAMPLE'],
            pair = config["pairs"]
        ),
    log:
        os.path.join(config["working_dir"], "logs/copy_f2_seq_data/{F2_sample}.log")
    shell:
        """
        cp {input[0]} {output[0]} ;
        cp {input[1]} {output[1]}
        """

rule bwa_mem2_mem:
    input:
        reads=expand(os.path.join(config["working_dir"], "fastqs/F2/{{F2_sample}}_{pair}.txt.gz"),
            F2_sample = F2_samples['SAMPLE'],
            pair = config["pairs"]
        ),
        idx=rules.bwa_mem2_index.output,
    output:
        os.path.join(config["working_dir"], "sams/F2/bwamem2/mapped/{F2_sample}.sam"),
    log:
        os.path.join(config["working_dir"], "logs/bwa_mem2_mem/{F2_sample}.log")
    params:
        index=lambda wildcards: config["ref_prefix"] + ".fasta",
        extra=r"-R '@RG\tID:{F2_sample}\tSM:{F2_sample}'",
    container:
        config["bwa-mem2"]
    resources:
        mem_mb = 15000
    threads:
        8
    shell:
        """
        bwa-mem2 mem \
            -t {threads} \
            {params.extra} \
            {params.index} \
            {input.reads} \
                > {output} \
                    2> {log}
        """
# One error for sample (LANE) 7171
# log file here: /hps/nobackup/birney/users/ian/somites/logs/bwa_mem2_mem/2604307_sample\=7171.err
# [mem_sam_pe_batch_post] paired reads have different names: "�0�����O&�����x�7�Z�HĻ�@�%Q�ü��", "ST-K00119:220:HKNVLBBXY:7:1101:1661:1138"
# That sample is hashed out of `config/F2_samples.txt` until the issue is resolved.
rule sort_sam_f2:
    input:
        os.path.join(config["working_dir"], "sams/F2/bwamem2/mapped/{F2_sample}.sam"),
    output:
        os.path.join(config["working_dir"], "bams/F2/bwamem2/sorted/{F2_sample}.bam"),
    log:
        os.path.join(config["working_dir"], "logs/sort_sam_f2/{F2_sample}.log")        
    params:
        sort_order="coordinate",
        extra=lambda wildcards: "VALIDATION_STRINGENCY=LENIENT TMP_DIR=" + config["tmp_dir"],
    container:
        config["picard"]
    resources:
        java_mem_mb = 1024,
        mem_mb = 20000
    shell:
        """
        picard SortSam \
            -Xmx{resources.java_mem_mb}M \
            {params.extra} \
            INPUT={input[0]} \
            OUTPUT={output[0]} \
            SORT_ORDER={params.sort_order} \
                2> {log}
        """
rule mark_duplicates_f2:
    input:
        os.path.join(config["working_dir"], "bams/F2/bwamem2/sorted/{F2_sample}.bam")
    output:
        bam=os.path.join(config["working_dir"], "bams/F2/bwamem2/marked/{F2_sample}.bam"),
        metrics=os.path.join(config["working_dir"], "bams/F2/bwamem2/marked/{F2_sample}.metrics.txt")
    log:
        os.path.join(config["working_dir"], "logs/mark_duplicates_f2/{F2_sample}.log") 
    params:
        extra = lambda wildcards: "REMOVE_DUPLICATES=true TMP_DIR=" + config["tmp_dir"]
    container:
        config["picard"]
    resources:
        mem_mb=1024
    shell:
        """
        picard MarkDuplicates \
            -Xmx{resources.mem_mb}M \
            {params.extra} \
            INPUT={input[0]} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metrics} \
                2> {log}
        """

rule samtools_index_f2:
    input:
        os.path.join(config["working_dir"], "bams/F2/bwamem2/marked/{F2_sample}.bam")
    output:
        os.path.join(config["working_dir"], "bams/F2/bwamem2/marked/{F2_sample}.bam.bai")
    log:
        os.path.join(config["working_dir"], "logs/samtools_index_f2/{F2_sample}.log") 
    container:
        config["samtools"]
    shell:
        """
        samtools index \
            {input[0]} \
            {output[0]} \
                2> {log}
        """
