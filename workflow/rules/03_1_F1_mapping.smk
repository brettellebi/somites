rule copy_F1_seq_data:
    input:
        lambda wildcards: F1_samples.loc[wildcards.F1_sample, ["fq1", "fq2"]].dropna().tolist(),
    output:
        expand(os.path.join(config["working_dir"], "fastqs/F1/{{F1_sample}}_{pair}.txt.gz"),
            F1_sample = F1_samples['SAMPLE'],
            pair = config["pairs"]
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/copy_F1_seq_data/{F1_sample}.log")
    resources:
        mem_mb = 100
    shell:
        """
        cp {input[0]} {output[0]} ;
        cp {input[1]} {output[1]}
        """

rule bwa_mem2_mem_F1:
    input:
        reads=expand(rules.copy_F1_seq_data.output,
            F1_sample = F1_samples['SAMPLE'],
            pair = config["pairs"]
        ),
        idx=rules.bwa_mem2_index.output,
    output:
        os.path.join(
            config["working_dir"],
            "sams/F1/bwamem2/mapped/{F1_sample}.sam"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/bwa_mem2_mem/{F1_sample}.log")
    params:
        index=lambda wildcards: config["ref_prefix"] + ".fasta",
        extra=r"-R '@RG\tID:{F1_sample}\tSM:{F1_sample}'",
    container:
        config["bwa-mem2"]
    resources:
        mem_mb = 10000
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

rule sort_sam_F1:
    input:
        os.path.join(
            config["working_dir"], 
            "sams/F1/bwamem2/mapped/{F1_sample}.sam"
        ),
    output:
        os.path.join(
            config["working_dir"], 
            "bams/F1/bwamem2/sorted/{F1_sample}.bam"
        ),
    log:
        os.path.join(
            config["working_dir"], 
            "logs/sort_sam_F1/{F1_sample}.log"
        )        
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

rule mark_duplicates_F1:
    input:
        os.path.join(
            config["working_dir"], 
            "bams/F1/bwamem2/sorted/{F1_sample}.bam"
        )
    output:
        bam=os.path.join(
            config["working_dir"], 
            "bams/F1/bwamem2/marked/{F1_sample}.bam"
        ),
        metrics=os.path.join(
            config["working_dir"], 
            "bams/F1/bwamem2/marked/{F1_sample}.metrics.txt"
        )
    log:
        os.path.join(
            config["working_dir"], 
            "logs/mark_duplicates_F1/{F1_sample}.log"
        ) 
    params:
        extra = lambda wildcards: "REMOVE_DUPLICATES=true TMP_DIR=" + config["tmp_dir"]
    container:
        config["picard"]
    resources:
        java_mem_mb = 1024,
        mem_mb=5000
    shell:
        """
        picard MarkDuplicates \
            -Xmx{resources.java_mem_mb}M \
            {params.extra} \
            INPUT={input[0]} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metrics} \
                2> {log}
        """

rule samtools_index_F1:
    input:
        os.path.join(
            config["working_dir"], 
            "bams/F1/bwamem2/marked/{F1_sample}.bam"
        )
    output:
        os.path.join(
            config["working_dir"], 
            "bams/F1/bwamem2/marked/{F1_sample}.bam.bai"
        )
    log:
        os.path.join(
            config["working_dir"], 
            "logs/samtools_index_F1/{F1_sample}.log"
        ) 
    container:
        config["samtools"]
    resources:
        mem_mb = 5000
    shell:
        """
        samtools index \
            {input[0]} \
            {output[0]} \
                2> {log}
        """
