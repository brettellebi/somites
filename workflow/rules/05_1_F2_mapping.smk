rule copy_F2_seq_data:
    input:
        lambda wildcards: F2_samples.loc[wildcards.F2_sample, ["fq1", "fq2"]].dropna().tolist(),
    output:
        expand(os.path.join(config["working_dir"], "fastqs/F2/{{F2_sample}}_{pair}.txt.gz"),
            F2_sample = F2_samples['SAMPLE'],
            pair = config["pairs"]
        ),
    resources:
        mem_mb = 1000
    log:
        os.path.join(config["working_dir"], "logs/copy_f2_seq_data/{F2_sample}.log")
    shell:
        """
        cp {input[0]} {output[0]} ;
        cp {input[1]} {output[1]}
        """

rule bwa_mem2_mem_F2:
    input:
        reads=expand(os.path.join(
            config["working_dir"],
            "fastqs/F2/{{F2_sample}}_{pair}.txt.gz"),
                F2_sample = F2_samples['SAMPLE'],
                pair = config["pairs"]
        ),
        idx=rules.bwa_mem2_index.output,
        ref = rules.get_genome.output,
    output:
        os.path.join(
            config["working_dir"],
            "sams/F2/{ref}/bwamem2/mapped/{F2_sample}.sam"
        ),
    log:
        os.path.join(
            config["working_dir"], 
            "logs/bwa_mem2_mem_F2/{ref}/{F2_sample}.log"
        ),
    params:
        extra=r"-R '@RG\tID:{F2_sample}\tSM:{F2_sample}'",
    container:
        config["bwa-mem2"]
    resources:
        mem_mb = 10000
    threads:
        1
    shell:
        """
        bwa-mem2 mem \
            -t {threads} \
            {params.extra} \
            {input.ref} \
            {input.reads} \
                > {output} \
                    2> {log}
        """

rule sort_sam_F2:
    input:
        rules.bwa_mem2_mem_F2.output,
    output:
        os.path.join(
            config["working_dir"],
            "bams/F2/bwamem2/sorted/{F2_sample}.bam"
        ),
    log:
        os.path.join(
            config["working_dir"], 
            "logs/sort_sam_F2/{F2_sample}.log"
        ),
    params:
        sort_order="coordinate",
        extra=lambda wildcards: "VALIDATION_STRINGENCY=LENIENT TMP_DIR=" + config["tmp_dir"],
    container:
        config["picard"]
    resources:
        java_mem_mb = 1024,
        mem_mb = 2000
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
rule mark_duplicates_F2:
    input:
        rules.sort_sam_F2.output,
    output:
        bam=os.path.join(
            config["working_dir"], 
            "bams/F2/{ref}/bwamem2/marked/{F2_sample}.bam"
        ),
        metrics=os.path.join(
            config["working_dir"], 
            "bams/F2/{ref}/bwamem2/marked/{F2_sample}.metrics.txt"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/mark_duplicates_F2/{ref}/{F2_sample}.log"
        ),
    params:
        extra = lambda wildcards: "REMOVE_DUPLICATES=true TMP_DIR=" + config["tmp_dir"]
    container:
        config["picard"]
    resources:
        java_mem_mb=1024,
        mem_mb=2000
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

rule samtools_index_F2:
    input:
        rules.mark_duplicates_F2.output,
    output:
        os.path.join(
            config["working_dir"],
            "bams/F2/{ref}/bwamem2/marked/{F2_sample}.bam.bai"
        ),
    log:
        os.path.join(
            config["working_dir"], 
            "logs/samtools_index_F2/{ref}/{F2_sample}.log"
        ),
    resources:
        mem_mb = 100
    container:
        config["samtools"]
    shell:
        """
        samtools index \
            {input[0]} \
            {output[0]} \
                2> {log}
        """
