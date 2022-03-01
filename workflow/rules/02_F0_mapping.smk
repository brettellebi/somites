rule map_reads_F0:
    input:
        target = config["ref_prefix"] + ".fasta",
        query = lambda wildcards: F0_samples.loc[(wildcards.F0_sample, wildcards.unit), ["fq1", "fq2"]].dropna().tolist(),
    output:
        os.path.join(
            config["working_dir"],
            "sams/F0/mapped/{F0_sample}-{unit}.sam"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/map_reads_F0/{F0_sample}/{unit}.log"
        ),
    params:
        extra="-ax sr",
    resources:
        mem_mb = 10000
    threads:
        8
    container:
        config["minimap2"]
    shell:
        """
        minimap2 \
            -t {threads} \
            {params.extra} \
            {input.target} \
            {input.query} \
                > {output[0]}
        """

rule replace_rg_F0:
    input:
        os.path.join(
            config["working_dir"],
            "sams/F0/mapped/{F0_sample}-{unit}.sam"
        ),
    output:
        os.path.join(
            config["working_dir"],
            "sams/F0/grouped/{F0_sample}-{unit}.sam"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/replace_rg_F0/{F0_sample}_{unit}.log"
        ),
    params:
        "RGLB=lib1 RGPL=ILLUMINA RGPU={unit} RGSM={F0_sample}"
    resources:
        mem_mb=1024
    container:
        config["picard"]
    shell:
        """
        picard AddOrReplaceReadGroups \
            -Xmx{resources.mem_mb}M \
            {params} \
            INPUT={input[0]} \
            OUTPUT={output[0]} \
                &> {log}
        """

rule sort_sam_F0:
    input:
        os.path.join(
            config["working_dir"],
            "sams/F0/grouped/{F0_sample}-{unit}.sam"
        ),
    output:
        os.path.join(
            config["working_dir"],
            "bams/F0/sorted/{F0_sample}-{unit}.bam"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/sort_sam_F0/{F0_sample}_{unit}.log"
        ),
    params:
        sort_order="coordinate",
        extra=lambda wildcards: "VALIDATION_STRINGENCY=LENIENT TMP_DIR=" + config["tmp_dir"]
    resources:
        java_mem_mb = 4096,
        mem_mb = 20000
    container:
        config["picard"]
    shell:
        """
        picard SortSam \
            -Xmx{resources.java_mem_mb}M \
            {params.extra} \
            INPUT={input[0]} \
            OUTPUT={output[0]} \
            SORT_ORDER={params.sort_order} \
                &> {log}
        """

rule mark_duplicates_F0:
    input:
        os.path.join(config["working_dir"], "bams/F0/sorted/{F0_sample}-{unit}.bam")
    output:
        bam=os.path.join(config["working_dir"], "bams/F0/marked/{F0_sample}-{unit}.bam"),
        metrics=os.path.join(config["working_dir"], "bams/F0/marked/{F0_sample}-{unit}.metrics.txt")
    log:
        os.path.join(config["working_dir"], "logs/mark_duplicates_F0/{F0_sample}_{unit}.log")
    params:
        lambda wildcards: "REMOVE_DUPLICATES=true TMP_DIR=" + config["tmp_dir"]
    resources:
        java_mem_mb=1024,
        mem_mb=10000
    container:
        config["picard"]
    shell:
        """
        picard MarkDuplicates \
            -Xmx{resources.java_mem_mb}M \
            {params} \
            INPUT={input[0]} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metrics} \
                &> {log}
        """

rule merge_bams_F0:
    input:
        expand(os.path.join(config["working_dir"], "bams/F0/marked/{{F0_sample}}-{unit}.bam"),
            unit = list(set(F0_samples['unit']))
        )
    output:
        os.path.join(config["working_dir"], "bams/F0/merged/{F0_sample}.bam")
    log:
        os.path.join(config["working_dir"], "logs/merge_bams_F0/{F0_sample}.log")
    params:
        extra = lambda wildcards: "VALIDATION_STRINGENCY=LENIENT TMP_DIR=" + config["tmp_dir"],
        in_files = lambda wildcards, input: " I=".join(input)
    resources:
        java_mem_mb=1024,
        mem_mb = 5000
    container:
        config["picard"]
    shell:
        """
        picard MergeSamFiles \
            -Xmx{resources.java_mem_mb}M \
            {params.extra} \
            INPUT={params.in_files} \
            OUTPUT={output} \
                &> {log}
        """

rule samtools_index_F0:
    input:
        os.path.join(config["working_dir"], "bams/F0/merged/{F0_sample}.bam")
    output:
        os.path.join(config["working_dir"], "bams/F0/merged/{F0_sample}.bam.bai")
    log:
        os.path.join(config["working_dir"], "logs/samtools_index_F0/{F0_sample}.log")
    container:
        config["samtools"]
    resources:
        mem_mb = 5000
    shell:
        """
        samtools index \
            {input[0]} \
            {output[0]} \
                &> {log}
        """
