rule map_reads:
    input:
        target = config["ref_prefix"] + ".fasta",
        query = get_fastq,
    output:
        os.path.join(config["working_dir"], "sams/F0/mapped/{sample}-{unit}.sam")
    params:
        extra="-ax sr"
    threads: 4
    wrapper:
        "0.74.0/bio/minimap2/aligner"

rule replace_rg:
    input:
        os.path.join(config["working_dir"], "sams/F0/mapped/{sample}-{unit}.sam")
    output:
        os.path.join(config["working_dir"], "sams/F0/grouped/{sample}-{unit}.sam")
    params:
        "RGLB=lib1 RGPL=ILLUMINA RGPU={unit} RGSM={sample}"
    resources:
        mem_mb=1024
    wrapper:
        "0.74.0/bio/picard/addorreplacereadgroups"

rule sort_sam:
    input:
        os.path.join(config["working_dir"], "sams/F0/grouped/{sample}-{unit}.sam")
    output:
        os.path.join(config["working_dir"], "bams/F0/sorted/{sample}-{unit}.bam")
    params:
        sort_order="coordinate",
        extra=lambda wildcards: "VALIDATION_STRINGENCY=LENIENT TMP_DIR=" + config["tmp_dir"] # optional: Extra arguments for picard.
    resources:
        mem_mb=1024
    wrapper:
        "0.74.0/bio/picard/sortsam"

rule mark_duplicates:
    input:
        os.path.join(config["working_dir"], "bams/F0/sorted/{sample}-{unit}.bam")
    output:
        bam=os.path.join(config["working_dir"], "bams/F0/marked/{sample}-{unit}.bam"),
        metrics=os.path.join(config["working_dir"], "bams/F0/marked/{sample}-{unit}.metrics.txt")
    params:
        lambda wildcards: "REMOVE_DUPLICATES=true TMP_DIR=" + config["tmp_dir"]
    resources:
        mem_mb=1024
    wrapper:
        "0.74.0/bio/picard/markduplicates"

rule merge_bams:
    input:
        expand(os.path.join(config["working_dir"], "bams/F0/marked/{{sample}}-{unit}.bam"),
            unit = list(set(F0_samples['unit']))
        )
    output:
        os.path.join(config["working_dir"], "bams/F0/merged/{sample}.bam")
    params:
        lambda wildcards: "VALIDATION_STRINGENCY=LENIENT TMP_DIR=" + config["tmp_dir"]
    resources:
        mem_mb=1024
    wrapper:
        "0.74.0/bio/picard/mergesamfiles"

rule samtools_index:
    input:
        os.path.join(config["working_dir"], "bams/F0/merged/{sample}.bam")
    output:
        os.path.join(config["working_dir"], "bams/F0/merged/{sample}.bam.bai")
    wrapper:
        "0.74.0/bio/samtools/index"
