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
#    log:
#        os.path.join(config["working_dir"], "log/replace_rg/{sample}_{unit}.log")
    params:
        "RGLB=lib1 RGPL=ILLUMINA RGPU={unit} RGSM={sample}"
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024
    wrapper:
        "0.74.0/bio/picard/addorreplacereadgroups"

rule sort_sam:
    input:
        os.path.join(config["working_dir"], "sams/F0/grouped/{sample}-{unit}.sam")
    output:
        os.path.join(config["working_dir"], "bams/F0/sorted/{sample}-{unit}.bam")
#    log:
#        os.path.join(config["working_dir"], "log/sort_sam/{sample}-{unit}.log")
    params:
        sort_order="coordinate",
        extra=lambda wildcards: "VALIDATION_STRINGENCY=LENIENT TMP_DIR=" + config["tmp_dir"] # optional: Extra arguments for picard.
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024
    wrapper:
        "0.74.0/bio/picard/sortsam"

rule mark_duplicates:
    input:
        os.path.join(config["working_dir"], "bams/F0/sorted/{sample}-{unit}.bam")
    # optional to specify a list of BAMs; this has the same effect
    # of marking duplicates on separate read groups for a sample
    # and then merging
    output:
        bam=os.path.join(config["working_dir"], "bams/F0/marked/{sample}-{unit}.bam"),
        metrics=os.path.join(config["working_dir"], "bams/F0/marked/{sample}-{unit}.metrics.txt")
#    log:
#        os.path.join(config["working_dir"], "log/mark_duplicates/{sample}-{unit}.log")
    params:
        lambda wildcards: "REMOVE_DUPLICATES=true TMP_DIR=" + config["tmp_dir"]
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024
    wrapper:
        "0.74.0/bio/picard/markduplicates"

rule merge_bams:
    input:
        expand(os.path.join(config["working_dir"], "bams/F0/marked/{{sample}}-{unit}.bam"),
            unit = list(set(samples['unit']))
        )
    output:
        os.path.join(config["working_dir"], "bams/F0/merged/{sample}.bam")
#    log:
#        os.path.join(config["working_dir"], "log/picard_mergesamfiles/{sample}.log")
    params:
        lambda wildcards: "VALIDATION_STRINGENCY=LENIENT TMP_DIR=" + config["tmp_dir"]
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024
    wrapper:
        "0.74.0/bio/picard/mergesamfiles"

rule samtools_index:
    input:
        os.path.join(config["working_dir"], "bams/F0/merged/{sample}.bam")
    output:
        os.path.join(config["working_dir"], "bams/F0/merged/{sample}.bam.bai")
#    params:
#        "" # optional params string
    wrapper:
        "0.74.0/bio/samtools/index"
