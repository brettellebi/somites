rule get_genome:
    output:
        config["ref_prefix"] + ".fasta"
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    wrapper:
        "0.74.0/bio/reference/ensembl-sequence"

rule genome_faidx:
    input:
        config["ref_prefix"] + ".fasta",
    output:
        config["ref_prefix"] + ".fasta.fai",
    wrapper:
        "0.74.0/bio/samtools/faidx"

rule genome_dict:
    input:
        config["ref_prefix"] + ".fasta"
    output:
        config["ref_prefix"] + ".dict",
    conda:
        "../envs/samtools_1.12.yaml"
    shell:
        "samtools dict {input} > {output}"

rule bwa_index:
    input:
        config["ref_prefix"] + ".fasta",
    output:
        multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    resources:
        mem_mb=369000,
    wrapper:
        "0.74.0/bio/bwa/index"