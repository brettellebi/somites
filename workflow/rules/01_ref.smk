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

#rule bwa_index:
#    input:
#        config["ref_prefix"] + ".fasta",
#    output:
#        multiext(config["ref_prefix"] + ".fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
#    log:
#        "logs/bwa_index.log",
#    resources:
#        mem_mb=369000,
#    wrapper:
#        "0.74.0/bio/bwa/index"

rule bwa_mem2_index:
    input:
        config["ref_prefix"] + ".fasta",
    output:
        multiext(config["ref_prefix"] + ".fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    params:
        prefix=lambda w, input: input[0]
    resources:
        mem_mb=369000,
    wrapper:
        "v0.75.0/bio/bwa-mem2/index"