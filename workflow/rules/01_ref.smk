rule get_genome:
    output:
        config["ref_prefix"] + ".fasta"
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
#    cache: True
    wrapper:
        "0.74.0/bio/reference/ensembl-sequence"

rule genome_faidx:
    input:
        config["ref_prefix"] + ".fasta",
    output:
        config["ref_prefix"] + ".fasta.fai",
#    cache: True
    wrapper:
        "0.74.0/bio/samtools/faidx"

rule genome_dict:
    input:
        config["ref_prefix"] + ".fasta"
    output:
        config["ref_prefix"] + ".dict",
    conda:
        "../envs/samtools_1.12.yaml"
#    cache: True
    shell:
        "samtools dict {input} > {output}"