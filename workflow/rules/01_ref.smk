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
    container:
        config["samtools"]
    shell:
        """
        samtools faidx {input[0]}  > {output[0]}
        """

rule genome_dict:
    input:
        config["ref_prefix"] + ".fasta"
    output:
        config["ref_prefix"] + ".dict",
    container:
        config["samtools"]
    shell:
        """
        samtools dict {input} > {output}
        """

rule bwa_mem2_index:
    input:
        config["ref_prefix"] + ".fasta",
    output:
        multiext(config["ref_prefix"] + ".fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    params:
        prefix=lambda w, input: input[0]
    resources:
        mem_mb=50000,
    container:
        config["bwa-mem2"]
    shell:
        """
        bwa-mem2 index {params.prefix} {input[0]}
        """