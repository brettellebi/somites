rule get_genome:
    output:
        config["ref_prefix"] + ".fasta"
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_genome/get_genome.log"
        ),
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    resources:
        mem_mb = 1000
    wrapper:
        "0.74.0/bio/reference/ensembl-sequence"

rule genome_faidx:
    input:
        config["ref_prefix"] + ".fasta",
    output:
        config["ref_prefix"] + ".fasta.fai",
    log:
        os.path.join(
            config["working_dir"],
            "logs/genome_faidx/genome_faidx.log"
        ),
    container:
        config["samtools"]
    resources:
        mem_mb = 1000
    shell:
        """
        samtools faidx {input[0]}  > {output[0]}
        """

rule genome_dict:
    input:
        config["ref_prefix"] + ".fasta"
    output:
        config["ref_prefix"] + ".dict",
    log:
        os.path.join(
            config["working_dir"],
            "logs/genome_dict/genome_dict.log"
        ),
    container:
        config["samtools"]
    resources:
        mem_mb = 5000
    shell:
        """
        samtools dict {input} > {output}
        """

rule bwa_mem2_index:
    input:
        config["ref_prefix"] + ".fasta",
    output:
        multiext(config["ref_prefix"] + ".fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    log:
        os.path.join(
            config["working_dir"],
            "logs/bwa_mem2_index/bwa_mem2_index.log"
        ),    
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