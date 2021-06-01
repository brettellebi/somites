rule get_genome:
    output:
        config["ref_path"]
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
#    cache: True
    wrapper:
        "0.74.0/bio/reference/ensembl-sequence"

checkpoint genome_faidx:
    input:
        config["ref_path"],
    output:
        config["ref_path"] + ".fai",
#    cache: True
    wrapper:
        "0.74.0/bio/samtools/faidx"

