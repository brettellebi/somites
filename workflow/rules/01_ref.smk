rule get_genome:
    output:
        os.path.join(
            config["ref_dir"],
            "{ref}.fasta"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_genome/{ref}.log"
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
        rules.get_genome.output,
    output:
        os.path.join(
            config["ref_dir"],
            "{ref}.fasta.fai"
        )
    log:
        os.path.join(
            config["working_dir"],
            "logs/genome_faidx/{ref}.log"
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
        rules.get_genome.output,
    output:
        os.path.join(
            config["ref_dir"],
            "{ref}.dict"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/genome_dict/{ref}.log"
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
        rules.get_genome.output,
    output:
        multiext(
            os.path.join(
                config["ref_dir"],
                "{ref}.fasta"),
            ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    log:
        os.path.join(
            config["working_dir"],
            "logs/bwa_mem2_index/{ref}.log"
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

rule get_chrom_lengths:
    input:
        rules.get_genome.output
    output:
        csv = "config/{ref}_chrom_lengths.csv"
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_chrom_lengths/{ref}.log"
        ),
    container:
        config["bash"]
    resources:
        mem_mb = 200
    shell:
        """
        grep "^>" {input} | cut -f1 -d' ' | sed 's/>//g' > tmp1.txt ;
        grep "^>" {input} | cut -f3 -d' ' | cut -f5 -d':' > tmp2.txt ;
        paste -d',' tmp1.txt tmp2.txt > {output.csv} ;
        rm tmp1.txt tmp2.txt \
            2> {log}
        """
    