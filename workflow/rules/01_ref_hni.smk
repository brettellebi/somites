rule get_genome_hni:
    output:
        os.path.join(
            config["ref_dir"],
            "{ref}.fasta"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_genome_hni/get_genome.log"
        ),
    params:
        species=config["ref_hni"]["species"],
        datatype="dna",
        build=config["ref_hni"]["build"],
        release=config["ref_hni"]["release"],
    resources:
        mem_mb = 100
    wrapper:
        "v1.3.2/bio/reference/ensembl-sequence"

rule genome_faidx_hni:
    input:
        rules.get_genome_hni.output,
    output:
        config["ref_prefix_hni"] + ".fasta.fai",
    log:
        os.path.join(
            config["working_dir"],
            "logs/genome_faidx_hni/genome_faidx_hni.log"
        ),
    container:
        config["samtools"]
    resources:
        mem_mb = 100
    shell:
        """
        samtools faidx {input[0]}  > {output[0]}
        """

rule genome_dict_hni:
    input:
        rules.get_genome_hni.output,
    output:
        config["ref_prefix_hni"] + ".dict",
    log:
        os.path.join(
            config["working_dir"],
            "logs/genome_dict_hni/genome_dict_hni.log"
        ),
    container:
        config["samtools"]
    resources:
        mem_mb = 100
    shell:
        """
        samtools dict {input} > {output}
        """

rule bwa_mem2_index_hni:
    input:
        config["ref_prefix_hni"] + ".fasta",
    output:
        multiext(config["ref_prefix_hni"] + ".fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    log:
        os.path.join(
            config["working_dir"],
            "logs/bwa_mem2_index_hni/bwa_mem2_index_hni.log"
        ),    
    params:
        prefix=lambda w, input: input[0]
    resources:
        mem_mb=20000,
    container:
        config["bwa-mem2"]
    shell:
        """
        bwa-mem2 index {params.prefix} {input[0]}
        """

rule get_chrom_lengths_hni:
    input:
        rules.get_genome_hni.output
    output:
        "config/hni_chrom_lengths.csv"
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_chrom_lengths_hni/all.log"
        ),
    container:
        config["bash"]
    resources:
        mem_mb = 200
    shell:
        """
        grep "^>" {input} | cut -f1 -d' ' | sed 's/>//g' > tmp1.txt ;
        grep "^>" {input} | cut -f3 -d' ' | cut -f5 -d':' > tmp2.txt ;
        paste -d',' tmp1.txt tmp2.txt > {output} ;
        rm tmp1.txt tmp2.txt \
            2> {log}
        """
    