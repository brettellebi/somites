rule haplotype_caller_F1:
    input:
        bam=rules.mark_duplicates_F1.output.bam,
        bai=rules.samtools_index_F1.output,
        ref=rules.get_genome.output,
        ref_index = rules.genome_faidx.output,
        ref_dict = rules.genome_dict.output,
    output:
        gvcf=os.path.join(
            config["working_dir"],
            "vcfs/F1/gvcfs/{F1_sample}/{contig}.g.vcf"
        ),
    log:
        os.path.join(
            config["working_dir"], 
            "logs/haplotype_caller/{F1_sample}_{contig}.log"
        ),
    params:
        extra = lambda wildcards: "-L " + wildcards.contig + " --tmp-dir " + config["tmp_dir"]
    resources:
        java_mem_gb=4,
        mem_mb=10000
    container:
        config["gatk"]
    shell:
        """
        gatk --java-options \"-Xmx{resources.java_mem_gb}G\" \
            HaplotypeCaller \
                {params.extra} \
                -R {input.ref} \
                -I {input.bam} \
                -ERC GVCF \
                -O {output.gvcf} \
                    > {log} 2>&1
        """

rule genotype_variants_F1:
    input:
        ref=rules.get_genome.output,
        gvcf=rules.haplotype_caller_F1.output.gvcf,
    output:
        os.path.join(
            config["working_dir"], 
            "vcfs/F1/genotyped/{F1_sample}/{contig}.vcf.gz"
        ),
    log:
        os.path.join(
            config["working_dir"], 
            "logs/genotype_variants/{F1_sample}/{contig}.log"
        ),
    resources:
        java_mem_gb=1,
        mem_mb = 5000
    container:
        config["gatk"]
    shell:
        """
        gatk --java-options \"-Xmx{resources.java_mem_gb}G\" \
            GenotypeGVCFs \
                -V {input.gvcf} \
                -R {input.ref} \
                -O {output[0]} \
                    > {log} 2>&1
        """

rule merge_variants_F1:
    input:
        vcfs=expand(os.path.join(
            config["working_dir"], 
            "vcfs/F1/genotyped/{{F1_sample}}/{contig}.vcf.gz"),
                contig=get_contigs()
        ),
    output:
        os.path.join(
            config["working_dir"],
            "vcfs/F1/final/{F1_sample}.vcf.gz"),
    log:
        os.path.join(
            config["working_dir"], 
            "logs/merge_variants/{F1_sample}.log"
        ),
    params:
        in_files = lambda wildcards, input: " ".join("INPUT={}".format(f) for f in input.vcfs)
    resources:
        mem_mb = 5000
    container:
        config["picard"]
    shell:
        """
        picard MergeVcfs \
            {params.in_files} \
            OUTPUT={output[0]} \
                &> {log}
        """
