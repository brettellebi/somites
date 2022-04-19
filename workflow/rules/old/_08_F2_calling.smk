rule haplotype_caller_F2:
    input:
        bam=os.path.join(config["working_dir"], "bams/F2/bwamem2/no_repeat_keep/{F2_sample}.bam"),
        bai=os.path.join(config["working_dir"], "bams/F2/bwamem2/no_repeat_keep/{F2_sample}.bam.bai"),
        ref=config["ref_prefix"] + ".fasta",
        ref_index = config["ref_prefix"] + ".fasta.fai",
        ref_dict = config["ref_prefix"] + ".dict"
    output:
        gvcf=os.path.join(config["working_dir"], "vcfs/F2/gvcfs/{F2_sample}/{contig}.g.vcf"),
    log:
        os.path.join(config["working_dir"], "logs/haplotype_caller_F2/{F2_sample}_{contig}.log")
    params:
        extra = lambda wildcards: "-L " + wildcards.contig + " --tmp-dir " + config["tmp_dir"]
    resources:
        java_mem_gb=4,
        mem_mb=5000
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

rule combine_calls_F2:
    input:
        ref=config["ref_prefix"] + ".fasta",
        gvcfs=expand(
            os.path.join(config["working_dir"], "vcfs/F2/gvcfs/{F2_sample}/{{contig}}.g.vcf"),
            F2_sample = F2_samples['SAMPLE']
        ),
    output:
        os.path.join(config["working_dir"], "vcfs/F2/combined/all.{contig}.g.vcf.gz"),
    log:
        os.path.join(config["working_dir"], "logs/combine_calls/{contig}.log")
    params:
        in_files = lambda wildcards, input: " -V ".join(input.gvcfs)
    resources:
        java_mem_gb=4,
        mem_mb=10000 
    container:
        config["gatk"]
    shell:
        """
        gatk --java-options \"-Xmx{resources.java_mem_gb}G\"\
            CombineGVCFs \
                -V {params.in_files} \
                -R {input.ref} \
                -O {output[0]} \
                    > {log} 2>&1
        """

rule genotype_variants_F2:
    input:
        ref=config["ref_prefix"] + ".fasta",
        gvcf=os.path.join(config["working_dir"], "vcfs/F2/combined/all.{contig}.g.vcf.gz"),
    output:
        os.path.join(config["working_dir"], "vcfs/F2/genotyped/all.{contig}.vcf.gz"),
    log:
        os.path.join(config["working_dir"], "logs/genotype_variants/{contig}.log")
    resources:
        java_mem_gb=1   
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

rule merge_variants_F2:
    input:
        vcfs=expand(
            os.path.join(config["working_dir"], "vcfs/F2/genotyped/all.{contig}.vcf.gz"),
            contig=get_contigs()
        ),
    output:
        os.path.join(config["data_store_dir"], "vcfs/F2/final/all.vcf.gz"),
    log:
        os.path.join(config["working_dir"], "logs/merge_variants/all.log")
    params:
        in_files = lambda wildcards, input: " ".join("INPUT={}".format(f) for f in input.vcfs)
    container:
        config["picard"]
    shell:
        """
        picard MergeVcfs \
            {params.in_files} \
            OUTPUT={output[0]} \
                &> {log}
        """
