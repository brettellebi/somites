rule haplotype_caller:
    input:
        bam=os.path.join(config["working_dir"], "bams/F0/merged/{sample}.bam"),
        bai=os.path.join(config["working_dir"], "bams/F0/merged/{sample}.bam.bai"),
        ref=config["ref_prefix"] + ".fasta",
        ref_index = config["ref_prefix"] + ".fasta.fai",
        ref_dict = config["ref_prefix"] + ".dict"
    output:
        gvcf=os.path.join(config["working_dir"], "vcfs/F0/gvcfs/{sample}/{contig}.g.vcf"),
    params:
        extra= lambda wildcards: "-L " + wildcards.contig + " --tmp-dir " + config["tmp_dir"] # optional
    resources:
        mem_mb=1024
    wrapper:
        "0.74.0/bio/gatk/haplotypecaller"

rule combine_calls:
    input:
        ref=config["ref_prefix"] + ".fasta",
        gvcfs=expand(
            os.path.join(config["working_dir"], "vcfs/F0/gvcfs/{sample}/{{contig}}.g.vcf"),
            sample=list(set(samples["sample"]))
        ),
    output:
        gvcf=os.path.join(config["working_dir"], "vcfs/F0/combined/all.{contig}.g.vcf.gz"),
    wrapper:
        "0.74.0/bio/gatk/combinegvcfs"

rule genotype_variants:
    input:
        ref=config["ref_prefix"] + ".fasta",
        gvcf=os.path.join(config["working_dir"], "vcfs/F0/combined/all.{contig}.g.vcf.gz"),
    output:
        vcf=os.path.join(config["working_dir"], "vcfs/F0/genotyped/all.{contig}.vcf.gz"),
    wrapper:
        "0.74.0/bio/gatk/genotypegvcfs"

rule merge_variants:
    input:
        vcfs=expand(
            os.path.join(config["working_dir"], "vcfs/F0/genotyped/all.{contig}.vcf.gz"),
            contig=get_contigs()
        ),
    output:
        vcf=os.path.join(config["working_dir"], "vcfs/F0/final/all.vcf.gz"),
    wrapper:
        "0.74.0/bio/picard/mergevcfs"