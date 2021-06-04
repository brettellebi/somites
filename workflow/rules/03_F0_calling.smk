rule haplotype_caller:
    input:
        # single or list of bam files
        bam=os.path.join(config["working_dir"], "bams/F0/merged/{sample}.bam"),
        bai=os.path.join(config["working_dir"], "bams/F0/merged/{sample}.bam.bai"),
        ref=config["ref_prefix"] + ".fasta",
        ref_index = config["ref_prefix"] + ".fasta.fai",
        ref_dict = config["ref_prefix"] + ".dict"
        # known="dbsnp.vcf"  # optional
    output:
        gvcf=os.path.join(config["working_dir"], "vcfs/F0/gvcfs/{sample}/{contig}.g.vcf"),
#   bam="{sample}.assemb_haplo.bam",
#    log:
#        os.path.join(config["working_dir"], "logs/haplotypecaller/{sample}.log")
    params:
        extra= lambda wildcards: "-L " + wildcards.contig + " --tmp-dir " + config["tmp_dir"] # optional
#        java_opts="", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
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
#    log:
#        os.path.join(config["working_dir"], "logs/combine_calls/{contig}.log"),
    wrapper:
        "0.74.0/bio/gatk/combinegvcfs"

rule genotype_variants:
    input:
        ref=config["ref_prefix"] + ".fasta",
        gvcf=os.path.join(config["working_dir"], "vcfs/F0/combined/all.{contig}.g.vcf.gz"),
    output:
        vcf=os.path.join(config["working_dir"], "vcfs/F0/genotyped/all.{contig}.vcf.gz"),
#    params:
#        extra=config["params"]["gatk"]["GenotypeGVCFs"],
#    log:
#        os.path.join(config["working_dir"], "logs/genotype_variants/{contig}.log"),
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
#    log:
#        os.path.join(config["working_dir"], "logs/merge_variants/F0.log"),
    wrapper:
        "0.74.0/bio/picard/mergevcfs"