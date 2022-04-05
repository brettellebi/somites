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

# Extract all homozygous sites per F0 sample
rule get_genos_F1:
    input:
        rules.merge_variants_F1.output,
    output:
        os.path.join(
            config["working_dir"],
            "raw_genos/F1/{F1_sample}.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_genos_F1/{F1_sample}.log"
        ),
    resources:
        mem_mb = 5000
    container:
        config["bcftools_1.14"]
    shell:
        """
        bcftools view \
            --regions 22 \
            --max-alleles 2 \
            --output-type u \
            {input} |\
        bcftools query \
            --format '%CHROM,%POS,%REF,%ALT,[%GT]\\n' \
            --output {output} \
                2> {log}
        """


###############
# Call F0s and F1 as trio
###############

# Create two vectors with GEN and SAMPLE
import numpy as np

F0_gens = np.repeat("F0", len(config["F0_lines"]))
F1_gens = np.repeat("F1", len(F1_samples['SAMPLE']))
ALL_GENS = [y for x in [F0_gens, F1_gens] for y in x]
ALL_SAMPLES = [y for x in [config["F0_lines"], F1_samples['SAMPLE']] for y in x]

rule combine_calls_F0_and_F1:
    input:
        ref=rules.get_genome.output,
        gvcfs=expand(
            os.path.join(
                config["working_dir"],
                "vcfs/{all_gens}/gvcfs/{all_samples}/{{contig}}.g.vcf"),
            zip,
            all_gens = ALL_GENS,
            all_samples = ALL_SAMPLES
        ),
    output:
        os.path.join(
            config["working_dir"],
            "vcfs/F0_and_F1/combined/all.{contig}.g.vcf.gz"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/combine_calls_F0_and_F1/{contig}.log"
        ),
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

rule genotype_variants_F0_and_F1:
    input:
        ref=rules.get_genome.output,
        gvcf=rules.combine_calls_F0_and_F1.output,
    output:
        os.path.join(
            config["working_dir"], 
            "vcfs/F0_and_F1/genotyped/all.{contig}.vcf.gz"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/genotype_variants_F0_and_F1/{contig}.log"
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

rule merge_variants_F0_and_F1:
    input:
        vcfs=expand(
            rules.genotype_variants_F0_and_F1.output,
            contig=get_contigs()
        ),
    output:
        os.path.join(
            config["working_dir"],
            "vcfs/F0_and_F1/final/all.vcf.gz"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/merge_variants_F0_and_F1/all.log"
        ),
    params:
        in_files = lambda wildcards, input: " ".join("INPUT={}".format(f) for f in input.vcfs)
    resources:
        mem_mb = 5000
    container:
        config["picard"],
    shell:
        """
        picard MergeVcfs \
            {params.in_files} \
            OUTPUT={output[0]} \
                &> {log}
        """

rule extract_trio_genos:
    input:
        rules.merge_variants_F0_and_F1.output
    output:
        os.path.join(
            config["working_dir"],
            "genos/F0_and_F1/final/all.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/extract_trio_genos/all.log"
        ),
    resources:
        mem_mb = 2000
    container:
        config["bcftools_1.14"]
    shell:
        """
        bcftools view \
            --max-alleles 2 \
            --output-type u \
            {input} |\
        bcftools query \
            --print-header \
            --format '%CHROM,%POS,%REF,%ALT[,%GT]\\n' \
            --output {output} \
                2> {log}
        """
