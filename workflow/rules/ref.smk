
######################
# Libraries
######################

import os

######################
# Rules
######################

rule get_genome:
    output:
        config["ref_path"]
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    wrapper:
        "0.74.0/bio/reference/ensembl-sequence"

checkpoint genome_faidx:
    input:
        config["ref_path"],
    output:
        config["ref_path"] + ".fai",
    cache: True
    wrapper:
        "0.74.0/bio/samtools/faidx"

rule align_to_reference:
    input:
        expand(os.path.join(config["fasta_dir"], "{{kaga_sample}}_{read}.fastq.gz"),
                read = READS)
    output:
        os.path.join(config["sam_dir"], "{kaga_sample}.sam")
    conda:
        "envs/minimap2_2.19.yaml"
    shell:
        """
        /minimap2-2.17_x64-linux/minimap2 -ax sr {config[hdrr_reference]} {input} > {output}
        """

rule add_kaga_read_groups:
    input:
        os.path.join(config["sam_dir"], "{kaga_sample}.sam")
    output:
        os.path.join(config["sam_dir"], "{kaga_sample}_withRGs.sam")
    params:
        id = lambda wildcards: wildcards.kaga_sample
    singularity:
        config["picard"]
    shell:
        """
        java -Xmx12g -jar /usr/picard/picard.jar AddOrReplaceReadGroups \
            I={input} \
            O={output} \
            RGID={params.id} \
            RGLB=lib1 \
            RGPL=ILLUMINA \
            RGPU=unit1 \
            RGSM=kaga
        """

rule sort_kaga_sam:
    input:
        os.path.join(config["sam_dir"], "{kaga_sample}_withRGs.sam")
    output:
        os.path.join(config["bam_dir"], "sorted/{kaga_sample}.bam")
    singularity:
        config["picard"]
    shell:
        """
        java -Xmx12g -jar /usr/picard/picard.jar SortSam \
            I={input} \
            O={output} \
            SORT_ORDER=coordinate \
            TMP_DIR={config[tmp_dir]}
        """

rule mark_kaga_duplicates:
    input:
        os.path.join(config["bam_dir"], "sorted/{kaga_sample}.bam")
    output:
        bam = os.path.join(config["bam_dir"], "marked/{kaga_sample}.bam"),
        metrics = os.path.join(config["bam_dir"], "marked/{kaga_sample}_metrics.txt")
    singularity:
        config["picard"]
    shell:
        """
        java -jar /usr/picard/picard.jar MarkDuplicates \
            I={input} \
            O={output.bam} \
            M={output.metrics} \
            TMP_DIR={config[tmp_dir]}
        """

rule kaga_merge_bams:
    input:
        expand(os.path.join(config["bam_dir"], "marked/{kaga_sample}.bam"),
                kaga_sample = KAGA_SAMPLES)
    output:
        os.path.join(config["bam_dir"], "merged/kaga.bam")
    params:
        files = lambda wildcards, input: " I=".join(input)
    singularity:
        config["picard"]
    shell:
        """
        java -jar /usr/picard/picard.jar MergeSamFiles \
            I={params.files} \
            O={output} \
            TMP_DIR={config[tmp_dir]}
        """

rule kaga_index_bam:
    input:
        os.path.join(config["bam_dir"], "merged/kaga.bam")
    output:
        os.path.join(config["bam_dir"], "merged/kaga.bai")
    singularity:
        config["picard"]
    shell:
        """
        java -jar /usr/picard/picard.jar BuildBamIndex \
            I={input}
        """

rule kaga_call_haplotypes:
    input:
        bam = os.path.join(config["bam_dir"], "merged/kaga.bam"),
        index = os.path.join(config["bam_dir"], "merged/kaga.bai")
    output:
        os.path.join(config["vcf_dir"], "gvcfs/kaga_{chr}.g.vcf.gz")
    params:
        chr = lambda wildcards: wildcards.chr
    singularity:
        config["gatk"]
    shell:
        """
        gatk --java-options "-Xmx12g" HaplotypeCaller  \
            -R {config[hdrr_reference]} \
            -I {input.bam} \
            -O {output} \
            -L {params.chr} \
            -ERC GVCF \
            --tmp-dir {config[tmp_dir]}
        """

rule kaga_genotype:
    input:
        os.path.join(config["vcf_dir"], "gvcfs/kaga_{chr}.g.vcf.gz")
    output:
        os.path.join(config["vcf_dir"], "final/kaga_{chr}.vcf.gz")
    singularity:
        config["gatk"]
    shell:
        """
        gatk --java-options "-Xmx12g" GenotypeGVCFs \
            -R {config[hdrr_reference]} \
            -V {input} \
            -O {output}
        """

rule kaga_index_vcf:
    input:
        os.path.join(config["vcf_dir"], "final/kaga_{chr}.vcf.gz")
    output:
        os.path.join(config["vcf_dir"], "final/kaga_{chr}.vcf.gz.tbi")
    singularity:
        config["bcftools"]
    shell:
        """
        bcftools index --tbi {input}
        """

rule kaga_merge_chr_vcfs:
    input:
        vcfs = expand(os.path.join(config["vcf_dir"], "final/kaga_{chr}.vcf.gz"),
                chr = CHRS_WITH_MT),
        inds = expand(os.path.join(config["vcf_dir"], "final/kaga_{chr}.vcf.gz.tbi"),
                chr = CHRS_WITH_MT)
    output:
        os.path.join(config["vcf_dir"], "merged/kaga.vcf.gz")
    params:
        files = lambda wildcards, input: " I=".join(input.vcfs)
    singularity:
        config["picard"]
    shell:
        """
        java -jar /usr/picard/picard.jar MergeVcfs \
            I={params.files} \
            O={output} \
            TMP_DIR={config[tmp_dir]}
        """