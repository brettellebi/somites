# Get genotypes for trio including both SNPs and indels
rule extract_trio_genos:
    input:
        rules.merge_variants_F0_and_F1.output
    output:
        os.path.join(
            config["working_dir"],
            "genos/F0_and_F1/{ref}/final/all.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/extract_trio_genos/{ref}.log"
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

# Extract only biallelic SNPs
rule extract_trio_snps:
    input:
        rules.merge_variants_F0_and_F1.output
    output:
        os.path.join(
            config["working_dir"],
            "genos/F0_and_F1/snps_with_AD/{ref}.txt"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/extract_trio_snps/{ref}.log"
        ),
    resources:
        mem_mb = 2000
    container:
        config["bcftools_1.14"]
    shell:
        """
        bcftools view \
            --min-alleles 2 \
            --max-alleles 2 \
            --types snps \
            --output-type u \
            {input} |\
        bcftools query \
            --print-header \
            --format '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT\\t%AD]\\n' \
            --output {output} \
                2> {log}
        """

# Get counts of all genotypes within bins
rule trio_gt_counts_in_bins:
    input:
        genos = rules.extract_trio_genos.output,
        chrom_lens = rules.get_chrom_lengths.output,
    output:
        os.path.join(
            config["working_dir"],
            "genos/F0_and_F1/{ref}/counts/{sample}/{bin_length}.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/trio_gt_counts_in_bins/{ref}/{sample}/{bin_length}.log"
        ),
    params:
        bin_length = "{bin_length}",
        sample = "{sample}"
    container:
        config["tidyverse_4.1.3"]
    resources:
        mem_mb = 6000
    script:
        "../scripts/trio_gt_counts_in_bins.R"

# Create Circos plots for homozygosity
rule circos_homozygosity:
    input:
        gt_counts = rules.trio_gt_counts_in_bins.output,
        chrom_lens = rules.get_chrom_lengths.output,
    output:
        plot = "book/plots/circos/trio_homo/{ref}/{bin_length}/{sample}.png",
    log:
        os.path.join(
            config["working_dir"],
            "logs/circos_homozygosity/{ref}/{bin_length}/{sample}.log"
        ),
    params:
        ref = "{ref}",
        bin_length = "{bin_length}",
        sample = "{sample}",
        palette = lambda wildcards: config["palette"][wildcards.sample]
    container:
        config["R_4.1.3"]
    resources:
        mem_mb = 5000
    script:
        "../scripts/circos_homozygosity.R"
    
# Extract target SNPs for F2 mapping
## We want SNPs that are homozygous-divergent between Cab and Kaga
## And are heterozygous in F1, with a decent allele depth
rule extract_homo_div_snps:
    input:
        genos = rules.extract_trio_snps.output,
    output:
        full = os.path.join(
            config["working_dir"],
            "genos/F0_and_F1/{ref}/homo_div/snps_all.csv"
        ),
        sites = os.path.join(
            config["working_dir"],
            "sites_files/F0_Cab_Kaga/{ref}/homo_divergent/F1_het_min_DP.txt"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/extract_homo_div_snps/{ref}.log"
        ),
    params:
        # set minimum allele depth
        min_ad = 5
    container:
        config["tidyverse_4.1.3"]
    resources:
        mem_mb = 20000
    script:
        "../scripts/extract_homo_div_snps.R"