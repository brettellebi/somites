# Get counts of genotypes within bins
rule trio_gt_counts_in_bins:
    input:
        genos = rules.extract_trio_genos.output,
        chrom_lens = rules.get_chrom_lengths.output,
    output:
        os.path.join(
            config["working_dir"],
            "genos/F0_and_F1/counts/{sample}/{bin_length}.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/trio_gt_counts_in_bins/{sample}/{bin_length}.log"
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

rule circos_homozygosity:
    input:
        gt_counts = rules.trio_gt_counts_in_bins.output,
        chrom_lens = rules.get_chrom_lengths.output,
    output:
        plot = "book/plots/circos/trio_homo/{bin_length}/{sample}.png",
    log:
        os.path.join(
            config["working_dir"],
            "logs/circos_homozygosity/{bin_length}/{sample}.log"
        ),
    params:
        bin_length = "{bin_length}",
        sample = "{sample}",
        palette = lambda wildcards: config["palette"][wildcards.sample]
    container:
        config["R_4.1.3"]
    resources:
        mem_mb = 5000
    script:
        "../scripts/circos_homozygosity.R"
    
