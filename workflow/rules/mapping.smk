#rule trim_reads_pe:
#    input:
#        unpack(get_fastq),
#    output:
#        r1=temp(os.path.join(config["working_dir"], "fastqs/trimmed/{sample}-{unit}.1.fastq.gz")),
#        r2=temp(os.path.join(config["working_dir"], "fastqs/trimmed/{sample}-{unit}.2.fastq.gz")),
#        r1_unpaired=temp(os.path.join(config["working_dir"], "fastqs/trimmed/{sample}-{unit}.1.unpaired.fastq.gz")),
#        r2_unpaired=temp(os.path.join(config["working_dir"], "fastqs/trimmed/{sample}-{unit}.2.unpaired.fastq.gz")),
#        trimlog=os.path.join(config["working_dir"], "fastqs/trimmed/{sample}-{unit}.trimlog.txt"),
#    params:
#        **config["params"]["trimmomatic"]["pe"],
#        extra=lambda w, output: "-trimlog {}".format(output.trimlog),
#    wrapper:
#        "0.74.0/bio/trimmomatic/pe"

rule map_reads:
    input:
        target = config["ref_path"],
        query = get_fastq,
    output:
        os.path.join(config["working_dir"], "sams/mapped/{sample}-{unit}.sam")
#    params:
#        extra="-ax sr"
    threads: config["params"]["minimap2"]["threads"]
#    wrapper:
#        "0.74.0/bio/minimap2/aligner"
    conda:
        "../envs/minimap2_2.19.yaml"
    shell:
        """
        minimap2 -ax sr {input.target} {input.query} > {output}
        """