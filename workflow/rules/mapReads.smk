SAMPLE=config["Sample name"]
ID=config["ID"]
fasta=config["fasta file"]
# identifying, orienting and removing adapter sequences and poly-A tails using Pychopper and cutadapt.
rule pychopper:
    output: 
        expand("results/processed_fastq/pychopper/pychopper_{sample}_{ident}/{sample}_enriched_{ident}_full_length_output.fq", sample=SAMPLE, ident=ID),
        expand("results/processed_fastq/pychopper/pychopper_{sample}_{ident}/{sample}_control_{ident}_full_length_output.fq", sample=SAMPLE, ident=ID)
    params:
        enriched=config["enriched fastq"],
        control=config["control fastq"],
        sample=SAMPLE,
        ident=ID
    conda:
        "../envs/env_read_mapping.yaml"
    shell:
        """
        	pychopper -r results/processed_fastq/pychopper/pychopper_{params.sample}_{params.ident}/report_enriched.pdf -S results/processed_fastq/pychopper/pychopper_{params.sample}_{params.ident}/statistics_enriched.tsv {params.enriched} {output[0]}

        	pychopper -r results/processed_fastq/pychopper/pychopper_{params.sample}_{params.ident}/report_control.pdf -S results/processed_fastq/pychopper/pychopper_{params.sample}_{params.ident}/statistics_control.tsv {params.control} {output[1]}
        """

rule cutadapt:
    input: "results/processed_fastq/pychopper/pychopper_{sample}_{ident}/{sample}_{cond}_{ident}_full_length_output.fq"
    output: temp("results/processed_fastq/cutadapt/{sample}_{cond}_{ident}_cutadapt_temp1.fq"),
            temp("results/processed_fastq/cutadapt/{sample}_{cond}_{ident}_cutadapt_temp2.fq"),
            temp("results/processed_fastq/cutadapt/{sample}_{cond}_{ident}_cutadapt_2.fq")
    conda:
        "../envs/env_read_mapping.yaml"
    shell:
        """
            cutadapt --poly-a -j 0 -o {output[0]} {input}
            ## reads that are completely trimmed (i.e., remaining sequence is empty) are not removed by default. If not removed, result in ZeroDivisionError for --max-aer argument
            chmod +x ./remove_empty_seqs.py
            ./remove_empty_seqs.py {output[0]} {output[1]}
            ## --max-aer recommended for long-read sequencing
            cutadapt -j 0 -o {output[2]} {output[1]} --max-aer 0.1
        """
rule cutadapt2:
    input: "results/processed_fastq/cutadapt/{sample}_{cond}_{ident}_cutadapt_2.fq"
    output: "results/processed_fastq/cutadapt/{sample}_{cond}_{ident}_cutadapt.fq"
    conda:
        "../envs/env_read_mapping.yaml"
    shell:
        """
            cutadapt -g TTTCTGTTGGTGCTGATATTGCTGGG -e 0.1 -j 0 -o {output} {input}
        """
# Mapping reads onto genome
rule minimap2:
    input: "results/processed_fastq/cutadapt/{sample}_{cond}_{ident}_cutadapt.fq"
    output: temp("results/alignments/{sample}_{cond}_{ident}_minimap.sam")
    conda:
        "../envs/env_read_mapping.yaml"
    params:
        fa=fasta
    shell:
        """
        minimap2 -ax map-ont -k14 -t 8 -p 0.99 {params.fa} {input} > {output}
        """
# creates index file for fasta file when needed
rule fastaIndex:
    output: expand("{fa}.fai",fa=fasta)
    conda:
        "../envs/env_read_mapping.yaml"
    params:
        fa=fasta
    shell:
        """
        export LC_ALL=C
        samtools faidx {params.fa}
        """
# Soft clipping and converting SAM file to sorted BAM file
rule clipping:
    input: "results/alignments/{sample}_{cond}_{ident}_minimap.sam", expand("{fa}.fai",fa=fasta)
    output: temp("results/alignments/{sample}_{cond}_{ident}_clipped.sam")
    conda:
        "../envs/env_read_mapping.yaml"
    params:
        fa=fasta
    shell:
        """
        export LC_ALL=C
        samclip --max 10 --ref  {params.fa} < {input[0]} > {output}
        """
# convert SAM files to BAM files
rule samToBam:
    input: "results/alignments/{sample}_{cond}_{ident}_clipped.sam"
    output: temp("results/alignments/{sample}_{cond}_{ident}.bam")
    conda:
        "../envs/env_read_mapping.yaml"
    shell:
        """
        export LC_ALL=C
        samtools view -S -b {input} > {output}
        """
# sort BAM files
rule sortBam:
    input: "results/alignments/{sample}_{cond}_{ident}.bam"
    output: "results/alignments/BAM_files_{sample}/{sample}_{cond}_{ident}.sorted.bam"
    conda:
        "../envs/env_read_mapping.yaml"
    shell:
        """
        export LC_ALL=C
        samtools sort {input} -o {output}
        """
#index BAM files
rule indexSortedBAM:
    input: "results/alignments/BAM_files_{sample}/{sample}_{cond}_{ident}.sorted.bam"
    output: "results/alignments/BAM_files_{sample}/{sample}_{cond}_{ident}.sorted.bam.bai"
    conda:
        "../envs/env_read_mapping.yaml"
    shell:
        """
        export LC_ALL=C
        samtools index {input}
        """
