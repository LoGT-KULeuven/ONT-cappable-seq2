# create separate genomefiles for bedtools from each nucleotide sequence in the fasta file.
# Checkpoint is the alternative to dynamic() in previous snakemake versions
# TODO: find cleaner way to do this
checkpoint getGenomeFiles:
    input: expand("{fa}.fai",fa=config["fasta file"])
    output: directory("results/transcript_boundaries/gfiles/")
    shell:
        """
        mkdir {output};
        while IFS='\t', read first second; do
          hash=$(echo -n $first | tr -d '"' | md5sum | awk NF=1)
          echo -e "$first\t$second" > {output}/"$hash"_genome.txt
        done < {input}
        """

# get the names of the genomefiles
# input function for rule aggregate, return paths to all files produced by the checkpoint 'somestep'
def aggregate_input(wildcards):
    checkpoint_output = checkpoints.getGenomeFiles.get(**wildcards).output[0]
    return expand("results/transcript_boundaries/gfiles/{i}_genome.txt", i=glob_wildcards(os.path.join(checkpoint_output, "{i}_genome.txt")).i)

rule PosEffRatios:
    input: 
        '{sample}_peak_calling/{sample}_enriched_{ident}.3end.plus.peaks.oracle.narrowPeak.counts.clustered.csv',
        'results/alignments/BAM_files_{sample}/{sample}_enriched_{ident}.sorted.bam',
        aggregate_input
    output: 
        'results/transcript_boundaries/TTS_{sample}/TTS_{sample}_{ident}/eff_ratios_{sample}_{ident}.plus.drop.coverage'
    params:
        t=config["TTS threshold"]
    conda:
        "../envs/env_annotation.yaml"
    shell:
        """
        export LC_ALL=C

        bedFile=$(bedtools bamtobed -i {input[1]} | sort -k1,1 -k2,2n)
        posBedFile=$(grep -w "+" <(echo "$bedFile"))
        [ ! -e {output} ] || rm {output}
        > {output}

        sort -t ',' -k15 -n -u {input[0]} | awk -F ',' 'NR>1 {{print $2, $15}}' | while read -r chr peak;
        do
            chrtr=$(echo $chr | tr -d '"')
            hash_file="results/transcript_boundaries/gfiles/"$(echo -n "$chrtr" | md5sum | awk NF=1)"_genome.txt"

            posBedFileTTS=$(awk -v TTS="$peak" -v chr="$chrtr" '$1 == chr && $2 < TTS-10' <(echo "$posBedFile"))
            genomeCovFile=$(bedtools genomecov -g $hash_file -i <(echo "$posBedFileTTS") -d)

            upstreamTTSaverageCoverage=$(awk -v TTS="$peak" '$2 < TTS && $2 >= TTS-20' <(echo "$genomeCovFile") | \
            awk -v TTS="$peak" '{{total += $3}} END {{print TTS, TTS-20, total/NR}}')
            downstreamTTSaverageCoverage=$(awk -v TTS="$peak" '$2 > TTS && $2 <= TTS+20' <(echo "$genomeCovFile") |\
            awk -v TTS="$peak" '{{total += $3}} END {{print TTS, TTS+20, total/NR}}')
            upstreamValue=$(awk '{{ print $3 }}' <(echo "$upstreamTTSaverageCoverage"))

            if (( $(echo $upstreamValue'>'0 | bc -l) ));
            then
                paste -d " " <(echo "$upstreamTTSaverageCoverage") <(echo "$downstreamTTSaverageCoverage") |\
                awk -v chr=$chr '{{print chr, $1, $2, $3, $4, $5, $6, $6/$3, 1-$6/$3}}' |\
                awk '{{ if($9 >= {params.t}) {{ print }} }}' >> {output}
            fi
        done
        """
rule NegEffRatios:
    input: 
        '{sample}_peak_calling/{sample}_enriched_{ident}.3end.minus.peaks.oracle.narrowPeak.counts.clustered.csv',
        'results/alignments/BAM_files_{sample}/{sample}_enriched_{ident}.sorted.bam',
        aggregate_input
    output: 
        'results/transcript_boundaries/TTS_{sample}/TTS_{sample}_{ident}/eff_ratios_{sample}_{ident}.minus.drop.coverage'
    params:
        t=config["TTS threshold"]
    conda:
        "../envs/env_annotation.yaml"
    shell:
        """
        export LC_ALL=C

        bedFile=$(bedtools bamtobed -i {input[1]} | sort -k1,1 -k2,2n)
        negBedFile=$(grep -w "-" <(echo "$bedFile"))
        [ ! -e {output} ] || rm {output}
        > {output}

        sort -t ',' -k15 -n -u {input[0]} | awk -F ',' 'NR>1 {{print $2, $15}}' | while read -r chr peak;
        do
            chrtr=$(echo $chr | tr -d '"')
            hash_file="results/transcript_boundaries/gfiles/"$(echo -n "$chrtr" | md5sum | awk NF=1)"_genome.txt"

            negBedFileTTS=$(awk -v TTS="$peak" -v chr="$chrtr" '$1 == chr && $2 < TTS+10' <(echo "$negBedFile"))
            genomeCovFile=$(bedtools genomecov -g "$hash_file" -i <(echo "$negBedFileTTS") -d)
            downstreamTTSaverageCoverage=$(awk -v TTS="$peak" '$2 < TTS && $2 >= TTS-20' <(echo "$genomeCovFile") | \
            awk -v TTS="$peak" '{{total += $3}} END {{print TTS, TTS-20, total/NR}}')
            upstreamTTSaverageCoverage=$(awk -v TTS="$peak" '$2 > TTS && $2 <= TTS+20' <(echo "$genomeCovFile") |\
            awk -v TTS="$peak" '{{total += $3}} END {{print TTS, TTS+20, total/NR}}')
            upstreamValue=$(awk '{{ print $3 }}' <(echo "$upstreamTTSaverageCoverage"))

            if (( $(echo $upstreamValue'>'0 | bc -l) ));
            then
                paste -d " " <(echo "$upstreamTTSaverageCoverage") <(echo "$downstreamTTSaverageCoverage") |\
                awk -v chr=$chr '{{print chr, $1, $2, $3, $4, $5, $6, $6/$3, 1-$6/$3}}' |\
                awk '{{ if($9 >= {params.t}) {{ print }} }}' >> {output}
            fi
        done
        """
rule PosCoverageDropToBed:
    input: 'results/transcript_boundaries/TTS_{sample}/TTS_{sample}_{ident}/eff_ratios_{sample}_{ident}.plus.drop.coverage'
    output: temp('results/transcript_boundaries/TTS_{sample}/TTS_{sample}_{ident}/TTS_{sample}_{ident}.plus.bed')
    params:
        up=config["TTS sequence extraction"]["upstream"],
        down=config["TTS sequence extraction"]["downstream"]
    shell:
        """
        awk -v up={params.up} -v down={params.down} -v FS='\t' -v OFS='\t' -F ' ' '{{print $1, $2 - up, $2 + down, "TTS_POS_" NR,0 , "+"}}' {input} | uniq > {output}
        sed -i 's/\"//g' {output}
        """
rule NegCoverageDropToBed:
    input: 'results/transcript_boundaries/TTS_{sample}/TTS_{sample}_{ident}/eff_ratios_{sample}_{ident}.minus.drop.coverage'
    output: temp('results/transcript_boundaries/TTS_{sample}/TTS_{sample}_{ident}/TTS_{sample}_{ident}.minus.bed')
    params:
        up=config["TTS sequence extraction"]["upstream"],
        down=config["TTS sequence extraction"]["downstream"]
    shell:
        """
        awk -v up={params.up} -v down={params.down} -v FS='\t' -v OFS='\t' -F ' ' '{{print $1, $2 - down, $2 + up, "TTS_NEG_" NR,0 , "-"}}' {input} | uniq > {output}
        sed -i 's/\"//g' {output}
        """
rule combinePosAndNegBedFilesTTS:
    input: 'results/transcript_boundaries/TTS_{sample}/TTS_{sample}_{ident}/TTS_{sample}_{ident}.plus.bed', 'results/transcript_boundaries/TTS_{sample}/TTS_{sample}_{ident}/TTS_{sample}_{ident}.minus.bed'
    output: 'results/transcript_boundaries/TTS_{sample}/TTS_{sample}_{ident}/TTS_{sample}_{ident}.bed'
    shell:
        """
        cat {input} > {output}
        """

rule correctNegativeValuesTTS:
    input:
        'results/transcript_boundaries/TTS_{sample}/TTS_{sample}_{ident}/TTS_{sample}_{ident}.bed'
    output:
        'results/transcript_boundaries/TTS_{sample}/TTS_{sample}_{ident}/TTS_{sample}_{ident}_corrected.bed'
    run:
        with open(input[0], 'r') as infile, open(output[0], 'w') as outfile:
            for line in infile:
                parts = line.split()
                if int(parts[1]) < 0:
                    parts[1] = '0'
                outfile.write('\t'.join(parts) + '\n')

rule extractSequencesTTS:
    input: 
        'results/transcript_boundaries/TTS_{sample}/TTS_{sample}_{ident}/TTS_{sample}_{ident}_corrected.bed'
    output:
        'results/transcript_boundaries/TTS_{sample}/TTS_{sample}_{ident}/TTS_seq_{sample}_{ident}.fa.out'
    params:
        fasta=config["fasta file"]
    conda:
        "../envs/env_annotation.yaml"
    shell:
        """
        bedtools getfasta -fi {params.fasta} -bed {input} -fo {output} -s -name 
        """
