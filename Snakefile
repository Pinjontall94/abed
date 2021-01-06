# Specify NCBI's forward and reverse fastq file naming scheme for snakemake
DIRECTION = ["1", "2"]

# Generate sample list by reading NCBI accession list file line-by-line
#TODO: Check if this works...
with open("SRR_Acc_List.txt") as f:
    SAMPLES = [line.rstrip("\n") for line in f]


configfile: "config.yaml"

rule all:
    input:
        expand("concat/{author}.fasta", author=config["AUTHOR"])

# Download fastqs from NCBI, reading from SRR_Acc_List.txt
rule srrMunch:
    input: "SRR_Acc_List.txt"
    output:
        expand("data/{sample}_{direction}.fastq", sample=SAMPLES, direction=DIRECTION)
    log: "logs/srrMunch/output.log"
    shell:
        """
        (while read line; do
            fasterq-dump -O data $line
        done < {input}) 2> {log}
        """

# Merge fastqs with bbmerge
rule mergeSeqs:
    input:
        fwd="data/{sample}_1.fastq",
        rev="data/{sample}_2.fastq"
    output:
        "merged/{sample}.fastq"
    log:
        "logs/mergeSeqs/{sample}.log"
    shell:
        """
        (bbmerge.sh in1={input.fwd} in2={input.rev} out={output}) 2> {log}
        """

# Convert fastqs to fastas with JGI reformat.sh
rule q2aReformat:
    input:
        "merged/{sample}.fastq"
    output:
        "fastas/{sample}.fasta"
    log:
        "logs/q2aReformat/{sample}.log"
    shell:
        """
        (reformat.sh in={input} out={output}) 2> {log}
        """

# Create group file for mothur, using the merged fastas
rule groupFormatter:
    input:
        "fastas/{sample}.fasta"
    output:
        "mothur_in/{author}.group"
    log:
        "logs/groupFormatter/{author}.log"
    shell:
        """
        (./scripts/groupFormatter.sh {wildcards.author} {output}) 2> {log}
        """

# Trim 5' primers (specified in config file) with cutadapt
rule fwdTrim:
    input:
        "fastas/{sample}.fasta"
    output:
        "fwd_trimmed/{sample}.fasta"
    params:
        fwd=config["primers"]["FWD"]
    log:
        "logs/fwdTrim/{sample}.fwd.log"
    shell:
        """
        (cutadapt -g {params.fwd} -o {output} {input} --discard-untrimmed)
            2> {log}
        """

# Trim 3' primers (specified in config file) with cutadapt
rule revTrim:
    input:
        "fwd_trimmed/{sample}.fasta"
    output:
        "rev_trimmed/{sample}.fasta"
    params:
        rev=config["primers"]["REV"]
    log:
        "logs/revTrim/{sample}.rev.log"
    shell:
        """
        (cutadapt -a {params.rev} -o {output} {input} --discard-untrimmed)
            2> {log}
        """

# Trim nucleotides off the 5' end with cutadapt
#   (specified by FWD offset in config file)
rule fwdOffset:
    input:
        "rev_trimmed/{sample}.fasta"
    output:
        "fwd_offset/{sample}.fasta"
    params:
        fwd=config["offset"]["FWD"]
    log:
        "logs/fwdOffset/{sample}.log"
    shell:
        """
        (cutadapt -u {params.fwd} -o {output} {input}) 2> {log}
        """

# Trim nucleotides off the 3' end with cutadapt
#   (specified by FWD offset in config file)
rule revOffset:
    input:
        "fwd_offset/{sample}.fasta"
    output:
        "rev_offset/{sample}.fasta"
    params:
        rev=config["offset"]["REV"]
    log:
        "logs/revOffset/{sample}.log"
    shell:
        """
        (cutadapt -u {params.rev} -o {output} {input}) 2> {log}
        """

# Screen fastas for PhiX contamination with bowtie2 and a specified PhiX db file
rule phixScreen:
    input:
        fasta="rev_offset/{sample}.fasta",
        db=expand("{phixdb}", phixdb=config["phixdb"])
    output:
        fasta="screened/{sample}.fasta",
        phix="screened/{sample}.merged.PhiX",
        bowtie="screened/{sample}.merged.bowtie"
    log:
        "logs/phixScreen/{sample}.log"
    shell:
        """
        (bowtie2 -f {input.fasta} -x {input.db}/PhiX_bowtie_db \
        -S {output.bowtie} \
        --un {output.fasta} \
        --al {output.phix}) 2> {log}
        """

# Generate list of accession numbers with PhiX contamination (if present)
#TODO: Figure out if this works as expected...
rule phixAccnos:
    input:
        "screened/{sample}.merged.bowtie"
    output:
        "PhiX_out/PhiX.accnos"
    log:
        "logs/phixAccnos/Phix.log"
    shell:
        """
        (grep -Eo "[A-Z]{3,6}[0-9]+\.[0-9]+" {input} | awk -F: '{ print $2 }' \
            >> {output}) 2> {log}
        """

# Concatenate fastas, label with the study's author name (from config file)
rule concat:
    input: expand("screened/{sample}.fasta", sample=SAMPLES)
    output:
        expand("concat/{author}.fasta", author=config["AUTHOR"])
    log: expand("logs/concat/{author}.log", author=config["AUTHOR"])
    shell:
        """
        (cat {input} >> {output}) 2> {log}
        """

def phix_contam(wildcards):
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with checkpoints.phixAccnos.get(sample=wildcards.sample).output[0].open() as f:
        if f.read().strip() == "a":
            return "PhiX_out/PhiX.accnos"
        else:
            return ""

# Screen fasta file with mothur
#rule mothurScreen:
#    input:
#        fasta="concat/{author}.fasta",
#        groups="mothur_in/{author}.groups"
#        phix_contam
#    output:
#        "mothur_out/{author}.good.fasta"
#    shell:
#        """
#        (./scripts/mothurScreen.sh {input.fasta} {input.groups} {input.phix} {output}) 2> {log}
#        """
#
# Split the screened fasta back into separate files with the group file
#rule groupSplit:
#    input:
#        "mothur_out/{author}.good.fasta"
#    output:
#        "split/{sample}.fasta"
#    shell:
#        "./scripts/groupSplit.sh"
#
# Relabel the fastas for compatibility with vsearch
#rule fastaHeaderrelabel:
#    input:
#        "split/{sample}.fasta"
#    output:
#        "barcoded/{author}.bar.fasta"
#    shell:
#        "./scripts/fastaHeaderrelabel.sh"

# Reconcatenate fastas to final author-labeled fasta file
#rule finalConcat:
#    input: expand("barcoded/{sample}.barcoded.fasta", sample=SAMPLES)
#    output: "final/{author}.fasta"
#    shell:
#        """
#        (cat {input} > {output}) 2> {log}
#        """
