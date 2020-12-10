#TODO: Find a way to extract the sample names line-wise from accession list
SAMPLES = ["SRR10007909", "SRR10007910"]
DIRECTION = ["1", "2"]

configfile: "config.yaml"

rule all:
    input:
        expand("concat/{author}.fasta", author=config["AUTHOR"])

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

rule mergeSeqs:
    input:
        fwd="data/{sample}_1.fastq",
        rev="data/{sample}_2.fastq"
    output:
        "merged/{sample}.fastq"
    log:
        "logs/mergeSeqs/{sample}.log"
    shell:
        "(bbmerge.sh in1={input.fwd} in2={input.rev} out={output}) 2> {log}"

rule q2aReformat:
    input:
        "merged/{sample}.fastq"
    output:
        "fastas/{sample}.fasta"
    log:
        "logs/q2aReformat/{sample}.log"
    shell:
        "(reformat.sh in={input} out={output}) 2> {log}"

rule groupFormatter:
    input:
        "fastas/{sample}.fasta"
    output:
        "mothur_in/{author}.group"
    log:
        "logs/groupFormatter/{author}.log"
    shell:
        "(./scripts/groupFormatter.sh {wildcards.author} {output}) 2> {log}"


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
        "(cutadapt -g {params.fwd} -o {output} {input} --discard-untrimmed) "
        "2> {log}"

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
        "(cutadapt -a {params.rev} -o {output} {input} --discard-untrimmed) "
        "2> {log}"

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
        "(cutadapt -u {params.fwd} -o {output} {input}) 2> {log}"

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
        "(cutadapt -u {params.rev} -o {output} {input}) 2> {log}"

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
        "(bowtie2 -f {input.fasta} -x {input.db}/PhiX_bowtie_db "
        "-S {output.bowtie} "
        "--un {output.fasta} "
        "--al {output.phix}) 2> {log}"

#TODO: Figure out if this works as expected...
rule phixAccnos:
    input:
        "screened/{sample}.merged.bowtie"
    output:
        "PhiX_out/PhiX.accnos"
    log:
        "logs/phixAccnos/Phix.log"
    shell:
        '(grep -Eo "[A-Z]{3,6}[0-9]+\.[0-9]+" {input} '
        "| awk -F: '{ print $2 }' >> {output}) 2> {log}"

rule concat:
    input: expand("screened/{sample}.fasta", sample=SAMPLES)
    output:
        expand("concat/{author}.fasta", author=config["AUTHOR"])
    log: expand("logs/concat/{author}.log", author=config["AUTHOR"])
    shell: "(cat {input} >> {output}) 2> {log}"

#rule mothurScreen:
#    input:
#        fasta="concat/{author}.fasta",
#        groups="mothur_in/{author}.groups"
#    output:
#        "mothur_out/{author}.good.fasta"
#    shell:
#        "./scripts/mothurScreen.sh {input.fasta} {input.groups} {output}"
#
#rule groupSplit:
#    input:
#        "mothur_out/{author}.good.fasta"
#    output:
#        "split/{sample}.fasta"
#    shell:
#        "./scripts/groupSplit.sh"
#
#rule fastaHeaderrelabel:
#    input:
#        "split/{sample}.fasta"
#    output:
#        "barcoded/{author}.barcoded.fasta"
#    shell:
#        "./scripts/fastaHeaderrelabel.sh"
