
configfile: "config.yaml"

rule all:
    input:
        expand("concat/{author}.fasta", author=config["author"])

#rule srrMunch:
#    input:
#        "SRR10007909"
#    output:
#        "fastqs/SRR10007909_1.fastq",
#        "fastqs/SRR10007909_2.fastq",
#        "fastqs/SRR10007909.fastq",
#    log:
#        "logs/srrMunch/SRR10007909.log"
#    shell:
#        "fasterq-dump -o {output} {input}) 2> {log}"

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
        "(./scripts/groupFormatter.sh {author} {output}) 2> {log}"


rule fwdTrim:
    input:
        "fastas/{sample}.fasta"
    output:
        temp("trimmed/{sample}.temp.fasta")
    params:
        forward=config["primers"]["fwd"]
    log:
        "logs/fwdTrim/{sample}.fwd.log"
    shell:
        "(cutadapt -g {params.forward} -o {output} {input} --discard-untrimmed) "
        "2> {log}"

rule revTrim:
    input:
        "trimmed/{sample}.temp.fasta"
    output:
        "trimmed/{sample}.fasta"
    params:
        reverse=config["primers"]["rev"]
    log:
        "logs/revTrim/{sample}.rev.log"
    shell:
        "(cutadapt -a {params.reverse} -o {output} {input} --discard-untrimmed) "
        "2> {log}"

rule fwdOffset:
    input:
        "trimmed/{sample}.fasta"
    output:
        temp("offset/{sample}.temp.fasta")
    params:
        forward=config["offset"]["fwd"]
    log:
        "logs/fwdOffset/{sample}.log"
    shell:
        "(cutadapt -u {params.forward} -o {output} {input}) 2> {log}"

rule revOffset:
    input:
        "offset/{sample}.temp.fasta"
    output:
        "offset/{sample}.fasta"
    params:
        reverse=config["offset"]["rev"]
    log:
        "logs/revOffset/{sample}.log"
    shell:
        "(cutadapt -u {params.reverse} -o {output} {input}) 2> {log}"

rule phixScreen:
    input:
        fasta="offset/{sample}.fasta",
        db=config["phixdb"]
    output:
        fasta="screened/{sample}.fasta",
        phix="screened/{sample}.merged.PhiX",
        bowtie="screened/{sample}.merged.bowtie"
    log:
        "logs/phixScreen/{sample}.log"
    shell:
        "(bowtie2 -f {input.fasta} -x {input.db} "
        "-S {output.bowtie} "
        "--un {output.fasta} "
        "--al {output.phix}) 2> {log}"

#TODO: Figure out if this works as expected...
rule phixAccnos:
    input:
        "screened/{sample}.merged.bowtie"
    output:
        "screened/PhiX.accnos"
    log:
        "logs/phixAccnos/Phix.log"
    shell:
        "(grep -Eo "[A-Z]{3,6}[0-9]+\.[0-9]+" {input} "
        "| awk -F: '{ print $2 }' >> {output}) 2> {log}"

rule concat:
    input: "screened/{sample}.fasta"
    output: "concat/{author}.fasta"
    log: "logs/concat/{author}.log"
    shell: "(cat {input} > {output}) 2> {log}"

rule mothurScreen:
    input:
        "concat/{author}.fasta",
        "mothur_in/{author}.groups"
    output:
        "mothur_out/{author}"
    shell:
        "./scripts/mothurScreen.sh"

rule groupSplit:
    input:
        "mothur_out/good.fasta"
    output:
        "split/{sample}.fasta"
    shell:
        "./scripts/groupSplit.sh"

#rule fastaHeaderrelabel:
#    input:
#        "split/{sample}.fasta"
#    output:
#        "barcoded/{author}.barcoded.fasta"
#    shell:
#        "./scripts/fastaHeaderrelabel.sh"
