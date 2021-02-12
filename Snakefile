# Specify NCBI's forward and reverse fastq file naming scheme for snakemake
DIRECTION = ["1", "2"]

# Generate sample list by reading NCBI accession list file line-by-line
#TODO: Check if this works...
with open("SRR_Acc_List.txt") as f:
    SAMPLES = [line.rstrip("\n") for line in f]


configfile: "config.yaml"

rule all:
    input:
        #expand("barcoded/{author}.fasta", author=config["AUTHOR"])
        "metaan.otus.final.readmap.table"

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


#TODO: Use input functions & config file to clean up the next four rules
#       (as well as letting the offset rules be optional)
#
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

# Create group file for mothur, using the merged fastas
# TODO: Create input function based on presence/absence of config file values
#           (primers, offsets if they exist)
rule groupFormatter:
    input:
        expand("rev_offset/{sample}.fasta", sample=SAMPLES)
    output:
        "mothur_in/{author}.groups"
    log:
        "logs/groupFormatter/{author}.log"
    shell:
        """
        (./scripts/groupFormatter.sh {wildcards.author} {output} {input}) 2> {log}
        """

# Screen fastas for PhiX contamination with bowtie2 and a specified PhiX db file
# TODO: Replace bowtie2 with BBDuk in kmer filtering mode (requires only
#           phix.fasta)
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
# NOTE: Used {sample}.merged.bowtie...shouldn't this be screening
#   {sample}.merged.PhiX files instead?
rule phixAccnos:
    input:
        expand("screened/{sample}.merged.PhiX", sample=SAMPLES)
    output:
        "PhiX_out/PhiX.accnos"
    log:
        "logs/phixAccnos/Phix.log"
    shell:
        """
        (./scripts/phixContam.sh {output} {input}) 2> {log}
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

# Screen fasta file with mothur
rule mothurScreen:
    input:
        fasta="concat/{author}.fasta",
        groups="mothur_in/{author}.groups",
        phix="PhiX_out/PhiX.accnos"
    output:
        "mothur_out/{author}.good.fasta"
    log: "logs/mothurScreen/{author}.log"
    shell:
        """
        (./scripts/mothurScreen.sh {input.fasta} {input.groups} {input.phix} {output}) 2> {log}
        """

# Split the screened fasta back into separate files with the group file
rule groupSplit:
    input:
        fasta=expand("mothur_out/{author}.good.fasta", author=config["AUTHOR"]),
        groups=expand("mothur_in/{author}.groups", author=config["AUTHOR"])
    output:
        "split/{author}.good.{sample}.fasta"
    log: "logs/groupSplit/{author}.good.{sample}.log"
    shell:
        "(./scripts/groupSplit.sh {input.fasta} {input.groups} split) 2> {log}"


# Relabel the fastas for compatibility with vsearch
rule fastaHeaderrelabel:
    input:
        expand("split/{author}.good.{sample}.fasta", sample=SAMPLES, author=config["AUTHOR"])
    output:
        "barcoded/{author}.fasta"
    log: "logs/fastaHeaderrelabel/{author}.log"
    shell:
        "(./scripts/fastaHeaderrelabel.sh {output} {input}) 2> {log}"

#TODO: Add vsearch pipeline from otuGen script

#rule otuGen:
#    input:
#        expand("barcoded/{author}.fasta", author=config["AUTHOR"])
#    output:
#        "metaan.otus.final.readmap.table"
#    log: "logs/otuGen/metaan.log"
#    shell:
#        """
#        (./scripts/otuGen -d silva_v138_mothur -f {input}) 2> {log}
#        """

# Dereplicate the barcoded fasta
#
# Poret-Peterson:
# "Singletons (OTUs represented 1 sequence) are discarded
# via the --minuniquesize command; setting this to 2 would
# generate OTUs represented by 2 or more sequences"
rule derep:
    input:
        expand("barcoded/{author}.fasta", author=config["AUTHOR"])
    output: "derep/{author}.fasta"
    params:
        unique_size=2
    log: "logs/derep/{author}.log"
    shell:
        """
        (vsearch --derep_fulllength {input} \
        --sizein --sizeout --minuniquesize {params.unique_size} \
        --output {output}) 2> {log}
        """


##   # Cluster seqs into OTUs
##   #
##   # Poret-Peterson:
##   #
##   #   "Cluster sequences into OTUs (operational taxonomic units) at a
##   #   97% cutoff (97% or more identity among sequences in an OTU).
##   #   Set threads to 40 in this example. This will need to be changed."
#rule clusterOTUs:
#    input:
#        expand("derep/{author}.fasta", author=config["AUTHOR"])
#    output:
#        clstr="clustered/{author}.otus.fasta",
#        uc="clustered/{author}.otus.uc"
#    params:
#        cutoff=0.97
#    log: "logs/clusterOTUs/{author}.log"
#    shell:
#        """
#        (vsearch --cluster_fast {input} --id 0.97 \
#        --centroids {output.clstr} --uc {output.uc} \
#        --relabel OTU_ --sizein --sizeout) 2> {log}
#        """
#
## Screen for chimeric sequences
#rule chimeraScreen:
#    input:
#        expand("clustered/{author}.otus.fasta", author=config["AUTHOR"])
#    output: "nonchimeras/{author}.fasta"
#    log: "logs/chimeraScreen/{author}.log"
#    shell:
#        """
#        (vsearch --uchime_denovo {input} --abskew 1.5\
#        --nonchimeras {output} --fasta_width 0) 2> {log}
#        """
#
#
#rule taxoClassify:
#    input:
#        nc=expand("nonchimeras/{author}.fasta", author=config["AUTHOR"]),
#        lambda wildcards: config["taxonomy"]["template"][wildcards.template],
#        lambda wildcards: config["taxonomy"]["tax"][wildcards.tax]
#    output:
#        pick="taxo/{author}.pick.fasta",
#        readmap="taxo/{author}.readmap.uc"
#    log: "logs/taxoClassify/{author}.log"
#    shell:
#        """
#        (./scripts/taxoClassify.sh {input.nc} {wildcards.template} {wildcards.tax} \
#        {output.pick} {output.readmap}) 2> {log}
#        """
##   # Uses mothur to taxonomically classify OTUs
##   #   [Note: May need to update silva version prior to running]
##   TAXO_TEMPLATE=$(find $TAXO_DB -maxdepth 1 -name "*.ng.fasta")
##   TAXO_NOMY=$(find $TAXO_DB -maxdepth 1 -name "*.tax")
##
##   pathTester mothur $MOTHUR_MODULE $MOTHUR_MOD_VER
##   TAXO_PARAMS+=("fasta=metaan.otus.nc.fasta,"\
##       "template=$TAXO_TEMPLATE,"\
##       "taxonomy=$TAXO_NOMY,"\
##       "cutoff=80", "output=metaan.otus.nc.good.wang.taxonomy")
##   LINEAGE_PARAMS+=("fasta=metaan.otus.nc.fasta,"\
##       "taxonomy=metaan.otus.nc.nr_v138.wang.taxonomy, taxon=unknown")
##   mothurBatch classify seqs "${TAXO_PARAMS[*]}" > taxo_batch.txt
##   mothurBatch remove lineage "${LINEAGE_PARAMS[*]}" >> taxo_batch.txt
##   mothur taxo_batch.txt
##   rm taxo_batch.txt
##}
##
##
##mapOTUs(){
##   # Maps OTUs against reads and converts them to an OTU table
##   vsearch --usearch_global metaan.barcoded.fasta\
##       -db metaan.otus.nc.pick.fasta\
##       --uc metaan.otus.nc.readmap.uc --id 0.97\
##       --strand plus --threads $THREADS
##}
##
##
##readmap2file(){
##   # Convert the readmap to an OTU table
##   # TODO: Reimplement this script so this can all be one file
##   create_otu_table_from_uc_file.py\
##       -i metaan.otus.nc.readmap.uc\
##       -o metaan.otus.final.readmap.table
