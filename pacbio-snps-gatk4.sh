#!/bin/bash

# This pipeline call snps from pacbio reads through gatk4.
# Follow these steps:
# minimap2 -> samtools sort -> MarkDuplicates -> gatk HaplotypeCaller


ref=$1
fastq=$2
sample=$3
thread=$4


function help {
        echo
        echo -e "The pipeline of pacbio reads call snps by gatk4."
        echo -e "This will execute follow these steps:"
        echo -e "   minimap2 -> samtools sort -> MarkDuplicates -> gatk HaplotypeCaller"
        echo
}

function usage {
        echo -e "Usage : $(basename $0) reference pacbio_fastq out_prefix [thread=1]"
        echo -e "Usage : $(basename $0) ref.fasta pacbio.reads.fastq pacbio 4"
}


if [[ -z $ref || -z $fastq || -z $sample ]];then
        usage;
        help;
        exit
fi

logfile=run.log
exec 2>$logfile

if [ ! -e $ref ];then
    echo -e "[ERROR] No such file of $ref";
    exit
fi
if [ ! -e $fastq ];then
    echo -e "[ERROR] No such file of $fastq";
    exit
fi

if [ -z $thread ];then
    thread=1
fi


create_fa_index.sh -f $ref -t $thread

# mapping
echo -e "[INFO] mapping..."
minimap2 --secondary=no -R "@RG\tID:${sample}\tSM:${sample}" -ax asm20 ${ref} ${fastq} | samtools view -bS > ${sample}.aln.bam

# sort bam
echo -e "[INFO] sorting bam..."

#java -jar $picard SortSam \
#    INPUT=${sample}.aln.bam \
#    OUTPUT=${sample}.sorted.bam \
#    SORT_ORDER='coordinate' \
#    VALIDATION_STRINGENCY=LENIENT \
#    MAX_RECORDS_IN_RAM=300000

samtools sort -@ ${thread} -o ${sample}.sorted.bam ${samtools}.aln.bam
if [ ! -s ${sample}.sorted.bam ];then
    echo -e "[ERROR] Failed to execute SortSam, please check the detail in $logfile";
    exit
fi

# Markduplicates, remove duplicates and create bam index
echo -e "[INFO] removing duplicates..."
java -jar $picard MarkDuplicates \
    INPUT=${sample}.sorted.bam \
    OUTPUT=${sample}.rmd.bam \
    METRICS_FILE=${sample}.rmd.metrics \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true #create bam index
if [ ! -s ${sample}.rmd.bam ];then
    echo -e "[ERROE] Failed to execute MarkDuplicates, please check the detail in $logfile";
    exit
fi


# call snps by gatk4 HaplotypeCaller
echo -e "[INFO] caculating snps by gatk4 HaplotypeCaller..."
gatk \
    HaplotypeCaller \
    --java-options "-Xms5g -Xmx20g -XX:ParallelGCThreads=$thread" \
    -R ${ref} \
    -I ${sample}.rmd.bam \
    -O ${sample}.vcf \
    --native-pair-hmm-threads ${thread} \
    -stand-call-conf 30
if [ -s ${sample}.vcf ];then
    echo -e "[INFO] Done, result is ${sample}.vcf"
else
    echo -e "[ERROE] Failed to execute gatk, please check the detail in $logfile"
fi
