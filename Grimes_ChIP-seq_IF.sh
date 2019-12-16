 #!/bin/bash

printf "====================\n  Grimes: ChIP-seq Pipeline\n====================\n"

if [[ $# -eq 0 ]]; then
	echo " grimes [prefix] [FASTQ file 1] [FASTQ file 2] [Adapter sequece]"
	echo 
	exit 0
fi
prefix=$1
fastq_1=$2
fastq_2=$3


######################
## Adaptor trimming ##
######################


trim_adaptors(){

	Sample=$1
	Adaptor_file=$2

	file1=/mnt/san2/lab5/ewilkie/Heidi/Run3/50bp_reads/Trimmed_reads/"$Sample"

	cutadapt=$(cutadapt -b "$Adaptor_file" -o "$Sample"_output "$file1")
	#echo $fastqc
	echo "$Sample"
	echo "$Adaptor_file"
	

}

trim_adaptors Input_50bp.fastq /mnt/san2/lab5/ewilkie/Heidi/Run3/50bp_reads/Adaptors/adaptors_input.fasta
trim_adaptors Twist_50bp.fastq /mnt/san2/lab5/ewilkie/Heidi/Run3/50bp_reads/Adaptors/adaptors_twist.fasta


###########QC1##################

fastqc(){

	Sample=$1

	file1=/mnt/san2/lab5/ewilkie/Heidi/Run3/50bp_reads/Cutadaptreads/"$Sample"

	fastqc=$(/mnt/san2/tools/FastQC/0.11.2/fastqc -o /mnt/san2/lab5/ewilkie/Heidi/Run3/50bp_reads/Fastqc_trimmed "$file1")
	echo $fastqc

	printf " fastqc for %s has been run"

}

fastqc Input_50bp_cutadapt_trimmed.fastq
fastqc Twist_50bp_cutadapt_trimmed.fastq

genome_dir="/home/nsalehin/processing/gencode/GRCm38"# chekc
Q="tools/qchipseq/2017-03-08/bin/Q"
bwa_index=""#check

echo "Aligning reads"
bwa mem -t ${THREAD_NO} -M ${bwa_index} ${fastq_1} ${fastq_2} > chip.sam #check mm10

echo "Prefiltering alignments"
SAMTOOLS="/tools/samtools/1.3/samtools"
 prefiltered="${prefix}.prededup.bam"
${SAMTOOLS} sort --threads 4 Aligned.out.sam | ${SAMTOOLS} view --threads 4 -F 1804 -f 2 -b -o ${prefiltered}

echo "Removing duplicates"
MARKDUP="/tools/picard-tools/1.119/MarkDuplicates.jar"
dedup="${prefix}.dedup.bam"
java -Xmx4G -jar ${MARKDUP} INPUT=${prefiltered} OUTPUT=${dedup} METRICS_FILE=REPORT.markdup.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true


echo "Removing orphaned reads"
cleaned="${prefix}.cleaned.bam"
${SAMTOOLS} view -@ 4 -F 1804 -f 2 -b ${dedup} | ${SAMTOOLS} sort -@ 4 -n -T preshift_sort -o ${cleaned}

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 --extsize=250 --nomodel -g mm -t bam -n test

echo "Thank your for flying with us. Have a nice stay!"

