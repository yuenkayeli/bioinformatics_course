#!/usr/bin/env bash

#Trimmomatic version 0.39
#must input .gz files ensure to convert

trimmomatic PE \
-threads 4 \
-phred33 \
$1 $2 \
-baseout /home/ubuntu/ngs/dnaseq_pipeline/data/trimmed_fastq/trimmed_data \
ILLUMINACLIP:/home/ubuntu/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
TRAILING:25 MINLEN:50

#FastQC v0.11.9
#on paired trimmed fastq data output to results/ fastqc_trimmed_reads
fastqc -t 4 -o /home/ubuntu/ngs/dnaseq_pipeline/results/fastqc_trimmed_reads \
/home/ubuntu/ngs/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_1P \
/home/ubuntu/ngs/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_2P

#fastqc on fastqc_untrimmed_reads
fastqc -t 4 -o /home/ubuntu/ngs/dnaseq_pipeline/results/fastqc_untrimmed_reads \
/home/ubuntu/ngs/dnaseq_pipeline/data/untrimmed_fastq/NGS0001.R1.fastq.gz \
/home/ubuntu/ngs/dnaseq_pipeline/data/untrimmed_fastq/NGS0001.R2.fastq.gz

#Alignment >>BWA MEM<<

#alignment with bwa trimmed data with read groups
bwa mem \
-t 4 -v 1 \
-R '@RG\tID:11V6WR1.111.D1375ACXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera-NGS0001\tDT:2021-04-06\tPU:11V6WR1' \
-I 250,50 \
/home/ubuntu/ngs/dnaseq_pipeline/data/reference/hg19.fa.gz \
/home/ubuntu/ngs/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_1P \
/home/ubuntu/ngs/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_2P > /home/ubuntu/ngs/dnaseq_pipeline/data/aligned_data/NGS0001.sam

#>>SAMTOOLS<< --version 1.7
# converts sam to bam > sorts bam > index bam
samtools view -h -b /home/ubuntu/ngs/dnaseq_pipeline/data/aligned_data/NGS0001.sam > /home/ubuntu/ngs/dnaseq_pipeline/data/aligned_data/NGS0001.bam # convert .sam to .bam
samtools sort /home/ubuntu/ngs/dnaseq_pipeline/data/aligned_data/NGS0001.bam > /home/ubuntu/ngs/dnaseq_pipeline/data/aligned_data/NGS0001_sorted.bam #sorts bam file in alphanum order
samtools index /home/ubuntu/ngs/dnaseq_pipeline/data/aligned_data/NGS0001_sorted.bam #creates .bai index file

#basic alignment post processing

#>>PICARD<< marks duplicates- version:2.25.1
cd /home/ubuntu/ngs/dnaseq_pipeline/data/aligned_data
picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt
mv *.txt ~/ngs/dnaseq_pipeline/results/  #send to results
samtools index NGS0001_sorted_marked.bam #index sorted_marked.bam

#>>SAMTOOLS filtered BAM on mapping quality and bitwise flags<< - version 1.7
samtools view -F 1796 -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam
samtools index NGS0001_sorted_filtered.bam # index filtered.bam

#variant calling

#>>FREEBAYES<< version-- v0.9.21
# call variants using freebayes

#perform variant calling
freebayes --bam /home/ubuntu/ngs/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered.bam \
--fasta-reference /home/ubuntu/ngs/dnaseq_pipeline/data/reference/hg19.fa \
--vcf /home/ubuntu/ngs/dnaseq_pipeline/results/NGS0001.vcf # perform variant calling and create .vcf
bgzip /home/ubuntu/ngs/dnaseq_pipeline/results/NGS0001.vcf #zips .vcf file
tabix -p vcf /home/ubuntu/ngs/dnaseq_pipeline/results/NGS0001.vcf.gz #index compressed vcf with tabix

#>>VCFFILTER<< 1.0.2 --version
#Hard filter variants
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
/home/ubuntu/ngs/dnaseq_pipeline/results/NGS0001.vcf.gz > /home/ubuntu/ngs/dnaseq_pipeline/results/NGS0001_filtered.vcf

#>>BEDTOOLS<< -- version v2.30.0
#align filtered.vcf to annotation.bed file to hone alignment
bedtools intersect -header -wa -a /home/ubuntu/ngs/dnaseq_pipeline/results/NGS0001_filtered.vcf \
-b /home/ubuntu/ngs/dnaseq_pipeline/data/annotation.bed > /home/ubuntu/ngs/dnaseq_pipeline/results/NGS0001_filtered_bed.vcf
bgzip /home/ubuntu/ngs/dnaseq_pipeline/results/NGS0001_filtered_bed.vcf # compress filtered .vcf file
tabix -p vcf /home/ubuntu/ngs/dnaseq_pipeline/results/NGS0001_filtered_bed.vcf.gz #index filtered zipped .vcf

#Annotation

#>>ANNOVAR<< -- annotate with annovar
#requires downloaded, updated and installed annovar and download databases for annotation
#Convert filtered_bed.vcf.gz to annovar input format- must stay in annovar directory to run

cd /home/ubuntu/annovar
./convert2annovar.pl -format vcf4 ~/ngs/dnaseq_pipeline/results/NGS0001_filtered_bed.vcf.gz > ~/ngs/dnaseq_pipeline/results/NGS0001_filtered_bed.avinput

#create cvs file
./table_annovar.pl ~/ngs/dnaseq_pipeline/results/NGS0001_filtered_bed.avinput humandb/ -buildver hg19 \
-out ~/ngs/dnaseq_pipeline/results/NGS0001_filtered_bed \
-remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro \
-operation g,g,f,f,f -otherinfo -nastring . -csvout

#END
