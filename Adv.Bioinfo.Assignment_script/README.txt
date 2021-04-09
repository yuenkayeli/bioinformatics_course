SETUP BEFORE RUNNING BASH PIPELINE.SH

1.##install programmes needed to run pipeline:
 install conda -- 4.10.0

$ wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-x86_64.sh
$ chmod +x ./Miniconda3-py39_4.9.2-Linux-x86_64.sh
$ bash ./Miniconda3-py39_4.9.2-Linux-x86_64.sh
$ source ~/.bashrc

install tools:

$ conda install samtools
$ conda install bwa
$ conda install freebayes
$ conda install picard
$ conda install bedtools
$ conda install trimmomatic
$ conda install fastqc
$ conda install vcflib
'''

2.##install ANNOVAR for annotation
 # download, unzip and install into /home/ubuntu
$ tar -zxvf annovar.latest.tar.gz 
# to use ANNOVAR command must be run while in directory
$ cd annovar

#download databases required
$ ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
$ ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
$ ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
$ ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
$ ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
$ ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/

#zip file can be removed
$ rm annovar.latest.tar.gz

'''
3.#Create directory tree before running bash script
 
>>directiory tree<<

|annovar
|miniconda3
  |-ngs
    |-README.txt
    |-dnaseq_pipeline
       |- data
       |   |- reference
       |   |    |---hg19.fa.gz
       |   |
       |   |- annotation.bed
       |   |
       |   |- untrimmed_fastq
       |   |    |---NGS0001.R1.fastq.qz (change to .gz)
       |   |    |---NGS0001.R2.fastq.qz (change to .gz)
       |   |
       |   |- trimmed_fastq
       |   |- aligned_data    
       |
       |- results
       |   |- fastqc_untrimmed_reads
       |   |- fastqc_trimmed_reads
       |              
       |- scripts  

4.# download files in correct folders following tree

$ cd ~/ngs/dnaseq_pipeline/data/untrimmed_fastq
$ wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
$ wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz

#convert fastq.qz to fastq.gz
$ mv NGS0001.R1.fastq.qz NGS0001.R1.fastq.gz
$ mv NGS0001.R2.fastq.qz NGS0001.R2.fastq.gz

total 1.5G
-rw-rw-r-- 1 ubuntu ubuntu 743M Feb 28  2017 NGS0001.R1.fastq.gz
-rw-rw-r-- 1 ubuntu ubuntu 759M Feb 28  2017 NGS0001.R2.fastq.gz

'''
#download annotation.bed file
$ cd ~/ngs/dnaseq_pipeline/data/aligned_data
$ wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed
-rw-rw-r-- 1 ubuntu ubuntu  11M Feb 27  2017 annotation.bed

'''
#download human genome reference
$ cd ~/ngs/dnaseq_pipeline/data/reference
$ wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
-rw-rw-r-- 1 ubuntu ubuntu 905M Aug 21  2018 hg19.fa.gz

'''
5.#creates index for reference genome hg19.fa.gz
bwa index ~/ngs/dnaseq_pipeline/data/reference/hg.19.fa.gz
#when completed there should be
#5 files (.amb .ann .bwt .pac .sa)
-rw-rw-r-- 1 ubuntu ubuntu 8.4K Apr  6 15:44 hg19.fa.gz.amb
-rw-rw-r-- 1 ubuntu ubuntu 4.0K Apr  6 15:44 hg19.fa.gz.ann
-rw-rw-r-- 1 ubuntu ubuntu 3.0G Apr  6 15:43 hg19.fa.gz.bwt
-rw-rw-r-- 1 ubuntu ubuntu 748M Apr  6 15:43 hg19.fa.gz.pac
-rw-rw-r-- 1 ubuntu ubuntu 1.5G Apr  6 15:58 hg19.fa.gz.sa

'''
6.#uncompress hg19.fa.gz file to make fasta file, then create index of .fa file
zcat /home/ubuntu/ngs/dnaseq_pipeline/data/reference/hg19.fa.gz > /home/ubuntu/ngs/dnaseq_pipeline/data/reference/hg19.fa
samtools faidx /home/ubuntu/ngs/dnaseq_pipeline/data/reference/hg19.fa
-rw-rw-r-- 1 ubuntu ubuntu 3.0G Apr  9 12:03 hg19.fa
-rw-rw-r-- 1 ubuntu ubuntu 3.5K Apr  7 12:59 hg19.fa.fai

''' 

7.#go to scripts and run bash pipeline.sh
 
$cd ~/ngs 
$ bash pipeline.sh ~/ngs/dnaseq_pipeline/data/untrimmed_fastq/NGS0001.R1.fastq.gz \
~/ngs/dnaseq_pipeline/data/untrimmed_fastq/NGS0001.R2.fastq.gz 

'''
CONTENTS
runs the following in order

1. Trimmomatic- trims low quality 3'reads and adapter sequences if present
	  OUTPUT: 4 files
                         (trimmed_data_1P & trimmed_data_2P)
                         (timmed_data_1U & trimmed_data_2P)
total 4.2G
-rw-rw-r-- 1 ubuntu ubuntu 2.1G Apr  9 15:37 trimmed_data_1P
-rw-rw-r-- 1 ubuntu ubuntu  58M Apr  9 15:37 trimmed_data_1U
-rw-rw-r-- 1 ubuntu ubuntu 2.1G Apr  9 15:37 trimmed_data_2P
-rw-rw-r-- 1 ubuntu ubuntu  12M Apr  9 15:37 trimmed_data_2U

2. FASTQC- performs fastqc on trimmed sequence and untrimmed sequences
        OUTPUT: 4 file (2x html and 2x.fastq.gz)
trimmed:
-rw-rw-r-- 1 ubuntu ubuntu 617K Apr  9 15:23 trimmed_data_1P_fastqc.html
-rw-rw-r-- 1 ubuntu ubuntu 345K Apr  9 15:23 trimmed_data_1P_fastqc.zip
-rw-rw-r-- 1 ubuntu ubuntu 619K Apr  9 15:23 trimmed_data_2P_fastqc.html
-rw-rw-r-- 1 ubuntu ubuntu 348K Apr  9 15:23 trimmed_data_2P_fastqc.zip

untrimmed:
total 2.0M
-rw-rw-r-- 1 ubuntu ubuntu 612K Apr  9 15:39 NGS0001.R1_fastqc.html
-rw-rw-r-- 1 ubuntu ubuntu 354K Apr  9 15:39 NGS0001.R1_fastqc.zip
-rw-rw-r-- 1 ubuntu ubuntu 615K Apr  9 15:39 NGS0001.R2_fastqc.html
-rw-rw-r-- 1 ubuntu ubuntu 367K Apr  9 15:39 NGS0001.R2_fastqc.zip

3. BWA MEM- aligns trimmed data with reference genome and adds read groups
	OUTPUT: NGS0001.sam
-rw-rw-r-- 1 ubuntu ubuntu 6.0G Apr  9 16:03 NGS0001.sam

   SAMTOOLS- converts .SAM to .BAM, then sorts .BAM, then index .BAM 
	NGS0001.sam > NGS0001.bam > NGS0001_sorted.bam > NGS0001_sorted.bam.bai
-rw-rw-r-- 1 ubuntu ubuntu 1.9G Apr  9 12:44 NGS0001.bam
-rw-rw-r-- 1 ubuntu ubuntu 1.5G Apr  9 13:45 NGS0001_sorted.bam
-rw-rw-r-- 1 ubuntu ubuntu 3.7M Apr  9 13:46 NGS0001_sorted.bam.bai

4. PICARD- marks duplicate sequences
	OUTPUT:2 files NGS0001_sorted_marked.bam and marked_dup_metrics.txt
-rw-rw-r-- 1 ubuntu ubuntu 1.5G Apr  9 13:52 NGS0001_sorted_marked.bam
-rw-rw-r-- 1 ubuntu ubuntu 3.7M Apr  9 13:53 NGS0001_sorted_marked.bam.bai
move NGS0001_sorted_marked.bam to /results and index:NGS0001_sorted_marked.bam.bai

5. SAMTOOLS-post-alignment filtering, filter sorted-marked.bam on mapping quality and bitwise flags
	OUTPUT: NGS0001_sorted_filtered.bam 
-rw-rw-r-- 1 ubuntu ubuntu 1.3G Apr  9 13:59 NGS0001_sorted_filtered.bam
-rw-rw-r-- 1 ubuntu ubuntu 3.7M Apr  9 13:59 NGS0001_sorted_filtered.bam.bai

6. FREEBAYES-variant calling
	OUTPUT:NGS0001.vcf file
 index compressed .vcf and compressed
-rw-rw-r-- 1 ubuntu ubuntu  11M Apr  9 14:32 NGS0001.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 532K Apr  9 14:32 NGS0001.vcf.gz.tbi

7. VCFfilter-Variant filtering for hard variants
	OUTPUT: NGS0001_filtered.vcf
-rw-rw-r-- 1 ubuntu ubuntu  19M Apr  9 14:37 NGS0001_filtered.vcf

  BEDTOOLS- align filted.vcf to .bed file to focus on desired genome regions
	OUTPUT: NGS0001_filtered_bed.vcf
-rw-rw-r-- 1 ubuntu ubuntu 7.1M Apr  9 14:49 NGS0001_filtered_bed.vcf
	compress:NGS0001_filtered_bed.vcf.gz
-rw-rw-r-- 1 ubuntu ubuntu 1.3M Apr  9 14:52 NGS0001_filtered_bed.vcf.gz
	index:NGS0001_filtered_bed.vcf.gz.tbi
-rw-rw-r-- 1 ubuntu ubuntu 122K Apr  9 14:52 NGS0001_filtered_bed.vcf.gz.tbi

8. ANNOVAR--Variant annotation with downloaded databases * must run in annovar directory
 OUTPUT: NGS0001_filtered_bed.avinput
-rw-rw-r-- 1 ubuntu ubuntu 724K Apr  9 15:04 NGS0001_filtered_bed.avinput 
	1 x CVS file downloadable and usable in excel
-rw-rw-r-- 1 ubuntu ubuntu 7.0M Apr  9 15:07 NGS0001_filtered_bed.hg19_multianno.csv

--END--	 

