% fastq.gz: zip raw reads      K38_ACAGTG_L001_R1_001.fastq.gz     1 GB  *
% 
% unzip [1 min]
% fastq: raw reads             K38_ACAGTG_L001_R1_001.fastq        3 GB
%
% ??? [1 min]
% filtered fastq                                                   1.5 GB 
% (reads with all phred score > 20) K38_filtered.fastq
%
% bowtie < fastq, indexed genome NewmapSAP (dos cmd) [30 min]
% SAM: read alignment          K38_NewmanSAP.sam                   1 GB
%
% samtools (dos cmd) < sam [2 min]
% BAM: binary alignment        K38_NewmanSAP.bam                   800 MB
% 
% samtool < bam [5 min]
% sorted bam                   K38_newmanSAP.sorted.bam            500 MB *
% 
% mpileup (cygwin) < sorted bam, NewmanSAP.txt [1 min]
% bcf: binary call file        tmp1                                500 MB
%
% bcftools view <bcf [1 min]
% raw vcf:                     K38_newmanSAP.vcf                   200 MB *
%
% bcftools view                
% variant vcf                  K38_newmanSAP.vcf                   2 MB *
%
%
%
%
