# This code is to run in the HPC to pull out unmapped reads
# from the SAM files to investigate and look at the cigar string

S=SRR11021460

E=/scratch/home/lfung/PFAS_Data_NF_ref_bowtie/results/results/bowtie_alignment/${S}_bowtie_Aligned.out.bam
L=/scratch/home/lfung/PFAS_Data_NF_ref_bowtie/results/results/unmapped_reads_local/bowtie_alignment/${S}_bowtie_local_Aligned.out.bam

#1
# names unmapped in each run
samtools view -f 4 $E | cut -f1 | sort -u > 1ref_end2end_unmapped.names
samtools view -f 4 $L | cut -f1 | sort -u > 1ref_local_unmapped.names
# getting all the unmapped reads in this sample for both end2end and local in ref


#2
# find which ones unmapped in end2end, but was mapped in local (rescued)
comm -23 1ref_end2end_unmapped.names 1ref_local_unmapped.names > 2rescued.names # only print uniq file 2

# view local alignments BAM for those rescued, with CIGAR
samtools view -N 2rescued.names $L > 2rescued.sam

# table: sample, dataset, read id, CIGAR, reference probe
# awk -v s=$S 'BEGIN{OFS="\t"}{print s,"PFAS",$1,$6,$3}' 2rescued.sam > ${S}_2rescued_table.tsv

# FASTA of the reads for BLAT
# awk '{print ">"$1"\n"$10}' 2rescued.sam > ${S}_2rescued.fa

#3
comm -12 1ref_end2end_unmapped.names 1ref_local_unmapped.names > 3unmapped_both.names

# get the SAM file for the FASTA sequence to BLAT (only for 100 sequences)
samtools view -N 3unmapped_both.names $L | head -n 100 > 3unmapped_both.sam

samtools view -N 3unmapped_both.names $L | head -n 100 | awk '{print ">"$1"\n"$10}' > ${S}_3unmapped_both.fa
