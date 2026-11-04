We are investigating the unmapped reads for SRR29882081 sample for Reference Dataset. We are trying to find reads that can be BLAT and reads that cannot be BLAT.
Objective: blatable=(chimeric), not blatable=(primer dimer of illumina sequence artifact and/or primer landing site )
- see how many of each there are and find percentage

BLAT by downloading BLAT to run locally because of internet api restriction they noted

`./blat hg38.2bit ref/end2end.fasta ./output_blat/end2end_output.psl`
`./blat hg38.2bit ref/local.fasta ./output_blat/local_output.psl`

You can download the 2bit file for your assembly of interest from the download server under the "Genome sequence files and select annotations (2bit, GTF, GC-content, etc)" link. For example, the hg38.2bit sequence can be found here. 

No. Command-line BLAT only aligns your query to the genome, so the PSL output contains genomic coordinates (target chr, start, end, strand, matches, mismatches, block sizes), not gene annotations.

To get gene names, intersect the hits with an annotation file:

Get a gene annotation for hg38 (e.g. GENCODE GTF, or hg38.ncbiRefSeq.gtf.gz from the same UCSC bigZips/genes folder).
Convert PSL to BED:

`pslToBed ./output_blat/end2end_output.psl ./intersect/end2end_output.bed`
`pslToBed ./output_blat/local_output.psl ./intersect/local_output.bed`

(UCSC utility; psl2bed from BEDOPS also works.)
3. Intersect:

`bedtools intersect -split -a ./intersect/end2end_output.bed -b hg38.ncbiRefSeq.gtf -wa -wb > ./intersect/end2end_hits_with_genes.txt`
`bedtools intersect -split -a ./intersect/local_output.bed -b hg38.ncbiRefSeq.gtf -wa -wb > ./intersect/local_hits_with_genes.txt`


4. Filter only exons and primary chromosomes
# awk -F'\t' '$15 == "exon"' ./intersect/end2end_hits_with_genes.txt > ./intersect/filtered/end2end_exon_hits.txt
# awk -F'\t' '$15 == "exon"' ./intersect/local_hits_with_genes.txt > ./intersect/filtered/local_exon_hits.txt

```
awk -F'\t' '$1 ~ /^chr([0-9]+|X|Y|M)$/ && $15 == "exon" {
  match($21, /gene_name "[^"]*"/); g = substr($21, RSTART+11, RLENGTH-12)
  print $4, $1, $2, $3, g
}' ./intersect/end2end_hits_with_genes.txt | sort -u > ./intersect/end2end_primary_exon_genes.txt
```
```
awk -F'\t' '$1 ~ /^chr([0-9]+|X|Y|M)$/ && $15 == "exon" {
  match($21, /gene_name "[^"]*"/); g = substr($21, RSTART+11, RLENGTH-12)
  print $4, $1, $2, $3, g
}' ./intersect/local_hits_with_genes.txt | sort -u > ./intersect/local_primary_exon_genes.txt
```


Find how many sequences for unmapped read share: CCTGTGGTCGTAGCATCAGCTAGTAGGACAGACACTCATGTA (Illumina P7)
- primer dimer
```
grep -B1 CCTGTGGTCGTAGCATCAGCTAGTAGGACAGACACTCATGTA ref/end2end.fasta
grep -c CCTGTGGTCGTAGCATCAGCTAGTAGGACAGACACTCATGTA ref/end2end.fasta

grep -B1 CCTGTGGTCGTAGCATCAGCTAGTAGGACAGACACTCATGTA ref/local.fasta
grep -c CCTGTGGTCGTAGCATCAGCTAGTAGGACAGACACTCATGTA ref/local.fasta
```


## Now we compare the names on the fasta file to those filtered ( see what excluded=primer dimer, and what included in both=mispriming DO )
## - then count to see %unBLATable and %chimeric
## - do for local and end2end
take ./ref/end2end.fasta and ./ref/local.fasta to just get the names of those reads
compare the names of those reads in uniq names in ./intersect/end2end_primary_exon_genes.txt and ./intersect/local_primary_exon_genes.txt

1.
get a list of the read names uniq on the first column and save in a txt as end2end_blatable_reads and local_blatable_reads
```
cut -d' ' -f1 intersect/end2end_primary_exon_genes.txt | sort -u > intersect/end2end_blatable_readsnames.txt
cut -d' ' -f1 intersect/local_primary_exon_genes.txt | sort -u > intersect/local_blatable_readsnames.txt
```
Check the counts to answer number of blatable reads for EXON only:
`wc -l end2end_blatable_readsnames.txt local_blatable_readsnames.txt`

2.
take those end2end.fasta and local.fasta to extract those radnames first into end2end_readnames and local_readnames in ref/
```
grep '>' ref/end2end.fasta | sed 's/^>//' | sort -u > ref/end2end_readnames.txt
grep '>' ref/local.fasta   | sed 's/^>//' | sort -u > ref/local_readnames.txt
```
compare those intersect/ readnames (blatable) with ref/ readnames (all) = get unblatable
```
mkdir -p comparison
comm -23 ref/end2end_readnames.txt intersect/end2end_blatable_readsnames.txt > comparison/end2end_not_blatable_reads.txt
comm -23 ref/local_readnames.txt   intersect/local_blatable_readsnames.txt   > comparison/local_not_blatable_reads.txt
```
check the counts to answer number of unblatable reads:
```
wc -l ref/end2end_readnames.txt intersect/end2end_blatable_readsnames.txt comparison/end2end_not_blatable_reads.txt
wc -l ref/local_readnames.txt   intersect/local_blatable_readsnames.txt   comparison/local_not_blatable_reads.txt
```


Total number of unmapped reads (end2end)./ref: 600806
Total number of unmapped reads (local)./ref: 447151
1.
Number of reads Blatable/chimeric-uniq in ./intersect (end2end):68870
Number of reads Blatable/chimeric-uniq in ./intersect (local):3807
2.
Number of reads unBlatable/primerdimer-uniq, not found in ./intersect/end2end_blatable_readsnames.txt (end2end):531936
Number of reads unBlatable/primerdimer-uniq, not found in ./intersect/local_blatable_readsnames.txt (local):443344



Final check with the original SRR29882081 sample fastq file with all the read:
- count number of reads in original file = 
```
grep -Gc '^@SRR.*' ref/SRR29882081_original.fastq
#or faster
grep -Gc '^@SRR' ref/SRR29882081_original.fastq
```
