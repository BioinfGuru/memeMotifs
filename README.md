# Motif discovery with [meme suite](http://meme-suite.org/index.html)

This is a work in progress.

This repository provides a useful workflow for many of the tools in meme suite including [MEME-ChIP](http://meme-suite.org/doc/meme-chip.html?man_type=web) for motif discovery in sequencing data e.g. RNA-seq, ChIP-seq and ATAC-seq. To use this tutorial you will need to install [Meme Suite](http://alternate.meme-suite.org/doc/install.html?man_type=web) and the [Meme motif databases](http://alternate.meme-suite.org/doc/install.html?man_type=web#motif_db.csv).

## Meme tools
| Name | Description |
|------|-------------|
| [Meme](http://alternate.meme-suite.org/doc/meme.html?man_type=web) | De novo motif discovery, highly specific, slow, no TF database used |
| [Dreme](http://alternate.meme-suite.org/doc/dreme.html?man_type=web) | De novo motif discovery, highly sensitive, fast, no TF database used |
| [TomTom](http://alternate.meme-suite.org/doc/tomtom.html?man_type=web) | Takes all the de novo motifs found by meme and dreme and looks in the provided TF database for similar motifs |
| [Centrimo](http://alternate.meme-suite.org/doc/centrimo.html?man_type=web) | Takes all motifs (denovo + from databases) and checks for their enrichment in your dataset at the peak centers, very useful for finding cofactors that bind with a ChiP-ed TFs |
| [Spamo](http://alternate.meme-suite.org/doc/spamo.html?man_type=web) | Uses a discovered (known or novel) motif as the "primary" motif, and each of the other discovered motifs (known or novel) as "secondary" motifs and reports the secondary motifs whose occurrences are enriched at particular distances relative to the primary motif's occurrences in the input sequences i.e. good for finding co-factors of your motif |
| [Fimo](http://alternate.meme-suite.org/doc/fimo.html?man_type=web) | Scans sequences to Find all Individual Motif Occurrences |
| [Meme-ChIP](http://alternate.meme-suite.org/doc/meme-chip.html?man_type=web) | Mmotif discovery on LARGE sets of nucleotide sequences e.g. ChIP-seq data |

## Prepare your reference fasta file
Download the latest masked reference for your species from [UCSC Golden Path](http://hgdownload.cse.ucsc.edu/downloads.html). For example, for mouse:
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFaMasked.tar.gz
tar xvzf chromFaMasked.tar.gz
cat ./chromFaMasked/* >mm10_masked.fa
```
## Prepare your sequences fasta file
Here I am using the narrowPeak file output from MACS2 as an example

```
# sort by q-value
sort -k9nr sample.narrowPeak >sample.sorted.narrowPeak
# select the top 1000 peaks
head -1000 sample.sorted.narrowPeak >sample.top1000.narrowPeaks
# create a bed file of 500bp regions centered on the peak summits
awk 'BEGIN{ OFS="\t";}{ midPos=$2+$10; print $1, midPos-250, midPos+250; }' sample.top1000.narrowPeaks >sample.regions.bed
# create fasta file
fastaFromBed -fi mm10_masked.fa -bed sample.regions.bed -fo sample.sequences.fa
```
## Prepare your background model
Here I have followed advice from the [MEME suite google forum](https://groups.google.com/forum/#!msg/meme-suite/yNascbE8Tig/rb27JMuZlwsJ;context-place=forum/meme-suite). The order of the Markov model is the number of preceding positions considered when calculating the character frequencies for the current position. A 0th order Markov model assumes that character frequencies at each position in the sequence are independent of the characters found in the previous positions. In many cases this is a reasonable assumption, but in other cases it may be an invalid assumption (CpG islands, for example). They recommend trying models up to order three (-m 3).Typically, you should not specify an order larger than 3 for DNA sequences, or larger than 2 for protein sequences. However, if your input sequences contain higher-order non-random effects that are getting in the way of motif finding, you can follow the following "rules of thumb": 

1 Use a background model at least four orders less than the shortest motifs you are looking for. So, if you want to find motifs as short as six, don't use a model higher than order two. 

2 For an accurate model of order N, the fasta-get-markov input sequences.fa file should have at least 10(4^(N+1)) DNA characters - the more the better.
* order 2 requires 640 characters
* order 3 requires 2560 characters 
* order 4 requires 10240 characters 
* order 5 requires 40960 characters

```
fasta-get-markov -m 2 -dna -nostatus -nosummary sample.sequences.fa background.model
```
### Important for ATAC-seq peaks:
For each peak pf 500bp the default is to use the centre 100 for motif searching and the remaining 400 for background. but for ATAC-seq i want to search the whole 500 bases. So I must provide a background. The max data size limit is 100000, calculated by ccut x nmeme (default: 100x600=60000). To increase max size use -meme-maxsize [600 x mean peak length] and -ccut 0, so no cutting occurs and all of each region is used. Assuming MACS2 was used to call peaks:

* calculate the mean peak length from narrow_peak.xls (column 4)
* get peak summit from narrow_peak.xls (column 5) or narrow_summits.bed (column 3)
* bed file region length = midPos+/- (1/2 mean peak length)
* create background file (fasta-get-markov -m 2)
* add to meme-chip command:		-ccut = 0 --meme-maxsize = [600 x mean peak length]

## Basic meme-chip command
```
meme-chip
-oc memechip_out /
-dna /
-bfile background.model / # or '-order N' if no background model
-norand / 
-meme-mod zoops /
-meme-nmotifs 10 /
-meme-minw 6 /
-meme-maxw 30 /
-meme-p 12 /
-spamo-skip / # disables spamo, takes too long, just run on individual motifs you are interested in from the results
-db meme_install/db/motifs/MOUSE/HOCOMOCOv10_MOUSE_mono_meme_format.meme /
-db meme_install/db/motifs/MOUSE/uniprobe_mouse.meme /
-db meme_install/db/motifs/JASPAR/JASPAR_CORE_2016_vertebrates.meme / 
-db meme_install/db/motifs/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme /
sample.sequences.fa
```
## Identifying co-factors from meme-chip results
* click "distribution" graph of top motif to display centrimo output, the list of enriched motifs is ranked by E-value
* the top motifs (novel and known) will be similar - hover over the id to see the PWM
* move down the list until you see a dissimilar motif, and select it
* copy the list of matching sequences
* use this list in go enrichment and pathway analysis
* copy the consensus sequence + make a note of the program (meme or dreme?) + alt ID
* go to main results html --> programs --> TomTom link (1st = meme, 2nd = dreme)
* search for 'consensus sequence' or 'alt ID' --> scroll down to "matches to 'consensus sequence'"
* NAME --> this protein is the co-factor --> click name for more information from the selected TF database --> UniProt ID

## Identifying co-regulators from meme-chip results
* the top hit is the most enriched motif
* the next best motif may be binding site of a co-regulator
* If novel (meme-dreme) --> show more --> are the other motifs in the cluster known (i.e. found by CentriMo)?
* In the cluster, what is the most common motif in "Known or Similar Motifs" column? - this is the potential co-regulator.

## Extract promoter peaks from a peak file
1 Create a file containing the coordinates of all promotors in your species.
* Download TSS coords from [biomart](https://www.ensembl.org/biomart/martview/cdbbc908631658bc94eaeee55643e20b) to tss.txt for each chromosome
* Remove the header
* Sort and adjust the TSS to include promoter regions:
```
sed -i -e 's/^/chr/' tss.txt
sort -k1,1 -k2,2n tss.txt > tss.sorted.txt
awk '{ print $1"\t"$2"\t"$2}' tss.sorted.txt >tss.sorted.bed
for i in tss.sorted.bed; do awk -F $'\t' 'BEGIN {OFS=FS}{{if ($4 > 0){$2 = $2 - 2000; $3 = $3 + 200;}; if ($4 < 0) {$2 = $2 - 200; $3 = $3 + 2000;}; if ($2 < 0) {$2 = 0}} print $0}' $i >$i.adj; done;
cut -f 1-3 tss.sorted.bed.adj > pro.bed
```

2 Filter the peak file
```
sed -i '/seqnames/d' peaks.bed
sort -u -k1,1 -k2,2n peaks.bed
awk '$9 >= 2 || $9 <= -2' peaks.bed > peaks.filtered.bed # removees peaks with cold change < 2
```
3 Extract promoter peaks
```
# replace "-wa -u" with  "-v" to get peaks not occuring in enhancer regions
bedtools intersect -a peaks.filtered.bed -b pro.bed -wa -u > promoter.peaks.bed
```
