#!/usr/bin/perl 
# Runs meme-chip on the top 1000 peaks of a single narrowPeak file or diffbind format
# Kenneth Condon Aug 2017
# Important: 	The key difference between looking at all peaks or just the 
#				top 1000 is the resulting E-values. Unlike in the "all peak" 
#				results, the "top1000" results give a clear drop off in E-value
#				after the first few peaks.

## Choosing oops/zoops/anr ###################################################

# This is where you tell MEME how you believe occurrences of the motifs are
# distributed among the sequences. Selecting the correct type of distribution
# improves the sensitivity and quality of the motif search.

## ZOOPS

# Zero or one occurrence per sequence MEME assumes that each sequence may contain
# at most one occurrence of each motif. This option is useful when you suspect
# that some motifs may be missing from some of the sequences. In that case, the
# motifs found will be more accurate than using the one occurrence per sequence
# option. This option takes more computer time than the one occurrence per
# sequence option (about twice as much) and is slightly less sensitive to weak
# motifs present in all of the sequences.

## OOPS

# One occurrence per sequence MEME assumes that each sequence in the dataset
# contains exactly one occurrence of each motif. This option is the fastest and
# most sensitive but the motifs returned by MEME may be "blurry" if any of the
# sequences is missing them.

## ANR

# Any number of repetitions MEME assumes each sequence may contain any number of
# non-overlapping occurrences of each motif. This option is useful when you
# suspect that motifs repeat multiple times within a single sequence. In that
# case, the motifs found will be much more accurate than using one of the other
# options. This option can also be used to discover repeats within a single
# sequence. This option takes much more computer time than the one occurrence per
# sequence option (about ten times as much) and is somewhat less sensitive to weak
# motifs which do not repeat within a single sequence than the other two options.


###############################################################################

use strict;
use warnings;
use Getopt::Long;

###############################################################################

my ($infile,$outdir,$mode, $help) = "" ;
GetOptions(
    'help' => \$help,
    'i=s' => \$infile,
    'o=s' => \$outdir,
    'm=s' => \$mode
  
) or die "\n**********  Incorrect usage!  ***********\nrun with -help option to see the useage\n\n"; 

sub useage { die(qq/
	USAGE : perl memechip.pl [options]
	ARGUMENTS : 
                    REQUIRED
                    -i -> input file (must be narrowPeak or diffbind format)
                    -o -> output directory
                    -m -> 	d (input diffbind format)
                                n (input MACS2 narrowPeaks format)

                    OPTIONAL
                    -help -> prints this message
                \n/);
}

if($help) { &useage ;}
if(!$infile || !$outdir || !$mode) { print "\n MISSING ARGUMENTS:\tPlease give all the required options\n" ; &useage ;}
if (($mode ne "d") and ($mode ne "n")) {print "\n Please input mode as 'd' or 'n'\n"; &useage;}

###############################################################################

if ($mode eq "n")
    {
        # prepare sequence.fa + background files
        system "sort -k9nr $infile >$outdir/sorted.peaks\n"; # -k9nr sorts by p-value, -k9nr sorts by q-value
        system "head -1000 $outdir/sorted.peaks >$outdir/top1000.peaks\n";
        system qq(awk 'BEGIN{ OFS="\\t";}{ midPos=\$2+\$10; print \$1, midPos-250, midPos+250; }' $outdir/top1000.peaks >$outdir/top1000.coords); # extract cordinates of peaks (length 500bps)
        system "fastaFromBed -fi /NGS/Software/meme_install/db/sequences/mm10_masked.fa -bed $outdir/top1000.coords -fo $outdir/top1000.sequences.fa"; # get the sequences of the peaks
        system "fasta-get-markov -m 2 -dna -nostatus -nosummary $outdir/top1000.sequences.fa $outdir/top1000.backgroundModel"; # build the background markov model
        system "rm $outdir/sorted.peaks $outdir/top1000.peaks $outdir/top1000.coords";
        
        # run meme-chip (runtime ~40min)
        system "meme-chip -oc $outdir/zoops -dna -bfile $outdir/top1000.backgroundModel -norand -meme-mod zoops -meme-nmotifs 10 -meme-minw 6 -meme-maxw 30 -meme-p 8 -spamo-skip -db /NGS/Software/meme_install/db/motifs/MOUSE/HOCOMOCOv10_MOUSE_mono_meme_format.meme -db /NGS/Software/meme_install/db/motifs/MOUSE/uniprobe_mouse.meme -db /NGS/Software/meme_install/db/motifs/JASPAR/JASPAR_CORE_2016_vertebrates.meme -db /NGS/Software/meme_install/db/motifs/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme $outdir/top1000.sequences.fa";
        # system "meme-chip -oc $outdir/oops -dna -bfile $outdir/top1000.backgroundModel -norand -meme-mod oops -meme-nmotifs 10 -meme-minw 6 -meme-maxw 30 -meme-p 8 -spamo-skip -db /NGS/Software/meme_install/db/motifs/MOUSE/HOCOMOCOv10_MOUSE_mono_meme_format.meme -db /NGS/Software/meme_install/db/motifs/MOUSE/uniprobe_mouse.meme -db /NGS/Software/meme_install/db/motifs/JASPAR/JASPAR_CORE_2016_vertebrates.meme -db /NGS/Software/meme_install/db/motifs/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme $outdir/top1000.sequences.fa";
        # system "meme-chip -oc $outdir/anr -dna -bfile $outdir/top1000.backgroundModel -norand -meme-mod anr -meme-nmotifs 10 -meme-minw 6 -meme-maxw 30 -meme-p 8 -spamo-skip -db /NGS/Software/meme_install/db/motifs/MOUSE/HOCOMOCOv10_MOUSE_mono_meme_format.meme -db /NGS/Software/meme_install/db/motifs/MOUSE/uniprobe_mouse.meme -db /NGS/Software/meme_install/db/motifs/JASPAR/JASPAR_CORE_2016_vertebrates.meme -db /NGS/Software/meme_install/db/motifs/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme $outdir/top1000.sequences.fa";
    }

elsif ($mode eq "d")
    {
        # choose a sorting function
        # system "sort -k10 $infile >$outdir/sorted.peaks"; # -k10 sorts file by p-value, -k11 sorts file by FDR
        system "sed 's/-//g' $infile >$outdir/noNegs.peaks; sort -k9 $outdir/noNegs.peaks >$outdir/sorted.peaks"; # sort by fold change
         
        # prepare sequence.fa + background files
        system "sed -i '/names/d' $outdir/sorted.peaks"; # remove header
        system "head -1000 $outdir/sorted.peaks >$outdir/top1000.diffbind.peaks\n"; # get top 1000 peaks
        system "cut -f1-3 $outdir/top1000.diffbind.peaks >$outdir/top1000.diffbind.coords"; # extract cordinates of peaks (length 500bps)
        system "fastaFromBed -fi /NGS/Software/meme_install/db/sequences/mm10_masked.fa -bed $outdir/top1000.diffbind.coords -fo $outdir/top1000.diffbind.sequences.fa"; # get the sequences of the peaks
        system "fasta-get-markov -m 2 -dna -nostatus -nosummary $outdir/top1000.diffbind.sequences.fa $outdir/top1000.diffbind.backgroundModel"; # build the background markov model
        system "rm $outdir/top1000.diffbind.peaks $outdir/top1000.diffbind.coords $outdir/sorted.peaks $outdir/noNegs.peaks";
        
        # run meme-chip (runtime ~40min)
        system "meme-chip -oc $outdir/zoops -dna -bfile $outdir/top1000.diffbind.backgroundModel -norand -meme-mod zoops -meme-nmotifs 10 -meme-minw 6 -meme-maxw 30 -meme-p 8 -spamo-skip -db /NGS/Software/meme_install/db/motifs/MOUSE/HOCOMOCOv10_MOUSE_mono_meme_format.meme -db /NGS/Software/meme_install/db/motifs/MOUSE/uniprobe_mouse.meme -db /NGS/Software/meme_install/db/motifs/JASPAR/JASPAR_CORE_2016_vertebrates.meme -db /NGS/Software/meme_install/db/motifs/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme $outdir/top1000.diffbind.sequences.fa";
        # system "meme-chip -oc $outdir/oops -dna -bfile $outdir/top1000.diffbind.backgroundModel -norand -meme-mod oops -meme-nmotifs 10 -meme-minw 6 -meme-maxw 30 -meme-p 8 -spamo-skip -db /NGS/Software/meme_install/db/motifs/MOUSE/HOCOMOCOv10_MOUSE_mono_meme_format.meme -db /NGS/Software/meme_install/db/motifs/MOUSE/uniprobe_mouse.meme -db /NGS/Software/meme_install/db/motifs/JASPAR/JASPAR_CORE_2016_vertebrates.meme -db /NGS/Software/meme_install/db/motifs/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme $outdir/top1000.diffbind.sequences.fa";
        # system "meme-chip -oc $outdir/anr -dna -bfile $outdir/top1000.diffbind.backgroundModel -norand -meme-mod anr -meme-nmotifs 10 -meme-minw 6 -meme-maxw 30 -meme-p 8 -spamo-skip -db /NGS/Software/meme_install/db/motifs/MOUSE/HOCOMOCOv10_MOUSE_mono_meme_format.meme -db /NGS/Software/meme_install/db/motifs/MOUSE/uniprobe_mouse.meme -db /NGS/Software/meme_install/db/motifs/JASPAR/JASPAR_CORE_2016_vertebrates.meme -db /NGS/Software/meme_install/db/motifs/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme $outdir/top1000.diffbind.sequences.fa";
     }
     
#####
exit;