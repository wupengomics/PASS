#!/usr/bin/perl

use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
=pod

=head1 DESCRIPTION

    Convert peptide spectal matches to alignment file.

=head1 USAGE 

    generateSAM [options] -m <PSM> -f <.gtf> -t <transcript.fa> -p <protein.fa> -o <.sam>

    Options:
    -m    Peptide spectral matches.
          - Original msms files require conversion format via preparePSM.pl.
    -f    File name of gene annotation, .gtf format.
    -t    File name of transcript sequences, .fa format.
    -p    File name of protein sequences, .fa format.
    -o    Output file, .sam format.
    -h    Help message.

=head1 AUTHOR

    Contact:     Peng Wu; wupeng1@ihcams.ac.cn
    Last update: 2018-10-24

=cut

## Parsing arguments from command line
($psm, $genesgtf, $transcript, $protein, $help);

GetOptions(
	'm:s' => \$psm,
    'f:s' => \$genesgtf,
    't:s' => \$transcript,
    'p:s' => \$protein,
    'o:s' => \$out,
    'h|help' => \$help
);

## Print usage  
pod2usage( { -verbose => 3, -output => \*STDERR } ) if ( $help );
( $psm and $genesgtf and $transcript and $protein) or pod2usage();


## Set default
$out ||= "out.sam";

open R_SCRIPT, '>', "generateSAM.R" or die $!;

print R_SCRIPT 
'library("proBAMr")
gtfFile="$genesgtf"
CDSfasta="$transcript"
pepfasta="$protein"
PrepareAnnotationGENCODE(gtfFile, CDSfasta, pepfasta,annotation_path="./", dbsnp=NULL,splice_matrix=FALSE, COSMIC=FALSE)

load("exon_anno.RData")
load("ids.RData")
load("proseq.RData")
load("procodingseq.RData")

passedPSM <- read.table("$psm", sep='\t', header=TRUE)
SAM <- PSMtab2SAM(passedPSM, XScolumn='mvh', exon, proteinseq,procodingseq)
write.table(SAM, file="$out", sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
';

`cat generateSAM.R | R --vanilla --slave`;

`rm generateSAM.R exon_anno.RData ids.RData proseq.RData procodingseq.RData`;