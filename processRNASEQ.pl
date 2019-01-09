#!/usr/bin/perl

use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
=pod

=head1 DESCRIPTION

    Align RNA-Seq reads to the reference genome and reconstruct transcripts.

=head1 USAGE 

    processRNASEQ [options] -g <genome> -f <.gtf> -r <.fastq>

    Options:
    -g    Genome bowtie2 index name.
    -f    Gene annotation file, .gtf format.
    -r    File names for sequencing reads, .fastq format.
          - Compressed files (.fastq.gz) are also supported.
          - Paired-end files separated by commas.
    -t    Path to tophat, eg. /home/user/bin/tophat
          - By default, we try to search tophat in system PATH.
    -c    Path to cufflinks, eg. /home/user/bin/cufflinks
          - By default, we try to search cufflinks in system PATH. 
    -p    Number of used threads. [Default: 12]
    -o    Output folder. [Default: ./RNASEQ_output]
    -h    Help message.

=head1 AUTHOR

    Contact:     Peng Wu; wupeng1@ihcams.ac.cn
    Last update: 2018-10-24

=cut

## Parsing arguments from command line
($genome, $genegtf, $reads, $tophat, $cufflinks, $threads, $outfolder);

GetOptions(
    'g:s' => \$genome,
    'f:s' => \$genegtf,
    'r:s' => \$reads,
    't:s' => \$tophat,
    'c:s' => \$cufflinks,
    'p:i' => \$threads,
    'o:s' => \$outfolder,
    'h|help' => \$help
);

## Print usage  
pod2usage( { -verbose => 2, -output => \*STDERR } ) if ( $help );
( $genome and  $genegtf and $reads ) or pod2usage();

## Set default
$tophat ||= `which tophat`;
chomp $tophat;
$tophat or pod2usage();

$cufflinks ||= `which cufflinks`;
chomp $cufflinks;
$cufflinks or pod2usage();

$threads ||= 12;
$outfolder ||= "RNASEQ_output";

if($reads=~/,/){
    @readsfile=split/,/,$reads;
    `$tophat -p $threads -G $genegtf -o $outfolder $genome $readsfile[0] $readsfile[1]`;
}else{
    `$tophat -p $threads -G $genegtf -o $outfolder $genome $reads`;
}

`$cufflinks -o $outfolder -p $threads -g $genegtf $outfolder/accepted_hits.bam`;