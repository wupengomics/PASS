#!/usr/bin/perl

use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
=pod

=head1 DESCRIPTION

    Sort annotation file.

=head1 USAGE 

    sortGTF [options] -f <.gtf>

    Options:
    -f    File name of gene annotation, .gtf format.
          - Recommend cufflinks to generate this file.  
    -h    Help message.

=head1 AUTHOR

    Contact:     Peng Wu; wupeng1@ihcams.ac.cn
    Last update: 2018-10-24

=cut

## Parsing arguments from command line
($gtffile);

GetOptions(
    'f:s' => \$gtffile,
    'h|help' => \$help
);

## Print usage  
pod2usage( { -verbose => 0, -output => \*STDERR } ) if ( $help );
( $gtffile ) or pod2usage();

#$gtffile="./test.gtf";
open(GTF,$gtffile);
$gtffile=~s/.gtf/.sorted.gtf/;
open(OUT,">$gtffile");

while(<GTF>){
	chomp;
	$line=$_;
	@lines=split/\t/;
	if(/\ttranscript\t/){
		if($strand eq "-"){
			@order=reverse @order;
			for(@order){
				print OUT "$_\n";
			}
		}
		@order=();
		print "$line\n";
		$strand=$lines[6];
	}
	if(/\texon\t/){
		if($lines[6] eq "+"){
			print OUT "$line\n";
		}else{
			push @order,$line;
		}
	}
}
if($strand eq "-"){
	@order=reverse @order;
	for(@order){
		print OUT "$_\n";
	}
}
@order=();

close GTF;
close OUT;
