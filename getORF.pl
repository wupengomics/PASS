#!/usr/bin/perl

use List::MoreUtils qw/uniq/;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
=pod

=head1 DESCRIPTION

    Get ORF sequences from gene annotation files and also prepare files for BAM convertion.

=head1 USAGE 

    getORF [options] -f <.gtf> -g <genome.fa>

    Options:
    -f    File name of gene annotation, .gtf format.
          - Recommend cufflinks to generate this file. 
    -g    Reference genome file name, fasta format.
    -h    Help message.

=head1 AUTHOR

    Contact:     Peng Wu; wupeng1@ihcams.ac.cn
    Last update: 2018-10-24

=cut

## Parsing arguments from command line
($gtffile, $genome, $help);

GetOptions(
    'f:s' => \$gtffile,
    'g:s' => \$genome,
    'h|help' => \$help
);

## Print usage  
pod2usage( { -verbose => 1, -output => \*STDERR } ) if ( $help );
( $genome and $gtffile ) or pod2usage();


## Set default
$aafile=$gtffile;
$ntfile=$gtffile;
$aafile=~s/\.gtf$/.protein.fa/;
$ntfile=~s/\.gtf$/.transcript.fa/;


$gffread ||= `which gffread`;
chomp $gffread;
$gffread or pod2usage();

$getorf ||= `which getorf`;
chomp $getorf;
$getorf or pod2usage();

`$gffread $gtffile -g $genome -w $ntfile`;
`$getorf -sequence $ntfile -find 1 -outseq $aafile`;


#use List::AllUtils qw(min max);
#$aafile="./cuff.fa.onlycuff.fa"; #format: generated by getorf, includes the CDS postion
#$ntfile="./exons.new1362.v2.fa"; #format: generated by gffread
#$gtffile="./transcripts.new1362.v2.sorted.gtf"; #format: generated by cufflinks -g, .gtf
#$gtffile="./test.gtf"; # strand +: exon order from small to big; strand -: exon order from big to small
open(AA,$aafile);#The sequence ID example: >CUFF.11.2_1 [63 - 218] gene=ENSMUSG00000033845.9
open(NT,$ntfile);#The sequence ID example: >CUFF.11.2 gene=ENSMUSG00000033845.9
open(GTF,$gtffile);

$aafile=~s/\.fa$/.out.fa/;
$ntfile=~s/\.fa$/.out.fa/;
$gtffile=~s/\.gtf$/.out.gtf/;

open(OUTAA,">$aafile"); #contain the ID which is same to the ID in transcript.fa 
open(OUTNT,">$ntfile"); #contain the ID and CDS:start-end, separated by "|"
open(OUTGTF,">$gtffile");

while(<AA>){
	chomp;
	if(/^>([^ ]*) \[(\d+) - (\d+)\] gene=(.*)/){
		print OUTAA ">$1|$4\n";
		$start{$1}=$2; #cds start in transcript sequence
		$end{$1}=$3; #cds end in transcript sequence
		@pidsplit=split/_/,$1;
		push @{$tid2pid{$pidsplit[0]}},$1;
	}else{
		print OUTAA "$_\n";
	}
}

$seq="";
while(<NT>){
	chomp;
	if(/^>([^ ]*) gene=(.*)/){
		if($seq ne ""){
			for(@{$tid2pid{$tid}}){
				
				print OUTNT ">$_|$gene|CDS:$start{$_}-$end{$_}\n";
				print OUTNT "$seq\n";
				$seq_cds=substr($seq,$start{$_}-1,$end{$_}-$start{$_}+1);
				$length{$_}=length $seq;
				
				
			}
		}
		
		$seq="";
		$tid=$1;
		$gene=$2;
	}else{
		$seq.="$_";
	}
}
if($seq ne ""){
	for(@{$tid2pid{$tid}}){
		
		print OUTNT ">$_|$gene|CDS:$start{$_}-$end{$_}\n";
		print OUTNT "$seq\n";
		$seq_cds=substr($seq,$start{$_}-1,$end{$_}-$start{$_}+1);
							
	}
}

while(<GTF>){
	chomp;
	$line=$_;
	@lines=split/\t/;

	if(/\ttranscript\t/){
		if(/gene_id "([^"]*)"; transcript_id "([^"]*)/){
			$transcript=$2;
			$gene=$1;
			$transcript2gene{$transcript}=$gene;
			for(@{$tid2pid{$transcript}}){
				$strand{$_}=$lines[6];
				print OUTGTF "$lines[0]\t$lines[1]\t$lines[2]\t$lines[3]\t$lines[4]\t$lines[5]\t$lines[6]\t$lines[7]\tgene_id \"$gene\"; transcript_id \"$_\"; gene_type \"protein_coding\"; gene_status \"NEW\"; gene_name \"$gene\"; transcript_type \"protein_coding\"; transcript_status \"NEW\"; transcript_name \"$transcript\"; level 0;\n";
			}
		}
		$exonLengthSum=0;
	}
	if(/\texon\t/){
		$exonLengthSum+=$lines[4]-$lines[3]+1;
		$exonLengthStart=$exonLengthSum-$lines[4]+$lines[3];
		if($lines[6] eq "+"){
			for(@{$tid2pid{$transcript}}){
				print OUTGTF "$lines[0]\t$lines[1]\t$lines[2]\t$lines[3]\t$lines[4]\t$lines[5]\t$lines[6]\t$lines[7]\tgene_id \"$gene\"; transcript_id \"$_\"; gene_type \"protein_coding\"; gene_status \"NEW\"; gene_name \"$gene\"; transcript_type \"protein_coding\"; transcript_status \"NEW\"; transcript_name \"$_\"; level 0;\n";
				$information1{$_}="$lines[0]\t$lines[1]\tCDS";
				$information1_2{$_}="$lines[5]\t$lines[6]";
				$information2{$_}="gene_id \"$gene\"; transcript_id \"$_\"; gene_type \"protein_coding\"; gene_status \"NEW\"; gene_name \"$gene\"; transcript_type \"protein_coding\"; transcript_status \"NEW\"; transcript_name \"$_\"; level 0;";
				push @{$position_t{$_}},$lines[3];
				push @{$position_t{$_}},$lines[4];
				push @{$position_r{$_}},$exonLengthStart;
				push @{$position_r{$_}},$exonLengthSum;
			}
		}
		if($lines[6] eq "-"){
			for(@{$tid2pid{$transcript}}){
				#print OUTGTF "$_\t",$length{$_},"\n";
				print OUTGTF "$lines[0]\t$lines[1]\t$lines[2]\t$lines[3]\t$lines[4]\t$lines[5]\t$lines[6]\t$lines[7]\tgene_id \"$gene\"; transcript_id \"$_\"; gene_type \"protein_coding\"; gene_status \"NEW\"; gene_name \"$gene\"; transcript_type \"protein_coding\"; transcript_status \"NEW\"; transcript_name \"$_\"; level 0;\n";
				$information1{$_}="$lines[0]\t$lines[1]\tCDS";
				$information1_2{$_}="$lines[5]\t$lines[6]";
				$information2{$_}="gene_id \"$gene\"; transcript_id \"$_\"; gene_type \"protein_coding\"; gene_status \"NEW\"; gene_name \"$gene\"; transcript_type \"protein_coding\"; transcript_status \"NEW\"; transcript_name \"$_\"; level 0;";
				push @{$position_t{$_}},$lines[4];
				push @{$position_t{$_}},$lines[3];
				push @{$position_r{$_}},$exonLengthStart;
				push @{$position_r{$_}},$exonLengthSum;
			}
		}
		
	}
}


$phase{0}=0;
$phase{1}=2;
$phase{2}=1;
for (keys %position_r){
	print "$_\n";
	@position2_r=@{$position_r{$_}};
	@position2_t=@{$position_t{$_}};
	print "@position2_r\n";
	print "@position2_t\n";
	if($strand{$_} eq "+"){
		for $i (1..scalar(@position2_r)/2){
			if($start{$_}>=$position2_r[2*$i-1-1] && $start{$_}<=$position2_r[2*$i-1]){
				$start_t=$start{$_}-$position2_r[2*$i-1-1]+$position2_t[2*$i-1-1];
				$start_i=$i;
			}
			if($end{$_}>=$position2_r[2*$i-1-1] && $end{$_}<=$position2_r[2*$i-1]){
				$end_t=$end{$_}-$position2_r[2*$i-1-1]+$position2_t[2*$i-1-1];
				$end_i=$i;
			}
		}
		if($start_i==$end_i){
			print OUTGTF "$information1{$_}\t$start_t\t$end_t\t$information1_2{$_}\t0\t$information2{$_}\n";
		}else{
			for $i ($start_i..$end_i){
				if($i==$start_i){
					$cds_length=0;
					
					print OUTGTF "$information1{$_}\t$start_t\t$position2_t[2*$i-1]\t$information1_2{$_}\t0\t$information2{$_}\n";
					$cds_length+=$position2_t[2*$i-1]-$start_t+1;
				}elsif($i==$end_i){
					
					print OUTGTF "$information1{$_}\t$position2_t[2*$i-1-1]\t$end_t\t$information1_2{$_}\t$phase{$cds_length%3}\t$information2{$_}\n";
					$cds_length+=$end_t-$position2_t[2*$i-1-1]+1;
				}else{
					print OUTGTF "$information1{$_}\t$position2_t[2*$i-1-1]\t$position2_t[2*$i-1]\t$information1_2{$_}\t$phase{$cds_length%3}\t$information2{$_}\n";
					$cds_length+=$position2_t[2*$i-1]-$position2_t[2*$i-1-1]+1;
				}
			}
		}
	}
	if($strand{$_} eq "-"){
		for $i (1..scalar(@position2_r)/2){
			if($start{$_}>=$position2_r[2*$i-1-1] && $start{$_}<=$position2_r[2*$i-1]){
				$start_t=$position2_r[2*$i-1-1]+$position2_t[2*$i-1-1]-$start{$_};
				$start_i=$i;
			}
			if($end{$_}>=$position2_r[2*$i-1-1] && $end{$_}<=$position2_r[2*$i-1]){
				$end_t=$position2_r[2*$i-1-1]+$position2_t[2*$i-1-1]-$end{$_};
				$end_i=$i;
			}
		}
		if($start_i==$end_i){
			print OUTGTF "$information1{$_}\t$end_t\t$start_t\t$information1_2{$_}\t0\t$information2{$_}\n";
		}else{
			for $i ($start_i..$end_i){
				if($i==$start_i){
					$cds_length=0;
					print OUTGTF "$information1{$_}\t$position2_t[2*$i-1]\t$start_t\t$information1_2{$_}\t0\t$information2{$_}\n";
					$cds_length+=$start_t-$position2_t[2*$i-1]+1;
				}elsif($i==$end_i){
					print OUTGTF "$information1{$_}\t$end_t\t$position2_t[2*$i-1-1]\t$information1_2{$_}\t$phase{$cds_length%3}\t$information2{$_}\n";
					$cds_length+=$position2_t[2*$i-1-1]-$end_t+1;
				}else{
					print OUTGTF "$information1{$_}\t$position2_t[2*$i-1]\t$position2_t[2*$i-1-1]\t$information1_2{$_}\t$phase{$cds_length%3}\t$information2{$_}\n";
					$cds_length+=$position2_t[2*$i-1-1]-$position2_t[2*$i-1];
				}
			}
		}
	}
}

close AA;
close NT;
close GTF;

close OUTAA;
close OUTNT;
close OUTGTF;
