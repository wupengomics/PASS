#!/usr/bin/perl
use List::MoreUtils qw/uniq/;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
=pod

=head1 DESCRIPTION

    Prepare PSM from MaxQuant msms.txt results

=head1 USAGE 

    preparePSM -m <msms.txt>

    Options:
    -m    File name of peptide spectral matches.
          - Recommend MaxQuant to generate this file. 
    -h    Help message.

=head1 AUTHOR

    Contact:     Peng Wu; wupeng1@ihcams.ac.cn
    Last update: 2018-10-24

=cut

## Parsing arguments from command line
($msms);

GetOptions(
    'm:s' => \$msms,
    'h|help' => \$help
);

## Print usage  
pod2usage( { -verbose => 0, -output => \*STDERR } ) if ( $help );
( $msms ) or pod2usage();




open(MSMS,$msms);
$msms=~s/\.\w+$//;
open(OUT,">$msms.psm.tab");

print OUT "spectrum\tspectrumNativeID\tassumed_charge\thit_rank\tpeptide\tnum_missed_cleavages\tmvh\tmodification\tNTT\n";


while(<MSMS>){
	chomp;
	@l=split/\t/;
	if(/^Raw file/){
		for $i (0..$#l){
			if($l[$i] eq "Raw file"){
				$rawfile=$i;
			}elsif($l[$i] eq "Scan number"){
				$scannumber=$i;
			}elsif($l[$i] eq "Sequence"){
				$sequence=$i;
			}elsif($l[$i] eq "Missed cleavages"){
				$missedcleavages=$i;
			}elsif($l[$i] eq "Modifications"){
				$modifications=$i;
			}elsif($l[$i] eq "Charge"){
				$charge=$i;
			}elsif($l[$i] eq "Score"){
				$score=$i;
			}
		}
	}else{
		print OUT "$l[$rawfile].$l[$scannumber].$l[$scannumber].$l[$charge]\t";
		print OUT "controllerType=0 controllerNumber=1 scan=$l[$scannumber]\t";
		print OUT "$l[$charge]\t";
		print OUT "1\t";
		print OUT "$l[$sequence]\t";
		print OUT "$l[$missedcleavages]\t";
		print OUT "$l[$score]\t";
		print OUT "$l[$modifications]\t";
		print OUT "2\n";
	}
}

close MSMS;
close OUT;