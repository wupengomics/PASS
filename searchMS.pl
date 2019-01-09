#!/usr/bin/perl

use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
=pod

=head1 DESCRIPTION

    Search MS file against protein sequence database.

=head1 USAGE 

    searchMS [options] -b <MaxQuantCmd.exe> -p <mqpar.xml>

    Options:
    -p    Parameter configuration file.
          - Recommend to preconfigure the mqpar.xml file in MaxQuant GUI.
    -b    Path to MaxQuantCmd.exe eg. /home/user/MaxQuant/bin/MaxQuantCmd.exe.
          - By default, we try to search MaxQuantCmd.exe in system PATH.
    -h    Help message.

=head1 AUTHOR

    Contact:     Peng Wu; wupeng1@ihcams.ac.cn
    Last update: 2018-10-24

=cut

## Parsing arguments from command line
($par);

GetOptions(
    'p:s' => \$par,
    'b:s' => \$exe,
    'h|help' => \$help
);

## Print usage  
pod2usage( { -verbose => 0, -output => \*STDERR } ) if ( $help );
( $par ) or pod2usage();

$exe ||= `which MaxQuantCmd.exe`;
chomp $exe;
$exe or pod2usage();

`mono $exe $par`;