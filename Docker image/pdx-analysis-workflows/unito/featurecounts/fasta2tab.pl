#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
$,="\t";
$\="\n";

my $usage = "fasta2tab [-s|--single-line] <FASTA >TAB
        -s write each sequence on a single line\n";

my $single_line = 0;
GetOptions (
        'single-line|s' => \$single_line
) or die($usage);

my $HEADER=undef;
my $first_line=undef;
while(<>){
        chomp;
        if(m/^>/){
                if ($single_line && defined($HEADER)) {
                        printf "\n";
                }
                $HEADER = $_;
                $HEADER =~ s/^>//;
                $first_line = 1;
        }else{
                die("HEADER not defined") if !defined($HEADER);
                if ($single_line) {
                        if ($first_line) {
                                printf "%s\t%s",$HEADER,$_;
                                $first_line=0;
                        }
                        else {
                                printf "%s",$_;
                        }
                }
                else {
                        print $HEADER,$_;
                }
        }
}
