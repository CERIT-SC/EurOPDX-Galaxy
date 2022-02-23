#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
$,="\t";
$\="\n";

$SIG{__WARN__} = sub {die @_};

my $usage = "$0 [-h|help] [-n|dry_run] '*someting' > matrix\nconcatenate the content of *someting files, each bloch identified by a fasta header '>P' where P is th pattern matched by * in the argument.\n\t-l|limit S\t\tspace separed list of filenames, limit the glob of '*someting' to those file, usefull to exclude some of the files that match the glob.\n";

my $help=0;
my $dry_run=0;
my $limit="";
my %limit=();
GetOptions (
        'h|help' => \$help,
        'n|dry_run' => \$dry_run,
        'l|limit=s' => \$limit
) or die($usage);

if($help){
        print $usage;
        exit(0);
}

if($limit){
        %limit =  map { $_ => 1 } split(/\s+/,$limit);
}

die("Unexpected argument number") if scalar(@ARGV)!=1;

my $pattern=shift(@ARGV);
my $regexp=$pattern;
$regexp=~s/\*/(.*)/g;

for(glob($pattern)){
        if($limit and not $limit{$_}){
                next
        }
        my @caputre = (m/$regexp/);
        my $prefix = join(";",@caputre);
        print ">$prefix";
        
        my $cat="cat";
        $cat="zcat" if m/\.gz$/;

        system("$cat '$_'") if not $dry_run;
}
