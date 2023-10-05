#process_fa.pl

#Converts fasta sequences to have the sequence on a single line

#!/usr/bin/perl -w

use strict;
use Getopt::Std;

my %OPTS;
getopts('i:o:nsm',\%OPTS);

my $file = $OPTS{"i"};
my $FQ_FILE;

open (FQ_FILE, "< $file") or die "Cannot open $file";

my $index=0;
my $num=0;
my $header;
my $seq='';
while (<FQ_FILE>) {
  chomp;
  my @words = split;
  my @first = split(//, $words[0]);
  if ($first[0] eq '>') {
    if ($index==0) {
      print "$_\n";
    } else {
      print "$seq\n";
      print "$_\n";
      $seq='';
    }
  } else {
    $seq=$seq.$_;
  }  
  $index=1;
}
print "$seq\n";

