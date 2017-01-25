# input from Reorder_pre.pl file <prefix>_sorted.txt
# input is 4 tab delimited columns contig1, contig2, position in contig1 and position in contig2

use strict;
use warnings;

#my $infileCon = $outfile;
my $infileCon = shift @ARGV;

#initialize and declare variables;

my $Con1Prev = "Contig1";
my $Con2Prev = "Contig2";
my $Con1;
my $Con2;
my $C1Pos;
my $C2Pos;
my $count = 1;
my $totalC1 = 0;
my $totalC2 = 0;
my $avgC1Pos;
my $avgC2Pos;
my $outfileCon = "linkers_sorted.txt";
my $ConLine;

open CONIN, '<', $infileCon or die;
open CONOUT, '>', $outfileCon or die;
while ($ConLine = <CONIN>) {
    chomp $ConLine;
    ($Con1, $Con2, $C1Pos, $C2Pos) = split( /\t/, $ConLine);
    if ($Con1 eq $Con1Prev && $Con2 eq $Con2Prev) {
        $totalC1 += $C1Pos;
        $totalC2 += $C2Pos;
        $count++;
    }else {
        $avgC1Pos = $totalC1 / $count; #calculate avg postion
        $avgC2Pos = $totalC2 / $count;
        #output to file
        print CONOUT "$Con1Prev\t$Con2Prev\t$count\t$avgC1Pos\t$avgC2Pos\n";
        #reset values to new contig pairs and begin counting again
        $Con1Prev = $Con1;
        $Con2Prev = $Con2;
        $totalC1 = $C1Pos;
        $totalC2 = $C2Pos;
        $count = 1;
    }
}

close CONIN;
close CONOUT;

#sort output

my @ConPairs;
my $removed_line;

open CONOUT, '<', $outfileCon or die;

while (<CONOUT>) {
    chomp;
    push @ConPairs, $_;
}

close CONOUT;

$removed_line = shift @ConPairs; #removes header line
print "$removed_line removed.\n";

#sort remaineder
my @sorted_ConPairs =   map     { $_->[0] }
                        sort    { $b->[1][0] <=> $a->[1][0] }
                        map  { [$_, [(split /\t/, $_)[2]]] } @ConPairs;

my $sortedConPairs = "pairs_sorted.txt";

open SORTEDCON, '>', $sortedConPairs or die;
print SORTEDCON"$_\n" for @sorted_ConPairs;
close SORTEDCON;

exit 0;

#calculate average is of format
#while (my $line = <>) {
#   $total +=$line;
#   $count ++=;
#}
#print "Average =", $total / $count, "\n";