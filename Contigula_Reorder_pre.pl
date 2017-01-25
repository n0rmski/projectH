use strict;
use warnings;

# extract columns 3, 4, 7, and 8 from linkers.pre file
# has linked contigs and positions of reads
# usage Reorder_pre.pl <infile>
# example perl Reorder_pre.pl linkers.pre

my $infile = shift @ARGV;

my $prefix = $infile;
my $pathname = quotemeta ('/home/paul/PH_Contigula/forLG/');
$prefix =~s/$pathname//; #remove pathname from file
$prefix =~s/\.pre//; #remove .pre from filename to make prefix
#print $prefix\n;
my $outfile = $prefix."_extraction.txt";
my @columns;

open IN, '<', $infile  or die "File not found";
open OUT, '>', $outfile or die;
while( <IN> ) {
    @columns = split /\t/, $_;
    print OUT "$columns[2]\t$columns[6]\t$columns[3]\t$columns[7]\n";
}
close IN;
close OUT;


print "Columns extracted\n";

# remove 'tom'

#my $noTom = $prefix."_extr_noTom.txt";

#open TOM, '<', $outfile or die;
#open NOTOM, '>', $noTom or die;

#while (<TOM>) {
#   tr/tom//d;
#   print NOTOM;
#}
#close TOM;
#close NOTOM;

#print "'tom' removed\n";

# reorder columns 1 and 2 so that smallest value is in column 1 also reorder 3
# and 4 to correspond

my $reordered = $prefix."_reordered.txt";
my $odd = $prefix."_odd.txt";
my $contig1;
my $contig2;
my $posCon1;
my $posCon2;
my $line;

#open REORDER, '<', $noTom or die;
open REORDERTXT, '<', $outfile or die;
open REORDERED, '>', $reordered or die;
open ODD, '>', $odd or die;

#while ($line = <REORDER>) {
while ($line = <REORDERTXT>) {
    chomp $line;
    ($contig1, $contig2, $posCon1, $posCon2) = split( /\t/, $line);
    if ($contig2 eq "="){
        print ODD "$contig1\t$contig2\t$posCon1\t$posCon2\n";
    }
    #elsif ($contig1 < $contig2) {
    elsif ($contig1 lt $contig2) {
        print REORDERED "$contig1\t$contig2\t$posCon1\t$posCon2\n";
    }
    #elsif ($contig1 > $contig2) {
    elsif ($contig1 gt $contig2) {
        print REORDERED "$contig2\t$contig1\t$posCon2\t$posCon1\n";
    }
    else {
        print ODD "$contig1\t$contig2\t$posCon1\t$posCon2\n";
    }
}
print "reordered\n";

# sort on values first in Column 1 and then in Column 2

my $sorted = $prefix."_sorted.txt";
my @linkers;
#my $sorted_row = $prefix."sorted_row.txt";

open UNSORTED, '<', $reordered or die;

while (<UNSORTED>) {
    chomp;
    push @linkers, $_;
}

close   UNSORTED;
my @sorted_linkers =    map     { $_->[0] }
                        sort    { $a->[1][0] cmp $b->[1][0]
                                            ||
                                  $a->[1][1] cmp $b->[1][1]}
                        map  { [$_, [(split /\t/, $_)[0,1,2,3]]] } @linkers;

open SORTED, ">", $sorted or die;
    print SORTED "$_\n" for @sorted_linkers;
close SORTED;

#while (<UNSORTED>) {
#   push @linkers, [split /\t/];
#   }

#close UNSORTED;
#my @sorted;

#my @sorted = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @linkers;
#    @sorted = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] ||$a->[2] <=> $b->[2] || $a->[3] <=> $b->[3] } @linkers;
#my @sorted = sort { $a->[0] cmp $b->[0] || $a->[1] cmp $b->[1] } @linkers;
#my @sorted = sort { $a->[0] cmp $b->[0] || $a->[1] cmp $b->[1] ||$a->[2] <=> $b->[2] || $a->[3] <=> $b->[3]} @linkers;

#open SORTED, '>', $sorted or die;
#foreach my $row ( @sorted ) {
#        print SORTED join ("\t", @{$row});
#    }

#close SORTED;

#open SORTED_ROW, '>', $sorted_row or die;
#foreach my $row ( @sorted ) {
#    print SORTED_ROW $row->[0], "\t", $row->[1], "\t", $row->[2], "\t", $row->[3], "\n";
#}
#close SORTED_ROW;

print "sorted\n";



exit 0;