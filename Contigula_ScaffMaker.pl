use strict;
use warnings;

# usage perl ScaffMaker.pl <input map list> <input pairs from ConDense/Consort> <output prefix> <number of linkers for positve link>
# example perl ScaffMaker.pl InList InPairs OutScaffPrefix
# add '> outfile name' to save log data in file other than screen/.out
# will also output list of contigs in scaffolds and diagnostics with prefix

my $infile = $ARGV[0];
my $pairfile = $ARGV[1];
my $outPrefix = $ARGV[2];
my $LowNumberLink = $ARGV[3];

# define variables
my $lengthfile = $outPrefix."_length.txt";
my $Con;
my $ConSeq;
my $ConLength;
my $ConIn;
my $prevCon; # retain identity of contig to be added to scaff, test matches
my $prevCon2; # for backtracking where more than one match

my $Scaffold = $outPrefix."_scaffolds.txt";
my $GoodPairs = $outPrefix."_good_pairs.txt";
my $BadPairsPos = $outPrefix."_bad_pairs_pos.txt";
my $BadPairsNum = $outPrefix."_bad_pairs_num.txt";
my $ScaffList = $outPrefix."_scaff_contig_list.txt";
my $ScaffCount = 1;
my $SeedContig = "true";
my $tempPairs = "pairs.temp";
my $ConLen;
my $PairLen;
my $usePairs = "usePairs.temp";
my $DistFromEnd = 500;
my $countGoodPair = 0;
my $UpperLimPair;
#my $LowNumberLink = 50;
my $CurrConOrient = "forward";
my $ConPair;
my $Co2oreient = "forward";
my $ConPos;
my $SkipppedLinker = "false";
my $ConCount = 1;
my $pair;
my $Inpairs;
#my $INremPairs;
my %ConLengthHash;
my %ConSeqHash;
my $Co1;
my $Co2;
my $numberPairs;
my $posCo1;
my $posCo2;
my $OUTbadPairsPos;
my $OUTbadPairsNum;
my $INmatchPairs;
my $INscaff;
my $ScaffNum = 1;
my $pairCount;
my $INusePairs;
my $Co2orient;
my $INusePair;
my $SkipPair;
my @Pairs_sorted;
my @Values_sorted;
my $testPair2;
my $LonePair;
my $PrevPairCount;
my $checkPair1;
my $checkPair2;
my $compPair1;
my $compPair2;
my $INcount;
my $INpairs;
#my $TMPpairs;
my $INConCount;
my $backMatch;
my $ConInLen;
my $INlength;
my $OUTlength;
my $INconLength;
my $INconSeq;
my $ConSeqIn;
my $OUTscaff;
#my $OUTgoodPairs;
my $OUTscaffList;
my $usePair;
my $revConSeq;
#my %Co1PosHash;
my %Co2PosHash;
my $ConLenIn;
my $CoLenIn;
my $Pair1;
my $Pair2;
my $PairNumber;
my $posPair1;
my $posPair2;
my $INNextPairs;
my $NextPair;
my $minPos;
my $maxPos;
my $UpperLimCon;
my $walkPair;
my $testPair;
my $Seed1;
my $Seed2;
my $TempPairfile = "pairs.tmp";
#my $OUTremPairs;
my $TempConFile = "contigs.tmp";
#my $INremContig;
#my $OUTremContigs;
my $OUTusePairs;
my $NextPairs;



#create length file from input map list (format is Con\tConSeq) output is file with Con\tConLength where ConLength is the character count of ConSeq


open ($INlength, '<', $infile) or die;
open ($OUTlength, '>>', $lengthfile) or die;

while ($ConInLen = <$INlength>){
    
    chomp $ConInLen;
    ($Con, $ConSeq) = split(/\t/, $ConInLen);
     #print "$Con read in\n";
        $ConLength = length($ConSeq);
        print $OUTlength "$Con\t$ConLength\n";
    
}


close $INlength;
close $OUTlength;

print "Length file created.\n";


# use subroutines for removing pairs from pair file, removing contigs from map file,
# searcing pairs file for match to current contig. Placed at end. all calls are subName ()


# BEGIN ASSEMBLY LOOP

# create hash for lookup of length data

open ($INconLength, '<', $lengthfile) or die;

    while ($ConLenIn = <$INconLength>)
    {
        chomp $ConLenIn;
        ($Con, $ConLength) = split(/\t/, $ConLenIn);
        $ConLengthHash{$Con} = $ConLength;
    }
close $INconLength;

# create hash for lookup of sequence data

open ($INconSeq, '<', $infile) or die;

while ($ConSeqIn = <$INconSeq>)
{
    chomp $ConSeqIn;
    ($Con, $ConSeq) = split( /\t/, $ConSeqIn);
    $ConSeqHash{$Con} = $ConSeq;
}
close $INconSeq;

# open output files and close at end
open ($OUTscaff, '>>', $Scaffold) or die;
open ($OUTscaffList, '>>', $ScaffList) or die;
#open ($OUTgoodPairs, '>>', $GoodPairs) or die;
#open ($OUTbadPairsPos, '>>', $BadPairsPos) or die;
#open ($OUTbadPairsNum, '>>', $BadPairsNum) or die;
 
# initialize outside of loop as conCount = 1 then count at end after removing used so loop is do while count greater than 0

# start with first contig in $infile then search $pairfile for all pairs linked
# add found pairs to temporary file
do {
    if ($SeedContig eq "true") {
        (open $INscaff, '<', $infile) or die;
        $ConIn = <$INscaff>;
        chomp $ConIn;
        ($Con, $ConSeq) = split(/\t/, $ConIn);
        close $INscaff;
    
        print $OUTscaff ">Scaff_$ScaffNum\n";
        print $OUTscaffList ">Scaff_$ScaffNum\n";
        print "$Con read in.\n";
    }
    else {
        print "New round of addition.\n";
    }
    
    
    MatchCurrentContig ();
    
    if ($pairCount == 0) {
        print $OUTscaff "$ConSeq\n";
        print $OUTscaffList "$Con\n";
        RemoveContig ();
        $ScaffNum++;
        $SeedContig = "true";
        print "$Con has no linkers.\n";
        }

    # test for goodness of pairs
    else {
        TestGoodPair ();
        # the following is only checked if there were potential linkers
    
        # if no good pairs add contig to file and print message to STDOUT
        if ($countGoodPair == 0){
                print $OUTscaff "$ConSeq\n";
                print $OUTscaffList "$Con\n";
                RemoveContig ();
                $ScaffNum++;
                $SeedContig = "true";
                print "$Con has no good linkers.\n";
                }
        # else if 1 good pair test for orientation and add to file
        elsif ($countGoodPair == 1){
                   print "$Con has 1 good linker.\n";
                    open ($INusePair, '<', $usePairs) or die;
                    while ($usePair = <$INusePair>) {
                        chomp $usePair;
                        ($Co1, $Co2, $numberPairs, $posCo1, $posCo2) = split(/\t/, $usePair);
                        if ($ConLen < 1000) {
                            print $OUTscaff "$ConSeq";
                            print $OUTscaff "NNNNNNNNNN";
                            print $OUTscaffList "$Con\t";
                        }
                        elsif ($posCo1 > 500) {
                            print $OUTscaff "$ConSeq";
                            print $OUTscaff "NNNNNNNNNN";
                            print $OUTscaffList "$Con\t";
                        }
                        else {
                            #need to reverse complement contig
                            $revConSeq = reverse $ConSeq;
                            $revConSeq =~ tr/ACGT/TGCA/;
                            print $OUTscaff "$revConSeq";
                            print $OUTscaff "NNNNNNNNNN";
                            print $OUTscaffList "rev$Con\t";
                            print "$Con reverse complemented.\n"
                        }
                        # test for orientation of Co2 and match later
                        if ($posCo2 < 500) {
                            $Co2orient = "forward";
                        }
                        else {
                            $Co2orient = "reverse";
                        }
                    }
                    close $INusePair;
                    RemoveContig ();
                    RemovePairs ();
                    $Con = $Co2;
                    $ConSeq = $ConSeqHash{$Con};
                    $ConLen = $ConLengthHash{$Con};
                    $SeedContig = "false";
                    $countGoodPair = 0;
                }
                else {
                    # if more than 1 good pair test whether overlap or bi-directional
                    # make hash of Co1 positions and hash of Co2 positions to compare
                    print "$Con has $countGoodPair linkers.\n";
                    open ($INusePairs, '<', $usePairs) or die;
                    
                    my %Co1PosHash;
                    #my %Co2PosHash;
                    
                    $PrevPairCount = $countGoodPair;
                    print "PrevPairCount is $PrevPairCount.\n";
                    while ($CoLenIn = <$INusePairs>)
                    {
                        chomp $CoLenIn;
                        ($Co1, $Co2, $numberPairs, $posCo1, $posCo2) = split(/\t/, $CoLenIn);
                        $Co1PosHash{$Co2} = $posCo1;
                        $Co2PosHash{$Co2} = $posCo2;
                        print "$Co2\t$Co1PosHash{$Co2}\t$Co2PosHash{$Co2}\n";
                    }
                    close $INusePairs;
                    
                    # test if values of Co1 are forward, reverse, or both
                    # if length of contig is <1000 assume forward orientation
                    
                    while ($PrevPairCount > 1) {
                        if ($ConLen < 1000) {
                            # assume contig under test is in forward orientation and append to file
                            
                            print $OUTscaff "$ConSeq";
                            print $OUTscaff "NNNNNNNNNN";
                            print $OUTscaffList "$Con\t";
                            RemoveContig ();
                            $SeedContig = "false";
                            
                            # order linkers based on position and test linkage
                            # @sorted = sort { $hash{$a} cmp $hash{$b} } keys %hash;
                            
                            @Pairs_sorted = sort { $Co1PosHash{$a} <=> $Co1PosHash{$b} } keys %Co1PosHash;
                            AddLinkers ();

                            
                        } #end of ConLength less than 1000 and pair count > 1
                        
                        else { # contig is greater than 1000 and pair count > 1
                            #as of now does not have the contig addition segment included, either move this up or copy contig addition
                            @Pairs_sorted = sort { $Co1PosHash{$a} <=> $Co1PosHash{$b} } keys %Co1PosHash;
                            # test that all sorted pairs on same side
                            #print "$Con for values sorted list. Hash pos. 1 is $Co1PosHash{$a} and pos. 2 is $Co1PosHash{$b}.\n";
                            #@Values_sorted = sort { $Co1PosHash{$a} <=> $Co1PosHash{$b} } values %Co1PosHash;
                            @Values_sorted = sort { $a <=> $b } values %Co1PosHash;
                            $minPos = shift @Values_sorted;
                            $maxPos = pop @Values_sorted;
                            print "Min position is $minPos.\n";
                            print "Max position is $maxPos.\n";
                            $ConLen = $ConLengthHash{$Con};
                            $UpperLimCon = $ConLen - $DistFromEnd;
                            
                            #print "Pairs array @Pairs_sorted\n";
                            #print "Values array @Values_sorted\n";

                          
                            if ($minPos < $DistFromEnd and $maxPos < $DistFromEnd){ # both at beginning, reverse contig
                                $revConSeq = reverse $ConSeq;
                                $revConSeq =~ tr/ACGT/TGCA/;
                                print $OUTscaff "$revConSeq";
                                print $OUTscaff "NNNNNNNNNN";
                                print $OUTscaffList "rev$Con\t";
                                print "$Con reverse complemented.\n";
                                RemoveContig ();
                                $SeedContig = "false";
                                AddLinkers ();
                                # set orientation of addition? replicate test above for adding more than one contig, make sub?
                                
                                
                                
                            }
                            elsif ($minPos > $UpperLimCon and $maxPos > $UpperLimCon) { # both at end, add contig
                                print $OUTscaff "$ConSeq";
                                print $OUTscaff "NNNNNNNNNN";
                                print $OUTscaffList "$Con\t";
                                RemoveContig ();
                                $SeedContig = "false";
                                AddLinkers ();
                                # as above, need to check pairs for consistency and add, end with new ending contig
                                
                            }
                            else { #contigs tested to be good so this will be where the min and max values are split and contig can extend two directions
                                # walk backward to get to beginning of scaffold
                                # take first contig in file, test for match, take low good match not current contig set to current contig, continue until no match
                                
                                
                                if ($SeedContig eq "false") {
                                    print "Interior contig $Con has bidirectional linkers. Restart with new seed contig.\n";
                                    print $OUTscaff "\n";
                                    print $OUTscaffList "next in line has bidirectional linkers.\n";
                                    $ScaffNum++;
                                    $SeedContig = "true";
                                }

                                else {  #($SeedContig eq "true")
                                    $backMatch = "true";
                                    $Con = $Pairs_sorted[0]; #get lowest linker
                                    #$Con = $Pairs_sorted[1];    # resests current contig and parameters
                                    $ConSeq = $ConSeqHash{$Con};
                                    $ConLen = $ConLengthHash{$Con};
                                    $ConPos  = $Co2PosHash{$Con};
                                
                                    do {
                                        
                                        MatchCurrentContig ();
                                        TestGoodPair (); # can perform together as always guaranteed one good match from previous
                                        
                                        if ($pairCount == 1) { # will always have at least match to previous as has passed good pair criteria
                                            $SeedContig = "true";
                                            $backMatch = "false";
                                        }
                                        else { # matches more than just previous, check if farther forward and reset to current contig
                                            
                                            open ($INusePairs, '<', $usePairs) or die;
                                            while ($walkPair = <$INusePairs>)
                                            {
                                                chomp $walkPair;
                                                ($Co1, $Co2, $numberPairs, $posCo1, $posCo2) = split(/\t/, $walkPair);
                                                $Co1PosHash{$Co2} = $posCo1;
                                                $Co2PosHash{$Co2} = $posCo2;
                                            }
                                            close $INusePairs;
                                            # have hash tables for contig positions keyed on pairing contig
                                            # can get array of sorted pairs...test remove prev pair and retain rest for testing need to confirm orientation
                                            # if pair == 2 then use non-prev pair as new $Con, else loop to get single $Con
                                            
                                            @Pairs_sorted = sort { $Co1PosHash{$a} <=> $Co1PosHash{$b} } keys %Co1PosHash;
                                            #@Values_sorted = sort { $Co1PosHash{$a} <=> $Co1PosHash{$b} } values %Co1PosHash;
                                            
                                            if ($pairCount == 2){ # if only two matches test for which is not previous con and set to current; both can be to previous if matched to two contigs test that pair matches neither before setting; will always have at least one match to prev
                                                $testPair = shift @Pairs_sorted;
                                                $testPair2 = shift @Pairs_sorted;
                                                if ($testPair ne $prevCon and $testPair ne $prevCon2) {
                                                    $prevCon = $Con;
                                                    $Con = $testPair;
                                                    $ConSeq = $ConSeqHash{$Con};
                                                    $ConLen = $ConLengthHash{$Con};
                                                    $ConPos  = $Co2PosHash{$Con};
                                                }
                                                elsif ($testPair2 ne $prevCon and $testPair2 ne $prevCon2) {
                                                    $prevCon = $Con;
                                                    $Con = $testPair2;
                                                    $ConSeq = $ConSeqHash{$Con};
                                                    $ConLen = $ConLengthHash{$Con};
                                                    $ConPos  = $Co2PosHash{$Con};
                                                }
                                                else {
                                                    $SeedContig = "true";
                                                    $backMatch = "false";
                                                }
                                            }
                                            else {
                                                # assumes that only other possibility is 3, based on our size distribution
                                                # @array = grep !/dog/, @array; to remove anything containg dog from the array
                                                # remove $prevCon from array also need to remove $prevCon2 where there was more than one testPair
                                                @Pairs_sorted = grep !/"$prevCon"/, @Pairs_sorted;
                                                @Pairs_sorted = grep !/"$prevCon2"/, @Pairs_sorted;
                                                $testPair = shift @Pairs_sorted;
                                                $testPair2 = shift @Pairs_sorted;
                                                $Seed1 = "false";
                                                $Seed2 = "false";
                                                
                                                
                                                $prevCon = $Con;
                                                $Con = $testPair;
                                                $ConSeq = $ConSeqHash{$Con};
                                                $ConLen = $ConLengthHash{$Con};
                                                $ConPos  = $Co2PosHash{$Con};
                                                
                                                MatchCurrentContig ();
                                                TestGoodPair ();
                                                
                                                if ($pairCount == 1){
                                                    $Seed1 = "true";
                                                }
                                                elsif ($pairCount == 2) {
                                                    open ($INusePairs, '<', $usePairs) or die;
                                                    while ($walkPair = <$INusePairs>)
                                                    {
                                                        chomp $walkPair;
                                                        ($Co1, $Co2, $numberPairs, $posCo1, $posCo2) = split(/\t/, $walkPair);
                                                        $Co1PosHash{$Co2} = $posCo1;
                                                        $Co2PosHash{$Co2} = $posCo2;
                                                    }
                                                    close $INusePairs;
                                                    @Pairs_sorted = sort { $Co1PosHash{$a} <=> $Co1PosHash{$b} } keys %Co1PosHash;
                                                    @Pairs_sorted = grep !/"$prevCon"/, @Pairs_sorted;
                                                    @Pairs_sorted = grep !/"$prevCon2"/, @Pairs_sorted;
                                                    $checkPair1 = shift @Pairs_sorted;
                                                    if (defined $checkPair1){
                                                        #continues
                                                    }
                                                    else {
                                                        $Seed1 = "true";
                                                    }
                                                }
                                                else {
                                                    # matches three
                                                    open ($INusePairs, '<', $usePairs) or die;
                                                    while ($walkPair = <$INusePairs>)
                                                    {
                                                        chomp $walkPair;
                                                        ($Co1, $Co2, $numberPairs, $posCo1, $posCo2) = split(/\t/, $walkPair);
                                                        $Co1PosHash{$Co2} = $posCo1;
                                                        $Co2PosHash{$Co2} = $posCo2;
                                                    }
                                                    close $INusePairs;
                                                    @Pairs_sorted = sort { $Co1PosHash{$a} <=> $Co1PosHash{$b} } keys %Co1PosHash;
                                                    @Pairs_sorted = grep !/"$prevCon"/, @Pairs_sorted;
                                                    @Pairs_sorted = grep !/"$prevCon2"/, @Pairs_sorted;
                                                    $checkPair1 = shift @Pairs_sorted;
                                                    $checkPair2 = shift @Pairs_sorted;
                                                    if (defined $checkPair2){
                                                        #continues
                                                    }
                                                    else {
                                                        $pairCount--; #logic follows one pair not two in following
                                                    }
                                                }
                                                # test second contig if present
                                                if (defined $testPair2) {
                                                    $PrevPairCount = $pairCount;
                                                    $Con = $testPair2;
                                                    $ConSeq = $ConSeqHash{$Con};
                                                    $ConLen = $ConLengthHash{$Con};
                                                    $ConPos  = $Co2PosHash{$Con};
                                                    
                                                    MatchCurrentContig ();
                                                    TestGoodPair ();
                                                    
                                                    if ($pairCount == 1){
                                                        $Seed2 = "true";
                                                    }
                                                    elsif ($pairCount == 2) {
                                                        open ($INusePairs, '<', $usePairs) or die;
                                                        while ($walkPair = <$INusePairs>)
                                                        {
                                                            chomp $walkPair;
                                                            ($Co1, $Co2, $numberPairs, $posCo1, $posCo2) = split(/\t/, $walkPair);
                                                            $Co1PosHash{$Co2} = $posCo1;
                                                            $Co2PosHash{$Co2} = $posCo2;
                                                        }
                                                        close $INusePairs;
                                                        @Pairs_sorted = sort { $Co1PosHash{$a} <=> $Co1PosHash{$b} } keys %Co1PosHash;
                                                        @Pairs_sorted = grep !/"$prevCon"/, @Pairs_sorted;
                                                        @Pairs_sorted = grep !/"$prevCon2"/, @Pairs_sorted;
                                                        $compPair1 = shift @Pairs_sorted;
                                                        if (defined $compPair1){
                                                            #continues
                                                        }
                                                        else {
                                                            $Seed2 = "true";
                                                        }
                                                    }
                                                    else {
                                                        # matches three
                                                        open ($INusePairs, '<', $usePairs) or die;
                                                        while ($walkPair = <$INusePairs>)
                                                        {
                                                            chomp $walkPair;
                                                            ($Co1, $Co2, $numberPairs, $posCo1, $posCo2) = split(/\t/, $walkPair);
                                                            $Co1PosHash{$Co2} = $posCo1;
                                                            $Co2PosHash{$Co2} = $posCo2;
                                                        }
                                                        close $INusePairs;
                                                        @Pairs_sorted = sort { $Co1PosHash{$a} <=> $Co1PosHash{$b} } keys %Co1PosHash;
                                                        @Pairs_sorted = grep !/"$prevCon"/, @Pairs_sorted;
                                                        @Pairs_sorted = grep !/"$prevCon2"/, @Pairs_sorted;
                                                        $compPair1 = shift @Pairs_sorted;
                                                        $compPair2 = shift @Pairs_sorted;
                                                        if (defined $compPair2){
                                                            #continues
                                                        }
                                                        else {
                                                            $pairCount--;
                                                        }
                                                        
                                                    }
                                                }
                                                else { # no $testPair2
                                                    $LonePair = "true";
                                                    $testPair2 = "empty";
                                                }
                                                
                                                # test combinations for moving forward
                                                if ($Seed1 eq "true" and $LonePair eq "true"){
                                                    $Con = $testPair;
                                                    $SeedContig = "true";
                                                    $backMatch = "false";
                                                }
                                                elsif ($PrevPairCount == 2 and $LonePair = "true"){
                                                    if (defined $checkPair1){
                                                        $Con = $checkPair1;
                                                    }
                                                    else {
                                                        $Con = $testPair;
                                                        $SeedContig = "true";
                                                        $backMatch = "false";
                                                    }
                                                }
                                                elsif ($PrevPairCount > 2 and $LonePair = "true"){
                                                    $Con = $testPair;
                                                }
                                                elsif ($Seed1 eq "true" and $Seed2 eq "true"){
                                                    $Con = $testPair;
                                                    $SeedContig = "true";
                                                    $backMatch = "false";
                                                }
                                                elsif ($Seed1 eq "true" and $pairCount == 2){
                                                    $Con = $compPair1;
                                                }
                                                elsif ($Seed1 eq "true" and $pairCount > 2){
                                                    $Con = $testPair2;
                                                    # SeedContig remains false and backMatch remains true
                                                }
                                                elsif ($Seed2 eq "true" and $PrevPairCount == 2){
                                                    $Con = $checkPair1;
                                                }
                                                elsif ($Seed2 eq "true" and $PrevPairCount > 2){
                                                    $Con = $testPair;
                                                    # SeedContig remains false and backMatch remains true
                                                }
                                                elsif ($pairCount == 2 and $PrevPairCount == 2) {
                                                    if ($checkPair1 eq $testPair2 and $compPair1 eq $testPair) {
                                                        $Con = $testPair;
                                                        $SeedContig = "true";
                                                        $backMatch = "false";
                                                    }
                                                    else {
                                                        $Con = $checkPair1;
                                                    }
                                                }
                                                elsif ($PrevPairCount > 2 and $pairCount == 2){
                                                    if ($checkPair1 eq $testPair2) {
                                                        $Con = $checkPair2;
                                                    }
                                                    else {
                                                        $Con = $checkPair1;
                                                    }
                                                }
                                                elsif ($pairCount > 2 and $PrevPairCount == 2){
                                                    if ($compPair1 eq $testPair){
                                                        $Con = $compPair2;
                                                    }
                                                    else {
                                                        $Con = $compPair1;
                                                    }
                                                }
                                                else {
                                                    if ($checkPair1 eq $testPair2) { #one pair in compPairs must have testPair, choose alternate as new $Con
                                                        $Con = $checkPair2;
                                                    }
                                                    elsif ($checkPair2 eq $testPair2) {
                                                        $Con = $checkPair1;
                                                    }
                                                    else { #either both match both in compPair or different, assume first checkPair is farthest
                                                        $Con = $checkPair1;
                                                    }
                                                }
                                                $prevCon = $testPair;
                                                $prevCon2 = $testPair2;
                                                
                                                #$minPos = shift @Values_sorted;
                                                #$maxPos = pop @Values_sorted;
                                                #$UpperLimCon = $ConLen - $DistFromEnd;
                                            }
                                        }
                                    }while ($backMatch eq "true'");
                                
                            } #end else bidirectional where if was SeedContig false
                                print "New seed contig $Con.\n";
                                }
                            
                                
                        }
                        $PrevPairCount--;
                        # with contig length greater than 1000 must consider orientation and possible bi-directional addition
                        # bidirectional addition not allowed for non-Seed contig, terminate addition and start again
                        # bidirectional addition for seed means it is not the end, 'walk' back to beginning contig and then proceed with addition
                        # modify the above loop to include orientations for the linker greater than zero condition
                        
                    } # end of while loop for previous pair
                    
                } # end of else more than one good linker pair to seed contig

    } #end else contig has potential linkers
    
# add in routine to count lines in infile
    open ($INConCount, '<', $infile) or die;
    $ConCount = 0;
    while (<$INConCount>) {
        $ConCount++;
    }
    close $INConCount;

}while ($ConCount > 0);
    # add final contig
    print "All contigs added.\n";
    
    #close files
    close $OUTscaff;
    close $OUTscaffList;
    #close $OUTgoodPairs;
    #close $OUTbadPairsPos;
    #close $OUTbadPairsNum;


sub MatchCurrentContig {
    print "Match $Con.\n";
    $pairCount = 0;
    open ($INpairs, '<', $pairfile) or die "can't open $pairfile.\n";
    open (my $TMPpairs, '>', $tempPairs) or die "can't open $tempPairs.\n";
    while ($pair = <$INpairs>)
    {
        chomp $pair;
        ($Co1, $Co2, $numberPairs, $posCo1, $posCo2) = split(/\t/, $pair);
        if ($Co1 eq $Con) {
            print "$Con is Co1.\n";
            print "$Co1\t$Co2\t$numberPairs\t$posCo1\t$posCo2\n";
            print $TMPpairs "$Co1\t$Co2\t$numberPairs\t$posCo1\t$posCo2\n";
            $pairCount++;
        }
        elsif ($Co2 eq $Con){
            print "$Con is Co2.\n";
            print "$Co2\t$Co1\t$numberPairs\t$posCo2\t$posCo1\n";
            print $TMPpairs "$Co2\t$Co1\t$numberPairs\t$posCo2\t$posCo1\n";
            $pairCount++;
        }
    }
    print "$Con has $pairCount pairs total.\n";
    return $pairCount;
    close $INpairs;
    close $TMPpairs;
}

# subroutine to remove pairs from file
sub RemovePairs {
    print "Remove $Con pairs.\n";
    open (my $INremPairs, '<', $pairfile) or die;
    open (my $OUTremPairs, '>', $TempPairfile) or die;
    while (<$INremPairs>) {
        #print $OUTremPairs $_ unless ($_ =~ /\b"$Con"\b/);
        print $OUTremPairs $_ unless /\b$Con\b/;
    }
    close $INremPairs;
    close $OUTremPairs;
    rename $TempPairfile,$pairfile;
    print "$Con had pairs removed.\n";
    return $Con;
}


# subroutine to remove current contig from file
sub RemoveContig {
    print "Remove contig $Con.\n";
    open (my $INremContig, '<', $infile) or die;
    open (my $OUTremContigs, '>', $TempConFile) or die;
    while (<$INremContig>) {
        #print $OUTremContigs $_ unless ($_ =~ /\b"$Con"\b/);
        print $OUTremContigs $_ unless /\b$Con\b/;
    }
    close $INremContig;
    close $OUTremContigs;
    rename $TempConFile,$infile;
    print "$Con removed.\n";
    return $Con;
}

#subroutine to test for good pairs
sub TestGoodPair {
    print "Test good pairs to $Con.\n";
    open ($INmatchPairs, '<', $tempPairs) or die "can't open tempPairs $tempPairs.\n";
    open (my $OUTusePairs, '>', $usePairs) or die "can't open usePairs $usePairs.\n";
    #open ($OUTgoodPairs, '>>', $GoodPairs) or die;
    open (my $OUTbadPairsPos, '>>', $BadPairsPos) or die;
    open (my $OUTbadPairsNum, '>>', $BadPairsNum) or die;
    $ConLen = $ConLengthHash{$Con};
    $UpperLimCon = $ConLen - $DistFromEnd;
    $countGoodPair = 0;
    #print "Contig length is $ConLen.\n";
    #print "Upper limit is $UpperLimCon.\n";
    
    while ($testPair = <$INmatchPairs>){
        #print "In while loop.\n";
        chomp $testPair;
        ($Co1, $Co2, $numberPairs, $posCo1, $posCo2) = split(/\t/, $testPair);
        $PairLen = $ConLengthHash{$Co2};
        $UpperLimPair = $PairLen - $DistFromEnd;
        print "Pair length is $PairLen.\n";
        print "Upper pair limit is $UpperLimPair.\n";
        # test if position in Contig is at end, or write to badPairs file
        # remove bad pair first
        if ($DistFromEnd <= $posCo1 and $posCo1 <= $UpperLimCon)
        {
            print $OUTbadPairsPos "$Co1\t$Co2\t$numberPairs\t$posCo1\t$posCo2\n";
        }
        elsif ($DistFromEnd <= $posCo2 and $posCo2 <= $UpperLimPair)
        {
            print $OUTbadPairsPos "$Co1\t$Co2\t$numberPairs\t$posCo1\t$posCo2\n";
        }
        elsif ($numberPairs < $LowNumberLink)
        {
            print $OUTbadPairsNum "$Co1\t$Co2\t$numberPairs\t$posCo1\t$posCo2\n";
        }
        else
        {
            print $OUTusePairs "$Co1\t$Co2\t$numberPairs\t$posCo1\t$posCo2\n";
            $countGoodPair++;
        }
        
    }
    close $INmatchPairs;
    close $OUTusePairs;
    close $OUTbadPairsPos;
    close $OUTbadPairsNum;
    print "Left with $countGoodPair pairs.\n";
    print "$Con is Con.\n";
    RemovePairs ();
    return $countGoodPair;
}

sub AddLinkers{
    $Con = $Pairs_sorted[0];    # resests current contig and parameters
    $ConSeq = $ConSeqHash{$Con};
    $ConLen = $ConLengthHash{$Con};
    $ConPair = $Pairs_sorted[1];
    $ConPos  = $Co2PosHash{$Con};
    
    if ($ConPos < 500) {
        $Co2orient = "forward";
    }
    else {
        $Co2orient = "reverse";
    } # end determine orientation of contig to be added
    
    MatchCurrentContig (); # test if first contig in line has link to second
    
    if ($pairCount == 0 and $ConLen < 1000 and $Co2orient eq "forward") { #no link
        print $OUTscaff "$ConSeq";
        print $OUTscaff "NNNNNNNNNN";
        print $OUTscaffList "$Con no linker to next\t";
        RemoveContig ();
        RemovePairs ();
        $Con = $ConPair;
        $ConSeq = $ConSeqHash{$Con};
        $ConLen = $ConLengthHash{$Con};
        $SeedContig = "false";
        print "$Con added with no internal link to next.\n"
        # update contig, subtract from PrevPairCount or at end
    }
    elsif ($pairCount == 0 and $ConLen < 1000 and $Co2orient eq "reverse") {
        $revConSeq = reverse $ConSeq;
        $revConSeq =~ tr/ACGT/TGCA/;
        print $OUTscaff "$revConSeq";
        print $OUTscaff "NNNNNNNNNN";
        print $OUTscaffList "rev$Con no linker to next\t";
        RemoveContig ();
        RemovePairs ();
        $Con = $ConPair;
        $ConSeq = $ConSeqHash{$Con};
        $ConLen = $ConLengthHash{$Con};
        $SeedContig = "false";
        print "$Con reverse complemented and added with no internal link to next.\n"
    }
    else {
        TestGoodPair (); # if first linker has a link test if to orignal second $ConPair
        
        if ($countGoodPair == 0 and $Co2orient eq "forward") { # no link
            print $OUTscaff "$ConSeq";
            print $OUTscaff "NNNNNNNNNN";
            print $OUTscaffList "$Con no linker to next\t";
            RemoveContig ();
            RemovePairs ();
            $Con = $ConPair;
            $ConSeq = $ConSeqHash{$Con};
            $ConLen = $ConLengthHash{$Con};
            $SeedContig = "false";
            print "$Con added with no internal link to next.\n"
        }
        elsif ($countGoodPair == 0 and $Co2orient eq "reverse"){
            $revConSeq = reverse $ConSeq;
            $revConSeq =~ tr/ACGT/TGCA/;
            print $OUTscaff "$revConSeq";
            print $OUTscaff "NNNNNNNNNN";
            print $OUTscaffList "rev$Con no linker to next\t";
            RemoveContig ();
            RemovePairs ();
            $Con = $ConPair;
            $ConSeq = $ConSeqHash{$Con};
            $ConLen = $ConLengthHash{$Con};
            $SeedContig = "false";
            print "$Con reverse complemented and added with no internal link to next.\n"
        }
        else {
            # ($countGoodPair >= 1)
            #  test if next in line in usePair is same as second contig from orig
            #  test if additional pair found to link past next
            
            open ($INNextPairs, '<', $usePairs) or die;
            $NextPair = <$INNextPairs>;
            chomp $NextPair;
            ($Pair1, $Pair2, $PairNumber, $posPair1, $posPair2) = split(/\t/, $NextPair);
            
            if ($Pair2 eq $ConPair and $Co2orient eq "forward") { # matches with original
                print $OUTscaff "$ConSeq";
                print $OUTscaff "NNNNNNNNNN";
                print $OUTscaffList "$Con\t";
                RemoveContig ();
                RemovePairs ();
                $Con = $ConPair;
                $ConSeq = $ConSeqHash{$Con};
                $ConLen = $ConLengthHash{$Con};
                $SeedContig = "false";
                print "$Con added with internal link to next.\n"
            }
            elsif ($Pair2 eq $ConPair and $Co2orient eq "reverse") {
                $revConSeq = reverse $ConSeq;
                $revConSeq =~ tr/ACGT/TGCA/;
                print $OUTscaff "$revConSeq";
                print $OUTscaff "NNNNNNNNNN";
                print $OUTscaffList "rev$Con\t";
                RemoveContig ();
                RemovePairs ();
                $Con = $ConPair;
                $ConSeq = $ConSeqHash{$Con};
                $ConLen = $ConLengthHash{$Con};
                $SeedContig = "false";
                print "$Con reverse complemented and added with internal link to next."
            }
            else { # should link to end of section
                if ($Co2orient eq "forward") { # doesn't match to original second, but does have match
                    print $OUTscaff "$ConSeq";
                    print $OUTscaff "NNNNNNNNNN";
                    print $OUTscaffList "$Con no internal link to next\t";
                }
                else {
                    $revConSeq = reverse $ConSeq;
                    $revConSeq =~ tr/ACGT/TGCA/;
                    print $OUTscaff "$revConSeq";
                    print $OUTscaff "NNNNNNNNNN";
                    print $OUTscaffList "rev$Con no internal link to next\t";
                }
                print "$Con linked to $Pair2 and skips $ConPair.\n";
                $SkipPair = $Pair2;
                RemoveContig ();
                $SeedContig = "false";
                
            
                # after adding first linker, reset current contig to second linker and check for match to skipped.
                RemovePairs ();
                $Con = $Pairs_sorted[1];    # resests current contig and parameters
                $ConSeq = $ConSeqHash{$Con};
                $ConLen = $ConLengthHash{$Con};
                $ConPos  = $Co2PosHash{$Con};
            
                if ($ConPos < 500) {
                    $Co2orient = "forward";
                }
                else {
                    $Co2orient = "reverse";
                }
            
                MatchCurrentContig (); # test for matches to second linker if skipped
            
                if ($pairCount == 0 and $ConLen < 1000 and $Co2orient eq "forward") {
                    print $OUTscaff "$ConSeq";
                    print $OUTscaff "NNNNNNNNNN";
                    print $OUTscaffList "$Con no linker to next\t";
                    RemoveContig ();
                    $SeedContig = "false";
                    print "$Con added with no internal link to next.\n";
                    $Con = $SkipPair;
                    $ConSeq = $ConSeqHash{$Con};
                    $ConLen = $ConLengthHash{$Con};
                    $ConPos = 0;
                    $PrevPairCount--;
                }
                elsif ($pairCount == 0 and $ConLen < 1000 and $Co2orient eq "reverse") {
                    $revConSeq = reverse $ConSeq;
                    $revConSeq =~ tr/ACGT/TGCA/;
                    print $OUTscaff "$revConSeq";
                    print $OUTscaff "NNNNNNNNNN";
                    print $OUTscaffList "rev$Con no linker to next\t";
                    RemoveContig ();
                    $SeedContig = "false";
                    print "$Con reverse complemented and added with no internal link to next.\n";
                    $Con = $SkipPair;
                    $ConSeq = $ConSeqHash{$Con};
                    $ConLen = $ConLengthHash{$Con};
                    $ConPos = 0;
                    $PrevPairCount--;
                }
                else {
                    TestGoodPair ();
                
                    if ($countGoodPair == 0 and $Co2orient eq "forward") {
                        print $OUTscaff "$ConSeq";
                        print $OUTscaff "NNNNNNNNNN";
                        print $OUTscaffList "$Con no linker to next\t";
                        RemoveContig ();
                        RemovePairs ();
                        $SeedContig = "false";
                        print "$Con added with no internal link to next.\n";
                        $Con = $SkipPair;
                        $ConSeq = $ConSeqHash{$Con};
                        $ConLen = $ConLengthHash{$Con};
                        $ConPos = 0;
                        $PrevPairCount--;
                    }
                    elsif ($countGoodPair == 0 and $Co2orient eq "reverse"){
                        $revConSeq = reverse $ConSeq;
                        $revConSeq =~ tr/ACGT/TGCA/;
                        print $OUTscaff "$revConSeq";
                        print $OUTscaff "NNNNNNNNNN";
                        print $OUTscaffList "rev$Con no linker to next\t";
                        RemoveContig ();
                        RemovePairs ();
                        $SeedContig = "false";
                        print "$Con reverse complemented and added with no internal link to next.\n";
                        $Con = $SkipPair;
                        $ConSeq = $ConSeqHash{$Con};
                        $ConLen = $ConLengthHash{$Con};
                        $ConPos = 0;
                        $PrevPairCount--;
                    }
                    elsif ($countGoodPair == 1) {
                        #  test if match to previous first (third linker) is same as match to second contig from orig
                        # if not equal add the pair under test and set test Con to new linker
                        #  set the 'third linker' to Con to test; this will end loop if 'hopscotch' continuing
                    
                        open ($INNextPairs, '<', $usePairs) or die;
                        $NextPairs = <$INNextPairs>;
                        chomp $NextPairs;
                        ($Pair1, $Pair2, $PairNumber, $posPair1, $posPair2) = split(/\t/, $NextPairs);
                    
                        if ($Pair2 ne $SkipPair and $Co2orient eq "forward") {
                            print $OUTscaff "$ConSeq";
                            print $OUTscaff "NNNNNNNNNN";
                            print $OUTscaffList "$Con\t";
                            RemoveContig ();
                            RemovePairs ();
                            $SeedContig = "false";
                            print "$Con added with no internal link to next.\n";
                            $Con = $SkipPair;
                            $ConSeq = $ConSeqHash{$Con};
                            $ConLen = $ConLengthHash{$Con};
                            $ConPos = 0;
                            $PrevPairCount--;
                        
                        }
                        elsif ($Pair2 ne $SkipPair and $Co2orient eq "reverse") {
                            $revConSeq = reverse $ConSeq;
                            $revConSeq =~ tr/ACGT/TGCA/;
                            print $OUTscaff "$revConSeq";
                            print $OUTscaff "NNNNNNNNNN";
                            print $OUTscaffList "rev$Con\t";
                            RemoveContig ();
                            RemovePairs ();
                            $SeedContig = "false";
                            print "$Con reverse complemented and added with no internal link to next.\n";
                            $Con = $SkipPair;
                            $ConSeq = $ConSeqHash{$Con};
                            $ConLen = $ConLengthHash{$Con};
                            $ConPos = 0;
                            $PrevPairCount--;
                        }
                        else {
                            # if the pair found in the first loop is also present here begin addition again.
                            print "Search for match to $Con.\n"
                        
                        }
                    } # end elsif count good pair is exactly 1 for add of second linker
                } # end else pair count is not equal to 0 for add of second linker
            } # end else ConPair neq Pair2 test for skipped pairs
        } # end of else count good pair is greater than or equal to 1 for addition of original linker
        # where pair count is 2 or more check all in same orientation should be able to do with hash table. if all in same orientation order then add and check for internal linkage. else if in both orientations implement stepwise check to find beginning of contig.
        # this is within loop where ConLen < 1000
        $PrevPairCount--; # decrease count, loop to while
        
    } # end else pair count not equal to 0 for original linker
    return $Con;
    
}


exit 0;