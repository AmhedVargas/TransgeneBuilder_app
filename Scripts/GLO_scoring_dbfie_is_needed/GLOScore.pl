#!/usr/bin/perl

#####Documentation######
$num_args = $#ARGV + 1;
if ($num_args != 1) {
  print "\nUsage: Score_Dans.pl sequences_fasta\n\n";
  print "\tSeq_fasta: Fasta file containing sequence to score. One per line\n";
  print "\tPlease note this program does not round values and so it differs slightly to WebVersion. Also sequence_lib_scores.db should be present in same directory and packages should be properly installed\n";

  exit;
}
############
#Load libraries

###Installed by getting tar.gz and modifying its config file. Also BerkeleyDB was downloaded from: http://www.linuxfromscratch.org/blfs/view/svn/server/db.html
use BerkeleyDB;
use Bio::Seq;
use Math::Random;
use IO::Handle;

##Reading and verifying input commands
my $Sqinps=$ARGV[0];

###If the following files are not present -> stop the program
if (!(-e $Sqinps."")){die $Sqinps." file not found!\n";}

##Initialize parameters
my $error = '';

# Get database
my $dbfilename = 'sequence_lib_scores.db';
my $sequence_lib = new BerkeleyDB::Btree
	-Filename => $dbfilename
	or die "Cannot open $dbfilename: $! $BerkeleyDB::Error\n" ;
        
# Get input
my $dnaseq = '';
my $seqtype = 'DNA';

### If a nucleotide sequence was entered, calculate its score ###

my ( $input_sequence_score, $input_lowest_score, $input_n_w_lowest_score );

my $results = {}; #Get a pointer to an empty array that will hold the results

my @input_coding_sequence;
##Reading and hashing genomic sequence
open(op,$Sqinps) or die "cannot open Ref file\n";
while($line=<op>){
chomp($line);

$results = {}; #Get a pointer to an empty array that will hold the results
@input_coding_sequence = unpack("(A3)*", $line);
( $results->{'input_sequence_score'}, $results->{'input_lowest_score'}, $results->{'input_n_w_lowest_score'} ) = score_sequence( \@input_coding_sequence, $sequence_lib );
 

print $line."\t".$results->{'input_sequence_score'}."\t".$results->{'input_lowest_score'}."\t".$results->{'input_n_w_lowest_score'}."\n";

}
close(op);




##########################################################################################
#Define functions

# Syntax: ( $sequence_score, $lowest_scoring_word, $numer_with_lowest_score ) = score_sequence( \@sequence, $pointer_to_BerkeleyDB_word_library )

sub score_sequence {
    my $CDSref = shift;
    my $CDS_length = @$CDSref;
    my $sequence_lib = shift;
    
    # Calculate the total score and the lowest word score for the sequence
    my $totalscore = 0;
    my $lowscore;
    my $words_w_lowscore;
    B: for ( my $b = 0; $b <= ($CDS_length); $b += 4 ) {
        my $wordscore = score_word($b, $CDSref, $sequence_lib);
        $totalscore += $wordscore;
        if ( $b == 0 || $wordscore < $lowscore) {
            $lowscore = $wordscore;
            $words_w_lowscore = 1;
        } elsif ($wordscore == $lowscore) {
            $words_w_lowscore++;
        }
    }
    
    # Normalize the score to the sequence length
    my $normscore = $totalscore / $CDS_length;
    
    # Return the results
    return ( $normscore, $lowscore, $words_w_lowscore );
}

###########################################################################################

#Syntax: $wordscore = score_word( $position, \@sequence, $pointer_to_BerkeleyDB_word_library )

sub score_word {
    my $position = shift;
    my $CDSref = shift;
    my $length = @$CDSref;
    my $sequence_lib = shift;

    #Score for word that overlaps by 3 nt on the 5' side
    my $L1score = 0;
    if ( $position >= 3 && $position <= $length - 1) {
        my $L1word = join('', @$CDSref[$position-3..$position]);
        $sequence_lib -> db_get($L1word, $L1score);
    }
    
    #Score for word that overlaps by 6 nt on the 5' side
    my $L2score = 0;
    if ( $position >= 2 && $position <= $length - 2 ) {
        my $L2word = join('', @$CDSref[$position-2..$position+1]);
        $sequence_lib -> db_get($L2word, $L2score);
    }
    
    #Score for word that overlaps by 9 nt on the 5' side
    my $L3score = 0;
    if ( $position >= 1 && $position <= $length - 3 ) {
        my $L3word = join('', @$CDSref[$position-1..$position+2]);
        $sequence_lib -> db_get($L3word, $L3score);
    }
    
    #Score for central word
    my $mainscore = 0;
    if ( $position <= $length - 4 ) {
        my $mainword = join('', @$CDSref[$position..$position+3]);
        $sequence_lib -> db_get($mainword, $mainscore);
    }
        
    #Score for word that overlaps by 9 nt on the 3' side
    my $R3score = 0;
    if ( $position <= $length - 5 ) {
        my $R3word = join('', @$CDSref[$position+1..$position+4]);
        $sequence_lib -> db_get($R3word, $R3score);
    }
    
    #Score for word that overlaps by 6 nt on the 3' side
    my $R2score = 0;
    if ( $position <= $length - 6 ) {
        my $R2word = join('', @$CDSref[$position+2..$position+5]);
        $sequence_lib -> db_get($R2word, $R2score);
    }

    #Score for word that overlaps by 3 nt on the 3' side
    my $R1score = 0;
    if ( $position <= $length - 7 ) {
        my $R1word = join('', @$CDSref[$position+3..$position+6]);
        $sequence_lib -> db_get($R1word, $R1score);
    }
    
    #Add up scores
    my $wordscore = 0.25*$L1score + 0.5*$L2score + 0.75*$L3score + $mainscore + 0.75*$R3score + 0.5*$R2score + 0.25*$R1score;
    
    return $wordscore;
}