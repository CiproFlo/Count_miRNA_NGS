#!/usr/bin/perl
use warnings;
use strict;
use 5.18.0;
use autodie;
#use Data::Dumper; $Data::Dumper::Sortkeys = 1;
system("clear");
use POSIX qw/strftime/;
$| = 1;												## Turn on autoflush

####### User input #######

my $path_to_file		  = "";
my $fastq_file            = $path_to_file."34a_ctxt_199a_ctxt.fq_output2.fastq";

my %mirnas = (										    ## All miRNAs are shortened at both end for 2 nt for searching
	"hsa-miR-34a-5p"  => "AGTGTCTTAGCTGG", 				## "TGGCAGTGTCTTAGCTGGTTGT", 
	"hsa-miR-34a-3p"  => "CAGCAAGTATACTG", 				## "CAATCAGCAAGTATACTGCCCT",
	"hsa-miR-199a-5p" => "GTGTTCAGACTACCT", 			## "CCCAGTGTTCAGACTACCTGTTC", 
	"hsa-miR-199a-3p" => "TAGTCTGCACATTG",				## "ACAGTAGTCTGCACATTGGTTA",
);						  

##########################


##### Expected input #####
#
# @title an optional description and used barcodes
# sequence line
# +optional repeat of title line
# quality line
#
# EXAMPLE:
#
# @NS500786:89:HCV7MBGX2:1:11101:24717:1089 1:N:0:  ATCACG oder TAGCCT
# TATAGTGGATCCGACCGTGGTGCCGTGATCACGGTATCGGATTAGGCCCATACTTATCGCTTTTCTACCTACGTCG
# +
# AAAAAEEEEEEEEAEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEE/EEEEEEEEEE</EEEEEAEAEEEEE/EEE
#
##########################



my ( $read_lines, $line_counter, $current_identifier, $current_sequence) = ( 0, 0, 0, 0 );

open (FASTQ_FILE, "<$fastq_file"); 
say "Read file $fastq_file";
while (<FASTQ_FILE>){
    chomp();
    if    ($line_counter == 0) { $current_identifier = $_; die "Lineissue $_" if ( (substr($_,0,1) ne "@") ); }
	elsif ($line_counter == 1) { $current_sequence   = $_; }
    elsif ($line_counter == 2) { 			  }
	elsif ($line_counter == 3) {  			  }
	else                       { say "OOPS!"; }
    ++$line_counter;
	if ( $line_counter > 3) {
		$line_counter = 0;
		search_mirnas($current_identifier, $current_sequence);
	}
	if ( $read_lines%4000000 == 0 ) {									 	## Info all 1 Mio sequecnes read
		print strftime "[%d.%m.%Y - %H:%M:%S]\t", localtime();
		printf "%d sequences analyzed.\n", ($read_lines/4) 
	}
	++$read_lines;
    #last if $read_lines == 100000;                                                    ## Remove in final run!!!!
}
close FASTQ_FILE;


##########################


sub search_mirnas {

	my ($identifier, $sequence) = @_;
	foreach my $key ( sort keys %mirnas ){
	my $current_mirna    = $mirnas{$key};
	my $current_mirna_rc = reverse_complement($current_mirna);
	use re::engine::TRE max_cost => 0;
	if ( $sequence =~ /.*?$current_mirna.*?/i
				   or
	     $sequence =~ /.*?$current_mirna_rc.*?/i ){
			 	use re::engine::TRE max_cost => 2;
				$sequence =~ s/AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG.*?$//g;
				if ( length($sequence) < 50 && $sequence =~ /[^N]/gi){
					my $output_file = $path_to_file.$key."_ctxt.fasta";
					open (OUTPUT, ">>$output_file"); 
					print OUTPUT ">".$identifier."\n".$sequence."\n";
					close OUTPUT;
				}
			
			 }	
	}
	
}


sub reverse_complement {
	my ($rev_seq) = @_;
	$rev_seq = reverse($rev_seq);			
    $rev_seq =~ tr/ATGCatgcNn/TACGtacgNn/;	
	return $rev_seq;
}


## END OF FILE ##