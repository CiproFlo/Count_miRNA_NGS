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
my $fastq_file            = $path_to_file."34a_2_199a_2.fq_output2.fastq";



#my $adaptor_sequence = "G"; 
#my $barcode_length = 1;
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



my %fastq;
my ( $read_counter, $line_counter, $current_identifier ) = ( 0, 0, 0 );

open (FASTQ_FILE, "<$fastq_file"); 
say "Read file $fastq_file";
my $i = 0;
while (<FASTQ_FILE>){
    chomp();
    if      ($line_counter == 0) { $current_identifier                               = $_; 
                                   $fastq{$current_identifier}{identifier}           = $_;
                                   $fastq{$current_identifier}{counter}              = $read_counter;
                                   $fastq{$current_identifier}{barcode}              = 0;
    } elsif ($line_counter == 1) { $fastq{$current_identifier}{sequence}             = $_;
    } elsif ($line_counter == 2) { $fastq{$current_identifier}{optional_information} = $_; 
    } elsif ($line_counter == 3) { $fastq{$current_identifier}{quality_score}        = $_;
    } else                       { say "OOPS!";}
    die "Lineissue $_" if ( ($line_counter == 0 ) && (substr($_,0,1) ne "@"));                     ## First, short quality control
    ++$line_counter;
    if ( $line_counter > 3) { $line_counter = 0; ++$read_counter; }
	if ( $i%4000000 == 0 ) {																		## Info all 1 Mio sequecnes read
		print strftime "[%d.%m.%Y - %H:%M:%S]\t", localtime();
		printf "%d sequences read.\n", ($i/4) 

	}
	++$i;
    #last if $read_counter == 100000;                                                                  ## Remove in final run!!!!
}
close FASTQ_FILE;

say "Search for miRNAs...";
foreach my $key ( sort keys %mirnas ){
	my $output_file = $path_to_file.$key.".fasta";
	open (OUTPUT, ">>$output_file"); 
    my $current_mirna = $mirnas{$key};
	my $current_mirna_rc = reverse_complement($current_mirna);
	say "Looking for $key";
	my $i = 0;
	foreach my $key_fastaq ( sort keys %fastq){
		use re::engine::TRE max_cost => 0 ;
		if ( $fastq{$key_fastaq}{sequence} =~ /.*?$current_mirna.*?/i
										   or
			 $fastq{$key_fastaq}{sequence} =~ /.*?$current_mirna_rc.*?/i ){
			 	use re::engine::TRE max_cost => 2 ;
				$fastq{$key_fastaq}{sequence} =~ s/AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG.*?$//g;
				print OUTPUT ">".$fastq{$key_fastaq}{identifier}."\n".$fastq{$key_fastaq}{sequence}."\n" if (length($fastq{$key_fastaq}{sequence})<50);
				#say "$fastq{$key_fastaq}{sequence}";
				
			 }
	    ++$i;
    	printf "%d %%\n", $i/keys(%fastq)*100 if ($i%10_000 == 0);
	}
	no re::engine::TRE;
	close OUTPUT;
}

#print Dumper \%fastq;


#open (OUTPUT_ALL, ">>$output_file_all"); 
#open (OUTPUT_SEQ, ">>$output_file_sequences"); 
#foreach my $key ( sort { $fastq{$a}{counter} <=> $fastq{$b}{counter} } keys %fastq ){              ## Keep original order
#    print OUTPUT_ALL $fastq{$key}{identifier}           ."\n". 
#                     $fastq{$key}{sequence}             ."\n".
#                     $fastq{$key}{optional_information} ."\n".
#                     $fastq{$key}{quality_score}        ."\n";
#    print OUTPUT_SEQ $fastq{$key}{identifier}           ."\n". 
#                     $fastq{$key}{sequence}             ."\n";
#}
#close (OUTPUT_ALL); 
#close (OUTPUT_SEQ);
#
#say "__END__";
#
#
######## Subroutines #######
#
#sub sequence_quality {
#
#    return 0;               ## Diese Subroutine wird momentan nicht gebraucht und steht oben wie hier als Platzhalter drin
#
#}
#
#
#sub detect_adapter_and_barcode {
#    my ( $sequence, $identifier, ) = @_;
#    if ( $sequence =~ /(.*?)$adaptor_sequence(.*)/i ){
#        if ( length ($2) > ($barcode_length-1) ){
#            my $barcode = substr($2, 0, $barcode_length);
#            my $library = exists($barcodes{$barcode}) ? $barcodes{$barcode} : 0;
#            $fastq{$identifier}{"barcode"} = $library;
#            say "Barcode nicht erkannt" if !$library;
#            if ( $library ){
#                my $offset = length($1) + length($adaptor_sequence) + $barcode_length - 1;
#                substr($fastq{$identifier}{"sequence"},      0, $offset) = "";
#                substr($fastq{$identifier}{"quality_score"}, 0, $offset) = "";            
#            }            
#            return $library;
#        } else { return 0; }
#    } else { return 0; }
#}
#

sub reverse_complement {
	my ($rev_seq) = @_;
	$rev_seq = reverse($rev_seq);			
    $rev_seq =~ tr/ATGCatgcNn/TACGtacgNn/;	
	return $rev_seq;
}


## END OF FILE ##