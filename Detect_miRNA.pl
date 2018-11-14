#!/usr/bin/perl
use warnings;
use strict;
use 5.18.0;
use autodie;
use Data::Dumper; $Data::Dumper::Sortkeys = 1;
system("clear");


####### User input #######

my $fastq_file            = "/Volumes/Macintosh\ HD2/TEST/Overlapp_Thea.fastq";
my $output_file_all       = "/Volumes/Macintosh\ HD2/TEST/Thea_cleaned_R1_all.fq";
my $output_file_sequences = "/Volumes/Macintosh\ HD2/TEST/Thea_cleaned_R1_seq.fq";


my $adaptor_sequence = "G"; 
my $barcode_length = 1;
my %barcodes = (
	C => "C", 
	A => "A", 
	G => "G", 
	T => "T",
);

##########################


##### Expected input #####
#
# @title an optional description
# sequence line
# +optional repeat of title line
# quality line
#
# EXAMPLE:
#
# @NS500786:89:HCV7MBGX2:1:11101:24717:1089 1:N:0:ATTACTCG+AGGCTATA
# TATAGTGGATCCGACCGTGGTGCCGTGATCACGGTATCGGATTAGGCCCATACTTATCGCTTTTCTACCTACGTCG
# +
# AAAAAEEEEEEEEAEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEE/EEEEEEEEEE</EEEEEAEAEEEEE/EEE
#
##########################


my %fastq;
my ( $read_counter, $line_counter, $current_identifier ) = ( 0, 0, "Ätschibätch!" );

my $who = "Penner";

open (FASTQ_FILE, "<$fastq_file"); 
say "Read file $fastq_file";
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
    #last if $read_counter == 10;                                                                  ## Remove in final run!!!!
}
close FASTQ_FILE;

say "Check for bad sequences and recognize adaptors & barcodes...";
my $i = 0;
foreach my $key ( sort keys %fastq ){
    if ( $fastq{$key}{sequence} =~ /[^ACTG]/gi          
                        or 
         sequence_quality( $fastq{$key}{quality_score} ) 
                        or 
         !detect_adapter_and_barcode($fastq{$key}{sequence} , $key) )
    {
        say "\n\nBad sequence:"; print Dumper \$fastq{$key};
        delete $fastq{$key};
    }
    ++$i;
    printf "%d %%\n", $i/keys(%fastq)*100 if ($i%10_000 == 0); 
}

open (OUTPUT_ALL, ">>$output_file_all"); 
open (OUTPUT_SEQ, ">>$output_file_sequences"); 
foreach my $key ( sort { $fastq{$a}{counter} <=> $fastq{$b}{counter} } keys %fastq ){              ## Keep original order
    print OUTPUT_ALL $fastq{$key}{identifier}           ."\n". 
                     $fastq{$key}{sequence}             ."\n".
                     $fastq{$key}{optional_information} ."\n".
                     $fastq{$key}{quality_score}        ."\n";
    print OUTPUT_SEQ $fastq{$key}{identifier}           ."\n". 
                     $fastq{$key}{sequence}             ."\n";
}
close (OUTPUT_ALL); 
close (OUTPUT_SEQ);

say "__END__";


####### Subroutines #######

sub sequence_quality {

    return 0;               ## Diese Subroutine wird momentan nicht gebraucht und steht oben wie hier als Platzhalter drin

}


sub detect_adapter_and_barcode {
    my ( $sequence, $identifier, ) = @_;
    if ( $sequence =~ /(.*?)$adaptor_sequence(.*)/i ){
        if ( length ($2) > ($barcode_length-1) ){
            my $barcode = substr($2, 0, $barcode_length);
            my $library = exists($barcodes{$barcode}) ? $barcodes{$barcode} : 0;
            $fastq{$identifier}{"barcode"} = $library;
            say "Barcode nicht erkannt" if !$library;
            if ( $library ){
                my $offset = length($1) + length($adaptor_sequence) + $barcode_length - 1;
                substr($fastq{$identifier}{"sequence"},      0, $offset) = "";
                substr($fastq{$identifier}{"quality_score"}, 0, $offset) = "";            
            }            
            return $library;
        } else { return 0; }
    } else { return 0; }
}


## END OF FILE ##