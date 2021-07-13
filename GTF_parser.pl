#!/bin/perl
# Written by Will Gammerdinger at IST Austria on June 25-July 2nd, 2021.
# The purpose of this script is to parse a GTF file from Transdecoder (and maybe others) by a FASTA file that only has a subset of the transcripts.
# Note: Part of this tool has an explicit assumption that the "gene" lines for each transcript variant for a given gene are identical. I don't know if this is always true, but it was true for my GTF.
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

# Make variables for GetOptions
my ( $output_file, $gtf_file, $transcriptome_input, $genome_input, $help) = ( "empty", "empty", "empty", "empty");

# Use GetOptions to take entries from the command-line and assign them to variables
GetOptions(
    "gtf_file=s"               => \$gtf_file,
    "output_file=s"            => \$output_file,
    "transcriptome_input=s"    => \$transcriptome_input,
    "genome_input=s"           => \$genome_input,
    ) or Usage ();

# If any of the input variables is undefined then kill the program and trigger this message to show the proper input format and which variables are assigned to currently
if ($output_file eq "empty" || $gtf_file eq "empty" || $transcriptome_input eq "empty" || $genome_input eq "empty"){
        die "\nERROR: The format should be perl GTF_parser.pl  --gtf_file=gtf_file.gtf --transcriptome_input=transcriptome_input.fasta --genome_input=genome_input.fasta --output_file=output_file.gtf\n\nYour current inputs are:\n\ngtf_file: $gtf_file\ntranscriptome_input: $transcriptome_input\ngenome_input: $genome_input\noutput_file: $output_file\n\n"
}

# Open genome FASTA file
open (my $GENOME_INPUT, "<$genome_input");

# Create a hash to hold the sizes of chromosomes
my %chromosome_size_hash;
# Create a variable to hold chromosome accession numbers
my $genome_accession;

# Read in each line of the FASTA genome file
while (my $line_1 = <$GENOME_INPUT>){
    # Remove the hidden \n character
    chomp $line_1;
    # If the line is a header line
    if ($line_1 =~ m/^>/){
        # Strip the > from the header line
        $genome_accession = substr($line_1, 1);
    }
    # If it is the sequence line
    else {
	# Assign the length of the sequence to a variable
	my $chromosome_size = length($line_1);
	# Add the chromosome accession as the key and chromosome size to the chromosome_size_hash
	$chromosome_size_hash{$genome_accession} = $chromosome_size;
    }
}

# Close the FASTA genome input
close $GENOME_INPUT;

# Open FASTA transcriptome file with subset of transcripts
open (my $TRANSCRIPTOME_INPUT, "<$transcriptome_input");

# Create a hash to hold the transcriptome accessions that are in the header of the FASTA header line so that we know that they exist.
my %transcriptome_accessions_exist_hash;

# Read in each line from the FASTA transcriptome file
while (my $line_2 = <$TRANSCRIPTOME_INPUT>){
    # Remove the hidden newline from the end of the line
    chomp $line_2;
    # IF a line matches the header line
    if ($line_2 =~ m/^>/){
	# Strip the > from the header line
	my $transcriptome_accession = substr($line_2, 1);
	# Enter the transcritome accession as the key in the trancript exist hash. The value of this hash has no real meaning.
	$transcriptome_accessions_exist_hash{$transcriptome_accession} = "true";
    }
    # If the line is not a FASTA header line, then 
    else {
	# Go to the next line
	next;
    }
}

# Close the FASTA transcriptome file
close $TRANSCRIPTOME_INPUT;

# Intialize a hash that will use the transcript accession ID as the key and the parent gene accession as the value
my %ID_hash;
# Initialize a hash to hold the gene accession as the key and the whole "gene" line from GTF file as the value. As I noted eariler, every transcript from this gene had an identical "gene" line, so I wasn't concerned about overwriting, but you should probably check that before running this script. 
my %gene_line_hash;
# Intialize a hash that will hold a boolean of where a transcript has been used before or not. 
my %transcript_boolean;
# Initialize a hash similar to the "gene" line hash, but this time the hash holds the transcript accession and the "transcript" line from the GTF file.
my %transcript_line_hash;
# Initialize a hash that uses the transcript ID as the key and the "stop_codon" line as the value. Initially, I was printing this line to the output, but in the script's final form it is more for checking that a stop codon exists for the gene.
my %stop_codon_line_hash;
# Create an array that holds transcript IDs where Transdecoder has annotated a stop codon that is less than 3 bases long. I don't know what the deal is with these, but we are going to remove them.
my @stop_codon_remove_array;
# Initialize a hash that uses the transcript ID as the key and the "start_codon" line as the value. Initially, I was printing this line to the output, but in the script's final form it is more for checking that a start codon exists for the gene.  
my %start_codon_line_hash;
# Create an array that holds transcript IDs where Transdecoder has annotated a stop codon that is less than 3 bases long. I don't know what the deal is with these, but we are going to remove them. 
my @start_codon_remove_array;
# SnpEff gives a warning in the "gene" line exists for multiple transcripts, so we are making a hash boolean with the gene_ID as the key and the boolean status as the value.
my %gene_used_boolean;


# Open the GTF file for the first time.
open (my $GTF_FILE_1, "<$gtf_file");

# Read in each line of the GTF file.
while (my $line_3 = <$GTF_FILE_1>){
    # My GTF file that I was reading had each GTF entry separated by new lines, but I wasn't sure if this was was always true as I have seen GFF3 files without this formating. So if the line equals new line, then:
    if ($line_3 eq "\n"){
	# Go to the next line
	next;
    }
    # Cut off the new line character from the line being read in.
    chomp $line_3;
    # Split the line on the tabs and place it into an array
    my @array_of_line_1 = split(/\t/, $line_3);
    # Split the info element in the line array by semi-colons and put it into an array
    my @array_of_info_line_1 = split(/\;/, $array_of_line_1[8]);
    # If the line is a "gene" line, then:
    if ($array_of_line_1[2] eq "gene" && $array_of_line_1[4] < $chromosome_size_hash{$array_of_line_1[0]}){
	# Remove the gene_ID label and surrounding quotations from the gene accession
	my $gene_id_1 = substr(substr($array_of_info_line_1[0], 9),0,-1);
	# Add the gene ID as the key and the line as a value to the gene line hash
	$gene_line_hash{$gene_id_1}=$line_3;
	# Add the gene_ID as the key
	$gene_used_boolean{$gene_id_1}=0;
    }
    # If the line is a "transcript" line, then:
    if ($array_of_line_1[2] eq "transcript"){
	# Remove the gene_ID label and surrounding quotations from the gene accession 
	my $gene_id_2 = substr(substr($array_of_info_line_1[0], 9),0,-1);
	# Remove the transcript_ID label and surrounding quotations from the transcript accession  
	my $transcript_id_1 = substr(substr($array_of_info_line_1[1], 16),0,-1);
	# Add the transcript ID as the key and the gene ID as the value to the ID_hash
	$ID_hash{$transcript_id_1}=$gene_id_2;
	# Add the transcript ID as the key and the transcript line as the value to the transcript line hash
	$transcript_line_hash{$transcript_id_1}=$line_3;
	# Add the transcript ID as the key and "0" as the value to the transcript boolean hash
	$transcript_boolean{$transcript_id_1}=0;
    }
    # If the line is a "stop_codon" line, then:
    if ($array_of_line_1[2] eq "stop_codon"){
	# Remove the transcript_ID label and surrounding quotations from the transcript accession  
        my $transcript_id_2 = substr(substr($array_of_info_line_1[1], 16),0,-1);
	# Input into the stop_codon_line_hash the transcript_ID as the key and the line as the value
        $stop_codon_line_hash{$transcript_id_2}=$line_3;
	# If the stop codon is shorter than 3 bases (the distances are inclusive so we need to make it less than 2); then
	# Sidenote: it is unclear to me why transdecoder does this, and it is a very small fraction of the genes, so this is just a cleaning step.
	if ($array_of_line_1[4] - $array_of_line_1[3] < 2){
	    # Add the transcript_ID to the stop_codon_remove_array
	    push(@stop_codon_remove_array, $transcript_id_2);
	}
    }
    # If the line is a "start_codon" line, then:  
    if ($array_of_line_1[2] eq "start_codon"){
	# Remove the transcript_ID label and surrounding quotations from the transcript accession 
	my $transcript_id_3 = substr(substr($array_of_info_line_1[1], 16),0,-1);
	# Input into the start_codon_line_hash the transcript_ID as the key and the line as the value   
	$start_codon_line_hash{$transcript_id_3}=$line_3;
	# If the start codon is shorter than 3 bases (the distances are inclusive so we need to make it less than 2); then
        # Sidenote: it is unclear to me why transdecoder does this, and it is a very small fraction of the genes, so this is just a cleaning step. 
	if ($array_of_line_1[4] - $array_of_line_1[3] < 2){
	    # Add the transcript_ID to the start_codon_remove_array     
            push(@start_codon_remove_array, $transcript_id_3);
        }
    }
}


# For each transcript with a stop codon less than 3 bases
foreach my $stop_transcript ( @stop_codon_remove_array ){
    # If the corresponding gene exists in the gene line hash; then 
    if ( exists $gene_line_hash{$ID_hash{$stop_transcript}} ){
	# Remove this gene from the gene line hash
	delete($gene_line_hash{$ID_hash{$stop_transcript}});
    }
}

# For each transcript with a start codon less than 3 bases  
foreach my $start_transcript ( @start_codon_remove_array ){
    # If the corresponding gene exists in the gene line hash; then  
    if ( exists $gene_line_hash{$ID_hash{$start_transcript}} ){
	# Remove this gene from the gene line hash
        delete($gene_line_hash{$ID_hash{$start_transcript}});
    }
}

# Close the GTF file
close $GTF_FILE_1;

# Open the GTF file again
open (my $GTF_FILE_2, "<$gtf_file");
# Open the output file
open (my $OUTPUT_FILE, ">$output_file");

# Initialize a boolean to identify the first entry. This is important because the GTF format I am trying to replicate has new lines between entries, but there is not an initial newline at the beginning of the GTF file.
my $first_entry_boolean = "true";

# Read in each line of the GTF file
while (my $line_4 = <$GTF_FILE_2>){
    # If the line is a new line, then:
    if ($line_4 eq "\n"){
	# Go to the next line
        next;
    }
    # Remove the hidden new line character from the end of the line
    chomp $line_4;
    # Split the line on tabs and assign it to an array
    my @array_of_line_2 = split(/\t/, $line_4);
    # Split the info element on semi-colons and assign it to an array.
    my @array_of_info_line_2 = split(/\;/, $array_of_line_2[8]);
    # If the line is the "gene" line, then:
    if ($array_of_line_2[2] eq "gene"){
	# Skip the line
	next;
    }
    # If the line is the "transcript" line then:
    elsif ($array_of_line_2[2] eq "transcript"){
	# Skip the line
	next;
    }
    # If the line is the "start_codon" line then:
    elsif ($array_of_line_2[2] eq "start_codon"){
        # Skip the line                                                                                                                                                                                  
        next;
    }
    # If the line is the "stop_codon" line then:   
    elsif ($array_of_line_2[2] eq "stop_codon"){
        # Skip the line                                                                                                                                                                                  
        next;
    }
    # If the line is the "5UTR" line then:   
    elsif ($array_of_line_2[2] eq "5UTR"){
        # Skip the line
        next;
    }
    # If the line is the "3UTR" line then:   
    elsif ($array_of_line_2[2] eq "3UTR"){
        # Skip the line
        next;
    }
    # If the line is anything else like an exon/CDS/UTR line, then:
    else {
	# Remove the transcript_ID label and surrounding quotations from the transcript accession
	my $transcript_id_4 = substr(substr($array_of_info_line_2[1], 16),0,-1);
	# If the transcript accesion exists in the transcript exists hash, then:
	if (exists $transcriptome_accessions_exist_hash{$transcript_id_4} && exists $start_codon_line_hash{$transcript_id_4} && exists $stop_codon_line_hash{$transcript_id_4} && exists $gene_line_hash{$ID_hash{$transcript_id_4}}){
	    # If the transcript boolean is "0", i.e. the script has not seen this transcript ID in the CDS before, then:
	    if ($transcript_boolean{$transcript_id_4} == 0){
		# If the first line boolean is true (i.e. this is the first entry into the GTF output file), then:
		if ($first_entry_boolean eq  "true"){
		    # Make the boolean "false" for the rest of the script
		    $first_entry_boolean =  "false"
		}
		# If the boolean is not "true" (i.e. it is not the first entry in the GTF output file), then
		else{
		    # Print a new line
		    print $OUTPUT_FILE "\n";
		}
		# If the gene has not been called yet (i.e. the gene boolean is still set to zero); then
		if ($gene_used_boolean{$ID_hash{$transcript_id_4}} == 0){ 
		    # Print the gene ID line corresponding to the gene parent of this transcript to the output file
		    print $OUTPUT_FILE "$gene_line_hash{$ID_hash{$transcript_id_4}}\n";
		    # Increase the gene boolean to 1
		    $gene_used_boolean{$ID_hash{$transcript_id_4}} = 1;   
		}
		# Print the transcript line corresponding to this transcript to the output file
		print $OUTPUT_FILE "$transcript_line_hash{$transcript_id_4}\n";
		# Change the transcript boolen from "0" to "1"
		$transcript_boolean{$transcript_id_4} = 1;
	    }
	    # Print the current line to the output file
	    print $OUTPUT_FILE "$line_4\n";
	}
    }
}

# Close the GTF file
close $GTF_FILE_2;
# close the output file
close $OUTPUT_FILE;



# If the Usage subroutine is triggered
sub Usage
{
    # Print an usage example
    printf STDERR "\nThe format should be perl GTF_parser.pl  --gtf_file=gtf_file.gtf --transcriptome_input=transcriptome_input.fasta --genome_input=genome_input.fasta --output_file=output_file.gtf\n\nYour current inputs are:\n\ngtf_file: $gtf_file\ntranscriptome_input: $transcriptome_input\ngenome_input: $genome_input\noutput_file: $output_file\n\n";
    # Exit GTF_parser.pl
    exit;
}
