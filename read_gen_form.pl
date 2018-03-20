#!/usr/bin/perl
use strict;
use warnings;
use CGI ":standard";
use CGI::Carp qw(fatalsToBrowser);
use LWP::Simple;	#module for accessing web data

my $cgi = new CGI;      						# read in parameters
my @attributes = $cgi->param('attributes');		# store multiple attributes in array
my $viruses = $cgi->param('viruses');			# get virus (biological virus)
my $checkALL = $cgi->param('CHECKALL');
my $email = $cgi->param('mailto');				# get email (optional)

my $genFile;						# virus genbank file
my $attributeRef = \@attributes;	# reference for attributes to be found in genbank file
my $baseRef;						# sequence of file to be sent in for base count
my $processedFile;					# file attributes segregated

sub getGenFile($);		# gets the web file for specified organism from genbank
sub parseFile($$);		# Parses the file for specified atribute fields and appends them to a string
sub countBases($);		# Counts the number of bases in the DNA Sequence of a given file
sub sendMail($$);		# Prints the genbank file into a text file and emails it as an attachment

#########################
#     Main Program      #
#########################

# Print the HTML Page
print $cgi->header( ); 
print $cgi->start_html( );

# Get the corresponding Genbank file from the web
$genFile = getGenFile($viruses);

# Process the file to include only selected, segregated parameters
if($checkALL eq "on"){													# If CheckALL is selected, send the unmodified file containing all data
	print $genFile;														# Print the unprocessed file to HTML
}else{																	# If CheckALL is not selected, send ONLY the selected attributes
	$processedFile = parseFile($genFile, $attributeRef);				# Process the results to include only selected attributes, segregated
	print $processedFile;												# Print the processed file to HTML
}

# If the user entered an email ending in @myseneca
# Send the processed file as an attachment
if($email =~ /\@myseneca.ca$/){
	if($checkALL eq "on"){										# If CheckALL is selected, send the unmodified file 
		sendMail($email, $genFile);
	}else{														# If CheckALL is not selected, send the processed
		sendMail($email, $processedFile);								
	}	
}
print $cgi->end_html( );

#########################
#     Subroutines       #
#########################

sub getGenFile($){
	# Gets the web file for specified organism from genbank
	my $getURL = shift;	# get organism to retrieve the genbank file for
	my $webData;	
		
	# get web file from the URL
	$webData = get($getURL);
	$webData = "<pre>" . $webData . "</pre>";
	return $webData;
}

sub parseFile($$){
	# parses the file for specified atribute fields and appends them to a string
	my $rawFile = shift;			# file to be parsed for attributes
	my $arrayRef = shift;			# reference to array of attributes to be checked
	my @attributes = @$arrayRef;	# attributes to be checked
	my @lines = ();					# lines extracted by match 
	my $fullDoc = "<pre>";			# string containing all lines extracted

	foreach my $att(@attributes){															# segregate by attribute
		if ($att eq "LOCUS"){ 																# for example, if the attribute is locus
			@lines = $rawFile =~ /(LOCUS.*?)DEFINITION/gs;									# extract the portion of text beginning with LOCUS and ending with DEFINITION
			$fullDoc .= "@lines";															# and add it to the string
		}
		if ($att eq "DEFINITION"){ 
			@lines = $rawFile =~ /(DEFINITION.*?)ACCESSION/gs;
			$fullDoc .= "@lines";
		}
		if ($att eq "ACCESSION"){ 
			@lines = $rawFile =~ /(ACCESSION.*?)VERSION/gs;
			$fullDoc .= "@lines";
		}
		if ($att eq "VERSION"){ 
			@lines = $rawFile =~ /(VERSION.*?)DBLINK/gs;
			$fullDoc .= "@lines";
		}
		if ($att eq "KEYWORDS"){ 
			@lines = $rawFile =~ /(KEYWORDS.*?)SOURCE/gs;
			$fullDoc .= "@lines";
		}
		if ($att eq "SOURCE"){ 
			@lines = $rawFile =~ /(SOURCE.*?)ORGANISM/gs;
			$fullDoc .= "@lines";
		}
		if ($att eq "ORGANISM"){ 
			@lines = $rawFile =~ /(ORGANISM.*?)REFERENCE/gs;

			$fullDoc .= "@lines";
		}
		if ($att eq "REFERENCE"){ 															# REFERENCE has a corner case
			@lines = $rawFile =~ /(REFERENCE.*?)(?=AUTHORS|CONSRTM)/gs;						# AUTHORS and CONSRTM can show up 
			$fullDoc .= "@lines";															# Therefore, check for both!
		}
		if ($att eq "AUTHORS"){ 
			@lines = $rawFile =~ /(AUTHORS.*?)TITLE/gs;
			$fullDoc .= "@lines";
		}
		if ($att eq "TITLE"){ 
			@lines = $rawFile =~ /(TITLE.*?)JOURNAL/gs;
			$fullDoc .= "@lines"; 
		}
		if ($att eq "JOURNAL"){ 															# JOURNAL also has corner cases
			@lines = $rawFile =~ /(JOURNAL.*?)(?=PUBMED|REFERENCE|REMARK|COMMENT)/gs;		# Therefore, check for all four possible endings		
			$fullDoc .= "@lines";
		}
		if ($att eq "MEDLINE"){ 
			@lines = $rawFile =~ /(PUBMED.*?)REFERENCE/gs;
			$fullDoc .= "@lines"; 
		}
		if ($att eq "FEATURES"){ 
			@lines = $rawFile =~ /(FEATURES.*?)ORIGIN/gs;
			$fullDoc .= "@lines"; 
		}
		if ($att eq "BASECOUNT"){															# BASECOUNT gets the same sequence of text ORIGIN does
			@lines = $rawFile =~ /(ORIGIN.*?)\/\//gs;										# calls function to count bases
			$baseRef = "@lines";															# countBases removes everything that is not "actg"
			$fullDoc .= countBases($baseRef) . "<br />";									# and returns the number of each base
		}
		if ($att eq "ORIGIN"){ 
			@lines = $rawFile =~ /(ORIGIN.*?)\/\//gs;
			$fullDoc .= "@lines"; 
		}
	}
	$fullDoc .= "</pre>";
	return $fullDoc;
}

sub countBases($){
	#removes all characters that are not "atcg" from ORIGIN string
	#and then counts the number of each "actg" base
	my $sequence = shift;
	my $char;
	my $numA = 0;
	my $numT = 0;
	my $numC = 0;
	my $numG = 0;

	# remove everything that is not "atcg" from the string 
	$sequence =~ s/\bORIGIN\b//;	# remove ORIGIN from the start of the string
	$sequence =~ s/[0-9]//g;		# remove all numerical characters from the string
	$sequence =~ s/\s+//g;			# remove all whitespace from string

	for(my $i = 0; $i < length($sequence); $i++){ 	#count the number of each base in the string
												
		$char = uc(substr($sequence, $i, 1));		#iterate through each character in the string
		if($char eq "A"){
			$numA = $numA + 1;
		}elsif($char eq "T"){
			$numT = $numT + 1;
		}elsif($char eq "C"){
			$numC = $numC + 1;
		}elsif($char eq "G"){
			$numG = $numG + 1;
		}
	}
	return "BASE COUNT      $numA A    $numC C     $numG G     $numT T";	
}

sub sendMail($$){
	#Prints the genbank file into a text file and emails it as an attachment
	my $address = shift;	# email address to send the file to
	my $fullText = shift;	# data to send
	my $textFile = "genBankResults.txt";			# text file the data will be written into a sent as an attachment

	$fullText =~ s/<.+?>//g;	#remove html tags

	#open file to full text into
	open(OUT, "> $textFile") ||	warn("Error opening file: $textFile... (must be a permission issue) $!\n");
	# print full text to file
	print OUT $fullText;
	# close output file
	close(OUT) || die "Could not close output file '$textFile' properly... $!";

	#send the file to the specified email
	system("cat $textFile| mail -s 'Your GenBank Results!' $address");
}