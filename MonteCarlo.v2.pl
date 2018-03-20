#BIF 712 Assignment 1
#Robert Drutis 

#!/usr/bin/perl

use strict;
use warnings;

# CGI example using perl CGI module
use CGI ":standard";

my $cgi = new CGI;      # read in parameters

print $cgi->header( );  # print the HTML header
print $cgi->end_html( );

#initialize variables
my $range = 4; #range 1-2 to generate molecules H, O, N
my $mol; #molecule generated
my $compForm; #compound formula
my $compBond = 0; #compound open bond value
my $compMass = 0; #compound molecular weight
my $carbs = 0; #when number of carbons reaches 3, make the bond a double bond
my @molCount = (0,0,0,0); #array to count the number of (H,O,N,C) in compound

my $amino;
my $matchFound = 0;
my $trials = 0; #number of random trials performed
my $testDiff;
my $minDiff; 
my $closeAmino;
my $closeMass = 0;
my $closeFormula;
my $j = 0;

#for verbose commandline argument
#this is for the purpose of running multiple tests 
#to calculate average outcomes without the screen
my $verbose = 0; #verbose is FALSE by default
my $arg;
foreach $arg(@ARGV){
	if($arg eq "-v"){
		#print "the user entered -v !\n";
		$verbose = 1;
	}
}

my %weights = (
	Ala => "89.0929",  
	Asp => "133.1024",
	Glu => "147.1289",
	Phe => "165.1887",
	Gly => "75.0664",
	His => "155.1542",
	Ile => "131.1724",
	Lys => "146.1870",
	Leu => "131.1724",
	Asn => "132.1176",
	Pro => "115.1301",
	Gln => "146.1441",
	Arg => "174.2004",
	Ser => "105.0923",
	Thr => "119.1188",
	Val => "117.1459",
	Trp => "204.2247",
	Tyr => "181.1881",
);

my %highscores = (
	Ala => "Nothing",  
	Asp => "Nothing",
	Glu => "Nothing",
	Phe => "Nothing",
	Gly => "Nothing",
	His => "Nothing",
	Ile => "Nothing",
	Lys => "Nothing",
	Leu => "Nothing",
	Asn => "Nothing",
	Pro => "Nothing",
	Gln => "Nothing",
	Arg => "Nothing",
	Ser => "Nothing",
	Thr => "Nothing",
	Val => "Nothing",
	Trp => "Nothing",
	Tyr => "Nothing",
);

my %highscoreInfo = (
	Ala => "Nothing",  
	Asp => "Nothing",
	Glu => "Nothing",
	Phe => "Nothing",
	Gly => "Nothing",
	His => "Nothing",
	Ile => "Nothing",
	Lys => "Nothing",
	Leu => "Nothing",
	Asn => "Nothing",
	Pro => "Nothing",
	Gln => "Nothing",
	Arg => "Nothing",
	Ser => "Nothing",
	Thr => "Nothing",
	Val => "Nothing",
	Trp => "Nothing",
	Tyr => "Nothing",
);



#calculate compound OBV
sub calcBond($$$){
	my $addMol = shift; #molecule being added
	my $valence = shift; #current compound OBV
	my $doubleBond = shift; #TRUE if double bond is being added
	

	#if molecule is hydrogen, add +1 to the OBV
	if($addMol == 1){ 
		$valence = $valence + 1;

	#if molecule is oxygen, nitrogen, or carbon
	}elsif($addMol > 1 && $doubleBond == 0){ 

		#determine how to add bond, check if positive or negative
		if($valence >= 0){ #if the OBV is positive
			$valence = $valence + (-$addMol); #subtract the valence from the OBV

		}elsif($valence < 0){ #if the OBV is negative
			$valence = $valence + (-$addMol + 2); #add 2 to the calculation to represent bond formation
		}

	#if the molecule being added is a third carbon
	}elsif($addMol > 1 && $doubleBond == 1){
		#subtract 2 if OBV is positive, leave unchanged if it is negative
		if($valence >= 0){ #if the OBV is positive
			$valence = $valence + (-2); #subtract the valence from the OBV
			#print "\nDEBUG: DOUBLE BOND ADDED\n";
		}else{
			#print "\nDEBUG: DOUBLE BOND ADDED\n";
		}
	}

	return($valence); #return the new bond value for the compound
}

#calculate compound molecular weight
sub calcWeight($$){
	my $addMol = shift; #molecule being added
	my $mass = shift; #current compound OBV

	#add the weight of the new molecule to the compound
	if($addMol == 1){ #if 1, add weight for hydrogen
		$mass = $mass + 1.0079;
	}elsif($addMol == 2){ #if 2, add weight for oxygen
		$mass = $mass + 15.9994;
	}elsif($addMol == 3){ #if 3, add weight for nitrogen
		$mass = $mass + 14.0067;
	}elsif($addMol == 4){ #if 3, add weight for carbon
		$mass = $mass + 12.0107;
	}

	return($mass); #return the new molecular weight of the compound
}

#calculate compound formula
sub calcFormula($$$$){
	my $hydrogens = shift; #number of hydrogens
	my $oxygens = shift; #number of oxygens
	my $nitrogens = shift; #number of nitrogens
	my $carbons = shift; #number of carbons

	my $chemForm = ""; #chemical formula of the compound

	#Add the molecule symbol and its amount to the formula
	#if only one of the molecule is present, add only its symbol

	#add number of carbons to formula
	if($carbons > 0){ 
		if($carbons == 1){ 
			$chemForm = $chemForm . "C";
		}else{
			$chemForm = $chemForm . "C$carbons";
		}
	}
	#add number of hydrogens to formula
	if($hydrogens > 0){ 
		if($hydrogens == 1){ 
			$chemForm = $chemForm . "H";
		}else{ 
			$chemForm = $chemForm . "H$hydrogens";
		}
	}
	#add number of nitrogens to formula
	if($nitrogens > 0){ 
		if($nitrogens == 1){ 
			$chemForm = $chemForm . "N";
		}else{
			$chemForm = $chemForm . "N$nitrogens";
		}
	}
	#add number of oxygens to formula
	if($oxygens > 0){ 
		if($oxygens == 1){ 
			$chemForm = $chemForm . "O";
		}else{
			$chemForm = $chemForm . "O$oxygens";
		}
	}

	return($chemForm);
}

#DEBUG: display molecule info
sub molInfo($){
	my $molecule = shift;
	if($molecule == 1){
		print "hydrogen\n";
	}elsif($molecule == 2){
		print "oxygen\n";
	}elsif($molecule == 3){
		print "nitrogen\n";
	}elsif($molecule == 4){
		print "carbon\n";
	}
}

#########################
#     Main Program      #
#########################

#Continue generating compounds until an amino acid is found
while($matchFound != 1){
	do{ 
		#generate a random molecule (H,O,N,C)
		$mol = 1+ int(rand($range));
		#print "$mol\n"; #DEBUG: display number generated
		#molInfo($mol); #DEBUG: display molecule

		$molCount[$mol-1]++; #iterate the count by 1
		#print "DEBUG: $molCount[$mol-1]\n";
		if($mol == 4){ #if molecule is a carbon, iterate the carbon counter for double bonds
			$carbs++;
			#print "DEBUG: Number of carbons: $carbs\n";
		}

		#add the new molecule to the compound
		#the default case is that the molecule being added is not a 3rd carbon
		#otherwise, if it is a third carbon, double bond it
		if($carbs != 3){ 
			$compBond = calcBond($mol, $compBond, 0);
			#print "OBV of compound is: $compBond\n"; #DEBUG
		}else{ 
			$compBond = calcBond($mol, $compBond, 1); #double bond the carbon
			#print "OBV of compound is: $compBond\n"; #DEBUG
			$carbs = 0; #reset carbon counter
		}	

		#calculate the weight of the compound
		$compMass = calcWeight($mol, $compMass);
		#print "Molecular weight of compound is: $compMass\n"; #DEBUG

		#calculate the formula of the compound
		$compForm = calcFormula($molCount[0], $molCount[1], $molCount[2], $molCount[3]);
		#print "Chemical formula of compound is: $compForm\n"; #DEBUG

		#if valence drops below -4, add 4 hydrogens and print out new compound
		if($compBond < -4){

			#add 3 hydrogens to the molecule
			for(my $i=0; $i<3; $i++){
				$compBond = calcBond(1, $compBond, 0);
				$compMass = calcWeight(1, $compMass);
				$molCount[0] = $molCount[0] + 1; #add 3 hydrogens to the molecule counter
			}
			#print "\nDEBUG: ADDED 3 HYDROGENS\n";
			#print "OBV of compound is: $compBond\n"; #DEBUG
			#print "Molecular weight of compound is: $compMass\n"; #DEBUG
			$compForm = calcFormula($molCount[0], $molCount[1], $molCount[2], $molCount[3]);
			#print "Chemical formula of compound is: $compForm\n"; #DEBUG
		}


		#check it the mass matches any of the known amino acid weights
		for $amino (keys %weights) { 
		
			#if a match is found, print out the information
			if($weights{$amino} eq ("" . $compMass)){ 
				print "\ncompound: $compForm Mol wt.: $compMass amino acid: $amino\n";;
				$matchFound = 1;
			}

		}
		#if the bond is stable, display message
		#will not display if "-v" was entered in the commandline
		if($compBond == 0 && $verbose != 1){
			print "\nStable molecule has been generated!\n";
		}

		#sleep 1;

	#continue generating compounds until a match is found, valence is stable, or weight exceeds Trp
	}while($compBond != 0 && $compMass < 204.2247 && $matchFound != 1); 

		$trials++;

		for $amino (keys %weights){
			if($matchFound != 1){
				if($highscores{$amino} eq "Nothing"){
					$highscores{$amino} = $compMass;
					$highscoreInfo{$amino} = "compound: $compForm Trial: $trials"
				}elsif(abs($weights{$amino} - $highscores{$amino}) > abs($weights{$amino} - $compMass)){
					$highscores{$amino} = $compMass;
					$highscoreInfo{$amino} = "compound: $compForm Trial: $trials"
				}
			}
		}

		#store the closest match so far out of all of the amino acids
		if($matchFound != 1){
			for $amino (keys %highscores){
				$testDiff = abs($highscores{$amino} - $weights{$amino});

				if($j == 0){
					$minDiff = $testDiff;
					$closeAmino = $amino;
					$closeMass = $highscores{$amino};
					$j = 1;
				}elsif($minDiff > $testDiff){
					$minDiff = $testDiff;
					$closeAmino = $amino;
					$closeMass = $highscores{$amino};
				}
			}
		}
		#if the compound doesn't match an amino acid, print out its mass and formula
		#will not display if "-v" was entered in the commandline
		if($matchFound !=1 && $verbose != 1){
			print "compound: $compForm Mol wt.: $compMass amino acid: N/A\n";
			print "so far the closest match was: $closeAmino Mol wt.: $closeMass\n";
			print "Trials: $trials\n";
		}

		#Display Highscores every 10 trials
		#will not display if "-v" was entered in the commandline
		if($matchFound !=1 && $trials == 10 && $verbose != 1){
			print "CURRENT HIGHSCORES:\n";
			for $amino (keys %highscores) {
				print "$amino ------> $highscores{$amino} $highscoreInfo{$amino}\n";
			}
		}

		#reset values for generating compound
		$compBond = 0;
		$compMass = 0;
		$compForm = "";
		$carbs = 0;
		@molCount = (0,0,0,0);

}

print "Trials: $trials\n";















#############################################################

#Display Final Highscores
#will not display if "-v" was entered in the commandline
if($verbose != 1){
	print "FINAL HIGHSCORES:\n";
	for $amino (keys %highscores) {
		print "$amino ------> $highscores{$amino} $highscoreInfo{$amino}\n";
	}
}





print "Program Complete!";
