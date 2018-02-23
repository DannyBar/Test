use strict;
use warnings;

my $Input = $ARGV[0]; #input must be a sorted BED file. 
my $window = $ARGV[1]; #sliding window for averaging
my $treshold = $ARGV[2]; #minimal number of modifications
my $filename = $Input;

my $old_chr = "none";
my @stack = ();
my $oldLine = "";

open(my $fh, '<:encoding(UTF-8)', $filename)
  or die "Could not open file '$filename' $!";

# add option to account for missed modifications later

while (my $row = <$fh>) { 
  if (($row ne $oldLine) or (1 == 1)) { #Unique reads (not used for now).
	$oldLine = $row;
  	my @array = split(/\s+/, $row);
	my $chr = $array[0];
        my $pos = $array[1]; 
        my $posE = $array[2]; 
	if ($posE-$pos == 1) { #Make sure modification feature is 1 bp.
		if ($chr ne $old_chr) {  #If new chromosome, restart.
			@stack = ();
			$old_chr = $chr;
		}
		if ($stack[0]+$window<=$pos) {	#old positions will be excluded with the next position added, mark values for center now.
			if (scalar(@stack)>=$treshold) { #Write value for mid position of stack 
				my $coordinate = int(($stack[0] + $stack[scalar(@stack)-1])/2); 
				print($chr,"\t",$coordinate,"\t",$coordinate+1,"\t",scalar(@stack),"\n"); # Print moving average value
				for (my $i=1; $i < scalar(@stack); $i++) { print($chr,"\t",$coordinate,"\t",$coordinate+1,"\t","0","\n");} #For R histograms
			}
		}
		
		push(@stack, $pos);
	
		until ($stack[0]+$window>$pos) {
			shift(@stack);
		}
	}
	else { print "Error at $row";}
  }
}


