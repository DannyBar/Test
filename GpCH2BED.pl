#takes an indexed fq file and a file with the genome seq of the same coordinates. 
#Returns locations of all methylated GpCHs (i.e. perl GpCH2BED.pl Libname > LibName.GpCH.bed)
use strict;
use warnings;

my $Input = $ARGV[0];
				
my $filename = "$Input.sam";
open(my $fh, '<:encoding(UTF-8)', $filename)
  or die "Could not open file '$filename' $!";

my $HG = "$Input.MappedReads.seq";
open(my $hg, '<:encoding(UTF-8)', $HG)
  or die "Could not open file '$HG' $!";

my $clipping = 12; # can change.

while (my $row = <$fh>) {
#Ignore head
  while ($row =~/^@/) { $row = <$fh>  }

	my $Coordinates; my $HGseq; my $HGseqCT; my $HGseqGA;

	my $HGrow = <$hg>; 		
	if ($HGrow =~ /^>/) {
		$Coordinates = $HGrow;
		$HGrow = <$hg>;
		$HGseq = $HGrow;
		($HGseqCT = $HGseq) =~ s/C/T/g; 
		($HGseqGA = $HGseq) =~ s/G/A/g; 
	}
	else {die "Error, no > at $HGrow"}

# Split line
  	my @array = split(/\t/, $row);
	my $rev = $array[1]; #Strand
        my $chr = $array[2]; #$chr =~ s/^\s+|\s+$//g; #Chromosome
        my $start = $array[3]; #$start =~ s/^\s+|\s+$//g; #Start position
	my $Q = $array[4]; #Quality
	my $CIGAR = $array[5]; #Map Quality
	my $MD = $array[12]; #String for mismatching positions
	my $seq = $array[9]; #Sequence

	until ( (index($Coordinates,$chr) != -1) and (index($Coordinates,$start) != -1)) { #Find mapped read in the sequencing file
		$row = <$fh> ;
		@array = split(/\t/, $row);
		$rev = $array[1]; #Strand
        	$chr = $array[2]; #$chr =~ s/^\s+|\s+$//g; #Chromosome
        	$start = $array[3]; #$start =~ s/^\s+|\s+$//g; #Start position
		$Q = $array[4]; #Quality
		$CIGAR = $array[5]; #Map Quality
		$MD = $array[12]; #String for mismatching positions
		$seq = $array[9]; #Sequence
	}
	#print ($Coordinates,"\t",$chr,"\t",$start,"\n");
	
        if (($Q >= 60)) {  #and ($CIGAR =~ m/151M/)
		my $fwd = uc $seq; (my $seqCT = $fwd) =~ s/C/T/g;
		my $bck = uc $seq; (my $seqGA  = $bck) =~ s/G/A/g;
		if (substr($fwd,$clipping,length($seq)-(2*$clipping)) !~ m/[ATC]C[ATC]/) { #requireing complete conversion of the strand area being used.
	  	while ($fwd =~ /GC[ATC]/g) { #find GCH 
		   	    my $x = $-[0];  
			    if (($x > $clipping) and ($x < length($seq) - $clipping) #not at start or end
				and (index($HGseqCT,substr($seqCT,$x-7,15)) != -1) #15bp match in genome
				and (substr($HGseq,index($HGseqCT,substr($seqCT,$x-7,15))+7,2) eq "GC") and (substr($seqCT,$x-7,15) !~ m/T{8,}|G{8,}|C{8,}|A{8,}/)) { #Genome has GC at that point
	        	                my $line = join("\t",$chr,$start + index($HGseqCT,substr($seqCT,$x-7,15)) + 8, $start + index($HGseqCT,substr($seqCT,$x-7,15)) + 9); 
			 	       	print "$line\n";#}	#print($x,"\t",substr($fwd,$x-4,11),"\t",substr($seqCT,$x-4,11),"\n"); 
					#my $z = substr($seqCT,$x-7,15); die "$z"
			    }	
	        }}

        	if ((substr($bck,$clipping,length($seq)-(2*$clipping)) !~ m/[ATG]G[ATG]/)) {
  		while ($bck =~ /[ATG]GC/g) {
		   	    my $x = $-[0];  
			    if (($x > $clipping) and ($x < length($seq) - $clipping) 
			        and (index($HGseqGA,substr($seqGA,$x-7,15)) != -1) 
                                and (substr($HGseq,index($HGseqGA,substr($seqGA,$x-7,15))+8,2) eq "GC") and (substr($seqGA,$x-7,15) !~ m/T{8,}|G{8,}|C{8,}|A{8,}/)) {
        	                	my $line = join("\t",$chr,$start + index($HGseqGA,substr($seqGA,$x-7,15)) + 8, $start + index($HGseqGA,substr($seqGA,$x-7,15)) + 9); 
					print "$line\n";#}	#print($x,"\t",substr($bck,$x-4,11),"\t",substr($seqGA,$x-4,11),"\n");
			    }
	        }}
	}

  
}
