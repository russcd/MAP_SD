use strict ; 
use warnings ; 

my $pollen = $ARGV[0] ; 
my $leaf = $ARGV[1] ; 

my %pollen_lengths ; 

my $total = 0 ; 

print STDERR "IMPORTING POLLEN DISTRIBUTION\n" ;
open POL, "<$pollen" ;
while (<POL>) { 
	chomp $_ ; 
	my $read = <POL> ; 
	my $trash = <POL> ; 
	my $qual = <POL> ;
	chomp $read ; 
	$pollen_lengths{length($read)} ++ ;
	$total ++ ; 
}
close POL ; 

print STDERR "SUBSAMPLING LEAF READS\n" ;

open LEAF, "<$leaf" ;
while (<LEAF>) { 
	my $read = <LEAF> ; 
	my $trash = <LEAF> ; 
	my $qual = <LEAF> ; 
	chomp $read ; 
	chomp $qual ; 

	my $length = length( $read ) ; 	
	foreach my $l ( reverse( 20..$length ) ) { 
		if ( exists( $pollen_lengths{$l} ) && $pollen_lengths{$l} > 0 ) { 
			print $_, substr( $read, 0, $l ), "\n", $trash, substr( $qual, 0, $l ), "\n" ; 
			$pollen_lengths{$l} -- ; 
			if ( $pollen_lengths{$l} == 0 ) { 
				delete( $pollen_lengths{$l} ) ; 
			}
			last ; 
		}
	}

	if ( keys %pollen_lengths == 0 ) { 
		last ;
	}
}
