### simulate distortion along a chromosome conditional on read sequencing depths and on the distribution of recombination events

use strict ;
use warnings ; 
use Math::Random ; 

### uniform error rate in simulation
my $error = 1e-2 ;

system ("mkdir BOOT") ;

my %somatic_depth ;
my %germline_depth ;
my %cm ;

while (<STDIN>) { 
	chomp $_ ; 
	my @split = split ( /\t/, $_ ) ; 
	$somatic_depth{$split[1]} = $split[3] + $split[4] ;
	$germline_depth{$split[1]} = $split[5] + $split[6] ;
    $cm{$split[1]} = $split[2] ;
}

foreach my $effect ( 0, 0.01 ) {
    foreach my $rep ( 1 ) {

        my $site = (keys %cm)[rand keys %cm];
        open OUT, ">BOOT/${rep}_${effect}_${site}.txt" ;
	
        foreach my $position ( sort {$a<=>$b} keys %somatic_depth ) {
		
            my $drawA = Math::Random::random_binomial( 1, $somatic_depth{$position}, 0.5 ) ;
            print OUT "scaffold_X\t$position\t", $cm{$position}, "\t", $drawA, "\t", $somatic_depth{$position} - $drawA, "\t" ;
            
            ### recombination rate
            my $rec = 1 - exp( ( -2 * abs( $cm{$site} - $cm{$position} )/100 )/ 2 )  ;
            
            ### draw before recombination
            $drawA = Math::Random::random_binomial( 1, $germline_depth{$position}, 0.5+$effect ) ;
            
            ### swap with recombination
            my $swapA = Math::Random::random_binomial( 1, $drawA, $rec ) ;
            my $swapa = Math::Random::random_binomial( 1, $germline_depth{$position} - $drawA, $rec ) ;
            $drawA = $drawA - $swapA + $swapa ;
            
            ### swap for error
            $swapA = Math::Random::random_binomial( 1, $drawA, $error ) ;
            $swapa = Math::Random::random_binomial( 1, $germline_depth{$position} - $drawA, $error ) ;
            $drawA = $drawA - $swapA + $swapa ;
            
            print OUT $drawA, "\t", $germline_depth{$position} - $drawA, "\n" ;
            
        }
        
        close OUT ; 
    }
}
