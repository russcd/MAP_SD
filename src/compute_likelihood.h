#ifndef __COMPUTE_LIKELIHOOD_H
#define __COMPUTE_LIKELIHOOD_H

double compute_likelihood ( double &error, double &rate, vector<site> &sites, int site, double k ) {

    /// record log likelihood
    double lnl = 0 ;
    
    /// iterate through all sites
    for ( int s = 0 ; s < sites.size() ; s ++ ) {
        
        double rec = ( 1 - exp(-2*abs(sites[s].position-site) * rate ) ) / 2 ;
        
        /// deal with the A's first
        double l = 0 ;
        l += ( 1 - error ) * ( 1 - rec ) * k ;
        l += ( 1 - error ) * rec * ( 1 - k ) ;
        l += error * ( 1 - rec ) * ( 1 - k ) ;
        l += error * rec * k ;
        lnl += log( l ) * sites[s].A ;
        
        // now deal with the a's
        l = 0 ;
        l += ( 1 - error ) * ( 1 - rec ) * ( 1 - k ) ;
        l += ( 1 - error ) * rec * k ;
        l += error * ( 1 - rec ) * k ;
        l += error * rec * ( 1 - k ) ;
        lnl += log( l ) * sites[s].a ;
    }
    
    return ( lnl ) ;
}



#endif
