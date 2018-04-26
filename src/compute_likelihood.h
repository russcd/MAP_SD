#ifndef __COMPUTE_LIKELIHOOD_H
#define __COMPUTE_LIKELIHOOD_H

double compute_likelihood ( double &error, vector<site> &sites, int site, double k ) {

    /// record log likelihood
    double lnl = 0 ;
    
    /// find position
    double cm_at_site = 0 ;
    int min = 10000000 ;
    for ( int i = 0 ; i < sites.size() ; i ++ ) {
        if ( abs( sites[i].position - site ) < min ) {
            min = abs( sites[i].position - site ) ;
            cm_at_site = sites[i].cm ;
        }
    }
    
    /// iterate through all sites
    for ( int s = 0 ; s < sites.size() ; s ++ ) {
        
        double rec = ( 1 - exp(-2*(abs(sites[s].cm-cm_at_site)/100) ) ) / 2 ;
        
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
