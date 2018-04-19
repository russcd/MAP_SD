#ifndef __GOLDEN_SEARCH_H
#define __GOLDEN_SEARCH_H

/// golden section search for single parameter optimization
/// not currently included, but may be useful to legacy versions 
double golden_search_site ( cmd_line &options, vector<site> &s, double &rate, int site ) {
    
    /// now do golden search until we reach tolerance threshhold and stop
    double phi = ( sqrt(5) - 1 ) / 2 ;
    
    /// parameter values to hold during search
    double low_bracket = 0.01 ;
    double high_bracket = 0.99 ;
    double param_low = high_bracket + phi * ( low_bracket - high_bracket ) ;
    double param_high = low_bracket + phi * ( high_bracket - low_bracket ) ;
    
    /// likelihood information for test points
    double lnl_low = 0 ;
    double lnl_high = 0 ;
    double lnl_diff = 1000 ;
    int iteration = 0 ;
    
    /// optimize somatic first
    while ( options.tolerance < lnl_diff ) {
        
        /// compute likelihood of low point
        if ( lnl_low == 0 ) {
            lnl_low = compute_likelihood( options.error, rate, s, site, param_low ) ;
        }
                
        /// compute likelihood of high point
        if ( lnl_high == 0 ) {
            lnl_high = compute_likelihood( options.error, rate, s, site, param_high ) ;
        }
                
        /// record dfiference
        lnl_diff = abs( lnl_low - lnl_high ) ;
        
        /// if true, we know that the maximum is between param_low and high_bracket
        if ( lnl_high >= lnl_low ) {
            low_bracket = param_low ;
            param_low = param_high ;
            param_high = low_bracket + ( high_bracket - low_bracket ) * phi ;
            lnl_low = lnl_high ;
            lnl_high = 0 ;
        }
        
        /// otherwise, the maximum is between low_bracket and param_high
        else {
            high_bracket = param_high ;
            param_high = param_low ;
            param_low = high_bracket + ( low_bracket - high_bracket ) * phi ;
            lnl_high = lnl_low ;
            lnl_low = 0 ;
        }
        
        /// update by it
        iteration ++ ;
    }

    return ( ( param_low + param_high ) / 2 ) ;
}

/// golden section search for position
double golden_search_position ( cmd_line &options, vector<site> &somatic, vector<site> &germline, double &rate ) {
    
    /// now do golden search until we reach tolerance threshhold and stop
    double phi = ( sqrt(5) - 1 ) / 2 ;
    
    /// parameter values to hold during search
    double low_bracket = somatic[0].position ;
    double high_bracket = somatic.back().position ;
    double param_low = high_bracket + phi * ( low_bracket - high_bracket ) ;
    double param_high = low_bracket + phi * ( high_bracket - low_bracket ) ;
    
    /// likelihood information for test points
    double lnl_low = 0 ;
    double lnl_high = 0 ;
    double lnl_diff = 1000 ;
    int iteration = 0 ;
    
    /// optimize somatic first
    while ( options.tolerance < lnl_diff ) {
        
        /// compute likelihood of low point
        if ( lnl_low == 0 ) {
            
            double somatic_k = golden_search_site( options, somatic, rate, param_low ) ;
            double germline_k = golden_search_site( options, germline, rate, param_low ) ;
            
            lnl_low = compute_likelihood( options.error, rate, germline, param_low, germline_k ) -compute_likelihood( options.error, rate, germline, param_low, somatic_k ) ;
        }
        
        /// compute likelihood of high point
        if ( lnl_high == 0 ) {
            
            double somatic_k = golden_search_site( options, somatic, rate, param_high ) ;
            double germline_k = golden_search_site( options, germline, rate, param_high ) ;
            
            lnl_high = compute_likelihood( options.error, rate, germline, param_high, germline_k ) - compute_likelihood( options.error, rate, germline, param_high, somatic_k ) ;
        }
        
        /// record dfiference
        lnl_diff = abs( lnl_low - lnl_high ) ;
        
        //cout << iteration << "\t" << low_bracket << "\t" << param_low << "\t" << param_high << "\t" << high_bracket << "\t" << lnl_low << "\t" << lnl_high << endl ;
    
        /// if true, we know that the maximum is between param_low and high_bracket
        if ( lnl_high >= lnl_low ) {
            low_bracket = param_low ;
            param_low = param_high ;
            param_high = low_bracket + ( high_bracket - low_bracket ) * phi ;
            lnl_low = lnl_high ;
            lnl_high = 0 ;
        }
        
        /// otherwise, the maximum is between low_bracket and param_high
        else {
            high_bracket = param_high ;
            param_high = param_low ;
            param_low = high_bracket + ( low_bracket - high_bracket ) * phi ;
            lnl_high = lnl_low ;
            lnl_low = 0 ;
        }
        
        /// update by it
        iteration ++ ;
    }
    
    return ( ( param_low + param_high ) / 2 ) ;
}


#endif
