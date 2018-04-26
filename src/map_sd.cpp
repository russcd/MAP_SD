/// standard headers
#include <string.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <map>
#include <utility>
#include <random>
using namespace std ;

/// our headers
#include "cmd_line.h"
#include "site.h"
#include "import_data.h" 
#include "compute_likelihood.h" 
#include "golden_search.h"
#include "window.h" 

/// main 
int main ( int argc, char **argv ) {

    //// read command line options
    cmd_line options ;
    options.read_cmd_line( argc, argv ) ;
    cout.precision(10) ;

    /// data objects
    map<string,vector<site> > somatic ;
    map<string,vector<site> > germline ;
    
    /// read data object
    import_data( options.data_file, germline, somatic ) ;
    
    /// output windows across chromosomes
    if ( options.window != 0 && options.bootstrap == 0 ) {
        
        for ( auto i = germline.begin() ; i != germline.end() ; ++ i ) {
            
            window( i->first, germline[i->first], somatic[i->first],options ) ;
        }
    }
    
    /// otherwise, optimize position of max difference in lnl on chromosome
    else if ( options.bootstrap > 0 ) {
        
        /// iterate through chromosomes
        for ( auto i = germline.begin() ; i != germline.end() ; ++ i ) {
            
            /// bootstrap by position
            vector<double> boostrap_site ;
            if ( options.bootstrap > 0 ) {
                
                /// random number devices
                default_random_engine generator ;
                generator.seed(time(NULL)) ;
                uniform_int_distribution<int> dist(0,germline[i->first].size() ) ;
            
                for ( int b = 0 ; b < options.bootstrap ; b ++ ) {
                    
                    vector<site> bootstrap_germline ;
                    vector<site> bootstrap_somatic ;
                    
                    for ( int l = 0 ; l < germline[i->first].size() ; l ++ ) {
                        
                        ///// draw the index
                        int index = dist(generator) ;
                        bootstrap_germline.push_back( germline[i->first][index] ) ;
                        bootstrap_somatic.push_back( somatic[i->first][index] ) ;
                    }
                    
                    sort( bootstrap_somatic.begin(), bootstrap_somatic.end() ) ;
                    sort( bootstrap_germline.begin(), bootstrap_germline.end() ) ; 
                    
                    /// find the site with the maximum differnce in lnl
                    boostrap_site.push_back( golden_search_position( options, bootstrap_somatic, bootstrap_germline ) ) ;
                    
                    /// compute k at each site
                    double somatic_k = golden_search_site( options, somatic[i->first], boostrap_site.back() ) ;
                    double germline_k = golden_search_site( options, germline[i->first],  boostrap_site.back() ) ;
                    
                    cout << i->first << "\t" << b << "\t" << boostrap_site.back() << "\t" << somatic_k << "\t" << germline_k << "\t" << compute_likelihood( options.error, bootstrap_germline, boostrap_site.back(), somatic_k ) << "\t" << compute_likelihood( options.error, bootstrap_germline, boostrap_site.back(), germline_k ) << endl ;
                }
            }
        }
    }
    
    //// otherwise just do a single fit
    else {
        
        for ( auto i = germline.begin() ; i != germline.end() ; ++ i ) {

            /// find the site with the maximum differnce in lnl
            double optimum_site = golden_search_position( options, somatic[i->first], germline[i->first] ) ;
            
            /// compute k at each site
            double somatic_k = golden_search_site( options, somatic[i->first], optimum_site ) ;
            double germline_k = golden_search_site( options, germline[i->first], optimum_site ) ;
            
            cout << i->first << "\t" << optimum_site << "\t" << somatic_k << "\t" << germline_k << "\t" << compute_likelihood( options.error, germline[i->first], optimum_site, somatic_k ) << "\t" << compute_likelihood( options.error, germline[i->first], optimum_site, germline_k ) << endl ;
        }
    }
    
    //// sucess
    return(0) ;
}
