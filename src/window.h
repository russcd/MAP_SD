#ifndef __WINDOW_H
#define __WINDOW_H

void window( string id, vector<site> &germline, vector<site> &somatic, cmd_line &options ) {
 
    for ( int s = options.window/2 ; s < 22e6 ; s += options.window ) {
            
        double somatic_k = golden_search_site( options, somatic, s ) ;
        
        double germline_k = golden_search_site( options, germline, s ) ;
        
        cout << id << "\t" << s << "\t" << somatic_k << "\t" << germline_k << "\t" ;
        cout << compute_likelihood( options.error, germline, s, somatic_k ) << "\t" ;
        cout << compute_likelihood( options.error, germline, s, germline_k ) << "\t" ;
        cout << compute_likelihood( options.error, germline, s, somatic_k ) - compute_likelihood( options.error, germline, s, germline_k ) << endl ;
    }
}

#endif
