#ifndef __IMPORT_DATA_H
#define __IMPORT_DATA_H

void import_data( string &input, map<string,vector<site> > &germline, map<string,vector<site> > &somatic ) {
 
    string chrom ;
    ifstream in ( input.c_str() ) ;
    while( !in.eof() ) {
        site new_site ;
        site germ_site ;
        in >> chrom ;
        in >> new_site.position >> new_site.A >> new_site.a >> germ_site.A >> germ_site.a ;
        
        germ_site.position = new_site.position ;
        
        somatic[chrom].push_back( new_site ) ;
        germline[chrom].push_back( germ_site ) ; 
    }
}

#endif
