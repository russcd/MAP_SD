#ifndef __IMPORT_RATES_H
#define __IMPORT_RATES_H

void import_rates( string &input, map<string,double> &rates ) {
 
    string chrom ;
    double rate ;
    ifstream in ( input.c_str() ) ;
    while( !in.eof() ) {
        in >> chrom >> rate ;
        rates[chrom] = rate ; 
    }
}

#endif
