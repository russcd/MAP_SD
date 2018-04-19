#ifndef __CMD_LINE_H
#define __CMD_LINE_H

/// command line information and global parameters
class cmd_line {
public:
    
    /// read file
    string data_file ;
    
    /// rates files
    string rate_file ;
    
    /// minimum size
    int min_size ;
    
    /// error rate per site
    double error ;
    
    /// tolerance during convergence
    double tolerance ;
    
    /// windows
    int window ;
    
    /// number of bootstraps
    int bootstrap ;
    
    /// read options
    void read_cmd_line ( int argc, char *argv[] ) ;
    
} ;

void cmd_line::read_cmd_line ( int argc, char *argv[] ) {
    
    /// set parameter defaults
    data_file = "sites_ind1.txt" ;
    rate_file = "rates.txt" ;
    min_size = 1000 ;
    error = 1e-2 ;
    tolerance = 1e-4 ;
    window = 0 ;
    bootstrap = 0 ;
    
    /// accept command line parameters
    for (int i=1; i<argc; i++) {

        if ( strcmp(argv[i],"-d") == 0 ) {
            data_file = argv[++i] ;
        }
        if ( strcmp(argv[i],"-r") == 0 ) {
            rate_file = argv[++i] ;
        }
        if ( strcmp(argv[i],"-m") == 0 ) {
            min_size = atoi(argv[++i]) ;
        }
        if ( strcmp(argv[i],"-w") == 0 ) {
            window = atoi(argv[++i]) ;
        }
        if ( strcmp(argv[i],"-b") == 0 ) {
            bootstrap = atoi(argv[++i]) ;
        }
        if ( strcmp(argv[i],"-e") == 0 ) {
            error = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"-t") == 0 ) {
            tolerance = atof(argv[++i]) ;
        }
    }
    
    return ;
}


#endif
