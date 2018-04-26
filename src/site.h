#ifndef __REPRODUCE_H
#define __REPRODUCE_H

class site {
public:
    
    int position ;
    float cm ; 
    int A ;
    int a ;
    
    bool operator < (const site& s) const
    {
        return (position < s.position);
    }
    
} ;

#endif
