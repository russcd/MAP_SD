#ifndef __REPRODUCE_H
#define __REPRODUCE_H

class site {
public:
    
    int position ;
    int A ;
    int a ;
    
    bool operator < (const site& s) const
    {
        return (position < s.position);
    }
    
} ;

#endif
