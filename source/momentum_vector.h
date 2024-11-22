#ifndef MOMENTUM_VECTOR_H
#define MOMENTUM_H

#include<bits/stdc++.h>
#include "fastGL_based_functions.h"
#include "functions_mom_based.h"

typedef complex<double> comp;

void flavor_based_momentum_vector(  std::vector<comp> &qvec,
                                    std::vector<comp> &weights, 
                                    comp En, 
                                    double mi, //mass of the spectator
                                    double number_of_points )
{
    comp kmin = 0.0; 
    comp kmax = pmom(En, 0.0, mi);
    line_maker_with_weights(qvec, weights, kmin, kmax, number_of_points); 
}


#endif