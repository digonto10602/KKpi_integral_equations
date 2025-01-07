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
    comp epsilon_for_kvec = 1.0e-10; 
    comp kmin = 0.0 + epsilon_for_kvec; 
    comp kmax = pmom(En, 0.0, mi) - epsilon_for_kvec;
    line_maker_with_weights(qvec, weights, kmin, kmax, number_of_points); 
}


#endif