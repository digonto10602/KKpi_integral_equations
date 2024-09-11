#ifndef FASTGLFUNCTION_H
#define FASTGLFUNCTION_H

#include<bits/stdc++.h>
#include "./fastgl_algorithm/fastgl.hpp"
#include "./fastgl_algorithm/fastgl_cpp_convert.h"



using namespace std;

typedef complex<double> comp;

void line_maker_with_weights(   vector<comp> &qvec, 
                                vector<comp> &weights,
                                comp c,
                                comp d,
                                double points   )
{
    double a = -1.0;
    double b = 1.0;
    comp change_of_variable_multiplier = ((d-c)/(b-a));
    vector<comp> temp_weights;
    vector<comp> temp_qvec;

    for(int i=1;i<=points; ++i)
    {
        fastgl::QuadPair p = fastgl::GLPair(points,i);
        temp_weights.push_back(change_of_variable_multiplier * p.weight);
        temp_qvec.push_back(change_of_variable_multiplier * ( p.x() - a ) + c);
    }
    //double delz = abs(zfinal-zinitial)/points;
    for(int i=0;i<temp_qvec.size();++i)
    {
        //double z = zinitial + i*delz;
        //comp x_of_z = (b-a)*z + a;
        weights.push_back(temp_weights[temp_qvec.size() - i - 1]);
        qvec.push_back(temp_qvec[temp_qvec.size() - i - 1]);
        //qvec.push_back(x_of_z);

    }
}




#endif