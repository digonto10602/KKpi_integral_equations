#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include<bits/stdc++.h>
#include<Eigen/Dense>
using namespace std;

typedef complex<double> comp;

void interpolator_ds_integraleq(    Eigen::VectorXcd &dsol,
                                    vector<comp> &sigvec,
                                    comp s,
                                    comp sigp,
                                    comp sigk,
                                    double a,
                                    double m,
                                    double eps,
                                    comp &result  )
{
    comp ii = {0.0,1.0};

    comp Gs = GS(s,sigp,sigk,m,eps);
    comp delsig;
    comp integralsum = {0.0,0.0};

    for(int i=0;i<sigvec.size();++i)
    {
        if(i==0)
        {
            delsig = sigvec[0];
        }
        else
        {
            delsig = sigvec[i] - sigvec[i-1];
        }

        integralsum = integralsum + delsig*kernel(s,sigp,sigvec[i],a,m,eps)*dsol(i);


    }

    comp one = {1.0,0.0};
    result = -one*Gs - one*integralsum;
}

#endif