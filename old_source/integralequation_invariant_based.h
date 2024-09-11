#ifndef INTEGRALEQUATION_INVARIANT_H
#define INTEGRALEQUATION_INVARIANT_H

#include<bits/stdc++.h>
#include<Eigen/Dense>
using namespace std;

typedef complex<double> comp;

void Gvec_maker(    Eigen::VectorXcd &Gvec,
                    comp s,
                    vector<comp> &sigmap,
                    comp sigmak,
                    double m,
                    double epsilon  )
{
    for(int i=0;i<sigmap.size();++i)
    {
        Gvec(i) = GS(s,sigmap[i],sigmak,m,epsilon);
    }
}

void Bmat_maker(    Eigen::MatrixXcd &Bmat,
                    comp s,
                    vector<comp> &sigp,
                    vector<comp> &sigk,
                    double a,
                    double m,
                    double epsilon  )
{
    comp ii = {0.0,1.0};
    comp pi = acos(-1.0);
    comp delsig = {0.0,0.0};

    for(int i=0;i<sigp.size();++i)
    {
        for(int j=0;j<sigk.size();++j)
        {
            if(j==0)
            {
                delsig = sigk[j];
            }
            else
            {
                delsig = sigk[j] - sigk[j-1];
            }

            //make sure you supply both the invariant vectors
            //as the same, otherwise this condition will fail
            if(j==i)
            {
                comp one = {1.0,0.0};
                Bmat(i,j) = one + delsig*kernel(s,sigp[i],sigk[j],a,m,epsilon);
            }
            else 
            {
                Bmat(i,j) = delsig*kernel(s,sigp[i],sigk[j],a,m,epsilon);
            }
        }
    }

}


#endif