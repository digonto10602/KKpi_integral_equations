#ifndef SA_METHOD_FUNCTIONS_H
#define SA_METHOD_FUNCTIONS_H

#include<bits/stdc++.h>
#include<Eigen/Dense>
#include "functions_momrep_based.h"
using namespace std;

typedef complex<double> comp;


comp delta_epsilon_func(    comp sigk,
                            double a,
                            double m,
                            double epsilon  )
{
    double pi = acos(-1.0);
    comp sb = sigmab(a,m);
    //comp sigk = sigma_p(s,k,m);
    double delsfnc = real(sigk - sb);
    if(delsfnc==0.0)
    {
        return 1.0/(pi*epsilon);
    }
    else 
    {
        return epsilon/(pi*((sigk - sb)*(sigk - sb) + epsilon*epsilon));
    }
}

comp delta_M2kfunc(     comp sigk,
                        double a,
                        double m,
                        double epsilon  )
{
    comp M2k = M2kfunc(a,sigk,m,epsilon);
    double pi = acos(-1.0);
    comp ii = {0.0,1.0};

    comp gsq = gfuncConst(a,0.0,m)*gfuncConst(a,0.0,m);

    return M2k - gsq*ii*pi*delta_epsilon_func(sigk,a,m,epsilon);
}

comp kernel_with_deltaM2k_pk_2eps(    comp s,
                        comp p,
                        comp k,
                        double a,
                        double m,
                        double epsilon,
                        double eps_for_m2k )
{
    double pi = acos(-1.0);
    comp picomp = (comp) pi;
    //cout<<"Gs:"<<GS(s,sigp,sigk,m,epsilon)<<endl;
    //cout<<"tau:"<<tau(s,sigk,m)<<endl;
    //cout<<"M2:"<<M2kfunc(a,sigk,m,epsilon)<<endl;
    comp sigk = sigma_p(s,k,m);
    
    return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*GS_pk(s,p,k,m,epsilon)*delta_M2kfunc(sigk,a,m,eps_for_m2k);
    
}

void deltaBmat_maker_momrep_2eps_with_weights(       Eigen::MatrixXcd &Bmat,
                                                comp s,
                                                vector<comp> &p,
                                                vector<comp> &k,
                                                vector<comp> &weights,
                                                double a,
                                                double m,
                                                double epsilon,
                                                double eps_for_m2k  )
{
    comp ii = {0.0,1.0};
    comp pi = acos(-1.0);
    comp delk = {0.0,0.0};

    for(int i=0;i<p.size();++i)
    {
        for(int j=0;j<k.size();++j)
        {
            if(j==0)
            {
                delk = k[1] - k[0];
            }
            else
            {
                delk = k[j] - k[j-1];
            }

            comp weight = weights[j];
            //make sure you supply both the invariant vectors
            //as the same, otherwise this condition will fail
            if(j==i)
            {
                comp one = {1.0,0.0};
                Bmat(i,j) = one + weight*kernel_with_deltaM2k_pk_2eps(s,p[i],k[j],a,m,epsilon,eps_for_m2k);
                //Bmat(i,j) =  weight*kernel_pk_2eps(s,p[i],k[j],a,m,epsilon,eps_for_m2k);
            
            }
            else 
            {
                Bmat(i,j) = weight*kernel_with_deltaM2k_pk_2eps(s,p[i],k[j],a,m,epsilon,eps_for_m2k);
            }
        }
    }

}


void interpolator_Kphib_integraleq_momrep_2eps_with_weights(    Eigen::VectorXcd &Kphib,
                                                                vector<comp> &qvec,
                                                                vector<comp> &weights,
                                                                comp s,
                                                                comp p,
                                                                comp k,
                                                                double a,
                                                                double m,
                                                                double eps,
                                                                double eps_for_m2k,
                                                                comp &result  )
{
    comp ii = {0.0,1.0};

    comp Gs = GS_pk(s,p,k,m,eps);
    comp gsq = gfuncConst(a,0.0,m)*gfuncConst(a,0.0,m);
    comp delq;
    comp integralsum = {0.0,0.0};

    for(int i=0;i<qvec.size();++i)
    {
        if(i==0)
        {
            delq = qvec[1] - qvec[0];
        }
        else
        {
            delq = qvec[i] - qvec[i-1];
        }

        integralsum = integralsum + weights[i]*kernel_with_deltaM2k_pk_2eps(s,p,qvec[i],a,m,eps,eps_for_m2k)*Kphib(i);


    }

    comp one = {1.0,0.0};
    result = -one*gsq*Gs - one*integralsum;
}





















#endif