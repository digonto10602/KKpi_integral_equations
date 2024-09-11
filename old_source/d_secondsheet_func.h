#ifndef D_SECONDSHEET_FUNC_H
#define D_SECONDSHEET_FUNC_H

#include<bits/stdc++.h>
#include<Eigen/Dense>
#include "functions_momrep_based.h"
using namespace std;

typedef complex<double> comp;

void cmat_maker_with_weights(   Eigen::MatrixXcd &Cmat,
                                Eigen::MatrixXcd Dmatsol,
                                vector<comp> qvec,
                                vector<comp> weights,
                                comp s, 
                                double m,
                                double a,
                                double epsilon 
                                )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);

    for(int i=0;i<qvec.size();++i)
    {
        for(int j=0;j<qvec.size();++j)
        {
            if(i==j)
            {
                comp mom_p = qvec[i];
                comp mom_k = qvec[j];
                comp weight_p = weights[i];
                comp weight_k = weights[j];
                comp sigp = sigma_p(s,mom_p,m);
                comp sigk = sigma_p(s,mom_k,m);
                double theta = 0.0;
                comp m2k_sigp = M2kfunc(a,sigp,m,epsilon);
                comp m2k_sigk = M2kfunc(a,sigk,m,epsilon);
                comp firstterm = ii/2.0 - rho_func(sigp, m, theta)*m2k_sigk;
                comp secondterm = -mom_p*mom_p*weight_p/(pow(2.0*pi,2.0)*omega_comp(mom_p,m))*(rho_func(sigp,m,theta)*Dmatsol(i,j) - imag(GS_pk(s,mom_p,mom_k,m,epsilon))*m2k_sigk);
                comp temp_third = {0.0,0.0};
                for(int k=0;k<qvec.size();++k)
                {
                    comp mom_kprime = qvec[k];
                    comp weight_kprime = weights[k];
                    temp_third       =   temp_third + mom_p*mom_p*weight_p/(pow(2.0*pi,2.0)*omega_comp(mom_p,m))
                                     *   mom_kprime*mom_kprime*weight_kprime/(pow(2.0*pi,2.0)*omega_comp(mom_kprime,m))
                                     *   imag(GS_pk(s,mom_p,mom_kprime,m,epsilon))*Dmatsol(k,j);
                }
                comp thirdterm = temp_third;

                Cmat(i,j) = firstterm + secondterm + thirdterm;
            }
            else 
            {
                comp mom_p = qvec[i];
                comp mom_k = qvec[j];
                comp weight_p = weights[i];
                comp weight_k = weights[j];
                comp sigp = sigma_p(s,mom_p,m);
                comp sigk = sigma_p(s,mom_k,m);
                double theta = 0.0;
                comp m2k_sigp = M2kfunc(a,sigp,m,epsilon);
                comp m2k_sigk = M2kfunc(a,sigk,m,epsilon);
                comp firstterm = ii/2.0 - rho_func(sigp, m, theta)*m2k_sigk;
                comp secondterm = -mom_p*mom_p*weight_p/(pow(2.0*pi,2.0)*omega_comp(mom_p,m))*(rho_func(sigp,m,theta)*Dmatsol(i,j) - imag(GS_pk(s,mom_p,mom_k,m,epsilon))*m2k_sigk);
                comp temp_third = {0.0,0.0};
                for(int k=0;k<qvec.size();++k)
                {
                    comp mom_kprime = qvec[k];
                    comp weight_kprime = weights[k];
                    temp_third       =   temp_third + mom_p*mom_p*weight_p/(pow(2.0*pi,2.0)*omega_comp(mom_p,m))
                                     *   mom_kprime*mom_kprime*weight_kprime/(pow(2.0*pi,2.0)*omega_comp(mom_kprime,m))
                                     *   imag(GS_pk(s,mom_p,mom_kprime,m,epsilon))*Dmatsol(k,j);
                }
                comp thirdterm = temp_third;

                Cmat(i,j) = secondterm + thirdterm;
            }
        }
    }
}

//assumming k = q always
void Bvec_Dsecsht_maker_with_weights(   Eigen::VectorXcd &Bvec,
                                        Eigen::VectorXcd dsol,
                                        vector<comp> qvec,
                                        vector<comp> weights, 
                                        comp s,
                                        comp q,
                                        double m,
                                        double a, 
                                        double epsilon )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);

    for(int i=0;i<qvec.size();++i)
    {
        comp mom_p = qvec[i];
        comp weight_p = weights[i];
        comp sigp = sigma_p(s,mom_p,m);
        double theta = 0.0;
        comp m2k_sigp = M2kfunc(a,sigp,m,epsilon);
        comp sigq = sigma_p(s,q,m);
        comp m2k_sigq = M2kfunc(a,sigq,m,epsilon);
        comp firstterm = ii/2.0*dsol(i);
        comp secondterm = conj(m2k_sigp)*rho_func(sigp,m,theta)*dsol(i);
        comp thirdterm = {0.0,0.0};
        for(int j=0;j<qvec.size();++j)
        {
            comp mom_k = qvec[j];
            comp weight_k = weights[j];
            comp temp_third = - mom_k*mom_k*weight_k/(pow(2.0*pi,2.0)*omega_comp(mom_p,m))*conj(m2k_sigp)*imag(GS_pk(s,mom_p,mom_k,m,epsilon))*dsol(j);
            thirdterm = thirdterm + temp_third;
        }

        comp forthterm = -conj(m2k_sigp)*imag(GS_pk(s,mom_p,q,m,epsilon))*m2k_sigq;

        Bvec(i) = firstterm + secondterm + thirdterm + forthterm;
    }

}



#endif