#include <bits/stdc++.h>
#include "integral_equations_2plus1_system.h"
#include "mphib_testing.h"

void test_M2()
{
    double pmin = 1e-10;
    double pmax = 1.0; 
    double ppoints = 5; 

    double delp = std::abs(pmax - pmin)/ppoints; 

    std::vector<comp> pvec; 
    for(int i=0; i<ppoints; ++i)
    {
        double p = pmin + i*delp; 

        pvec.push_back(p); 
    }

    double E = 3.1; 
    double a = -2.0; 
    double r = 0.0; 
    double m1 = 1.0;
    double m2 = 0.5; 

    double eps = 0.0; 
    double total_P = 0.0; 

    std::vector<comp> weights; 
    Eigen::MatrixXcd kernmat;

    kernel_2plus1_system_ERE(kernmat, E, pvec, pvec, weights, weights, m1, m2, eps, eps, eps, total_P, a, r, a, r); 

}

void test_q2k_sigk()
{
    double En = 3.1; 
    double pmin = -1.5;
    double pmax = 1.5;
    double ppoints = 20;
    double mi = 1.0;
    double total_P = 0.0; 

    double delp = std::abs(pmax - pmin)/ppoints; 

    std::vector<comp> repvec;
    std::vector<comp> impvec; 

    for(int i=0;i<ppoints; ++i)
    {
        double p = pmin + i*delp; 

        repvec.push_back(p);
        impvec.push_back(ii*p); 
    }

    for(int i=0; i<repvec.size(); ++i)
    {
        comp sig1 = sigma(En, repvec[i], mi, total_P);
        comp sig2 = sigma(En, impvec[i], mi, total_P); 

        std::cout<<repvec[i]<<'\t'<<sig1<<'\t'<<impvec[i]<<'\t'<<sig2<<std::endl; 
    }
}

void test_sigk_q2k()
{
    double En = 3.1;
    double mi = 1.0; 

    double re_sigmin = 2.1;
    double im_sigmin = -1.0;
    double im_sigmax = 1.0;
    double sig_points = 20;
    double del_sig = std::abs(im_sigmin - im_sigmax)/sig_points; 

    for(int i=0; i<sig_points; ++i)
    {
        double im_sig = im_sigmin + i*del_sig; 
        comp sigk = re_sigmin + ii*im_sig; 

        comp p = pmom(En, sigk, mi);

        std::cout<<"sigk:"<<sigk<<'\t'<<"p="<<p<<std::endl;  
    }
}

void M2k_Plotter_KK()
{
    // atmpi = 0.06906
    // atmK = 0.09698
    // scattering_length_1_piK = 4.04;// - 0.2; //total uncertainty 0.05 stat 0.15 systematic 
    // scattering_length_2_KK = 4.07;// - 0.07; //total uncertainty 0.07 stat 
    // atinv = 5.666 GeV 
    // anisotropy xi = 3.444

    //double 
}

void test_kinematic_variables()
{
    double m1 = 1.0;
    double m2 = 0.5;

    double mi = m1;
    double mj = m1;
    double mk = m2; 

    double a01 = 2.0; 
    double a02 = 2.0; 

    double En = 2.6;

    comp kmax1 = pmom(En, 0.0, m1);
    comp kmax2 = pmom(En, 0.0, m2);
    comp sigb1 = sigma_b_plus(a01, mj, mk); 
    comp sigb2 = sigma_b_minus(a01, mj, mk); 

    comp qb = qb_i(En, sigb1, mi);

    std::cout<<"kmax1 = "<<kmax1<<std::endl; 
    std::cout<<"kmax2 = "<<kmax2<<std::endl; 
    std::cout<<"sigb1 = "<<sigb1<<std::endl; 
    std::cout<<"sigb2 = "<<sigb2<<std::endl; 
    std::cout<<"qb = "<<qb<<std::endl; 


}

void test_dpqb_building()
{
    double m1 = 1.0;
    double m2 = 0.5; 
    double a0_1 = 2.0; 
    double a0_2 = 2.0; 
    double r0_1 = 0.0; 
    double r0_2 = 0.0; 

    double En = 2.4;
    double total_P = 0.0; 
    double r = 0; 
    int number_of_points = 3;  

    double eps_for_m2k = 0.0;
    double eps_for_ope = 0.0; 
    double eps_for_cutoff = 0.0; 

    Eigen::MatrixXcd dpqbmat; 
    Eigen::MatrixXcd dqqmat; 
    std::vector<comp> pvec_for_m1m2; 
    std::vector<comp> weights_for_pvec_for_m1m2; 
    std::vector<comp> kvec_for_m1m1; 
    std::vector<comp> weights_for_kvec_for_m1m1; 
    comp qb; 

    dpqb_solver_ERE(dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points); 
    dqq_interpolator(dqqmat, dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points);

}
int main()
{
    //test_M2(); 
    //test_q2k_sigk();
    //test_sigk_q2k();
    //test_kinematic_variables();
    test_dpqb_building();
    
    return 0;
}