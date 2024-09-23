#include <bits/stdc++.h>
#include "integral_equations_2plus1_system.cpp"

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

    kernel_2plus1_system_ERE(E, pvec, pvec, weights, weights, m1, m2, eps, eps, eps, total_P, a, r); 

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

int main()
{
    //test_M2(); 
    //test_q2k_sigk();
    test_sigk_q2k();
    
    return 0;
}