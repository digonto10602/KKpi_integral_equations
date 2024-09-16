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

    kernel_2plus1_system_ERE(E, pvec, pvec, m1, m2, eps, eps, eps, total_P, a, r); 

}

int main()
{
    test_M2(); 
    
    return 0;
}