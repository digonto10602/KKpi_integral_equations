#ifndef FUNCTIONS_MOM_BASED_H
#define FUNCTIONS_MOM_BASED_H

#include<bits/stdc++.h>
#include<cmath>
#include<Eigen/Dense>

typedef std::complex<double> comp;

comp ii = {0.0,1.0};

/* Important Debug Functions */


/*---------------------------*/

comp mysqrt(    comp x  )
{
    comp ii = {0.0,1.0}; 
    return ii*std::sqrt(-x); 
}

comp kallentriangle(    comp x,
                        comp y,
                        comp z  )
{
    return x*x + y*y + z*z - 2.0*(x*y + y*z + z*x);
}

comp omega_func(    comp p, 
                    double m    )
{
    return std::sqrt(p*p + m*m);
}

comp sigma( comp En,
            comp spec_p,
            double mi,
            comp total_P    )
{
    comp A = En - omega_func(spec_p,mi);
    comp B = total_P - spec_p;

    return A*A - B*B;
}

/* This sigma takes the spectator momenta as a vector */
comp sigma_pvec_based(  comp En, 
                        std::vector<comp> p,
                        double mi, 
                        std::vector<comp> total_P   )
{
    comp px = p[0];
    comp py = p[1];
    comp pz = p[2];

    comp spec_p = std::sqrt(px*px + py*py + pz*pz);

    comp Px = total_P[0];
    comp Py = total_P[1];
    comp Pz = total_P[2];

    comp Pminusp_x = Px - px; 
    comp Pminusp_y = Py - py; 
    comp Pminusp_z = Pz - pz; 

    comp Pminusp_sq = Pminusp_x*Pminusp_x + Pminusp_y*Pminusp_y + Pminusp_z*Pminusp_z; 
    comp A = En - omega_func(spec_p,mi);

    return A*A - Pminusp_sq; 
}

comp Jfunc( comp z  )
{
    if(std::real(z)<=0.0)
    {
        return 0.0;
    }
    else if(std::real(z)>0.0 && std::real(z)<1.0)
    {
        comp A = -1.0/z;
        comp B = std::exp(-1.0/(1.0-z));
        return std::exp(A*B);
    }
    else
    {
        return 1.0;
    }
    
    

}

comp cutoff_function_1( comp sigma_i,
                        double mj, 
                        double mk, 
                        double epsilon_h   )
{

    if(mj==mk && epsilon_h==0.0)
    {
        comp Z = sigma_i/(4.0*mj*mj);
        return Jfunc(Z); 
    }
    else 
    {
        //std::cout<<"went to else"<<std::endl; 
        comp Z = (comp) (1.0 + epsilon_h)*( sigma_i - (comp) std::abs(mj*mj - mk*mk) )/( (mj + mk)*(mj + mk) - std::abs(mj*mj - mk*mk) );
        return Jfunc(Z);
        
    }

    //return 1.0; 

}

//The ij components of the G function depends on how 
//we define mi, mj, and mk. The momentum p corresponds 
//to mi and spectator momentum k corresponds to mj. The 
//exchanged particle has mass mk. 
comp GS_pk( comp En, 
            comp p, 
            comp k, 
            double mi,  
            double mj, 
            double mk,
            double eps, //This is the epsilon for OPE
            double epsilon_h //This epsilon is for cutoff function 
        )
{
    comp ii = {0.0,1.0};

    comp sigp = sigma(En, p, mi, 0.0);
    comp sigk = sigma(En, k, mj, 0.0); 

    comp Hp = cutoff_function_1(sigp, mj, mk, epsilon_h);
    comp Hk = cutoff_function_1(sigk, mi, mk, epsilon_h); 

    comp omgp = omega_func(p, mi);
    comp omgk = omega_func(k, mj);

    comp A = En - omgp - omgk; 
    comp zpk = A*A - p*p - k*k - mk*mk; 

    comp pk = p*k; 

    comp num = zpk + ii*eps - 2.0*pk;
    comp denom = zpk + ii*eps + 2.0*pk;

    comp result = - Hp*Hk/(4.0*pk) * std::log(num/denom); 
    //comp result = - 1.0/(4.0*pk) * std::log(num/denom); 

    return result; 

}

comp pmom(  comp En,
            comp sig_i, 
            double mi   )
{
    comp kalln = kallentriangle(En*En, mi*mi, sig_i);
    comp sqrtkalln = mysqrt(kalln);//std::sqrt(kalln); 

    return sqrtkalln/(2.0*En); 
}


comp q2k(   comp sig_i, 
            double mj, 
            double mk   )
{
    comp klntrngle = kallentriangle(sig_i, mj*mj, mk*mk);

    comp sqrtklntrngle = mysqrt(klntrngle); 

    comp denom = 2.0*std::sqrt(sig_i);

    return sqrtklntrngle/denom;  
}


//This shows the spectator momentum for a given flavor for 
//which the two-body subsystem forms a bound-state.
comp qb_i(  comp En, 
            comp sigb, 
            double mi   )
{
    comp klntrngle = kallentriangle(En*En, mi*mi, sigb);

    comp sqrtklntrngle = mysqrt(klntrngle);

    comp denom = 2.0*En; 

    return sqrtklntrngle/denom; 
}

//This shows the mass squared of the two-body bound-state
//depending on the scattering length a^{i}_0,
//the masses of the particles of the two-body 
//subsystem, mj and mk 
comp sigma_b_plus(  double scattering_length_a0i,
                    double mj, 
                    double mk )
{
    comp a0i = scattering_length_a0i; 
    comp firstterm = - 2.0/(a0i*a0i);
    comp secondterm = mj*mj + mk*mk; 
    comp midterm1 = mj*mj - 1.0/(a0i*a0i); 
    comp midterm2 = mk*mk - 1.0/(a0i*a0i); 
    comp thirdterm = 2.0*std::sqrt( midterm1*midterm2 );
    comp result = firstterm + secondterm + thirdterm;

    return result; 
}

comp sigma_b_minus( double scattering_length_a0i,
                    double mj, 
                    double mk )
{
    comp a0i = scattering_length_a0i; 
    comp firstterm = - 2.0/(a0i*a0i);
    comp secondterm = mj*mj + mk*mk; 
    comp midterm1 = mj*mj - 1.0/(a0i*a0i); 
    comp midterm2 = mk*mk - 1.0/(a0i*a0i); 
    comp thirdterm = 2.0*std::sqrt( midterm1*midterm2 );

    comp result = firstterm + secondterm - thirdterm;

    return result; 
}

//This is the kinematic phase-space for the 
//bound-state spectator system 
comp rhophib(   comp qb, 
                comp En )
{
    double pi = std::acos(-1.0); 
    return qb/(8.0*pi*En);
}


//This is the two-body residue for the bound-state
//for the degenerate mass case
comp gfunc( comp sigb, 
            double scattering_length_a0 )
{
    double pi = std::acos(-1.0); 
    return 8.0*std::sqrt(2.0*pi*std::sqrt(sigb)/scattering_length_a0);
}

//This is the two-body residue for 2+1 system
comp gfunc_i(   double eta_i, 
                comp sigb, 
                double mj, 
                double mk   )
{
    comp ii = {0.0, 1.0}; 
    double pi = std::acos(-1.0); 
    comp lamb = kallentriangle(sigb, mj*mj, mk*mk);

    comp termA = 32.0*pi*sigb*sigb*std::sqrt(lamb);
    comp termB = sigb - mj*mj + mk*mk; 
    comp termC = sigb + mj*mj - mk*mk; 

    comp res = std::sqrt(-ii*termA/(eta_i*termB*termC));

    if(std::abs(sigb.real())<1.0e-10 && std::abs(sigb.imag())<1.0e-10)
    {
        res = 0.0; 
        std::cout<<"g residue set to 0 because sigb=0"<<std::endl; 
    }
    
    return res; 
}

//Here we define the energy dependent epsilon 
//defined in https://arxiv.org/pdf/2010.09820
double energy_dependent_epsilon(    double eta, 
                                    comp En, 
                                    comp qb, //relative spectator momentum for two-body boundstate
                                    comp sigb, //mass square of the two-body bound-state
                                    comp kmax, //maximum spectator momentum
                                    double mi, //mass of the spectator 
                                    double number_of_points )
{
    double pi = std::acos(-1.0); 

    comp termA = (eta*kmax)/(2.0*pi*number_of_points);
    comp termB = (4.0*qb*En*En)/(En*En + mi*mi - sigb);

    return ((double)real(termA*termB)); 
}

//This is M2k_tilde, where the symmetry 
//factor is not multiplied
comp M2k_tilde_ERE( comp E, 
                    comp k, //spectator momentum
                    comp total_P,  
                    double a, //This is the scattering length
                    double r, //This is the effective range 
                    double mi, //mass of the spectator  
                    double mj, 
                    double mk, 
                    double eps  )
{
    comp ii = {0.0, 1.0}; 
    double pi = std::acos(-1.0); 

    comp sigk = sigma(E, k, mi, total_P);
    comp q = q2k(sigk + ii*eps, mj, mk); 

    comp num = 8.0*pi*sqrt(sigk);
    comp denom = (-1.0/a + r*q*q/2.0 - ii*q); 

    return num/denom; 

}

comp M2k_ERE(   double eta_i,
                comp E, 
                comp k, //spectator momentum
                comp total_P,  
                double a, //This is the scattering length
                double r, //This is the effective range 
                double mi, //mass of the spectator  
                double mj, 
                double mk, 
                double eps  )
{
    comp ii = {0.0, 1.0}; 
    double pi = std::acos(-1.0); 

    comp sigk = sigma(E, k, mi, total_P);
    comp q = q2k(sigk + ii*eps, mj, mk); 

    comp num = 8.0*pi*sqrt(sigk);
    comp denom = eta_i*(-1.0/a + r*q*q/2.0 - ii*q); 

    return num/denom; 

}

//This M2k is sigk based, used mostly for plotting
comp M2k_ERE_s2k(   double eta_i, 
                    comp E, 
                    comp sigk, //energy squared of the two-body system
                    comp total_P,  
                    double a, //This is the scattering length
                    double r, //This is the effective range 
                    double mi, //mass of the spectator  
                    double mj, 
                    double mk, 
                    double eps  )
{
    comp ii = {0.0, 1.0}; 
    double pi = std::acos(-1.0); 

    //comp sigk = sigma(E, k, mi, total_P);
    comp q = q2k(sigk + ii*eps, mj, mk); 

    comp num = 8.0*pi*std::sqrt(sigk + ii*eps);
    comp denom = eta_i*(-1.0/a + r*q*q/2.0 - ii*q); 

    return num/denom; 

}


#endif 

