#ifndef FUNCTIONS_MOMREP_BASED
#define FUNCTIONS_MOMREP_BASED

#include<bits/stdc++.h>
#include<Eigen/Dense>
using namespace std;

typedef complex<double> comp;

comp ii = {0.0,1.0};

comp mysqrt(    comp x  )
{
    comp ii = {0.0,1.0};
    return ii*sqrt(-x);
}

comp mysqrt1(    comp x  )
{
    comp ii = {0.0,1.0};
    return sqrt(-x);
}

comp kallentriangle(    comp x,
                        comp y,
                        comp z  )
{
    return x*x + y*y + z*z - 2.0*(x*y + y*z + z*x);
}

comp zppprime(   comp s,
                comp sigmap,
                comp sigmapprime,
                double m,
                double epsilon   )
{
    comp num = 2.0*s*(sigmap + ii*epsilon) - (s + sigmap - m*m)*(s + m*m - sigmapprime);
    comp denom = mysqrt(kallentriangle(s,sigmapprime,m*m))*mysqrt(kallentriangle(s,sigmap,m*m));

    return num/denom;
}

/* This are discarded because they are wrong.*/
/*
comp eps_energy_factor_plus(    comp s,
                                comp sigk,
                                double m    )
{
    return (1.0/(2.0*sigk) + (1.0/(2.0*sqrt(sigk)))*sqrt(kallentriangle(s,sigk,m*m))/(sqrt(sigk - 4.0*m*m)));
}

comp eps_energy_factor_minus(    comp s,
                                comp sigk,
                                double m    )
{
    return (1.0/(2.0*sigk) - (1.0/(2.0*sqrt(sigk)))*sqrt(kallentriangle(s,sigk,m*m))/(sqrt(sigk - 4.0*m*m)));
}*/

comp eps_energy_factor_plus(    comp s,
                                comp sigk,
                                double m    )
{
    return (1.0/(2.0*sigk)*(s + sigk - m*m) + (1.0/(2.0*sqrt(sigk)))*sqrt(kallentriangle(s,sigk,m*m))/(sqrt(sigk - 4.0*m*m)));
}

comp eps_energy_factor_minus(    comp s,
                                comp sigk,
                                double m    )
{
    return (1.0/(2.0*sigk)*(s + sigk - m*m) - (1.0/(2.0*sqrt(sigk)))*sqrt(kallentriangle(s,sigk,m*m))/(sqrt(sigk - 4.0*m*m)));
}


double Jfunc(   double x    )
{
    if(x<=0.0) return 0.0;
    else if(x>0.0 && x<1.0 ) return exp((-1.0/x)*exp(-1.0/(1.0-x)));
    else return 1.0;
}

comp Jfunc_comp( comp x  )
{
    if(real(x)<=0.0) return 0.0;
    else if(real(x)>0.0 && real(x)<1.0) return exp((-1.0/x)*exp(-1.0/(1.0-x)));
    else return 1.0;
}

double Hfunc(   comp sigmap,
                comp sigmapprime,
                double m    )
{
    return Jfunc(real(sigmap/(4.0*m*m)))*Jfunc(real(sigmapprime/(4.0*m*m)));
}

comp Hfunc_comp(    comp sigmap,
                    comp sigmapprime,
                    double m    )
{
    return Jfunc_comp(sigmap/(4.0*m*m))*Jfunc_comp(sigmapprime/(4.0*m*m));
}

comp GS(    comp s,
            comp sigmap,
            comp sigmapprime,
            double m,
            double epsilon      )
{
    comp p = mysqrt(kallentriangle(s,sigmap,m*m))/(2.0*mysqrt(s));
    comp pprime = mysqrt(kallentriangle(s,sigmapprime,m*m))/(2.0*mysqrt(s));
    //cout<<p<<'\t'<<pprime<<endl;
    comp zpp = zppprime(s,sigmap,sigmapprime,m,epsilon);

    return -((comp)Hfunc(sigmap,sigmapprime,m))/(4.0*p*pprime)*(log( (zpp - 1.0)/(zpp + 1.0) ));

}

//-------------MomRep------------------------------//

comp omega_comp(    comp p,
                    double m    )
{
    return sqrt(p*p + m*m);
}

comp sigma_p(   comp s,
                comp p,
                double m    )
{
    return (sqrt(s) - omega_comp(p,m))*(sqrt(s) - omega_comp(p,m)) - p*p;
}

comp M2kbranchcut_right_momrep_plus(     comp s,
                                        double m    )
{
    return + sqrt((s - 9.0*m*m)*(s - m*m))/(2.0*sqrt(s));
}

comp M2kbranchcut_right_momrep_minus(     comp s,
                                        double m    )
{
    return - sqrt((s - 9.0*m*m)*(s - m*m))/(2.0*sqrt(s));
}

comp M2kbranchcut_left_momrep_plus(     comp s,
                                        double m    )
{
    return + (s - m*m)/(2.0*sqrt(s));
}

comp M2kbranchcut_left_momrep_minus(     comp s,
                                        double m    )
{
    return - (s - m*m)/(2.0*sqrt(s));
}



comp M2kbranchcut_right_momrep_plus_eps(     comp s,
                                             double m,
                                             double eps)
{
    comp ii = {0.0,1.0};
    return + sqrt(pow(s + ii*eps-3.0*m*m,2.0) - 4.0*s*m*m)/(2.0*sqrt(s));
}
comp M2kbranchcut_right_momrep_minus_eps(     comp s,
                                             double m,
                                             double eps)
{
    comp ii = {0.0,1.0};
    return - sqrt(pow(s + ii*eps-3.0*m*m,2.0) - 4.0*s*m*m)/(2.0*sqrt(s));
}



comp GS_pk( comp s,
            comp p,
            comp k,
            double m,
            double eps  )
{
    comp ii = {0.0,1.0};
    comp sigp = sigma_p(s,p,m);
    comp sigk = sigma_p(s,k,m);
    comp alphapk = (sqrt(s)-omega_comp(p,m)-omega_comp(k,m))*(sqrt(s)-omega_comp(p,m)-omega_comp(k,m)) - p*p - k*k - m*m;
    comp num = alphapk - 2.0*p*k + ii*eps;
    comp denom = alphapk + 2.0*p*k + ii*eps;

    //smooth cutoff
    //return -Hfunc_comp(sigp,sigk,m)/(4.0*p*k)*(log(num/denom));
    
    //hard cutoff
    return -1.0/(4.0*p*k)*(log(num/denom));
    

    //above threshold
    //return -Hfunc(sigp,sigk,m)/(4.0*p*k)*(log(num/denom));
    
}

//This GS uses 2 indices and 2 tags, if the indices are between the 
//tags then the OPE is evaluated on the unphysical sheet by adding the 
//discontinuity with the original OPE 
comp GS_pk_withtags(    comp s,
                        comp p,
                        comp k,
                        double m,
                        double eps,
                        int index1,
                        int index2,
                        int tag1,
                        int tag2  )
{
    comp ii = {0.0,1.0};
    comp sigp = sigma_p(s,p,m);
    comp sigk = sigma_p(s,k,m);
    comp alphapk = (sqrt(s)-omega_comp(p,m)-omega_comp(k,m))*(sqrt(s)-omega_comp(p,m)-omega_comp(k,m)) - p*p - k*k - m*m;
    comp num = alphapk - 2.0*p*k + ii*eps;
    comp denom = alphapk + 2.0*p*k + ii*eps;
    double pi = acos(-1.0);

    //if(abs(p)<1.0e-10 || abs(k)<1.0e-10)
    //return Hfunc_comp(sigp,sigk,m)/(alphapk + ii*eps);
    //return -1.0/(alphapk + ii*eps);
    
    //else
    {
        if((index1>=tag1 && index1<=tag2) || (index2>=tag1 && index2<=tag2))
        //smooth cutoff
        //return -(Hfunc_comp(sigp,sigk,m)/(4.0*p*k))*(log(num/denom)) -(Hfunc_comp(sigp,sigk,m)/(4.0*p*k))*(2.0*pi*ii)  ;
        //hard cutoff
        return -(1.0/(4.0*p*k))*(log(num/denom)) -(1.0/(4.0*p*k))*(2.0*pi*ii)  ;
        
        else 
        //smooth cutoff
        //return -(Hfunc_comp(sigp,sigk,m)/(4.0*p*k))*(log(num/denom));
        
        //hard cutoff 
        return -(1.0/(4.0*p*k))*(log(num/denom));
    
    }

    
}

comp GS_pk_secondsheet(    comp s,
                        comp p,
                        comp k,
                        double m,
                        double eps )
{
    comp ii = {0.0,1.0};
    comp sigp = sigma_p(s,p,m);
    comp sigk = sigma_p(s,k,m);
    comp alphapk = (sqrt(s)-omega_comp(p,m)-omega_comp(k,m))*(sqrt(s)-omega_comp(p,m)-omega_comp(k,m)) - p*p - k*k - m*m;
    comp num = alphapk - 2.0*p*k + ii*eps;
    comp denom = alphapk + 2.0*p*k + ii*eps;
    double pi = acos(-1.0);

    //if(abs(p)<1.0e-10 || abs(k)<1.0e-10)
    //return Hfunc_comp(sigp,sigk,m)/(alphapk + ii*eps);
    //return -1.0/(alphapk + ii*eps);
    
    //smooth cutoff
    //return -(Hfunc_comp(sigp,sigk,m)/(4.0*p*k))*(log(num/denom)) -(Hfunc_comp(sigp,sigk,m)/(4.0*p*k))*(2.0*pi*ii)  ;
    //hard cutoff 
    return -(1.0/(4.0*p*k))*(log(num/denom)) -(1.0/(4.0*p*k))*(2.0*pi*ii)  ;
        
    
}


//------------------------------------------------------------//

comp sigmab(    double a,
                double m    )
{
    return 4.0*(m*m - 1.0/(a*a));
}

comp tau(   comp s,
            comp sigk,
            double m    )
{
    double pi = acos(-1.0);
    return sqrt(kallentriangle(s,sigk,m*m))/(8.0*pi*s);
}

comp tau1(  comp s,
            comp sigk,
            double m    )
{
    double pi = acos(-1.0);
    return (sqrt(-sigk + (sqrt(s) - m)*(sqrt(s) - m))*sqrt(-sigk + (sqrt(s) + m)*(sqrt(s) + m)))/(8.0*pi*s);
}

comp pmom(  comp s,
            comp sigk,
            double m    )
{
    return sqrt(kallentriangle(s,sigk,m*m))/(2.0*sqrt(s));
}

comp M2kfunc(   double a,
                comp sigk,
                double m,
                double epsilon    )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    //cout<<"sigk:"<<sigk<<endl;
    comp num = 16.0*pi*sqrt(sigk);
    
    comp denom = -1.0/a - ii*mysqrt((sigk+ii*epsilon)/4.0 - m*m);

    //cout<<"num:"<<num<<endl;
    //cout<<"denom:"<<denom<<endl;

    return num/denom;
}

comp M2kfunc_2mysqrt(   double a,
                comp sigk,
                double m,
                double epsilon    )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    //cout<<"sigk:"<<sigk<<endl;
    comp num = 16.0*pi*sqrt(sigk);
    
    comp denom = -1.0/a - ii*sqrt((sigk+ii*epsilon)/4.0 - m*m);

    //cout<<"num:"<<num<<endl;
    //cout<<"denom:"<<denom<<endl;

    return num/denom;
}

comp M2kfunc_secondsheet(   double a,
                            comp sigk,
                            double m, 
                            double epsilon  )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    //cout<<"sigk:"<<sigk<<endl;
    comp num = 16.0*pi*sqrt(sigk);
    
    comp denom = -1.0/a - ii*mysqrt((sigk+ii*epsilon)/4.0 - m*m);

    comp m2k = num/denom; 

    comp rho = sqrt(sigk - 4.0*m*m)/(32.0*pi*sqrt(sigk));

    return m2k/(1.0 + 2.0*ii*rho*m2k);
}

comp M2kfunc_1( double a,
                comp sigk,
                double m,
                double epsilon  )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    comp sigb = sigmab(a,m);

    comp A = -64.0*M_PI*sqrt(sigk);
    comp B = (1.0/a - ii*mysqrt(sigk/4.0 - m*m));
    comp C = 1.0/(sigk - sigb + ii*epsilon);

    return A*B*C;
}

comp M2kfunc_2(   double a,
                comp sigk,
                double m,
                double epsilon    )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    //cout<<"sigk:"<<sigk<<endl;
    comp num = 16.0*pi*sqrt(sigk);
    comp denom = -1.0/a - sqrt(-(sigk)/4.0 + m*m) + ii*epsilon;

    return num/denom;
}

comp arg_theta(    comp z,
                    double theta    )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    return arg(z*exp(-ii*(theta + pi))) + theta + pi;
}

comp sqrt_theta(    comp z,
                    double theta    )
{
    double pi = acos(-1.0);
    comp ii = {0.0,1.0};
    return sqrt(abs(z))*exp(ii/2.0 * arg_theta(z,theta));
}

comp m2k_cut_tracer(    comp sigk,
                        double m,
                        double theta ) 
{
    comp ii = {0.0,1.0};
    if(real(sigk)/(m*m)>=4.0)
    {
        double l = abs(real(sigk) - 4.0*m*m);
        double r = l/cos(theta);

        double yval = r*sin(theta);
        double xval = r*cos(theta);

        return 4.0 + xval + ii*yval;
    }
    else 
    return 0.0;
}


comp M2kfunc_rotatedcut(    double a,
                            comp sigk,
                            double m,
                            double epsilon,
                            double theta    )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    //cout<<"sigk:"<<sigk<<endl;
    comp num = 16.0*pi*sqrt(sigk);
    comp argument = (sigk+ii*epsilon)/4.0 - m*m;
    comp denom = -1.0/a - ii*sqrt_theta(argument, theta);

    return num/denom;
}

comp M2kfunc_secondsheet_rotatedcut(    double a,
                            comp sigk,
                            double m,
                            double epsilon,
                            double theta    )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    //cout<<"sigk:"<<sigk<<endl;
    comp num = 16.0*pi*sqrt(sigk);
    comp argument = (sigk+ii*epsilon)/4.0 - m*m;
    comp denom = -1.0/a - ii*sqrt_theta(argument, theta);

    comp m2k =  num/denom;

    comp rho = sqrt_theta(sigk - 4.0*m*m,theta)/(32.0*pi*sqrt(sigk));

    return m2k/(1.0 + 2.0*ii*rho*m2k);
}

//this was built for using inside D_second sheet testing
comp rho_func(  comp sigk,
                double m,
                double theta )
{
    double pi = acos(-1.0);
    comp rho = sqrt_theta(sigk - 4.0*m*m,theta)/(32.0*pi*sqrt(sigk));

    return rho;
}

comp kernel(    comp s,
                comp sigp,
                comp sigk,
                double a,
                double m,
                double epsilon )
{
    double pi = acos(-1.0);
    comp picomp = (comp) pi;
    //cout<<"Gs:"<<GS(s,sigp,sigk,m,epsilon)<<endl;
    //cout<<"tau:"<<tau(s,sigk,m)<<endl;
    //cout<<"M2:"<<M2kfunc(a,sigk,m,epsilon)<<endl;
    return ((comp)1.0/(2.0*picomp))*GS(s,sigp,sigk,m,epsilon)*tau1(s,sigk,m)*M2kfunc(a,sigk,m,epsilon);
    //return GS(s,sigp,sigk,m,epsilon);
}

comp kernel_1(    comp s,
                comp sigp,
                comp sigk,
                double a,
                double m,
                double epsilon )
{
    double pi = acos(-1.0);
    comp picomp = (comp) pi;
    //cout<<"Gs:"<<GS(s,sigp,sigk,m,epsilon)<<endl;
    //cout<<"tau:"<<tau(s,sigk,m)<<endl;
    //cout<<"M2:"<<M2kfunc(a,sigk,m,epsilon)<<endl;
    //return GS(s,sigp,sigk,m,epsilon)*tau(s,sigk,m)*M2kfunc_1(a,sigk,m,epsilon);
    return ((comp)1.0/(2.0*picomp))*GS(s,sigp,sigk,m,epsilon)*tau(s,sigk,m)*M2kfunc_1(a,sigk,m,epsilon);
    //return GS(s,sigp,sigk,m,epsilon);
}

//---------------------MomRep-----------------------//

comp kernel_pk(    comp s,
                comp p,
                comp k,
                double a,
                double m,
                double epsilon )
{
    double pi = acos(-1.0);
    comp picomp = (comp) pi;
    //cout<<"Gs:"<<GS(s,sigp,sigk,m,epsilon)<<endl;
    //cout<<"tau:"<<tau(s,sigk,m)<<endl;
    //cout<<"M2:"<<M2kfunc(a,sigk,m,epsilon)<<endl;
    comp sigk = sigma_p(s,k,m);
    
    //return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*GS_pk(s,p,k,m,epsilon)*M2kfunc(a,sigk,m,epsilon);
    return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*GS_pk(s,p,k,m,epsilon)*M2kfunc(a,sigk,m,0.0);
    
    
    //return ((comp)1.0/(2.0*picomp))*GS(s,sigp,sigk,m,epsilon)*tau1(s,sigk,m)*M2kfunc(a,sigk,m,epsilon);
    //return GS(s,sigp,sigk,m,epsilon);
}

comp kernel_pk_2eps(    comp s,
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
    
    //return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*GS_pk(s,p,k,m,epsilon)*M2kfunc(a,sigk,m,epsilon);
    return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*GS_pk(s,p,k,m,epsilon)*M2kfunc(a,sigk,m,eps_for_m2k);
    
    //return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*M2kfunc(a,sigk,m,eps_for_m2k);
    
    
    //return ((comp)1.0/(2.0*picomp))*GS(s,sigp,sigk,m,epsilon)*tau1(s,sigk,m)*M2kfunc(a,sigk,m,epsilon);
    //return GS(s,sigp,sigk,m,epsilon);
}

comp kernel_pk_2eps_1(    comp s,
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
    
    if(real(s)<9.0*m*m)
    return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*GS_pk(s,p,k,m,epsilon)*M2kfunc(a,sigk,m,eps_for_m2k);
    else 
    return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*GS_pk(s,p,k,m,epsilon)*M2kfunc_2mysqrt(a,sigk,m,eps_for_m2k);
}

comp kernel_pk_2eps_for_resonances(    comp s,
                        comp p,
                        comp k,
                        double a,
                        double m,
                        double epsilon,
                        double eps_for_m2k,
                        double theta,
                        int sheet_tag )
{
    double pi = acos(-1.0);
    comp picomp = (comp) pi;
    //cout<<"Gs:"<<GS(s,sigp,sigk,m,epsilon)<<endl;
    //cout<<"tau:"<<tau(s,sigk,m)<<endl;
    //cout<<"M2:"<<M2kfunc(a,sigk,m,epsilon)<<endl;
    comp sigk = sigma_p(s,k,m);
    
    if(sheet_tag==0)
    return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*GS_pk(s,p,k,m,epsilon)*M2kfunc_rotatedcut(a,sigk,m,eps_for_m2k,theta);
    //return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*GS_pk(s,p,k,m,epsilon);
    
    else 
    return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*GS_pk(s,p,k,m,epsilon)*M2kfunc_secondsheet_rotatedcut(a,sigk,m,eps_for_m2k,theta);
    //return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*GS_pk(s,p,k,m,epsilon);
    

}

comp kernel_pk_2eps_rotated_m2k_for_resonances(    comp s,
                        comp p,
                        comp k,
                        double a,
                        double m,
                        double epsilon,
                        double eps_for_m2k,
                        double theta )
{
    double pi = acos(-1.0);
    comp picomp = (comp) pi;
    //cout<<"Gs:"<<GS(s,sigp,sigk,m,epsilon)<<endl;
    //cout<<"tau:"<<tau(s,sigk,m)<<endl;
    //cout<<"M2:"<<M2kfunc(a,sigk,m,epsilon)<<endl;
    comp sigk = sigma_p(s,k,m);
    
    return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*GS_pk(s,p,k,m,epsilon)*M2kfunc_rotatedcut(a,sigk,m,eps_for_m2k,theta);
    

}



comp kernel_pk_2eps_withtags(       comp s,
                                    comp p,
                                    comp k,
                                    double a,
                                    double m,
                                    double epsilon,
                                    double eps_for_m2k,
                                    int index1,
                                    int index2,
                                    int tag1,
                                    int tag2  )
{
    double pi = acos(-1.0);
    comp picomp = (comp) pi;
    //cout<<"Gs:"<<GS(s,sigp,sigk,m,epsilon)<<endl;
    //cout<<"tau:"<<tau(s,sigk,m)<<endl;
    //cout<<"M2:"<<M2kfunc(a,sigk,m,epsilon)<<endl;
    comp sigk = sigma_p(s,k,m);
    
    //return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*GS_pk(s,p,k,m,epsilon)*M2kfunc(a,sigk,m,epsilon);
    return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*GS_pk_withtags(s,p,k,m,epsilon,index1,index2,tag1,tag2)*M2kfunc(a,sigk,m,eps_for_m2k);
    
    
    //return ((comp)1.0/(2.0*picomp))*GS(s,sigp,sigk,m,epsilon)*tau1(s,sigk,m)*M2kfunc(a,sigk,m,epsilon);
    //return GS(s,sigp,sigk,m,epsilon);
}

comp kernel_pk_2eps_using_Gvec(     Eigen::VectorXcd Gvec,
                                    int Gvec_index,
                                    comp s,
                                    comp p,
                                    comp k,
                                    double a,
                                    double m,
                                    double epsilon,
                                    double eps_for_m2k,
                                    int index1,
                                    int index2,
                                    int tag1,
                                    int tag2  )
{
    //momentum p here is redundant, it uses the outer momenta of Gvec
    //in most cases which is set to the two body bound state spectator 
    //momenta, q 

    double pi = acos(-1.0);
    comp picomp = (comp) pi;
    //cout<<"Gs:"<<GS(s,sigp,sigk,m,epsilon)<<endl;
    //cout<<"tau:"<<tau(s,sigk,m)<<endl;
    //cout<<"M2:"<<M2kfunc(a,sigk,m,epsilon)<<endl;
    comp sigk = sigma_p(s,k,m);
    
    //return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*GS_pk(s,p,k,m,epsilon)*M2kfunc(a,sigk,m,epsilon);
    //return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*GS_pk_withtags(s,p,k,m,epsilon,index1,index2,tag1,tag2)*M2kfunc(a,sigk,m,eps_for_m2k);
    return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*Gvec[Gvec_index]*M2kfunc(a,sigk,m,eps_for_m2k);
    
    
    //return ((comp)1.0/(2.0*picomp))*GS(s,sigp,sigk,m,epsilon)*tau1(s,sigk,m)*M2kfunc(a,sigk,m,epsilon);
    //return GS(s,sigp,sigk,m,epsilon);
}


//-------------------------------------------------//


comp sigmax(    comp s,
                double m    )
{
    return pow(sqrt(s) - m,2.0);
}

comp sigpplus(  comp s,
                comp sigk,
                double m    )
{
    return (1/2.0)*(s - sigk + 3.0*m*m) + (1.0/(2.0*sqrt(sigk))*sqrt(sigk - 4.0*m*m)*mysqrt(kallentriangle(s,sigk,m*m)));
}

comp sigpminus(  comp s,
                comp sigk,
                double m    )
{
    return (1/2.0)*(s - sigk + 3.0*m*m) - (1.0/(2.0*sqrt(sigk))*sqrt(sigk - 4.0*m*m)*mysqrt(kallentriangle(s,sigk,m*m)));
}

comp sigc1( comp s,
            comp sigk,
            double m    )
{
    return (m*m - s)*(m*m - sigk + s )/(m*m - sigk - s);
}

comp sigc2( comp s,
            comp sigk,
            double m    )
{
    return (m*m - s)*(2.0*m*m - sigk)/sigk;
}

comp sigpplus_witheps(  comp s,
                        comp sigk,
                        double m,
                        double eps  )
{
    comp ii = {0.0,1.0};
    return sigpplus(s,sigk,m) - ii*eps*eps_energy_factor_minus(s,sigk,m);
    
}

comp sigpminus_witheps(  comp s,
                        comp sigk,
                        double m,
                        double eps  )
{
    comp ii = {0.0,1.0};
    return sigpminus(s,sigk,m) - ii*eps*eps_energy_factor_plus(s,sigk,m);
    
}


//-----------------------MomRep---------------------//
comp q_plus(    comp s,
                double a,
                double m,
                double eps  )
{
    comp sigb = sigmab(a,m);
    comp sigppls = sigpplus_witheps(s,sigb,m,eps);

    return sqrt(kallentriangle(s,sigppls,m*m))/(2.0*sqrt(s));
}

comp q_minus(    comp s,
                double a,
                double m,
                double eps  )
{
    comp sigb = sigmab(a,m);
    comp sigpmns = sigpminus_witheps(s,sigb,m,eps);

    return sqrt(kallentriangle(s,sigpmns,m*m))/(2.0*sqrt(s));
}

comp q_c1(    comp s,
                double a,
                double m,
                double eps  )
{
    comp sigb = sigmab(a,m);
    comp sigpc1 = sigc1(s,sigb,m);

    return sqrt(kallentriangle(s,sigpc1,m*m))/(2.0*sqrt(s));
}

comp q_c2(    comp s,
                double a,
                double m,
                double eps  )
{
    comp sigb = sigmab(a,m);
    comp sigpc2 = sigc2(s,sigb,m);

    return sqrt(kallentriangle(s,sigpc2,m*m))/(2.0*sqrt(s));
}

comp betaxqs(   double x,
                comp s,
                comp q,
                double m      )
{
    return (omega_comp(q,m) - sqrt(s))*(omega_comp(q,m) - sqrt(s)) - x*x*q*q;
}

comp betaxqs_comp(   comp x,
                comp s,
                comp q,
                double m      )
{
    //return (omega_comp(q,m) - sqrt(s))*(omega_comp(q,m) - sqrt(s)) - x*x*q*q;
    return pow(omega_comp(q,m) - sqrt(s),2.0) - pow(x*q,2.0);

}

comp sqbeta0qs(     comp s,
                    comp q,
                    double m    )
{
    return omega_comp(q,m) - sqrt(s);
}


comp pcut_plus(     double x,
                    comp s,
                    comp q,
                    double m,
                    double eps  )
{
    comp ii = {0.0,1.0};
    comp firstterm = -q*x*betaxqs(1.0,s,q,m);
    comp secondterm = sqbeta0qs(s,q,m);
    comp thirdtermsq = (betaxqs(1.0,s,q,m)+ii*eps)*(betaxqs(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs(x,s,q,m);
    comp thirdterm;
    
    if(real(thirdtermsq)>=0.0 && imag(thirdtermsq)>=0.0)
    {
        thirdterm = mysqrt(abs(real(thirdtermsq)) + ii*abs(imag(thirdtermsq)));
    }
    else if(real(thirdtermsq)>=0.0 && imag(thirdtermsq)<0.0)
    {
        thirdterm = sqrt(abs(real(thirdtermsq)) - ii*abs(imag(thirdtermsq)));
    }
    else if(real(thirdtermsq)<0.0 && imag(thirdtermsq)>=0.0)
    {
        thirdterm = ii*sqrt(abs(real(thirdtermsq)) - ii*abs(imag(thirdtermsq)));
    }
    else if(real(thirdtermsq)<0.0 && imag(thirdtermsq)<0.0)
    {
        thirdterm = ii*sqrt(abs(real(thirdtermsq)) + ii*abs(imag(thirdtermsq)));
    }
    comp forthterm = 2.0*betaxqs(x,s,q,m);

    return (firstterm + secondterm*thirdterm)/forthterm;
}


comp pcut_plus_comp(     comp x,
                    comp s,
                    comp q,
                    double m,
                    double eps  )
{
    comp ii = {0.0,1.0};
    if(imag(x)==0.0)
    {
        x = real(x) + ii*1.0e-15;
    }
    comp firstterm = -q*x*betaxqs_comp(1.0,s,q,m);
    //cout<<"firstterm:"<<firstterm<<endl;
    comp secondterm = sqbeta0qs(s,q,m);
    //cout<<"secondterm:"<<secondterm<<endl;
    comp thirdtermsq = (betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs_comp(x,s,q,m);
    comp thirdterm;
    if(imag(thirdtermsq)>0.0)
    {
        if(real(thirdtermsq)<=0.0)
        {
            thirdterm = sqrt((betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs_comp(x,s,q,m));
        }
        else 
        {
            thirdterm = sqrt((betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs_comp(x,s,q,m));
        }
    }
    else 
    {
        if(real(thirdtermsq)<=0.0)
        {
            thirdterm = sqrt((betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs_comp(x,s,q,m));
        }
        else 
        {
            thirdterm = sqrt((betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs_comp(x,s,q,m));
        }
    }
    //cout<<"p_thirdtermsq:"<<(betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs_comp(x,s,q,m)<<endl;
    //cout<<"thirdterm"<<thirdterm<<endl;
    //cout<<"3.1part:"<<(betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps)<<endl;
    //cout<<"3.2part:"<<4.0*m*m*betaxqs_comp(x,s,q,m)<<endl;
    comp forthterm = 2.0*betaxqs_comp(x,s,q,m);
    //cout<<"forthterm:"<<forthterm<<endl;
    //cout<<"res:"<<(firstterm + secondterm*thirdterm)/forthterm<<endl;
    return (firstterm + secondterm*thirdterm)/forthterm;
}

comp pcut_minus(     double x,
                    comp s,
                    comp q,
                    double m,
                    double eps  )
{
    comp ii = {0.0,1.0};
    comp firstterm = -q*x*betaxqs(1.0,s,q,m);
    comp secondterm = sqbeta0qs(s,q,m);
    comp thirdtermsq = (betaxqs(1.0,s,q,m)+ii*eps)*(betaxqs(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs(x,s,q,m);
    comp thirdterm;
    
    if(real(thirdtermsq)>=0.0 && imag(thirdtermsq)>=0.0)
    {
        thirdterm = mysqrt(abs(real(thirdtermsq)) + ii*abs(imag(thirdtermsq)));
    }
    else if(real(thirdtermsq)>=0.0 && imag(thirdtermsq)<0.0)
    {
        thirdterm = sqrt(abs(real(thirdtermsq)) - ii*abs(imag(thirdtermsq)));
    }
    else if(real(thirdtermsq)<0.0 && imag(thirdtermsq)>=0.0)
    {
        thirdterm = ii*sqrt(abs(real(thirdtermsq)) - ii*abs(imag(thirdtermsq)));
    }
    else if(real(thirdtermsq)<0.0 && imag(thirdtermsq)<0.0)
    {
        thirdterm = ii*sqrt(abs(real(thirdtermsq)) + ii*abs(imag(thirdtermsq)));
    }
    comp forthterm = 2.0*betaxqs(x,s,q,m);

    return (firstterm - secondterm*thirdterm)/forthterm;
}

comp pcut_plus_withprinter(     double x,
                    comp s,
                    comp q,
                    double m,
                    double eps,
                    int n  )
{
    comp ii = {0.0,1.0};
    comp firstterm = -q*x*betaxqs(1.0,s,q,m);
    comp secondterm = sqbeta0qs(s,q,m);
    comp thirdtermsq = (betaxqs(1.0,s,q,m)+ii*eps)*(betaxqs(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs(x,s,q,m);
    comp thirdterm;
    /*if(imag(thirdtermsq)>0.0)
    {
        if(real(thirdtermsq)<=0.0)
        {
            thirdterm = sqrt((betaxqs(1.0,s,q,m)+ii*eps)*(betaxqs(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs(x,s,q,m));
        }
        else 
        {
            thirdterm = sqrt((betaxqs(1.0,s,q,m)+ii*eps)*(betaxqs(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs(x,s,q,m));
        }
    }
    else 
    {
        if(real(thirdtermsq)<=0.0)
        {
            thirdterm = mysqrt((betaxqs(1.0,s,q,m)+ii*eps)*(betaxqs(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs(x,s,q,m));
        }
        else 
        {
            thirdterm = sqrt((betaxqs(1.0,s,q,m)+ii*eps)*(betaxqs(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs(x,s,q,m));
        }
    }*/
    if(real(thirdtermsq)>=0.0 && imag(thirdtermsq)>=0.0)
    {
        thirdterm = mysqrt(abs(real(thirdtermsq)) + ii*abs(imag(thirdtermsq)));
    }
    else if(real(thirdtermsq)>=0.0 && imag(thirdtermsq)<0.0)
    {
        thirdterm = sqrt(abs(real(thirdtermsq)) - ii*abs(imag(thirdtermsq)));
    }
    else if(real(thirdtermsq)<0.0 && imag(thirdtermsq)>=0.0)
    {
        thirdterm = ii*sqrt(abs(real(thirdtermsq)) - ii*abs(imag(thirdtermsq)));
    }
    else if(real(thirdtermsq)<0.0 && imag(thirdtermsq)<0.0)
    {
        thirdterm = ii*sqrt(abs(real(thirdtermsq)) + ii*abs(imag(thirdtermsq)));
    }
    comp forthterm = 2.0*betaxqs(x,s,q,m);

    return (firstterm + secondterm*thirdterm)/forthterm;
}


comp pcut_minus_withprinter(     double x,
                    comp s,
                    comp q,
                    double m,
                    double eps,
                    int n  )
{
    comp ii = {0.0,1.0};
    comp firstterm = -q*x*betaxqs(1.0,s,q,m);
    comp secondterm = sqbeta0qs(s,q,m);
    comp thirdtermsq = (betaxqs(1.0,s,q,m)+ii*eps)*(betaxqs(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs(x,s,q,m);
    comp thirdterm;
    /*if(imag(thirdtermsq)>0.0)
    {
        if(real(thirdtermsq)<=0.0)
        {
            thirdterm = sqrt((betaxqs(1.0,s,q,m)+ii*eps)*(betaxqs(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs(x,s,q,m));
        }
        else 
        {
            thirdterm = sqrt((betaxqs(1.0,s,q,m)+ii*eps)*(betaxqs(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs(x,s,q,m));
        }
    }
    else 
    {
        if(real(thirdtermsq)<=0.0)
        {
            thirdterm = mysqrt((betaxqs(1.0,s,q,m)+ii*eps)*(betaxqs(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs(x,s,q,m));
        }
        else 
        {
            thirdterm = sqrt((betaxqs(1.0,s,q,m)+ii*eps)*(betaxqs(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs(x,s,q,m));
        }
    }*/
    if(n==0)
    {
        cout<<"m_thirdtermsq:"<<(betaxqs(1.0,s,q,m)+ii*eps)*(betaxqs(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs(x,s,q,m)<<endl;
    }
    
    if(real(thirdtermsq)>=0.0 && imag(thirdtermsq)>=0.0)
    {
        thirdterm = mysqrt(abs(real(thirdtermsq)) + ii*abs(imag(thirdtermsq)));
    }
    else if(real(thirdtermsq)>=0.0 && imag(thirdtermsq)<0.0)
    {
        thirdterm = sqrt(abs(real(thirdtermsq)) - ii*abs(imag(thirdtermsq)));
    }
    else if(real(thirdtermsq)<0.0 && imag(thirdtermsq)>=0.0)
    {
        thirdterm = ii*sqrt(abs(real(thirdtermsq)) - ii*abs(imag(thirdtermsq)));
    }
    else if(real(thirdtermsq)<0.0 && imag(thirdtermsq)<0.0)
    {
        thirdterm = ii*sqrt(abs(real(thirdtermsq)) + ii*abs(imag(thirdtermsq)));
    }
    comp forthterm = 2.0*betaxqs(x,s,q,m);

    return (firstterm - secondterm*thirdterm)/forthterm;
}

comp pcut_minus_comp(     comp x,
                    comp s,
                    comp q,
                    double m,
                    double eps  )
{
    comp ii = {0.0,1.0};
    if(imag(x)==0.0)
    {
        x = real(x) + ii*1.0e-15;
    }
    comp firstterm = -q*x*betaxqs_comp(1.0,s,q,m);
    //cout<<"firstterm:"<<firstterm<<endl;
    comp secondterm = sqbeta0qs(s,q,m);
    comp thirdtermsq = (betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs_comp(x,s,q,m);
    comp thirdterm;
    if(imag(thirdtermsq)>0.0)
    {
        if(real(thirdtermsq)<=0.0)
        {
            thirdterm = sqrt((betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs_comp(x,s,q,m));
        }
        else 
        {
            thirdterm = sqrt((betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs_comp(x,s,q,m));
        }
    }
    else 
    {
        if(real(thirdtermsq)<=0.0)
        {
            thirdterm = sqrt((betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs_comp(x,s,q,m));
        }
        else 
        {
            thirdterm = sqrt((betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs_comp(x,s,q,m));
        }
    }
    //cout<<"m_thirdtermsq:"<<(betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs_comp(x,s,q,m)<<endl;
    //cout<<"thirdterm"<<thirdterm<<endl;
    //cout<<"3.1part:"<<(betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps)<<endl;
    //cout<<"3.2part:"<<4.0*m*m*betaxqs_comp(x,s,q,m)<<endl;
    comp forthterm = 2.0*betaxqs_comp(x,s,q,m);
    //cout<<"forthterm:"<<forthterm<<endl;

    //cout<<"res:"<<(firstterm - secondterm*thirdterm)/forthterm<<endl;    
    return (firstterm - secondterm*thirdterm)/forthterm;
}

comp x_for_qc1_plus(   comp s,
                       comp q,
                       double a,
                       double m,
                       double eps    )
{
    comp ii = {0.0,1.0};
    comp qc1 = q_c1(s,a,m,eps);
    comp firstterm = q*qc1*betaxqs(1.0,s,q,m) + ii*eps*q*qc1;
    comp secondterm = betaxqs(0.0,s,q,m)*qc1*qc1*q*q*m*m + betaxqs(0.0,s,q,m)*q*q*qc1*qc1*qc1*qc1;
    comp thirdterm = 2.0*q*q*qc1*qc1;

    return (firstterm + 2.0*sqrt(secondterm))/thirdterm;
}

comp x_for_qc1_minus(   comp s,
                       comp q,
                       double a,
                       double m,
                       double eps    )
{
    comp ii = {0.0,1.0};
    comp qc1 = q_c1(s,a,m,eps);
    comp firstterm = q*qc1*betaxqs(1.0,s,q,m) + ii*eps*q*qc1;
    comp secondterm = betaxqs(0.0,s,q,m)*qc1*qc1*q*q*m*m + betaxqs(0.0,s,q,m)*q*q*qc1*qc1*qc1*qc1;
    comp thirdterm = 2.0*q*q*qc1*qc1;

    return (firstterm - 2.0*sqrt(secondterm))/thirdterm;
}

comp x_for_pk(  comp s,
                comp p,
                comp k,
                double m,
                double eps  )
{
    comp ii = {0.0,1.0};
    comp A = pow(sqrt(s) - omega_comp(p,m) - omega_comp(k,m),2.0);

    return (A - p*p - k*k - m*m + ii*eps)/(2.0*p*k);
}


//--------------------------------------------//


comp radius(    comp s,
                comp sigk,
                double m    )
{
    return (sigc2(s,sigk,m) - sigc1(s,sigk,m))/2.0;
}

comp s3lim( comp sigk,
            double m    )
{
    return pow(m*m - sigk,2.0)/(m*m);
}

comp gfuncConst(	double a,
					double r,
					double m	)
{
    double pi = acos(-1.0);
	double k2k = 1.0/a + r/(2.0*a*a);
	double sb = 4.0*(-k2k*k2k + m*m);
	double A = 128.0 * pi;
	comp B = (comp) sb;
	comp C1 = sqrt(B);
	double D = k2k;
	double E = r*k2k;
	double F = 1.0 - E;
	double G = A*D/F;
	comp H = (comp) G;
	comp I = H*C1;
	comp J = sqrt(I);

	return J;
}

comp rhophib(   comp s,
                double a,
                double m    )
{
    double pi = acos(-1.0);
    comp sb = sigmab(a,m);
    comp q = pmom(s,sb,m);

    return q/(8.0*pi*sqrt(s));
}



void list_of_points()
{
    double a = 16.0;
    double m = 1.0;


    for(double s3=1.0;s3<=9.0;s3=s3+0.01)
    {
        comp sigb = sigmab(a,m);
        cout<<"s:"<<s3<<'\n'<<"smax:"<<sigmax(s3,m)<<'\n'
            <<"sigb:"<<sigb<<'\n'
            <<"sigplus:"<<sigpplus(s3,sigb,m)<<'\n'
            <<"sigmin:"<<sigpminus(s3,sigb,m)<<'\n'
            <<"sigc1:"<<sigc1(s3,sigb,m)<<'\n'
            <<"sigc2:"<<sigc2(s3,sigb,m)<<'\n'
            <<"rad:"<<radius(s3,sigb,m)<<'\n';
        cout<<"slim:"<<real(s3lim(sigb,m))<<endl;
        cout<<"----------------------------------"<<endl;
        
    }
}

comp phibthreshold( double a,
                    double m    )
{
    comp sigb = sigmab(a,m);

    return pow(sqrt(sigb) + m,2.0);
}

comp Gleftbranchpoint(  double a,
                        double m    )
{
    comp sigb = sigmab(a,m);

    return pow(m*m - sigb,2.0)/(m*m);
}

comp Grightbranchpoint( double a,
                        double m    )
{
    comp sigb = sigmab(a,m);

    return m*m + 2.0*sigb;
}

/*comp eps_energy_factor_plus(    comp s,
                                comp sigk,
                                double m    )
{
    return (1.0/(2.0*sigk) + (1.0/(2.0*sqrt(sigk)))*sqrt(kallentriangle(s,sigk,m*m))/(sqrt(sigk - 4.0*m*m)));
}

comp eps_energy_factor_minus(    comp s,
                                comp sigk,
                                double m    )
{
    return (1.0/(2.0*sigk) - (1.0/(2.0*sqrt(sigk)))*sqrt(kallentriangle(s,sigk,m*m))/(sqrt(sigk - 4.0*m*m)));
}*/

comp eps_abovethreshold_energy_factor(  comp s,
                                        double a,
                                        double m    )
{
    return (s + m*m - sigmab(a,m))/(4.0*pmom(s,sigmab(a,m),m)*s);
}

double eps_above_threshold( double eta,
                            comp s,
                            double a,
                            double m,
                            double N    )
{
    double sigmx =real(sigmax(s,m));
    double pi = acos(-1.0);
    return real(eta/(2.0*pi*N/sigmx*eps_abovethreshold_energy_factor(s,a,m)));
}

double eps_momrep_above_threshold(  double eta,
                                    comp s,
                                    double a,
                                    double m,
                                    double N    )
{
    double pi = acos(-1.0);
    comp kmax = pmom(s,0.0,m);
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb,m);

    comp energy_dependent_term = (s + m*m - sigb)/(4.0*q*s);

    return real(eta*kmax/(2.0*pi*energy_dependent_term*N));

}

/*void opeplotcontour()
{
    ofstream fout;

    double a = 16.0;
    double m = 1.0;
    double s = 8.65;
    double eps = 1.0e-3;
    comp scomp = (comp) s;
    comp sigmapprime = sigmab(a,m);
    double sigpRinitial = 3.55;
    double sigpRfinal = 4.0;
    double delsigpR = abs(sigpRinitial-sigpRfinal)/1000.0;
    double sigpIinitial = -0.1;
    double sigpIfinal = 0.1;
    double delsigpI = abs(sigpIinitial - sigpIfinal)/1000.0;
    

     

    string filename = "OPEcontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename.c_str());
    for(double sigmapR=sigpRinitial;sigmapR<=sigpRfinal;sigmapR=sigmapR + delsigpR)
    {
        
        for(double sigmapI=sigpIinitial;sigmapI<=sigpIfinal;sigmapI=sigmapI + delsigpI)
        //double someeps = 1.0e-3;
        {
            comp sigmapcomp = (comp) sigmapR + ii*sigmapI;

            comp ope = GS(scomp,sigmapcomp,sigmapprime,m,eps);

            fout<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
            cout<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
        }
        fout<<endl;
    }   

    fout.close();

    string filename1 = "threshOPEcontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename1.c_str());

    fout<<real(sigpplus(s,sigmab(a,m),m))<<'\t'<<real(sigpminus(s,sigmab(a,m),m))<<'\t'<<real(sigc1(s,sigmab(a,m),m))<<'\t'<<real(sigc2(s,sigmab(a,m),m))<<'\t'<<pow(sqrt(s)-m,2.0)<<'\t'<<0<<'\t'<<0<<endl;

    fout.close();
}
*/

/*void kernelplotcontour()
{
    ofstream fout;

    double a = 16.0;
    double m = 1.0;
    double s = 8.65;
    double eps = 1.0e-3;
    comp scomp = (comp) s;
    comp sigmapprime = sigmab(a,m);
    double sigpRinitial = 3.55;
    double sigpRfinal = 5.0;
    double delsigpR = abs(sigpRinitial-sigpRfinal)/2000.0;
    double sigpIinitial = -0.1;
    double sigpIfinal = 0.1;
    double delsigpI = abs(sigpIinitial - sigpIfinal)/2000.0;
    

     

    string filename = "kernelcontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + "_1.dat";

    fout.open(filename.c_str());
    for(double sigmapR=sigpRinitial;sigmapR<=sigpRfinal;sigmapR=sigmapR + delsigpR)
    {
        
        for(double sigmapI=sigpIinitial;sigmapI<=sigpIfinal;sigmapI=sigmapI + delsigpI)
        //double someeps = 1.0e-3;
        {
            comp sigmapcomp = (comp) sigmapR + ii*sigmapI;

            //comp ope = kernel(scomp,sigmapcomp,sigmapprime,a,m,eps);
            comp ope = kernel(scomp,sigmapprime,sigmapcomp,a,m,eps);

            fout<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
            cout<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
        }
        fout<<endl;
    }   

    fout.close();

    string filename1 = "threshkernelcontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename1.c_str());

    fout<<real(sigpplus(s,sigmab(a,m),m))<<'\t'<<real(sigpminus(s,sigmab(a,m),m))<<'\t'<<real(sigc1(s,sigmab(a,m),m))<<'\t'<<real(sigc2(s,sigmab(a,m),m))<<'\t'<<pow(sqrt(s)-m,2.0)<<'\t'<<0<<'\t'<<0<<endl;

    fout.close();
}
*/

/*void mommaker()
{
    double s = 8.65;
    double m = 1.0;
    double a = 16.0;
    comp sigb = sigmab(a,m);
    vector<comp> sigvec;
    double smalldel = 1.0e-5;
    double sigpinitial = 0.0;
    double sigpfinal = pow(sqrt(s) - m,2.0);
    double delsigp = abs(sigpfinal - sigpinitial)/1000.0;

    comp sigpp = sigpplus(s,sigb,m);
    comp sigpm = sigpminus(s,sigb,m);
    comp sigc11 = sigc1(s,sigb,m);
    comp sigc22 = sigc2(s,sigb,m);
    ofstream fout;
    string filename = "sigvec.dat";

    fout.open(filename.c_str());
    for(double sigp=sigpinitial;sigp<=real(sigc11)-smalldel;sigp=sigp+delsigp)
    {
        fout<<sigp<<'\t'<<0<<'\t'<<0<<endl;
    }

    fout<<endl;

    for(double sigp=0;sigp<=real(radius(s,sigb,m));sigp=sigp+delsigp)
    {

        fout<<real(sigc11)-smalldel<<'\t'<<sigp<<'\t'<<0<<endl;
    }

    fout<<endl;

    for(double sigp=0;sigp<=real(sigc22+smalldel)-real(sigc11-smalldel);sigp=sigp+delsigp)
    {

        fout<<real(sigc11)-smalldel + sigp<<'\t'<<real(radius(s,sigb,m))<<'\t'<<0<<endl;
        fout<<endl;
    }

    fout<<endl;

    for(double sigp=0;sigp<=real(radius(s,sigb,m));sigp=sigp+delsigp)
    {

        fout<<real(sigc22+smalldel)<<'\t'<<real(radius(s,sigb,m)- sigp)<<'\t'<<0<<endl;
    }

    fout<<endl;

    for(double sigp=0;sigp<=real(sigc22-sigmax(s,m));sigp=sigp+delsigp)
    {

        fout<<real(sigc22+smalldel - sigp)<<'\t'<<0<<'\t'<<0<<endl;
        fout<<endl;
    }

    fout.close();

}
*/

/*void opeplot()
{
    ofstream fout;

    double a = 16.0;
    double m = 1.0;
    double s = 8.3;
    double eps = 1.0e-3;
    comp scomp = (comp) s;
    comp sigmapprime = sigmab(a,m);
     

    string filename = "OPE_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename.c_str());
    for(double sigmap=3.3;sigmap<=3.999;sigmap=sigmap + abs(3.3-3.999)/1000.0)
    {
        double someeps = 1.0e-7;
        comp sigmapcomp = (comp) sigmap + ii*someeps;

        comp ope = GS(scomp,sigmapcomp,sigmapprime,m,eps);

        fout<<sigmap<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
    }

    fout.close();

    
}
*/



#endif 