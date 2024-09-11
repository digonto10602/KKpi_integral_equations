#include<bits/stdc++.h>
#include<omp.h>
#include "fastGL_based_functions.h"
//#include "ope_functions.h"
using namespace std;

typedef complex<double> comp;

comp ii = {0.0,1.0};

comp mysqrt(    comp x  )
{
    return ii*sqrt(-x);
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

comp zppprime_1( comp s,
                comp sigmap,
                comp sigmapprime,
                double m,
                double epsilon  )
{
    comp ii = {0.0,1.0};
    comp num = (sigmapprime - m*m)*(sigmap + s - m*m) + s*(sigmap - s + m*m + 2.0*ii*epsilon);
    comp denom = mysqrt(kallentriangle(s,sigmapprime,m*m))*mysqrt(kallentriangle(s,sigmap,m*m));

    return num/denom; 
}

comp zppprime_complex_s(    comp s,
                            comp sigmap,
                            comp sigmapprime,
                            double m,
                            double eps1,
                            double eps2 )
{
    comp ii = {0.0,1.0};
    comp zpk = zppprime_1(s,sigmap,sigmapprime,m,eps1);
    comp eps_energyfactor_1 = (sigmap + sigmapprime - 2.0*s + 2.0*ii*eps1)/(mysqrt(kallentriangle(s,sigmapprime,m*m))*mysqrt(kallentriangle(s,sigmap,m*m)));
    comp eps_energyfactor_2 = zpk*((s-sigmap-m*m)/kallentriangle(s,sigmap,m*m) + (s-sigmapprime-m*m)/kallentriangle(s,sigmapprime,m*m));
    comp eps_factor = ii*eps2*(eps_energyfactor_1 - eps_energyfactor_2);

    return zpk + eps_factor;
}

comp zpk_energyfactor(  comp s,
                        comp sigmap,
                        comp sigmapprime,
                        double m,
                        double eps1 )
{
    comp ii = {0.0,1.0};
    comp zpk = zppprime_1(s,sigmap,sigmapprime,m,eps1);
    comp eps_energyfactor_1 = (sigmap + sigmapprime - 2.0*s + 2.0*ii*eps1)/(sqrt(kallentriangle(s,sigmapprime,m*m))*sqrt(kallentriangle(s,sigmap,m*m)));
    comp eps_energyfactor_2 = zpk*((s-sigmap-m*m)/kallentriangle(s,sigmap,m*m) + (s-sigmapprime-m*m)/kallentriangle(s,sigmapprime,m*m));
    comp eps_factor = (eps_energyfactor_1 - eps_energyfactor_2);

    return eps_factor;
}

double Jfunc(   double x    )
{
    if(x<=0) return 0.0;
    else if(x>0 && x<1 ) return exp((-1.0/x)*exp(-1.0/(1.0-x)));
    else return 1.0;
}

double Hfunc(   comp sigmap,
                comp sigmapprime,
                double m    )
{
    return Jfunc(real(sigmap/(4.0*m*m)))*Jfunc(real(sigmapprime/(4.0*m*m)));
}

comp Jfunc_comp( comp x  )
{
    if(real(x)<=0.0) return 0.0;
    else if(real(x)>0.0 && real(x)<1.0) return exp((-1.0/x)*exp(-1.0/(1.0-x)));
    else return 1.0;
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
    comp p = mysqrt(kallentriangle(s,sigmap,m*m))/(2.0*sqrt(s));
    comp pprime = mysqrt(kallentriangle(s,sigmapprime,m*m))/(2.0*sqrt(s));

    comp zpp = zppprime(s,sigmap,sigmapprime,m,epsilon);

    return -Hfunc(sigmap,sigmapprime,m)/(4.0*p*pprime)*(log( (zpp - 1.0)/(zpp + 1.0) ));
    //return -1.0/(4.0*p*pprime)*(log( (zpp - 1.0)/(zpp + 1.0) ));

}

comp GS_1(  comp s,
            comp sigmap,
            comp sigmapprime,
            double m,
            double epsilon      )
{
    comp p = mysqrt(kallentriangle(s,sigmap,m*m))/(2.0*sqrt(s));
    comp pprime = mysqrt(kallentriangle(s,sigmapprime,m*m))/(2.0*sqrt(s));

    comp zpp = zppprime_1(s,sigmap,sigmapprime,m,epsilon);

    return -Hfunc(sigmap,sigmapprime,m)/(4.0*p*pprime)*(log( (zpp - 1.0)/(zpp + 1.0) ));
    //return -1.0/(4.0*p*pprime)*(log( (zpp - 1.0)/(zpp + 1.0) ));

}


// ------------------Momentum Representation--------------------- //
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

    //smooth cut off
    //return -Hfunc_comp(sigp,sigk,m)/(4.0*p*k)*(log(num/denom));
    
    //hard cut off
    return -1.0/(4.0*p*k)*(log(num/denom));
    
    //above threshold
    //return -Hfunc(sigp,sigk,m)/(4.0*p*k)*(log(num/denom));
    
}

// ----------------------------------------------------------------- //

comp delta_epsilon( comp sigq,
                    comp sigpoint,
                    double eps      )
{
    double pi = acos(-1.0);

    if(real(sigq)-real(sigpoint)>1.0e-10)
    return eps/(pi*((sigq-sigpoint)*(sigq-sigpoint) + eps*eps));
    else
    return 1.0/(pi*eps);
}

comp deltaGS_1( comp s,
                comp sigp,
                comp sigk,
                comp sigpoint,
                double m,
                double eps  )
{
    double pi = acos(-1.0);
    double sigkreal = real(sigk);
    double sigpointreal = real(sigpoint);
    /*if(abs(sigkreal-sigpointreal)<1.0e-4)
    {
    cout<<"sigk:"<<sigk<<'\t'<<"sigpoint:"<<sigpoint<<endl;
    cout<<"gs1sigpsigk:"<<GS_1(s,sigp,sigk,m,eps)<<endl;
    cout<<"gs1sigpsigpoint:"<<GS_1(s,sigp,sigpoint,m,eps)<<endl;
    }*/
    cout<<exp(-abs(sigk - sigpoint))<<endl;
    return GS_1(s,sigp,sigk,m,eps) - GS_1(s,sigp,sigpoint,m,eps)*exp(-abs((sigk - sigpoint)*10.0));//*eps*pi*delta_epsilon(sigk,sigpoint,eps);
}

comp onebypk_complex_s( double sreal,
                        comp sigp,
                        comp sigk,
                        double m,
                        double eps2)
{
    comp ii = {0.0,1.0};
    comp p_with_sreal = mysqrt(kallentriangle(sreal,sigp,m*m))/(2.0*sqrt((comp) sreal));
    comp k_with_sreal = sqrt(kallentriangle(sreal,sigk,m*m))/(2.0*sqrt((comp) sreal));
    comp onebypk = 1.0/(p_with_sreal*k_with_sreal);
    comp energy_factor_1 = sreal*(m*m-sreal+sigp)/(sqrt(kallentriangle(sreal,sigk,m*m))*sqrt(kallentriangle(sreal,sigp,m*m))*kallentriangle(sreal,sigp,m*m));
    comp energy_factor_2 = (pow(sigk-m*m,2.0)-sreal*(sigk+m*m))/(sqrt(kallentriangle(sreal,sigp,m*m))*sqrt(kallentriangle(sreal,sigk,m*m))*kallentriangle(sreal,sigk,m*m));

    return onebypk + 4.0*ii*eps2*(energy_factor_1+energy_factor_2);
}

comp GS_complex_s(  comp s,
                    comp sigmap,
                    comp sigmapprime,
                    double m,
                    double eps1,
                    double eps2      )
{
    comp p = mysqrt(kallentriangle(s,sigmap,m*m))/(2.0*sqrt(s));
    comp pprime = mysqrt(kallentriangle(s,sigmapprime,m*m))/(2.0*sqrt(s));

    comp zpp = zppprime_complex_s(s,sigmap,sigmapprime,m,eps1,eps2);

    return -(Hfunc(sigmap,sigmapprime,m)/(4.0))*onebypk_complex_s(real(s),sigmap,sigmapprime,m,eps2)*(log( (zpp - 1.0)/(zpp + 1.0) ));
    //return -1.0/(4.0*p*pprime)*(log( (zpp - 1.0)/(zpp + 1.0) ));

}



comp sigmab(    double a,
                double m    )
{
    return 4.0*(m*m - 1.0/(a*a));
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

comp tau2(  comp s,
            comp sigk,
            double m    )
{
    double pi = acos(-1.0);
    return (sqrt(-sigk + (sqrt(s) - m)*(sqrt(s) - m))*sqrt(-sigk + (mysqrt(s) + m)*(mysqrt(s) + m)))/(8.0*pi*s);
}


comp pmom(  comp s,
            comp sigk,
            double m    )
{
    return sqrt(kallentriangle(s,sigk,m*m))/(2.0*sqrt(s));
}

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

    return num/denom;
}



comp kernel(    comp s,
                comp sigp,
                comp sigk,
                double a,
                double m,
                double epsilon )
{
    double pi = acos(-1.0);
    //cout<<"Gs:"<<GS(s,sigp,sigk,m,epsilon)<<endl;
    //cout<<"tau:"<<tau(s,sigk,m)<<endl;
    //cout<<"M2:"<<M2kfunc(a,sigk,m,epsilon)<<endl;
    //return (1.0/(2.0*pi))*GS(s,sigp,sigk,m,epsilon)*tau(s,sigk,m)*M2kfunc(a,sigk,m,epsilon);
    return GS(s,sigp,sigk,m,epsilon);
}




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

comp sigpplus_witheps(  comp s,
                        comp sigk,
                        double m,
                        double eps  )
{
    comp ii = {0.0,1.0};
    //cout<<"s:"<<s<<'\t'<<"sigp:"<<sigpplus(s,sigk,m)<<'\t'<<"en term:"<<eps*eps_energy_factor_minus(s,sigk,m)<<endl;
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

void sigma_vector_maker_linear( vector<comp> &sigvec,
                                comp min,
                                comp max,
                                double points,
                                double eps,
                                int ud   )
{
    double delsig = abs(max - min)/points;
    comp ii = {0.0,1.0};
    if(ud==1)
    {
        for(int i=1;i<=(int)points;++i)
        {
            comp sigma = min + (comp)i*delsig + ii*eps ;
            //if(sigma<=max)
            //{
                sigvec.push_back(sigma);
            //}
            //else break;
        } 
    }
    else if(ud==2)
    {
        for(int i=1;i<=(int)points;++i)
        {
            comp sigma = min + (comp)i*ii*delsig + ii*eps;
            //if(sigma<=max)
            //{
                sigvec.push_back(sigma);
            //}
            //else break;
        }

    }
    else if(ud==3)
    {
        for(int i=1;i<=(int)points;++i)
        {
            comp sigma = min - (comp)i*ii*delsig + ii*eps;
            //if(sigma<=max)
            //{
                sigvec.push_back(sigma);
            //}
            //else break;
        }

    }
    else if(ud==4)
    {
        for(int i=1;i<=(int)points;++i)
        {
            //cout<<"im here"<<endl;
            comp sigma = min - (comp)i*delsig + ii*eps;
            //if(sigma<=max)
            //{
                sigvec.push_back(sigma);
            //}
            //else break;
        }

    }
}

void sigma_vector_maker_2(  vector<comp> &sigvec,
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double point1,
                            double point2,
                            double eps   )
{
    comp ii = {0.0,1.0};
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);

    if(real(rad)>0.0)
    {
    

    cout<<"sigmin:"<<min<<endl;
    cout<<"sigmax:"<<max<<endl;
    cout<<"rad:"<<rad<<endl;
    cout<<"sigmac1:"<<sigmac1<<endl;
    cout<<"sigmac2:"<<sigmac2<<endl;

    sigma_vector_maker_linear(sigvec,min,sigmac1-100*eps,point1,eps,1); // the last 1 is for making vec on real axis
    sigma_vector_maker_linear(sigvec,sigmac1-100*eps,sigmac1+rad,point2,eps,2);
    sigma_vector_maker_linear(sigvec,sigmac1-100*eps+ii*(rad+100*eps),sigmac1+2.0*(rad+100*eps)+ii*(rad+100*eps),2.0*point2,eps,1);
    sigma_vector_maker_linear(sigvec,sigmac1+2.0*(rad+100*eps)+ii*(rad+100*eps),sigmac2 + 2.0*(100*eps) ,point2,eps,3);
    sigma_vector_maker_linear(sigvec,sigmac2+100*eps,max,point2,eps,4);

    cout<<"sigvec created with size = "<<sigvec.size()<<endl;
    }
    else
    {
        double length = real(max - min);
        comp delsig = length/(point1 + 5.0*point2) ;

        for(double tempsig=real(min);tempsig<=real(max);tempsig=tempsig+real(delsig))
        {
            comp sigma = tempsig + ii*eps;

            sigvec.push_back(sigma);
        }
    }
    

}


void sigma_vector_maker(    vector<comp> &sigvec,
                            comp s,
                            comp sigmin,
                            comp sigmax,
                            double a,
                            double m,
                            double points,
                            double eps   )
{
    comp ii = {0.0,1.0};
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);

    if(real(rad)>0.0)
    {
        //double eps = 1.0e-5;
        double length = + real(sigmac1 - eps) 
                        - real(sigmin)
                        + real(sigmac2 + eps)
                        - real(sigmax)
                        + 4.0*real(rad + eps);
        comp delsig = length/points;

        for(double tempsig=real(sigmin);tempsig<real(sigmac1 - eps);tempsig=tempsig+real(delsig))
        {
            comp sigma = tempsig + ii*eps;

            sigvec.push_back(sigma);
        }
        for(double tempsig=0.0;tempsig<real(rad+eps);tempsig=tempsig+real(delsig))
        {
            comp sigma = sigmac1 - eps + ii*tempsig + ii*eps;

            sigvec.push_back(sigma);
        }
        for(double tempsig=0.0;tempsig<2.0*real(rad+eps);tempsig=tempsig+real(delsig))
        {
            comp sigma = sigmac1 - eps + ii*(rad+eps) + tempsig + ii*eps;

            sigvec.push_back(sigma);
        }
        for(double tempsig=0.0;tempsig<real(rad+eps);tempsig=tempsig+real(delsig))
        {
            comp sigma = sigmac1 - eps + ii*(rad+eps) - ii*tempsig + 2.0*(rad+eps) + ii*eps;

            sigvec.push_back(sigma);
        }
        for(double tempsig=real(sigmac2 + eps);tempsig>=real(sigmax);tempsig=tempsig-real(delsig))
        {
            comp sigma = tempsig + ii*eps;

            sigvec.push_back(sigma);
        }
    }
    else
    {
        //double eps = 1.0e-5;
        double length = real(sigmax - sigmin);
        comp delsig = length/points;

        for(double tempsig=real(sigmin);tempsig<=real(sigmax);tempsig=tempsig+real(delsig))
        {
            comp sigma = tempsig + ii*eps;

            sigvec.push_back(sigma);
        }
    }

    cout<<"sigvec created with size = "<<sigvec.size()<<endl;
}

//We will make a linear sigma vector. It will have a    //
//starting point, we will provide with the delsig and   //
//number of points. The delsig will have a plus or minus//
//sign with it which tells it to go up or down and left //
//right to form the contour in the complex plane. It    //
//will also have a int,'1' sets it to form the contour  //
//along the real axis and '2' sets it to form the       //
//contour in the imaginary axis.                        //

void sigma_vector_maker_linear_1(   vector<comp> &sigvec,
                                    comp startingpoint,
                                    comp delsig,
                                    double points,
                                    int a               )
{
    comp ii = {0.0,1.0};
    if(a==1)
    {
        for(int i=1;i<=(int)points;++i)
        {
            comp tempsig = startingpoint + ((comp)i)*delsig;
            sigvec.push_back(tempsig);
        }
    }
    else if(a==2)
    {
        for(int i=1;i<=(int)points;++i)
        {
            comp tempsig = startingpoint + ii*((comp)i)*delsig;
            sigvec.push_back(tempsig);
        }
    }

}

//this one makes the contour for each s, we give it an  //
//epsratio, this determines how large the box is around //
//circular cut. The created sigvec do exclude the       //
//minimum value                                         //
void sigma_vector_maker_5(  vector<comp> &sigvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double points1,
                            double points2,
                            double epsratio     )
{
    comp ii = {0.0,1.0};
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp delsig = {0.0,0.0};

    cout<<"sigc1:"<<sigmac1<<endl;
    cout<<"sigc2:"<<sigmac2<<endl;
    cout<<"epsratio:"<<epsratio<<endl;

    if(real(rad)>0.0)
    {
        /*  this is the initial straight line along the real axis */
        comp startingpoint = min + ii*eps;
        delsig = abs((sigmac1-epsratio*eps)-min)/points1;
        cout<<"first c1:"<<sigmac1-epsratio*eps<<endl;
        cout<<"min:"<<min<<endl;
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,points1,1);
        
        /*  this is the first contour up in the imaginary axis   */
        startingpoint = sigvec[sigvec.size()-1];
        delsig = abs(rad+((comp)epsratio*eps))/points2;
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,points2,2);

        /*  this is again parallel to the real axis, to the right, this 
            goes for 2*radius and for 2*points2                         */
        startingpoint = sigvec[sigvec.size()-1];
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,2.0*points2,1);

        /*  the contour now comes down, so we flip the sign of the delsig */
        startingpoint = sigvec[sigvec.size()-1];
        sigma_vector_maker_linear_1(sigvec,startingpoint,-delsig,points2,2);

        /*  this is the end part, the contour goes to the left to max */
        startingpoint = sigvec[sigvec.size()-1];
        delsig = abs(startingpoint - max)/points2;
        sigma_vector_maker_linear_1(sigvec,startingpoint,-delsig,points2,1);

        /* print the size of the created sigma_vector */
        cout<<"sigvec created with size = "<<sigvec.size()<<endl;

    }
    else
    {
        double length = real(max - min);
        comp delsig = length/(points1 + 5.0*points2) ;

        int totpoints = points1 + 5*points2;
        for(int i=1;i<=totpoints;++i)
        {
            comp sigma = min + ii*eps + ((comp)i)*delsig;

            sigvec.push_back(sigma);
        }
    }

}


void sigma_vector_maker_6(  vector<comp> &sigvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double points1,
                            double points2,
                            double epsratio     )
{
    comp ii = {0.0,1.0};
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp delsig = {0.0,0.0};

    if(real(rad)>0.0)
    {
        /*  this is the initial straight line along the real axis */
        comp startingpoint = min + ii*eps;
        delsig = abs((sigmac1-epsratio*eps)-min)/points1;
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,points1,1);
        
        /*  this is the first contour up in the imaginary axis   */
        startingpoint = sigvec[sigvec.size()-1];
        delsig = abs(rad+((comp)epsratio*eps))/points2;
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,points2,2);

        /*  this is again parallel to the real axis, to the right, this 
            goes for 2*radius and for 2*points2                         */
        startingpoint = sigvec[sigvec.size()-1];
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,2.0*points2,1);

        /*  the contour now comes down, so we flip the sign of the delsig */
        startingpoint = sigvec[sigvec.size()-1];
        sigma_vector_maker_linear_1(sigvec,startingpoint,-delsig,points2,2);

        /*  this is the end part, the contour goes to the left to max */
        startingpoint = sigvec[sigvec.size()-1];
        delsig = abs(real(startingpoint) - max)/(2.0*points2);
        sigma_vector_maker_linear_1(sigvec,startingpoint,-delsig,2.0*points2,1);

        /* print the size of the created sigma_vector */
        cout<<"sigvec created with size = "<<sigvec.size()<<endl;

    }
    else
    {
        double length = real(max - min);
        comp delsig = length/(points1 + 5.0*points2) ;

        int totpoints = points1 + 5*points2;
        for(int i=1;i<=totpoints;++i)
        {
            comp sigma = min + ii*eps + ((comp)i)*delsig;

            sigvec.push_back(sigma);
        }
    }

}


void sigma_vector_maker_8(  vector<comp> &sigvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double points1,
                            double points2,
                            double epsratio     )
{
    comp ii = {0.0,1.0};
    comp rad = abs(radius(s,sigmab(a,m),m));
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp sigb = sigmab(a,m);
    comp delsig = {0.0,0.0};

    if(real(rad)>0.0)
    {
        /*  this is the initial straight line along the real axis */
        comp startingpoint = min + ii*eps;
        delsig = abs((sigmac1-epsratio*eps)-min)/points1;
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,points1,1);
        
        /*  this is the first contour up in the imaginary axis   */
        startingpoint = sigvec[sigvec.size()-1];
        delsig = abs(rad+((comp)epsratio*eps))/points2;
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,points2,2);

        /*  this is again parallel to the real axis, to the right, this 
            goes for 2*radius and for 2*points2                         */
        startingpoint = sigvec[sigvec.size()-1];
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,2.0*points2,1);

        /*  the contour now comes down, so we flip the sign of the delsig */
        startingpoint = sigvec[sigvec.size()-1];
        delsig = abs(rad + abs(imag(sigpplus_witheps(s,sigb,m,0.0)))+ ((comp)epsratio*eps + abs(imag(max) - imag(sigpplus_witheps(s,sigb,m,0.0)) )/2.0 ) )/points2;
        sigma_vector_maker_linear_1(sigvec,startingpoint,-delsig,points2,2);

        cout<<"imSigP:"<<imag(sigpplus_witheps(s,sigb,m,0.0))<<endl;
        cout<<"imSigMax:"<<imag(max)<<endl;
        cout<<"diff:"<<abs(imag(sigpplus_witheps(s,sigb,m,0.0)) - imag(max))/2.0<<endl;
        cout<<"what it should be:"<<imag(sigpplus_witheps(s,sigb,m,0.0)) - abs(imag(sigpplus_witheps(s,sigb,m,0.0)) - imag(max))/2.0<<endl;
        cout<<"startpoint:"<<sigvec[sigvec.size()-1]<<endl;
        /*  this is the end part, the contour goes to the left to max */
        startingpoint = sigvec[sigvec.size()-1];
        delsig = abs(real(startingpoint) - real(max))/(6.0*points2);
        sigma_vector_maker_linear_1(sigvec,startingpoint,-delsig,6.0*points2,1);

        /* print the size of the created sigma_vector */
        cout<<"sigvec created with size = "<<sigvec.size()<<endl;

    }
    else
    {
        double length = real(max - min);
        comp delsig = length/(points1 + 10.0*points2) ;

        int totpoints = points1 + 10*points2;
        for(int i=1;i<=totpoints;++i)
        {
            comp sigma = min + ii*eps + ((comp)i)*delsig;

            sigvec.push_back(sigma);
        }
    }

}

comp q_plus(    comp s,
                double a,
                double m,
                double eps  )
{
    comp sigb = sigmab(a,m);
    comp sigppls = sigpplus_witheps(s,sigb,m,eps);

    //cout<<"s:"<<s<<'\t'<<"ktriangle:"<<kallentriangle(s,sigppls,m*m)<<'\t'<<"sqktriangle:"<<sqrt(kallentriangle(s,sigppls,m*m))<<'\t'<<mysqrt(kallentriangle(s,sigppls,m*m))<<endl;
    //cout<<"breakdown:"<<s*s+sigppls*sigppls+m*m*m*m<<'\t'<<2.0*(s*sigppls+sigppls*m*m+m*m*s)<<endl;
    return sqrt(kallentriangle(s,sigppls,m*m))/(2.0*sqrt(s));
}

comp q_plus_k(    comp s,
                comp k,
                double m,
                double eps  )
{
    //comp sigb = sigmab(a,m);
    comp sigk = sigma_p(s,k,m);
    comp sigppls = sigpplus_witheps(s,sigk,m,eps);

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

comp q_minus_k(    comp s,
                comp k,
                double m,
                double eps  )
{
    //comp sigb = sigmab(a,m);
    comp sigk = sigma_p(s,k,m);
    comp sigpmns = sigpminus_witheps(s,sigk,m,eps);

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
    comp thirdterm = sqrt(thirdtermsq);
    
    /*
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
    */
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
    comp thirdterm = sqrt(thirdtermsq);
    
    /*if(real(thirdtermsq)>=0.0 && imag(thirdtermsq)>=0.0)
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
    */
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

void pcut_fixer(    vector<comp> &pcutp,
                    vector<comp> &pcutm  )
{
    int size = pcutp.size();
    for(int i=0;i<size;++i)
    {
        if(i==0)
        {
            comp checkval_pcutp = pcutp[i];
            comp checkprevval_pcutp = pcutp[i+1];
            comp checkval_pcutm = pcutm[i];

            double checkdiff_with_pcutp = abs(checkval_pcutp - checkprevval_pcutp);
            double checkdiff_with_pcutm = abs(checkval_pcutm - checkprevval_pcutp);

            if(checkdiff_with_pcutm<checkdiff_with_pcutp)
            {
                comp temp;
                temp = pcutp[i];
                pcutp[i] = pcutm[i];
                pcutm[i] = temp;
            }
            else continue;
        }
        comp checkval_pcutp = pcutp[i];
        comp checkprevval_pcutp = pcutp[i-1];
        comp checkval_pcutm = pcutm[i];

        double checkdiff_with_pcutp = abs(checkval_pcutp - checkprevval_pcutp);
        double checkdiff_with_pcutm = abs(checkval_pcutm - checkprevval_pcutp);

        if(checkdiff_with_pcutm<checkdiff_with_pcutp)
        {
            comp temp;
            temp = pcutp[i];
            pcutp[i] = pcutm[i];
            pcutm[i] = temp;
        }
        else continue;
    }
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

void mom_vector_maker_linear_1(     vector<comp> &qvec,
                                    comp startingpoint,
                                    comp delsig,
                                    double points,
                                    int a               )
{
    comp ii = {0.0,1.0};
    if(a==1)
    {
        for(int i=1;i<(int)points+1;++i)
        {
            comp tempsig = startingpoint + ((comp)i)*delsig;
            qvec.push_back(tempsig);
        }
    }
    else if(a==2)
    {
        for(int i=1;i<(int)points+1;++i)
        {
            comp tempsig = startingpoint + ii*((comp)i)*delsig;
            qvec.push_back(tempsig);
        }
    }

}


void mom_vector_maker_1(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double points1,
                            double multiplier    )
{
    comp ii = {0.0,1.0};
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    comp delq = {0.0,0.0};

    if(real(rad)>0.0)
    {
        /*  this is the initial straight line along the negative imaginary axis */
        comp startingpoint = min;
        delq = abs((abs(imag(qpls)) + abs(diffqplsqmns) + multiplier*eps)-min)/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,-delq,points1,2);
        
        /*  this is the contour going right in the imaginary axis   */
        startingpoint = qvec[qvec.size()-1];
        delq = abs((abs(real(qc1)) + eps) - real(startingpoint))/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,delq,points1,1);

        /*  This contour now goes up in the imaginary axis                       */
        startingpoint = qvec[qvec.size()-1];
        delq = abs(imag(startingpoint))/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,delq,points1,2);

        /*  This contour now goes right in the real axis                       */
        startingpoint = qvec[qvec.size()-1];
        delq = abs(real(max) - real(startingpoint))/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,delq,points1,1);

        
        /* print the size of the created sigma_vector */
        cout<<"sigvec created with size = "<<qvec.size()<<endl;

    }
    else
    {
        double length = real(max - min);
        comp delsig = length/(4.0*points1) ;

        int totpoints = 4.0*points1;
        for(int i=1;i<=totpoints;++i)
        {
            comp qmom = min + ii*eps + ((comp)i)*delsig;

            qvec.push_back(qmom);
        }
    }

}

void xy_coordinates_m(    double x1,
                        double y1,
                        double x2,
                        double y2,
                        double r,
                        double &x4,
                        double &y4   )
{
    double x3,y3;
    x3 = (x1 + x2)/2.0;
    y3 = (y1 + y2)/2.0;

    double m1 = (y2 - y1)/(x2 - x1);
    double m2 = -1.0/m1;

    if(m1==0.0||abs(m1)<1.0e-16)
    {
        m2 = INT_MAX;
    }
    
    
    x4 = x3 - r/sqrt(1.0 + m2*m2);
    y4 = y3 - m2*(x3 - x4);

    //cout<<"for m: m1:"<<m1<<'\t'<<" m2:"<<m2<<endl;
    //cout<<"       x3:"<<x3<<'\t'<<" y3:"<<y3<<endl;

}

void xy_coordinates_p(    double x1,
                        double y1,
                        double x2,
                        double y2,
                        double r,
                        double &x4,
                        double &y4   )
{
    double x3,y3;
    x3 = (x1 + x2)/2.0;
    y3 = (y1 + y2)/2.0;

    double m1 = (y2 - y1)/(x2 - x1);
    double m2 = -1.0/m1;

    if(m1==0.0||abs(m1)<1.0e-8)
    {
        m2 = INT_MAX;
    }
    
    
    x4 = x3 + r/sqrt(1.0 + m2*m2);
    y4 = y3 - m2*(x3 - x4);

    //cout<<"for p: m1:"<<m1<<'\t'<<" m2:"<<m2<<endl;
    //cout<<"       x3:"<<x3<<'\t'<<" y3:"<<y3<<endl;

}

int m1_checker( double x1,
                double y1,
                double x2,
                double y2   )
{
    double m1 = (y2 - y1)/(x2 - x1);

    if(m1>=0.0) return 1;
    else return -1;
}

void mom_vector_maker_along_pminus_xmtoxp(  vector<comp> &qvec,
                                            comp s,
                                            double a,
                                            double m,
                                            double eps,
                                            double eps_for_m2k,
                                            double r,
                                            double points   )
{
    comp ii = {0.0,1.0};
    double xinitial = -1.0;
    double xfinal = 1.0;
    comp qc1 = q_c1(s,a,m,eps);
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp x1 = x_for_pk(s,qc1,q,m,eps);
    xinitial = real(x1);
    double delx1 = abs(xinitial-xfinal)/(points + 1.0);
    double delx2 = abs(xinitial-xfinal)/(points);

    vector<comp> tempvec;
    
    for(int i=0;i<(points+1);++i)
    //for(double x=xinitial;x<=xfinal;x=x+delx1)
    {
        double x = xinitial + i*delx1;
        comp pcut = pcut_minus_comp(x,s,q,m,eps);

        tempvec.push_back(pcut);
    }

    for(int i=1;i<tempvec.size();++i)
    {
        double x4,y4;
        double x1 = (double)real(tempvec[i-1]);
        double y1 = (double)imag(tempvec[i-1]);
        double x2 = (double)real(tempvec[i]);
        double y2 = (double)imag(tempvec[i]);

        int m1val = m1_checker(x1,y1,x2,y2);
        if(m1val==1)
        {
            xy_coordinates_m(x1,y1,x2,y2,r,x4,y4);
        }
        else
        {
            xy_coordinates_p(x1,y1,x2,y2,r,x4,y4);
        }
        comp qpoint = x4 + ii*y4;
        //cout<<"from up to down: x:"<<x4<<'\t'<<" y4:"<<y4<<endl;
        qvec.push_back(qpoint);

    }
}

void mom_vector_maker_along_pminus_xptoxm(  vector<comp> &qvec,
                                            comp s,
                                            double a,
                                            double m,
                                            double eps,
                                            double eps_for_m2k,
                                            double r,
                                            double points   )
{
    comp ii = {0.0,1.0};
    double xinitial = -1.0;
    double xfinal = 1.0;
    comp qc1 = q_c1(s,a,m,eps);
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp x1 = x_for_pk(s,qc1,q,m,eps);
    xinitial = real(x1);
    double delx1 = abs(xinitial-xfinal)/(points + 1.0);
    double delx2 = abs(xinitial-xfinal)/(points);

    vector<comp> tempvec;
    
    for(int i=0;i<points+1;++i)
    //for(double x=xfinal;x>xinitial;x=x-delx1)
    {
        double x = xfinal - i*delx1;
        comp pcut = pcut_minus_comp(x,s,q,m,eps);

        tempvec.push_back(pcut);
    }

    for(int i=1;i<tempvec.size();++i)
    {
        double x4,y4;
        double x1 = (double)real(tempvec[i-1]);
        double y1 = (double)imag(tempvec[i-1]);
        double x2 = (double)real(tempvec[i]);
        double y2 = (double)imag(tempvec[i]);

        int m1val = m1_checker(x1,y1,x2,y2);
        if(m1val==1)
        {
            xy_coordinates_p(x1,y1,x2,y2,r,x4,y4);
        }
        else
        {
            xy_coordinates_m(x1,y1,x2,y2,r,x4,y4);
        }

        comp qpoint = x4 + ii*y4;

        //cout<<"from down to up: x4:"<<x4<<'\t'<<" y4:"<<y4<<endl;
        //cout<<qpoint<<endl;
        qvec.push_back(qpoint);

    }
}

void half_circle_drawer_momrep(             vector<comp> &qvec,
                                            comp s,
                                            double a,
                                            double m,
                                            double eps,
                                            double eps_for_m2k,
                                            double r,
                                            double points     )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    double thetainitial = 0.0+0.0000001;
    double thetafinal = pi-0.0000001;
    double deltheta = abs(thetainitial - thetafinal)/(points);

    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp pcutm = pcut_minus(1.0,s,q,m,eps);

    

    for(int i=1;i<points+1;++i)
    //for(double theta=thetainitial+deltheta;theta<thetafinal;theta=theta+deltheta)
    {
        double theta = thetainitial + i*deltheta;
        comp qpoint = pcutm + r*exp(ii*theta);
        qvec.push_back(qpoint);
    }
}


void half_circle_drawer_momrep_1(           vector<comp> &qvec,
                                            comp s,
                                            double a,
                                            double m,
                                            double eps,
                                            double eps_for_m2k,
                                            double r,
                                            double points     )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    

    comp sigb = sigmab(a,m);
    comp qc1 = q_c1(s,a,m,eps);    
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp pcutm = pcut_minus(1.0,s,q,m,eps);


    double xinitial = -1.0;
    double xfinal = 1.0;
    comp x1angle = x_for_pk(s,qc1,q,m,eps);
    xinitial = real(x1angle);
    double delx1 = abs(xinitial-xfinal)/(points + 1.0);
    double delx2 = abs(xinitial-xfinal)/(points);

    comp firstval = qvec[qvec.size()-1];
    double firstangle = atan((imag(firstval) - imag(pcutm))/(real(firstval)-real(pcutm)));

    double thetainitial = abs(firstangle) + 0.00001;
    double thetafinal = pi + abs(firstangle) - 0.00001;
    double deltheta = abs(thetainitial - thetafinal)/(points);

    for(int i=1;i<points+1;++i)
    //for(double theta=thetainitial+deltheta;theta<thetafinal;theta=theta+deltheta)
    {
        double theta = thetainitial + i*deltheta;
        comp qpoint = pcutm + r*exp(ii*theta);
        qvec.push_back(qpoint);
    }
}



void mom_vector_maker_2(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points1,
                            double r    )
{
    comp ii = {0.0,1.0};
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    comp delq = {0.0,0.0};

    if(real(rad)>0.05)
    {
        /*  this is the initial straight line along the positive axis */
        comp startingpoint = min;
        delq = abs((abs(real(qc1)) - 1.5*r)-min)/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,+delq,points1,1);
        /*  this is the contour going along pcut down   */
        startingpoint = qvec[qvec.size()-1];
        //delq = abs((abs(real(qc1)) + eps) - real(startingpoint))/points1;
        //mom_vector_maker_linear_1(qvec,startingpoint,delq,points1,1);
        mom_vector_maker_along_pminus_xmtoxp(qvec,s,a,m,eps,eps_for_m2k,r,points1);
        
        
        
        /*  circle around the p+                */
        half_circle_drawer_momrep_1(qvec,s,a,m,eps,eps_for_m2k,r,points1/2.0);
        /*  This contour now goes along the pcut up             */
        startingpoint = qvec[qvec.size()-1];
        delq = abs(imag(startingpoint))/points1;
        mom_vector_maker_along_pminus_xptoxm(qvec,s,a,m,eps,eps_for_m2k,r,points1);
        
        /*  This contour now goes right in the real axis                       */
        startingpoint = qvec[qvec.size()-1];
        double retemp = real(startingpoint);
        double imtemp = imag(startingpoint);
        startingpoint = retemp;
        delq = abs(real(max) - real(startingpoint))/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,delq,points1/2.0,1);
        //cout<<"qvec after fifth part:"<<qvec.size()-tmp1<<endl;
        
        
        /* print the size of the created sigma_vector */
        cout<<"qvec created with size = "<<qvec.size()<<endl;

    }
    else if(real(rad)>0.0 && real(rad)<0.05)
    {
        qvec.push_back(min);
        mom_vector_maker_along_pminus_xptoxm(qvec,s,a,m,eps,eps_for_m2k,r,2.0*points1-1.0);
        //cout<<"here"<<endl;
        /*  This contour now goes right in the real axis                       */
        comp startingpoint = qvec[qvec.size()-1];
        double retemp = real(startingpoint);
        double imtemp = imag(startingpoint);
        startingpoint = retemp;
        delq = abs(real(max) - real(startingpoint))/(2.0*points1);
        mom_vector_maker_linear_1(qvec,startingpoint,delq,2.0*points1,1);

        
        /* print the size of the created sigma_vector */
        cout<<"qvec created with size = "<<qvec.size()<<endl;

    }
    else
    {
        double length = real(max - min);
        comp delsig = length/(4.0*points1) ;

        int totpoints = 4.0*points1;
        for(int i=0;i<totpoints;++i)
        {
            comp qmom = min + ii*eps + ((comp)i)*delsig;

            qvec.push_back(qmom);
        }
        cout<<"qvec created with size = "<<qvec.size()<<endl;
    }

}


void check_which_pcut(  comp s,
                        comp q,
                        double m,
                        double eps,
                        int &flag   )
{
    //int flag;

    double xinitial = -1.0; //these are the 'x' values of the pcut
    double xfinal = 1.0;
    double points = 250.0;
    double delx = abs(xinitial - xfinal)/points;

    double tempxval1, tempxval2; //these are x axis values for the two cuts
    double xval1, xval2;
    int startcount = 0;
    for(double x=xinitial;x<=xfinal;x=x+delx)
    {
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        if(startcount==0)
        {
            xval1 = (double) real(pcutp);
            xval2 = (double) real(pcutm);
            startcount = 1;
        }

        tempxval1 = (double) real(pcutp);
        tempxval2 = (double) real(pcutm);

        if(tempxval1<xval1) xval1 = tempxval1;
        if(tempxval2<xval2) xval2 = tempxval2;
    }

    if(xval1<xval2) flag = 0;
    else flag = 1;

}

//this function tells which pcutvec is to the right 
//and which one is to the left. flag=0 means pcutpvec is 
//to the right and flag=1 means pcutmvec is to the right
void check_pcut_leftright(  vector<comp> pcutpvec,
                            vector<comp> pcutmvec,
                            int &flag   )
{
    //int flag;

    

    double tempxval1, tempxval2; //
    double xval1, xval2;
    int startcount = 0;
    //for(double x=xinitial;x<=xfinal;x=x+delx)
    for(int i=0;i<pcutpvec.size();++i)
    {
        comp pcutp = pcutpvec[i];//pcut_plus(x,s,q,m,eps);
        comp pcutm = pcutmvec[i];//pcut_minus(x,s,q,m,eps);

        if(startcount==0)
        {
            xval1 = (double) real(pcutp);
            xval2 = (double) real(pcutm);
            startcount = 1;
        }

        tempxval1 = (double) real(pcutp);
        tempxval2 = (double) real(pcutm);

        if(tempxval1>xval1) xval1 = tempxval1;
        if(tempxval2>xval2) xval2 = tempxval2;
    }

    if(xval1<xval2) flag = 1;
    else flag = 0;

}



//this function gives out the lowest y-values of 
//pcutpvec and pcutmvec

void lowest_yval_pcuts( vector<comp> &pcutpvec,
                        vector<comp> &pcutmvec,
                        double &yvalp,
                        double &yvalm,
                        int &lowest_p_index,
                        int &lowest_m_index       )
{
    double tempyvalp, tempyvalm;
    for(int i=0;i<pcutpvec.size();++i)
    {
        if(i==0)
        {
            yvalp = imag(pcutpvec[i]);
            yvalm = imag(pcutmvec[i]);
            lowest_p_index = i;
            lowest_m_index = i;
            //continue;
        }

        tempyvalp = imag(pcutpvec[i]);
        tempyvalm = imag(pcutmvec[i]);

        if(tempyvalp<yvalp)
        {
            yvalp = tempyvalp;
            lowest_p_index = i;
        }
        if(tempyvalm<yvalm)
        {
            yvalm = tempyvalm;
            lowest_m_index = i;
        }
    }
    //cout<<"yvalp = "<<yva
}

void highest_yval_pcuts( vector<comp> &pcutpvec,
                        vector<comp> &pcutmvec,
                        double &yvalp,
                        double &yvalm,
                        int &lowest_p_index,
                        int &lowest_m_index       )
{
    double tempyvalp, tempyvalm;
    for(int i=0;i<pcutpvec.size();++i)
    {
        if(i==0)
        {
            yvalp = imag(pcutpvec[i]);
            yvalm = imag(pcutmvec[i]);
            lowest_p_index = i;
            lowest_m_index = i;
            //continue;
        }

        tempyvalp = imag(pcutpvec[i]);
        tempyvalm = imag(pcutmvec[i]);

        if(tempyvalp>yvalp)
        {
            yvalp = tempyvalp;
            lowest_p_index = i;
        }
        if(tempyvalm>yvalm)
        {
            yvalm = tempyvalm;
            lowest_m_index = i;
        }
    }
    //cout<<"yvalp = "<<yva
}


void highest_xval_pcuts(    vector<comp> &pcutpvec,
                            vector<comp> &pcutmvec,
                            double &xvalp,
                            double &xvalm,
                            int &highest_p_index,
                            int &highest_m_index       )
{
    double tempxvalp, tempxvalm;
    for(int i=0;i<pcutpvec.size();++i)
    {
        if(i==0)
        {
            xvalp = real(pcutpvec[i]);
            xvalm = real(pcutmvec[i]);
            highest_p_index = i;
            highest_m_index = i;
            //continue;
        }

        tempxvalp = real(pcutpvec[i]);
        tempxvalm = real(pcutmvec[i]);

        if(tempxvalp>xvalp)
        {
            xvalp = tempxvalp;
            highest_p_index = i;
        }
        if(tempxvalm>xvalm)
        {
            xvalm = tempxvalm;
            highest_m_index = i;
        }
    }
    //cout<<"yvalp = "<<yva
}


int sign_checker(   double a    )
{
    if(a>=0.0) return 1;
    else return -1;
}

void pcutvec_realaxis_crossing_checker( vector<comp> pcutpvec,
                                        vector<comp> pcutmvec,
                                        int &flag   )
{
    double some_index1;
    double some_index2;
    double tempyvalp, tempyvalm, yvalp, yvalm;
    int somecountp = 0;
    int somecountm = 0;
    for(int i=0;i<pcutpvec.size();++i)
    {
        comp pcutp = pcutpvec[i];
        comp pcutm = pcutmvec[i];

        if(real(pcutp)>0.0)
        {
            if(somecountp==0)
            {
                yvalp = imag(pcutp);
                somecountp = 1;
                
            }

            tempyvalp = imag(pcutp);
            if(tempyvalp<yvalp) yvalp = tempyvalp;

            
        }

        if(real(pcutm)>0.0)
        {
            if(somecountm==0)
            {
                yvalp = imag(pcutm);
                somecountm = 1;
            }

            tempyvalm = imag(pcutm);
            if(tempyvalm<yvalm) yvalm = tempyvalm;

        }
    }

    if(yvalp*yvalm>0.0) flag = 0; //flag=0 means deform the contour
                                  //flag=1 means take a straight line through the real axis;
    else flag = 1;


}

void pcutvec_realaxis_crossing_checker_1(   vector<comp> pcutpvec,
                                            vector<comp> pcutmvec,
                                            int &flag   )
{
    double some_index1;
    double some_index2;
    double tempyvalp, tempyvalm, yvalp, yvalm;
    int somecountp = 0;
    int somecountm = 0;
    int prevpsign;
    int prevmsign;
    int signchanged_for_pcutp = 0;
    int signchanged_for_pcutm = 0;
    for(int i=0;i<pcutpvec.size();++i)
    {
        comp pcutp = pcutpvec[i];
        comp pcutm = pcutmvec[i];

        if(real(pcutp)>0.0)
        {
            int signp = sign_checker(imag(pcutp));

            if(somecountp==0)
            {
                prevpsign = signp;
                somecountp = 1;
                //continue;
                
            }
            else 
            {
                if(prevpsign==signp)
                {
                    signchanged_for_pcutp=0; 
                } 
                else
                {
                    signchanged_for_pcutp=1;
                    break;
                }
                prevpsign = signp;
            }
            //tempyvalp = imag(pcutp);
            //if(tempyvalp<yvalp) yvalp = tempyvalp;
            
            
        }

        if(real(pcutm)>0.0)
        {
            int signm = sign_checker(imag(pcutm));

            if(somecountm==0)
            {
                prevmsign = signm;
                somecountm = 1;
                //continue;
                
            }
            else 
            {

            //signchanged_for_pcut = 0 means sign didn't change
            //signchanged_for_pcut = 1 means sign did change 
                if(prevmsign==signm)
                {
                    signchanged_for_pcutm=0; 
                } 
                else
                {
                    signchanged_for_pcutm=1;
                    break;
                }
                prevmsign = signm;
            }

            
            //tempyvalp = imag(pcutp);
            //if(tempyvalp<yvalp) yvalp = tempyvalp;

        }
    }

    if(signchanged_for_pcutp==0 && signchanged_for_pcutm==0) flag = 1;  
    else if(signchanged_for_pcutp==1 && signchanged_for_pcutm==0) flag = 0;
    else if(signchanged_for_pcutp==0 && signchanged_for_pcutm==1) flag = 0;
    else if(signchanged_for_pcutp==1 && signchanged_for_pcutm==1) flag = 0;
    //flag=0 means deform the contour
    //flag=1 means take a straight line through the real axis;

}


void check_which_pcut_updown(   comp s,
                                comp q,
                                double m,
                                double eps,
                                int &flag   )
{
    //int flag;

    double xinitial = -1.0; //these are the 'x' values of the pcut
    double xfinal = 1.0;
    double points = 250.0;
    double delx = abs(xinitial - xfinal)/points;

    comp pcutp_ini = pcut_plus(xinitial,s,q,m,eps);
    comp pcutm_ini = pcut_minus(xinitial,s,q,m,eps);
    comp pcutp_fin = pcut_plus(xfinal,s,q,m,eps);
    comp pcutm_fin = pcut_minus(xfinal,s,q,m,eps);

    comp select_pcutp, select_pcutm;
    if(imag(pcutp_ini)>imag(pcutp_fin)) select_pcutp = pcutp_ini;
    else select_pcutp = pcutp_fin;

    if(imag(pcutm_ini)>imag(pcutm_fin)) select_pcutm = pcutm_ini;
    else select_pcutm = pcutm_fin;

    if(imag(select_pcutp)>imag(select_pcutm)) flag = 1;
    else flag = 0;

    //if flag = 0 then pcutm is above and pcutp is in the bottom
    //if flag = 1 then pcutp is above and pcutm is in the bottom


}

void mom_vector_maker_3(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    comp delq = {0.0,0.0};

    int flag; 
    check_which_pcut(s,q,m,eps,flag);
    if(flag==0)
    {
        double xinitial = -1.0;
        double xfinal = 1.0;
        double xpoints = 250;
        double delx = abs(xinitial-xfinal)/xpoints;
        double xval1,yval1;
        int somecount=0;

        for(double x=xinitial;x<=xfinal;x=x+delx)
        {
            comp pcutm = pcut_minus(x,s,q,m,eps);
            if(somecount==0)
            {
                xval1=(double) real(pcutm);
                yval1=(double) imag(pcutm);
                somecount=1;
            }
            double tempxval1 = (double) real(pcutm);
            double tempyval1 = (double) imag(pcutm);
            if(tempxval1<xval1) xval1 = tempxval1;
            if(tempyval1<yval1) yval1 = tempyval1;

        }

        if(xval1<=1.0e-8 && yval1<=1.0e-3)
        {
            comp startingpoint = min;
            comp delq = abs(startingpoint - r)/(points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint) - yval1 + r)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,1);

        }
        else if(xval1>1.0e-8 && yval1<=1.0e-3)
        {
            comp startingpoint = min;
            delq = abs(imag(startingpoint) - yval1 + r)/(points + points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points + points/10.0,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,1);
        }
        else if(yval1>1.0e-3)
        {
            comp startingpoint = min;
            double totpoints = 4.0*points + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);
        }
    }
    else 
    {
        double xinitial = -1.0;
        double xfinal = 1.0;
        double xpoints = 250;
        double delx = abs(xinitial-xfinal)/xpoints;
        double xval1,yval1;
        int somecount=0;

        for(double x=xinitial;x<=xfinal;x=x+delx)
        {
            comp pcutp = pcut_plus(x,s,q,m,eps); 
            if(somecount==0)
            {
                xval1=(double) real(pcutp);
                yval1=(double) imag(pcutp);
                somecount=1;
            }
            double tempxval1 = (double) real(pcutp);
            double tempyval1 = (double) imag(pcutp);
            if(tempxval1<xval1) xval1 = tempxval1;
            if(tempyval1<yval1) yval1 = tempyval1;

        }

        if(xval1<=1.0e-8 && yval1<=1.0e-3)
        {
            comp startingpoint = min;
            comp delq = abs(startingpoint - r)/(points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint) - yval1 + r)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,1);

        }
        else if(xval1>1.0e-8 && yval1<=1.0e-3)
        {
            comp startingpoint = min;
            delq = abs(imag(startingpoint) - yval1 + r)/(points + points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points + points/10.0,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,1);
        }
        else if(yval1>1.0e-3)
        {
            comp startingpoint = min;
            double totpoints = 4.0*points + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);
        }
    }
}



void closest_pcuts_bottom(  comp s,
                            comp q,
                            double m,
                            double eps,
                            comp &pcutp_close,
                            comp &pcutm_close,
                            double &x_pcutp,
                            double &x_pcutm )
{
    comp ii = {0.0,1.0};
    
    double initialx = 0.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    double tempdiff = 1.0e+16;
    double tempxp = 0.0;
    double tempxm = 0.0;

    for(int i=0;i<pcutpvec.size();++i)
    {
        for(int j=0;j<pcutmvec.size();++j)
        {
            double diff = abs(pcutpvec[i] - pcutmvec[j]);

            //cout<<"xp = "<<xvec[i]<<'\t'<<"xm = "<<xvec[j]<<endl;

            if(diff<=tempdiff)
            {
                //cout<<"diff minimum then previous"<<endl;
                //cout<<"diff = "<<diff<<endl;
                //if(imag(pcutpvec[i])<0.0 && imag(pcutmvec[j])<0.0)
                {
                    tempdiff = diff;
                    tempxp = xvec[i];
                    tempxm = xvec[j];
                }
            }
            
        }
    }

    pcutp_close = pcut_plus(tempxp,s,q,m,eps);
    pcutm_close = pcut_minus(tempxm,s,q,m,eps);
    x_pcutp = tempxp;
    x_pcutm = tempxm;
    cout<<"diff found in pcuts = "<<abs(pcutp_close-pcutm_close)<<endl;
    
}


//this function takes the pcutpvec, pcutmvec
//and finds the closest gap between the two 
//cuts and returns the closest pcutp1 and pcutm1,
//along with the second pair of these values 
void closest_pcuts(             vector<double> xvec,
                                vector<comp> pcutpvec,
                                vector<comp> pcutmvec,
                                comp s,
                                comp q,
                                double m,
                                double eps,
                                int &first_i,            //i is for pcutp index
                                int &first_j,            //j is for pcutm index
                                int &second_i,
                                int &second_j      )
{
    comp ii = {0.0,1.0};
    
    

    double tempdiff = 1.0e+16;
    double tempxp = 0.0;
    double tempxm = 0.0;
    
    for(int i=0;i<pcutpvec.size();++i)
    {
        for(int j=0;j<pcutmvec.size();++j)
        {
            double diff = abs(pcutpvec[i] - pcutmvec[j]);

            //cout<<"xp = "<<xvec[i]<<'\t'<<"xm = "<<xvec[j]<<endl;

            if(diff<=tempdiff)
            {
                //cout<<"diff minimum then previous"<<endl;
                //cout<<"diff = "<<diff<<endl;
                //if(imag(pcutpvec[i])<0.0 && imag(pcutmvec[j])<0.0)
                {
                    tempdiff = diff;
                    tempxp = xvec[i];
                    tempxm = xvec[j];
                    first_i = i;
                    first_j = j;
                }
            }
            
        }
    }

    comp pcutp_close1, pcutm_close1;
    pcutp_close1 = pcutpvec[first_i];//pcut_plus(tempxp,s,q,m,eps);
    pcutm_close1 = pcutmvec[first_j];//pcut_minus(tempxm,s,q,m,eps);
    
    //cout<<"first_i = "<<first_i<<endl;
    //cout<<"first_j = "<<first_j<<endl;

    //cout<<"pcutp_close1 = "<<pcutp_close1<<endl;
    //cout<<"pcutm_close1 = "<<pcutm_close1<<endl;
    
    tempdiff = 1.0e+16;

    for(int i=0;i<pcutpvec.size();++i)
    {
        for(int j=0;j<pcutmvec.size();++j)
        {
            double diff = abs(pcutpvec[i] - pcutmvec[j]);

            //cout<<"xp = "<<xvec[i]<<'\t'<<"xm = "<<xvec[j]<<endl;
            if(i==first_i && j==first_j)
            {
                
            }
            else
            {
                if(diff<=tempdiff)
                {
                    //cout<<"diff minimum then previous"<<endl;
                    //cout<<"diff = "<<diff<<endl;
                    //if(imag(pcutpvec[i])<0.0 && imag(pcutmvec[j])<0.0)
                
                    tempdiff = diff;
                    tempxp = xvec[i];
                    tempxm = xvec[j];
                    second_i = i;
                    second_j = j;
                
                }
            }
            
        }
    }

    comp pcutp_close2, pcutm_close2;
    pcutp_close2 = pcutpvec[second_i];//pcut_plus(tempxp,s,q,m,eps);
    pcutm_close2 = pcutmvec[second_j];//pcut_minus(tempxm,s,q,m,eps);

    //cout<<"second_i = "<<second_i<<endl;
    //cout<<"second_j = "<<second_j<<endl;

    //cout<<"pcutp_close2 = "<<pcutp_close2<<endl;
    //cout<<"pcutm_close2 = "<<pcutm_close2<<endl;
    
    cout<<"diff found in pcuts1 = "<<abs(pcutp_close1-pcutm_close1)<<endl;
    cout<<"diff found in pcuts2 = "<<abs(pcutp_close2-pcutm_close2)<<endl;
    
}


void closest_pcuts_bottom_withpcutvec(  vector<comp> &pcutpvec,
                                        vector<comp> &pcutmvec,
                                        vector<double> &xvec,
                                        comp s,
                                        comp q,
                                        double m,
                                        double eps,
                                        comp &pcutp_close,
                                        comp &pcutm_close,
                                        double &x_pcutp,
                                        double &x_pcutm )
{
    comp ii = {0.0,1.0};
    
    double initialx = 0.0;
    double finalx = 1.0;
    double size = 50000.0;
    double delx = abs(initialx-finalx)/size;

    int startingindex = 0;
    //check where x starts becoming positive 

    /*for(int i=0;i<xvec.size();++i)
    {
        double xval = xvec[i];
        if(xval>=0.0)
        {
            startingindex = i;
            break;
        }
    }*/

    
    double tempdiff = 1.0e+16;
    double tempxp = 0.0;
    double tempxm = 0.0;

    int savedindexi = 0;
    int savedindexj = 0;

    for(int i=startingindex;i<pcutpvec.size();++i)
    {
        for(int j=startingindex;j<pcutmvec.size();++j)
        {
            double diff = abs(pcutpvec[i] - pcutmvec[j]);

            //cout<<"xp = "<<xvec[i]<<'\t'<<"xm = "<<xvec[j]<<endl;

            if(diff<=tempdiff)
            {
                //cout<<"diff minimum then previous"<<endl;
                //cout<<"diff = "<<diff<<endl;
                //if(imag(pcutpvec[i])<0.0 && imag(pcutmvec[j])<0.0)
                {
                    tempdiff = diff;
                    tempxp = xvec[i];
                    tempxm = xvec[j];
                    savedindexi = i;
                    savedindexj = j;
                }
            }
            
        }
    }

    pcutp_close = pcutpvec[savedindexi];//pcut_plus(tempxp,s,q,m,eps);
    pcutm_close = pcutmvec[savedindexj];//pcut_minus(tempxm,s,q,m,eps);
    x_pcutp = tempxp;
    x_pcutm = tempxm;
    cout<<"diff found in pcuts = "<<abs(pcutp_close-pcutm_close)<<endl;
    
}



void slanted_line_maker(    vector<comp> &qvec,
                            comp startingpoint,
                            comp endingpoint,
                            double offset,
                            double points       )
{
    comp ii = {0.0,1.0};
    double x1 = real(startingpoint);
    double y1 = imag(startingpoint);
    double x2 = real(endingpoint);
    double y2 = imag(endingpoint);
    double m = (y2 - y1)/(x2 - x1);
    m = - abs(m);
    //cout<<"m:"<<m<<endl;
    double length = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
    double r = length + 2.0*offset;
    double theta = atan(m);

    double delr = r/points;

    for(int i=1;i<points+1;++i)
    {
        double rval = i*delr;

        comp qval = startingpoint + rval*(cos(theta) + ii*sin(theta));

        qvec.push_back(qval);
    }

    
}

//this is with arbitrary sign of m
void slanted_line_maker_1(    vector<comp> &qvec,
                            comp startingpoint,
                            comp endingpoint,
                            double offset,
                            double points       )
{
    comp ii = {0.0,1.0};
    double x1 = real(startingpoint);
    double y1 = imag(startingpoint);
    double x2 = real(endingpoint);
    double y2 = imag(endingpoint);
    double m = (y2 - y1)/(x2 - x1);
    cout<<"m found = "<<m<<endl;
    //m = - abs(m);
    //cout<<"m:"<<m<<endl;
    double length = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
    double r = length + 2.0*offset;
    double theta = atan(m);

    double delr = r/points;

    for(int i=1;i<points+1;++i)
    {
        double rval = i*delr;

        comp qval = startingpoint + rval*(cos(theta) + ii*sin(theta));

        qvec.push_back(qval);
    }

    
}

void slanted_line_maker_2(    vector<comp> &qvec,
                            comp startingpoint,
                            comp endingpoint,
                            double offset,
                            double points,
                            double thetaoffset       )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    double x1 = real(startingpoint);
    double y1 = imag(startingpoint);
    double x2 = real(endingpoint);
    double y2 = imag(endingpoint);
    double m = (y2 - y1)/(x2 - x1);
    cout<<"m found = "<<m<<endl;
    //m = - abs(m);
    //cout<<"m:"<<m<<endl;
    double length = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
    double r = length + 2.0*offset;
    double theta = atan(m) + thetaoffset;

    double delr = r/points;

    for(int i=1;i<points+1;++i)
    {
        double rval = i*delr;

        comp qval = startingpoint + rval*(cos(theta) + ii*sin(theta));

        qvec.push_back(qval);
    }

    
}

void line_maker(    vector<comp> &qvec, 
                    comp a,
                    comp b,
                    double points   )
{
    double zinitial = 0.0;
    double zfinal = 1.0;
    double delz = abs(zfinal-zinitial)/points;
    for(int i=1;i<points+1;++i)
    {
        double z = zinitial + i*delz;
        comp x_of_z = (b-a)*z + a;

        qvec.push_back(x_of_z);
    }
}


//This was being used previously that creates the momentum 
//vector, that goes to left->slanted through pcuts->right
//->up->right
void mom_vector_maker_4(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    comp delq = {0.0,0.0};

    comp pcutp_close;
    comp pcutm_close;
    double x_pcutp;
    double x_pcutm;
    closest_pcuts_bottom(s,q,m,eps,pcutp_close,pcutm_close,x_pcutp,x_pcutm);


    int flag; 
    check_which_pcut(s,q,m,eps,flag);
    if(flag==0)
    {
        double xinitial = -1.0;
        double xfinal = 1.0;
        double xpoints = 250;
        double delx = abs(xinitial-xfinal)/xpoints;
        double xval1,yval1;
        int somecount=0;



        for(double x=xinitial;x<=xfinal;x=x+delx)
        {
            comp pcutm = pcut_minus(x,s,q,m,eps);
            if(somecount==0)
            {
                xval1=(double) real(pcutm);
                yval1=(double) imag(pcutm);
                somecount=1;
            }
            double tempxval1 = (double) real(pcutm);
            double tempyval1 = (double) imag(pcutm);
            if(tempxval1<xval1) xval1 = tempxval1;
            if(tempyval1<yval1) yval1 = tempyval1;

        }



        if(xval1<=1.0e-8 && yval1<=1.0e-5)
        {
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close + pcutm_close)/2.0;
            slanted_line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,r,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max - startingpoint)/(3.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,3.0*points/10.0,1);

        }
        else if(xval1>1.0e-8 && yval1<=1.0e-5)
        {
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close + pcutm_close)/2.0;
            slanted_line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,r,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max - startingpoint)/(3.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,3.0*points/10.0,1);
        }
        else if(yval1>=0.0)
        {
            comp startingpoint = min;
            double totpoints = points;// + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);
        }
    }
    else 
    {
        double xinitial = -1.0;
        double xfinal = 1.0;
        double xpoints = 250;
        double delx = abs(xinitial-xfinal)/xpoints;
        double xval1,yval1;
        int somecount=0;

        for(double x=xinitial;x<=xfinal;x=x+delx)
        {
            comp pcutp = pcut_plus(x,s,q,m,eps); 
            if(somecount==0)
            {
                xval1=(double) real(pcutp);
                yval1=(double) imag(pcutp);
                somecount=1;
            }
            double tempxval1 = (double) real(pcutp);
            double tempyval1 = (double) imag(pcutp);
            if(tempxval1<xval1) xval1 = tempxval1;
            if(tempyval1<yval1) yval1 = tempyval1;

        }

        if(xval1<=1.0e-8 && yval1<=1.0e-5)
        {
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close + pcutm_close)/2.0;
            slanted_line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,r,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max - startingpoint)/(3.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,3.0*points/10.0,1);

        }
        else if(xval1>1.0e-8 && yval1<=1.0e-5)
        {
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close + pcutm_close)/2.0;
            slanted_line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,r,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max - startingpoint)/(3.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,3.0*points/10.0,1);
        }
        else if(yval1>=0.0)
        {
            comp startingpoint = min;
            double totpoints = points;// + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);
        }
    }
}


//this section does not work
void mom_vector_maker_5(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    comp delq = {0.0,0.0};

    comp pcutp_close;
    comp pcutm_close;
    double x_pcutp;
    double x_pcutm;
    //closest_pcuts_bottom(s,q,m,eps,pcutp_close,pcutm_close,x_pcutp,x_pcutm);


    int flag; 
    check_which_pcut(s,q,m,eps,flag);
    if(flag==0)
    {
        double xinitial = -1.0;
        double xfinal = 1.0;
        double xpoints = 250;
        double delx = abs(xinitial-xfinal)/xpoints;
        double xval1,yval1;
        int somecount=0;



        for(double x=xinitial;x<=xfinal;x=x+delx)
        {
            comp pcutm = pcut_minus(x,s,q,m,eps);
            if(somecount==0)
            {
                xval1=(double) real(pcutm);
                yval1=(double) imag(pcutm);
                somecount=1;
            }
            double tempxval1 = (double) real(pcutm);
            double tempyval1 = (double) imag(pcutm);
            if(tempxval1<xval1) xval1 = tempxval1;
            if(tempyval1<yval1) yval1 = tempyval1;

        }

        comp pcut_endpoint1 = xval1 + ii*yval1;
        comp pcut_endpoint2 = pcut_minus(-1.0,s,q,m,eps);

        comp startingpoint = min;
        comp delq = abs(startingpoint - rad/2.0)/(points/10.0);
        mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/10.0,1);

        startingpoint = qvec[qvec.size()-1];
        comp endingpoint_for_slanted_maker = (pcut_endpoint1 + pcut_endpoint2)/2.0;
        slanted_line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,0.0,2.0*points/10.0);//offset has been set to zero





        if(xval1<=1.0e-8 && yval1<=1.0e-5)
        {
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close + pcutm_close)/2.0;
            slanted_line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,r,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max - startingpoint)/(3.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,3.0*points/10.0,1);

        }
        else if(xval1>1.0e-8 && yval1<=1.0e-5)
        {
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close + pcutm_close)/2.0;
            slanted_line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,r,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max - startingpoint)/(3.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,3.0*points/10.0,1);
        }
        else if(yval1>=0.0)
        {
            comp startingpoint = min;
            double totpoints = points;// + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);
        }
    }
    else 
    {
        double xinitial = -1.0;
        double xfinal = 1.0;
        double xpoints = 250;
        double delx = abs(xinitial-xfinal)/xpoints;
        double xval1,yval1;
        int somecount=0;

        for(double x=xinitial;x<=xfinal;x=x+delx)
        {
            comp pcutp = pcut_plus(x,s,q,m,eps); 
            if(somecount==0)
            {
                xval1=(double) real(pcutp);
                yval1=(double) imag(pcutp);
                somecount=1;
            }
            double tempxval1 = (double) real(pcutp);
            double tempyval1 = (double) imag(pcutp);
            if(tempxval1<xval1) xval1 = tempxval1;
            if(tempyval1<yval1) yval1 = tempyval1;

        }

        if(xval1<=1.0e-8 && yval1<=1.0e-5)
        {
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close + pcutm_close)/2.0;
            slanted_line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,r,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max - startingpoint)/(3.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,3.0*points/10.0,1);

        }
        else if(xval1>1.0e-8 && yval1<=1.0e-5)
        {
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close + pcutm_close)/2.0;
            slanted_line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,r,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max - startingpoint)/(3.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,3.0*points/10.0,1);
        }
        else if(yval1>=0.0)
        {
            comp startingpoint = min;
            double totpoints = points;// + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);
        }
    }
}


//this one is used for checking the anomalous behavious
//in Mphib for Im(s)>0 
void mom_vector_maker_6(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double r,
                            int &tag1,
                            int &tag2     )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    comp delq = {0.0,0.0};

    comp pcutp_close;
    comp pcutm_close;
    double x_pcutp;
    double x_pcutm;
    //closest_pcuts_bottom(s,q,m,eps,pcutp_close,pcutm_close,x_pcutp,x_pcutm);


    int flag; 
    if(real(s)<8.925)
    check_which_pcut(s,q,m,eps,flag);
    else
    check_which_pcut_updown(s,q,m,eps,flag);
    //flag = 0;
    if(flag==0)
    {
        double xval = 0;

        comp pcutm = pcut_minus(xval+0.65,s,q,m,eps);
        comp pcutp = pcut_plus(xval-0.45,s,q,m,eps);

        comp startingpoint = min;
        double thetaoffset = acos(-1.0);

        comp endingpoint_for_slanted_maker = pcutp;

        line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
        //slanted_line_maker_2(qvec,startingpoint,endingpoint_for_slanted_maker,0.0,points,0.0);//offset set to zero
        tag1 = qvec.size() - 1;
        
        startingpoint = qvec[qvec.size()-1];
        endingpoint_for_slanted_maker = pcutm;
        line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
        //slanted_line_maker_2(qvec,startingpoint,endingpoint_for_slanted_maker,0.0,points,thetaoffset);//offset set to zero
        tag2 = qvec.size() - 1;

        startingpoint = qvec[qvec.size()-1];
        endingpoint_for_slanted_maker = max;
        line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
        //slanted_line_maker_2(qvec,startingpoint,endingpoint_for_slanted_maker,0.0,points,0.0);//offset set to zero

    }
    else 
    {
        double xval = 0;

        comp pcutm = pcut_minus(xval+0.65,s,q,m,eps);
        comp pcutp = pcut_plus(xval-0.45,s,q,m,eps);

        comp startingpoint = min;
        double thetaoffset = acos(-1.0);

        comp endingpoint_for_slanted_maker = pcutp;
        line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
        //slanted_line_maker_2(qvec,startingpoint,endingpoint_for_slanted_maker,0.0,points,0.0);//offset set to zero
        tag1 = qvec.size() - 1;
        
        startingpoint = qvec[qvec.size()-1];
        endingpoint_for_slanted_maker = pcutm;
        line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
        //slanted_line_maker_2(qvec,startingpoint,endingpoint_for_slanted_maker,0.0,points,thetaoffset);//offset set to zero
        tag2 = qvec.size() - 1;

        startingpoint = qvec[qvec.size()-1];
        endingpoint_for_slanted_maker = max;
        line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
        //slanted_line_maker_2(qvec,startingpoint,endingpoint_for_slanted_maker,0.0,points,0.0);//offset set to zero

    }
}

//same as 6 but a change in pcutp and pcutm
void mom_vector_maker_61(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double r,
                            int &tag1,
                            int &tag2     )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    comp delq = {0.0,0.0};

    comp pcutp_close;
    comp pcutm_close;
    double x_pcutp;
    double x_pcutm;
    //closest_pcuts_bottom(s,q,m,eps,pcutp_close,pcutm_close,x_pcutp,x_pcutm);


    int flag; 
    if(real(s)<8.925)
    check_which_pcut(s,q,m,eps,flag);
    else
    check_which_pcut_updown(s,q,m,eps,flag);
    //flag = 0;
    if(flag==0)
    {
        double xval = 0;

        comp pcutm = pcut_minus(xval-0.75,s,q,m,eps);
        comp pcutp = pcut_plus(xval+0.65,s,q,m,eps);

        comp startingpoint = min;
        double thetaoffset = acos(-1.0);

        comp endingpoint_for_slanted_maker = pcutp;

        line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
        //slanted_line_maker_2(qvec,startingpoint,endingpoint_for_slanted_maker,0.0,points,0.0);//offset set to zero
        tag1 = qvec.size() - 1;
        
        startingpoint = qvec[qvec.size()-1];
        endingpoint_for_slanted_maker = pcutm;
        line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
        //slanted_line_maker_2(qvec,startingpoint,endingpoint_for_slanted_maker,0.0,points,thetaoffset);//offset set to zero
        tag2 = qvec.size() - 1;

        startingpoint = qvec[qvec.size()-1];
        endingpoint_for_slanted_maker = max;
        line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
        //slanted_line_maker_2(qvec,startingpoint,endingpoint_for_slanted_maker,0.0,points,0.0);//offset set to zero

    }
    else 
    {
        double xval = 0;

        comp pcutm = pcut_minus(xval-0.75,s,q,m,eps);
        comp pcutp = pcut_plus(xval+0.65,s,q,m,eps);

        comp startingpoint = min;
        double thetaoffset = acos(-1.0);

        comp endingpoint_for_slanted_maker = pcutp;
        line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
        //slanted_line_maker_2(qvec,startingpoint,endingpoint_for_slanted_maker,0.0,points,0.0);//offset set to zero
        tag1 = qvec.size() - 1;
        
        startingpoint = qvec[qvec.size()-1];
        endingpoint_for_slanted_maker = pcutm;
        line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
        //slanted_line_maker_2(qvec,startingpoint,endingpoint_for_slanted_maker,0.0,points,thetaoffset);//offset set to zero
        tag2 = qvec.size() - 1;

        startingpoint = qvec[qvec.size()-1];
        endingpoint_for_slanted_maker = max;
        line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
        //slanted_line_maker_2(qvec,startingpoint,endingpoint_for_slanted_maker,0.0,points,0.0);//offset set to zero

    }
}


//This is the same vector maker as maker_4 but 
//only difference is the which pcut is to the 
//right has been changed, along with a pcut fixer 
//has been used and slanted line makers have been 
//changed with just line_maker which is a better 
//line maker maker
void mom_vector_maker_41(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    double qvec_r = 0.01;

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=0 means deform the contour first backward and then all the other steps
    //problem_index_flag=1 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);
    
    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,1);

            startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            startingpoint = qvec[qvec.size()-1];
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[lowest_p_index] - ii*0.01;
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[lowest_m_index] - ii*0.01;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            startingpoint = qvec[qvec.size()-1];
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
        }
        else
        {
            cout<<"deforming the contour downward first"<<endl;
            //comp startingpoint = min;
            //comp delq = abs(startingpoint - rad/2.0)/(points);
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,1);

            //startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            //comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
            }
            else
            {
                endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);

            startingpoint = qvec[qvec.size()-1];
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);
        }

    }
    else
    {
        cout<<"taking the straight line"<<endl;
        //cout<<"prob here with axis_flag="<<axis_flag<<endl;
        comp startingpoint = min;
        double totpoints = 4*points;// + points/10.0;
        delq = abs(real(min) - real(max))/(totpoints);
        mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);
    }



    
}


//this is the same as maker_41, only change that 
//has been included is that if the midpoint of the 
//difference of the cuts is lower than the lowest y
//value of the cut to the right, then the contour 
//goes straight to the other endpoint 'max'. Also
//if difference between the x axis and the lowest 
//yvalue of the right cut is less than 1.0e-2, we will
//still deform the contour downward.
void mom_vector_maker_42(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    double qvec_r = 0.01;

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=0 means deform the contour first backward and then all the other steps
    //problem_index_flag=1 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);
    

    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,1);

            startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            startingpoint = qvec[qvec.size()-1];
            double lowest_index_yval = 0.0;
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[lowest_p_index] - ii*0.01;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[lowest_m_index] - ii*0.01;
            }

            if(lowest_index_yval<=imag(startingpoint))
            {
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

                startingpoint = qvec[qvec.size()-1];
                endingpoint_for_slanted_maker = max;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
            }
            else
            {
                startingpoint = qvec[qvec.size()-1];
                endingpoint_for_slanted_maker = max;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);
            }
        }
        else
        {
            cout<<"deforming the contour downward first"<<endl;
            //comp startingpoint = min;
            //comp delq = abs(startingpoint - rad/2.0)/(points);
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,1);

            //startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            //comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
            }
            else
            {
                endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);

            startingpoint = qvec[qvec.size()-1];
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            lowest_index_yval = imag(pcutpvec[lowest_p_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
        }
        else
        {
            lowest_index_yval = imag(pcutmvec[lowest_m_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
        }


        if(abs(lowest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour downward first"<<endl;
            comp startingpoint = min;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);

            endingpoint_for_slanted_maker = max;
            startingpoint = qvec[qvec.size()-1];

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);


        }
        else 
        {
            cout<<"taking the straight line"<<endl;
            //cout<<"prob here with axis_flag="<<axis_flag<<endl;
            comp startingpoint = min;
            double totpoints = 4*points;// + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);

        }
    }



    
}

//this is the same as maker_42, with added
//qvec that goes to right, to qmax 
//again 
void mom_vector_maker_43(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            startingpoint = qvec[qvec.size()-1];

            cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[lowest_p_index] - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[lowest_m_index] - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            if(lowest_index_yval<=imag(startingpoint))
            {
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

                startingpoint = qvec[qvec.size()-1];
                cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points); 

                startingpoint = qvec[qvec.size()-1];
                cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = max;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
            }
            else
            {
                startingpoint = qvec[qvec.size()-1];
                cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
                
                startingpoint = qvec[qvec.size()-1];
                cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = max;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);
            }
        }
        else
        {
            cout<<"deforming the contour downward first"<<endl;
            //comp startingpoint = min;
            //comp delq = abs(startingpoint - rad/2.0)/(points);
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,1);

            //startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            //comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);


            startingpoint = qvec[qvec.size()-1];
            cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            startingpoint = qvec[qvec.size()-1];
            cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            lowest_index_yval = imag(pcutpvec[lowest_p_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
            highest_index_xval = real(pcutpvec[highest_p_index]);
        }
        else
        {
            lowest_index_yval = imag(pcutmvec[lowest_m_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
            highest_index_xval = real(pcutpvec[highest_m_index]);
        }


        if(abs(lowest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour downward first"<<endl;
            comp startingpoint = min;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);

            startingpoint = qvec[qvec.size()-1];
            cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            
            startingpoint = qvec[qvec.size()-1];
            cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);


        }
        else 
        {
            cout<<"taking the straight line"<<endl;
            //cout<<"prob here with axis_flag="<<axis_flag<<endl;
            comp startingpoint = min;
            double totpoints = 4*points;// + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);

        }
    }



    
}

//this is the same as maker_42, with added
//qvec that goes to right, to qmax 
//again, gives the weight vector back,
//piecewise linear line maker used 
void mom_vector_maker_43_with_weights(      vector<comp> &qvec, 
                                            vector<comp> &weights,
                                            comp s,
                                            comp min,
                                            comp max,
                                            double a,
                                            double m,
                                            double eps,
                                            double eps_for_m2k,
                                            double points,
                                            double qvec_r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    comp qc2 = q_c2(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    
    cout<<"axis_flag = "<<axis_flag<<endl;
    if(axis_flag==0)
    cout<<"cut to the right crosses real line"<<endl;
    else
    cout<<"cut to the right doesnot cross real line"<<endl;

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);
            line_maker_with_weights(qvec,weights,startingpoint,-rad/2.0,points/5.0);
            
            cout<<"turningpoint 0 = "<<qvec[0]<<endl;
            startingpoint = qvec[qvec.size()-1];
            cout<<"turningpoint 1 = "<<startingpoint<<endl;
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            cout<<"turningpoint 2 = "<<startingpoint<<endl;
            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[lowest_p_index] - ii*0.1*abs(imag(pcutpvec[lowest_p_index]));

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[lowest_m_index] - ii*0.1*abs(imag(pcutmvec[lowest_m_index]));

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            if(lowest_index_yval<=imag(startingpoint))
            {
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 3 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0); 

                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 4 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = max;
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);
                cout<<"end point 5 = "<<qvec[qvec.size() - 1]<<endl;
            }
            else
            {
                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 3 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
                
                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 4 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = max;
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
                cout<<"end point 5 = "<<qvec[qvec.size() - 1]<<endl;
            }
        }
        else
        {
            cout<<"deforming the contour downward first"<<endl;
            //comp startingpoint = min;
            //comp delq = abs(startingpoint - rad/2.0)/(points);
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,1);

            //startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            //comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.35*abs(imag(qc2));//abs(imag(pcutpvec[lowest_p_index]));
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.35*abs(imag(qc2));//*abs(imag(pcutmvec[lowest_m_index]));
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,2.0*points/5.0);


            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            lowest_index_yval = imag(pcutpvec[lowest_p_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.35*abs(imag(qc2));//*abs(imag(pcutpvec[lowest_p_index]));
            highest_index_xval = real(pcutpvec[highest_p_index]);
        }
        else
        {
            lowest_index_yval = imag(pcutmvec[lowest_m_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.35*abs(imag(qc2));//*abs(imag(pcutmvec[lowest_m_index]));
            highest_index_xval = real(pcutpvec[highest_m_index]);
        }


        if(abs(lowest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour downward first"<<endl;
            comp startingpoint = min;

            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,2.0*points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;

            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);


        }
        else 
        {
            cout<<"taking the straight line"<<endl;
            //cout<<"prob here with axis_flag="<<axis_flag<<endl;
            comp startingpoint = min;
            comp endingpoint = max;
            double totpoints = points;// + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            //mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint,points);
            

        }
    }



    
}


void mom_vector_maker_43_with_weights_with_seba_imsneg(      vector<comp> &qvec, 
                                            vector<comp> &weights,
                                            comp s,
                                            comp min,
                                            comp max,
                                            double a,
                                            double m,
                                            double eps,
                                            double eps_for_m2k,
                                            double points,
                                            double qvec_r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    comp qc2 = q_c2(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    
    cout<<"axis_flag = "<<axis_flag<<endl;
    if(axis_flag==0)
    cout<<"cut to the right crosses real line"<<endl;
    else
    cout<<"cut to the right doesnot cross real line"<<endl;

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            if(abs(imag(s))<0.10)
            {
                cout<<"deforming the contour backward first 1 "<<endl;
                comp q = pmom(s,sigmab(a,m),m);
                comp qc2 = q_c2(s,a,m,eps);
                comp pplus = pcut_plus(-1.0,s,q,m,eps);
                comp pminus = pcut_plus(+1.0,s,q,m,eps);
                cout<<"pplus = "<<pplus<<'\t'<<"pminus = "<<pminus<<endl;
                cout<<"pplus-pminus = "<<0.5*real(pplus - pminus)<<endl;
                comp qc1 = q_c1(s,a,m,eps);
                double rad = real(qc1);

                comp firstnode = min;
                comp secondnode = -2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2*qc2/qc1);
                line_maker_with_weights(qvec,weights,firstnode,secondnode,2.0*points/10.0);

                //secondnode - thirdnode
                comp thirdnode = 0.5*real(pplus - pminus) - abs(qc2)*ii;
                line_maker_with_weights(qvec,weights,secondnode,thirdnode,2.0*points/10.0);

                //thirdnode - fourthnode 
                comp fourthnode = (2.0/3.0)*(1.0 - 1.68*ii)*abs(qc2);//(2.0/3.0)*(1.0 - 2.0*ii)*abs(qc2);
                line_maker_with_weights(qvec,weights,thirdnode,fourthnode,2.0*points/10.0);

                //fourthnode - fifthnode
                comp fifthnode = (3.0/2.0)*abs(qc2);
                line_maker_with_weights(qvec,weights,fourthnode,fifthnode,2.0*points/10.0);

                //fifthnode - sixthnode 
                comp sixthnode = max;
                line_maker_with_weights(qvec,weights,fifthnode,sixthnode,2.0*points/10.0);
                //////////////////////////////////////////////////////
            }
            else 
            {
                cout<<"deforming contour background first 2"<<endl;
                comp startingpoint = min;
                comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
                //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);
                line_maker_with_weights(qvec,weights,startingpoint,-rad/2.0,points/5.0);
            
                cout<<"turningpoint 0 = "<<qvec[0]<<endl;
                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 1 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                //delq = abs(imag(startingpoint) - yval1 + r)/points;
                //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
                comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 2 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                double lowest_index_yval = 0.0;
            
                if(flag_pcut==0)
                {
                    lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                    endingpoint_for_slanted_maker = pcutpvec[lowest_p_index] - ii*0.1*abs(imag(pcutpvec[lowest_p_index]));

                    highest_index_xval = real(pcutpvec[highest_p_index]);
                    //cout<<"highest x val = "<<highest_index_xval<<endl;
                }
                else
                {
                    lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                    endingpoint_for_slanted_maker = pcutmvec[lowest_m_index] - ii*0.1*abs(imag(pcutmvec[lowest_m_index]));

                    highest_index_xval = real(pcutmvec[highest_m_index]);
                    //cout<<"highest x val = "<<highest_index_xval<<endl;
                }

                if(lowest_index_yval<=imag(startingpoint))
                {
                    line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

                    startingpoint = qvec[qvec.size()-1];
                    cout<<"turningpoint 3 = "<<startingpoint<<endl;
                    //cout<<"last point after step = "<<startingpoint<<endl;
                    endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                    line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0); 

                    startingpoint = qvec[qvec.size()-1];
                    cout<<"turningpoint 4 = "<<startingpoint<<endl;
                    //cout<<"last point after step = "<<startingpoint<<endl;
                    endingpoint_for_slanted_maker = max;
                    line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);
                    cout<<"end point 5 = "<<qvec[qvec.size() - 1]<<endl;
                }
                else
                {
                    startingpoint = qvec[qvec.size()-1];
                    cout<<"turningpoint 3 = "<<startingpoint<<endl;
                    //cout<<"last point after step = "<<startingpoint<<endl;
                    endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                    line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
                
                    startingpoint = qvec[qvec.size()-1];
                    cout<<"turningpoint 4 = "<<startingpoint<<endl;
                    //cout<<"last point after step = "<<startingpoint<<endl;
                    endingpoint_for_slanted_maker = max;
                    line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
                    cout<<"end point 5 = "<<qvec[qvec.size() - 1]<<endl;
                }
            }
            
        }
        else
        {
            if(abs(imag(s))>0.0230)
            {    
                cout<<"deforming the contour downward first"<<endl;
            

                comp startingpoint = min;//qvec[qvec.size()-1];
                comp endingpoint_for_slanted_maker;
                if(flag_pcut==0)
                {
                    endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.35*abs(imag(qc2));//abs(imag(pcutpvec[lowest_p_index]));
                    highest_index_xval = real(pcutpvec[highest_p_index]);
                }
                else
                {
                    endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.35*abs(imag(qc2));//*abs(imag(pcutmvec[lowest_m_index]));
                    highest_index_xval = real(pcutpvec[highest_m_index]);
                }
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,2.0*points/5.0);


                startingpoint = qvec[qvec.size()-1];
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);

                startingpoint = qvec[qvec.size()-1];
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = max;
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
            }
            else
            {
                cout<<"deforming the contour backward first 1 "<<endl;
                comp q = pmom(s,sigmab(a,m),m);
                comp qc2 = q_c2(s,a,m,eps);
                comp pplus = pcut_plus(-1.0,s,q,m,eps);
                comp pminus = pcut_plus(+1.0,s,q,m,eps);
                cout<<"pplus = "<<pplus<<'\t'<<"pminus = "<<pminus<<endl;
                cout<<"pplus-pminus = "<<0.5*real(pplus - pminus)<<endl;
                comp qc1 = q_c1(s,a,m,eps);
                double rad = real(qc1);

                comp firstnode = min;
                comp secondnode = -2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2*qc2/qc1);
                line_maker_with_weights(qvec,weights,firstnode,secondnode,2.0*points/10.0);

                //secondnode - thirdnode
                comp thirdnode = 0.5*real(pplus - pminus) - abs(qc2)*ii;
                line_maker_with_weights(qvec,weights,secondnode,thirdnode,2.0*points/10.0);

                //thirdnode - fourthnode 
                comp fourthnode = (2.0/3.0)*(1.0 - 1.68*ii)*abs(qc2);//(2.0/3.0)*(1.0 - 2.0*ii)*abs(qc2);
                line_maker_with_weights(qvec,weights,thirdnode,fourthnode,2.0*points/10.0);

                //fourthnode - fifthnode
                comp fifthnode = (3.0/2.0)*abs(qc2);
                line_maker_with_weights(qvec,weights,fourthnode,fifthnode,2.0*points/10.0);

                //fifthnode - sixthnode 
                comp sixthnode = max;
                line_maker_with_weights(qvec,weights,fifthnode,sixthnode,2.0*points/10.0);
                //////////////////////////////////////////////////////
            }
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            lowest_index_yval = imag(pcutpvec[lowest_p_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.35*abs(imag(qc2));//*abs(imag(pcutpvec[lowest_p_index]));
            highest_index_xval = real(pcutpvec[highest_p_index]);
        }
        else
        {
            lowest_index_yval = imag(pcutmvec[lowest_m_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.35*abs(imag(qc2));//*abs(imag(pcutmvec[lowest_m_index]));
            highest_index_xval = real(pcutpvec[highest_m_index]);
        }


        if(abs(lowest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour downward first"<<endl;
            comp startingpoint = min;

            if(imag(endingpoint_for_slanted_maker)<0.0)
            {
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,2.0*points/5.0);

                startingpoint = qvec[qvec.size()-1];
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);

            
                startingpoint = qvec[qvec.size()-1];
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = max;

                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
            }
            else 
            {
                endingpoint_for_slanted_maker = max;
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points);
            }


        }
        else 
        {
            cout<<"taking the straight line"<<endl;
            //cout<<"prob here with axis_flag="<<axis_flag<<endl;
            comp startingpoint = min;
            comp endingpoint = max;
            double totpoints = points;// + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            //mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint,points);
            

        }
    }



    
}


void mom_vector_maker_43_1(   vector<comp> &qvec,
                              int points   )
{
    
    comp point1 = {1.0e-7,1.0e-7};
    comp point2 = {-0.0619805,  -0.0619805};
    comp point3 = {0.00950383,  -0.13148};
    comp point4 = {0.105184,  -0.170925};
    comp point5 = {0.197221,  0};
    comp point6 = {1.2958,  -0.00475807};


    line_maker(qvec,point1,point2,points/5);
    line_maker(qvec,point2,point3,points/5);
    line_maker(qvec,point3,point4,points/5);
    line_maker(qvec,point4,point5,points/5);
    line_maker(qvec,point5,point6,points/5);

    
}


//This one has tags for taking the contour 
//under the second sheet of OPE
void mom_vector_maker_44(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            startingpoint = qvec[qvec.size()-1];

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);


            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);


            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points);
        
        }
        else 
        {
            cout<<"taking the straight line"<<endl;
            //cout<<"prob here with axis_flag="<<axis_flag<<endl;
            comp startingpoint = min;
            double totpoints = 4*points;// + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);

        }
    }



    
}


//This one has tags for taking the contour 
//under the second sheet of OPE, 
//with 2 tags that tells where ope should be taken
//to the unphysical sheet, the cut still goes through
//the second sheet of ope even when the ope cuts 
//are not crossing the real q axis.
void mom_vector_maker_46(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[6000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[6000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        

        }
    }



    
}


//This one has tags for taking the contour 
//under the second sheet of OPE, 
//with 2 tags that tells where ope should be taken
//to the unphysical sheet, the cut still goes through
//the second sheet of ope even when the ope cuts 
//are not crossing the real q axis. Also this 
//one goes through the pcuts only once 
void mom_vector_maker_47(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
        

        }
    }



    
}

void mom_vector_maker_47_imspos(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
        

        }
    }



    
}


void mom_vector_maker_47_imspos1em5(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        

        }
    }



    
}



void mom_vector_maker_47_2_A1(  int somepoint,  
                                vector<comp> &qvec,   
                                comp s,
                                comp min,
                                comp max,
                                double a,
                                double m,
                                double eps,
                                double eps_for_m2k,
                                double points,
                                double qvec_r,
                                int tag1,
                                int tag2    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    //int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index - somepoint];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index - somepoint];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        

        }
    }



    
}

//Same as contour_47 but we can take 
//the first line to left as much as we want
void mom_vector_maker_47_1( double somepoint,
                            vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r,
                            int &tag1,
                            int &tag2    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    //int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - somepoint*rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = somepoint*pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = somepoint*pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = somepoint*pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = somepoint*pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = somepoint*pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = somepoint*pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
        

        }
    }



    
}


//Same as contour_47 but we can take 
//the line that goes through the second OPE 
//sheet to bend to left as much as we want
void mom_vector_maker_47_2( double somepoint,
                            vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r,
                            int &tag1,
                            int &tag2    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    //int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                //endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;
                endingpoint_for_slanted_maker = pcutmvec[8000];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                //endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;
                endingpoint_for_slanted_maker = pcutpvec[8000];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[9700];;//pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[9700];//pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            comp end_temp = endingpoint_for_slanted_maker;
            startingpoint = qvec[qvec.size()-1];
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = -somepoint*1.5*pcutpvec[1000];
            }
            else 
            {
                endingpoint_for_slanted_maker = -somepoint*1.5*pcutmvec[1000];
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            endingpoint_for_slanted_maker = end_temp;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            //startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker = end_temp;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);



            /////////////////////////////////////////
            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
        

        }
    }



    
}

void mom_vector_maker_47_2_1( double somepoint,
                            vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r,
                            int &tag1,
                            int &tag2    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    //int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                //endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;
                endingpoint_for_slanted_maker = pcutmvec[8000];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                //endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;
                endingpoint_for_slanted_maker = pcutpvec[8000];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            comp temp_real = real(startingpoint);

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[9700];;//pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[9700];//pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            comp temp_imag = imag(endingpoint_for_slanted_maker) + somepoint;
            comp end_temp = endingpoint_for_slanted_maker;
            endingpoint_for_slanted_maker = temp_real + ii*temp_imag;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            endingpoint_for_slanted_maker = end_temp;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            //startingpoint = qvec[qvec.size()-1];
            //if(flag_pcut==0)
            //{
            //    endingpoint_for_slanted_maker = -somepoint*1.5*pcutpvec[1000];
            //}
            //else 
            //{
            //    endingpoint_for_slanted_maker = -somepoint*1.5*pcutmvec[1000];
            //}
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            //startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker = end_temp;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            //startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker = end_temp;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);



            /////////////////////////////////////////
            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
        

        }
    }



    
}

void mom_vector_maker_47_2_2( double somepoint,
                            vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r,
                            int &tag1,
                            int &tag2    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    //int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                //endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;
                endingpoint_for_slanted_maker = pcutmvec[8000];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                //endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;
                endingpoint_for_slanted_maker = pcutpvec[8000];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            double temp_real = real(startingpoint) - somepoint;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[9700];;//pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[9700];//pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            double temp_imag = imag(endingpoint_for_slanted_maker);
            comp end_temp = endingpoint_for_slanted_maker;
            endingpoint_for_slanted_maker = temp_real + ii*temp_imag;
        
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            endingpoint_for_slanted_maker = end_temp;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            //startingpoint = qvec[qvec.size()-1];
            //if(flag_pcut==0)
            //{
            //    endingpoint_for_slanted_maker = -somepoint*1.5*pcutpvec[1000];
            //}
            //else 
            //{
            //    endingpoint_for_slanted_maker = -somepoint*1.5*pcutmvec[1000];
            //}
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            //startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker = end_temp;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            //startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker = end_temp;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);



            /////////////////////////////////////////
            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
        

        }
    }



    
}


void mom_vector_maker_47_2_3( double somepoint,
                            vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r,
                            int &tag1,
                            int &tag2    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    //int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size();
            delq = abs(startingpoint - rad/4.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            
            double temp_real = real(startingpoint);// - somepoint;
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            //(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
            comp endingpoint_for_slanted_maker;
            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[9700];;//pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[9700];//pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            double temp_imag = imag(endingpoint_for_slanted_maker);
            comp end_temp = endingpoint_for_slanted_maker;
            endingpoint_for_slanted_maker = temp_real + ii*temp_imag;
        
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            endingpoint_for_slanted_maker = end_temp;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            /////////////////////////////////////////
            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
        

        }
    }



    
}

void mom_vector_maker_47_2_3_1( double somepoint,
                            vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r,
                            int &tag1,
                            int &tag2    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    //int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {   
            comp endingpoint_for_slanted_maker;
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            if(flag_pcut==0)
            {
                //lowest_index_yval = imag(pcutmvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutmvec[pcutmvec.size()/2];;//pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                //highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                //lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutpvec[pcutpvec.size()/2];//pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                //highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            } 
            //comp delq = abs(startingpoint - rad)/(points/5.0);
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size();
            delq = abs(startingpoint - rad/4.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            
            double temp_real = real(startingpoint);// - somepoint;
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            //(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
            
            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[9590];;//pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[9590];//pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            double temp_imag = imag(endingpoint_for_slanted_maker);
            comp end_temp = endingpoint_for_slanted_maker;
            endingpoint_for_slanted_maker = temp_real + ii*temp_imag;
        
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            endingpoint_for_slanted_maker = end_temp;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            /////////////////////////////////////////
            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
        

        }
    }



    
}




//Same as contour_47 but we can take 
//the OFFSET_R to be as long as we want
void mom_vector_maker_47_3( double somepoint,
                            vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r,
                            int &tag1,
                            int &tag2    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    //int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            //comp end_temp = endingpoint_for_slanted_maker;
            //startingpoint = qvec[qvec.size()-1];
            //if(flag_pcut==0)
            //{
            //    endingpoint_for_slanted_maker = -somepoint*1.5*pcutpvec[500];
            //}
            //else 
            //{
            //    endingpoint_for_slanted_maker = -somepoint*1.5*pcutmvec[500];
            //}
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            //startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker = end_temp;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            //startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker = end_temp;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);



            /////////////////////////////////////////
            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + somepoint*qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
        

        }
    }



    
}


void mom_vector_maker_opposite_1(   vector<comp> &qvec,   
                                    comp s,
                                    comp min,
                                    comp max,
                                    double a,
                                    double m,
                                    double eps,
                                    double points1    )
{
    comp ii = {0.0,1.0};
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qmns = q_minus(s,a,m,eps);
    comp qc1 = q_c1(s,a,m,eps);
    comp delq = {0.0,0.0};

    if(real(rad)>0.0)
    {
        /*  this is the initial straight line along the negative imaginary axis */
        comp startingpoint = min;
        delq = abs((abs(imag(qmns)) + eps)-min)/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,delq,points1,2);
        
        /*  this is the contour going right in the imaginary axis   */
        startingpoint = qvec[qvec.size()-1];
        delq = abs((abs(real(qc1)) + eps) - real(startingpoint))/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,delq,points1,1);

        /*  This contour now goes up in the imaginary axis                       */
        startingpoint = qvec[qvec.size()-1];
        delq = abs(imag(startingpoint))/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,-delq,points1,2);

        /*  This contour now goes right in the real axis                       */
        startingpoint = qvec[qvec.size()-1];
        delq = abs(real(max) - real(startingpoint))/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,delq,points1,1);

        
        /* print the size of the created sigma_vector */
        cout<<"sigvec created with size = "<<qvec.size()<<endl;

    }
    else
    {
        double length = real(max - min);
        comp delsig = length/(4.0*points1) ;

        int totpoints = 4.0*points1;
        for(int i=1;i<=totpoints;++i)
        {
            comp qmom = min + ii*eps + ((comp)i)*delsig;

            qvec.push_back(qmom);
        }
    }

}

comp phibthreshold( double a,
                    double m    )
{
    return pow(sqrt(sigmab(a,m)) + m,2.0);
}





void opeplot()
{
    ofstream fout;

    double a = 16.0;
    double m = 1.0;
    double s = 8.3;
    double eps = 1.0e-10;
    comp scomp = (comp) s;
    comp sigmapprime = sigmab(a,m);
     

    string filename = "OPE_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename.c_str());
    for(double sigmap=3.3;sigmap<=3.999;sigmap=sigmap + abs(3.3-3.999)/5000.0)
    {
        double someeps = 1.0e-6;
        //comp someeps = eps*eps_energy_factor_plus(s,sigmapprime,m);
        comp sigmapcomp = (comp) sigmap + ii*someeps;

        comp ope = GS_1(scomp,sigmapcomp,sigmapprime,m,eps);

        fout<<sigmap<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
        cout<<"sigk:"<<sigmap<<'\t'<<"eps:"<<someeps<<'\t'<<"ope:"<<real(ope)<<'\t'<<imag(ope)<<endl;
    }
    cout<<setprecision(15)<<sigc1(s,sigmab(a,m),m)<<'\t'
        <<sigpplus_witheps(s,sigmapprime,m,0.0)<<'\t'
        <<sigc2(s,sigmapprime,m)<<'\t'
        <<sigpminus_witheps(s,sigmapprime,m,0.0)<<endl;

    fout.close();

    
}



void opeplotcontour()
{
    ofstream fout;

    double a = 16.0;
    double m = 1.0;
    double s = 8.72;
    //comp s = 8.72 + ii*0.01
    double eps = 1.0e-7;
    comp scomp = (comp) s;
    comp sigmapprime = sigmab(a,m);
    double sigpRinitial = 3.6;
    double sigpRfinal = 4.0;
    double delsigpR = abs(sigpRinitial-sigpRfinal)/250.0;
    double sigpIinitial = -0.1;
    double sigpIfinal = 0.1;
    double delsigpI = abs(sigpIinitial - sigpIfinal)/250.0;
    

     

    string filename = "OPEcontour_a_" + to_string((int)a) + "_s+i_" + to_string(real(s)) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    int size = 250;
    int isize = size;
    int jsize = size;

    fout.open(filename.c_str());
    for(double sigmapR=sigpRinitial;sigmapR<sigpRfinal;sigmapR=sigmapR + delsigpR)
    {
        
        for(double sigmapI=sigpIinitial;sigmapI<sigpIfinal;sigmapI=sigmapI + delsigpI)
        //double someeps = 1.0e-3;
        {
            comp sigmapcomp = (comp) sigmapR + ii*sigmapI;

            comp ope = GS(scomp,sigmapcomp,sigmapprime,m,eps);

            fout<<setprecision(16)<<sigmapR<<'\t'<<setprecision(16)<<sigmapI<<'\t'<<setprecision(16)<<real(ope)<<'\t'<<setprecision(16)<<imag(ope)<<endl;
            cout<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
            //printf("%.16f\t %.16f\t %.16f\t %.16f\n",sigmapR,sigmapI,real(ope),imag(ope));
            
        }
        fout<<endl;
    }   

    fout.close();

    string filename1 = "threshOPEcontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename1.c_str());

    fout<<real(sigpplus(s,sigmab(a,m),m))<<'\t'<<real(sigpminus(s,sigmab(a,m),m))<<'\t'<<real(sigc1(s,sigmab(a,m),m))<<'\t'<<real(sigc2(s,sigmab(a,m),m))<<'\t'<<pow(sqrt(s)-m,2.0)<<'\t'<<0<<'\t'<<0<<endl;

    fout.close();

    string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename2.c_str());

    vector<comp> sigvec;
    //sigma_vector_maker(sigvec,s,0.0,sigmax(s,m),a,m,1000.0,1.0e-2);
    for(int i=0;i<sigvec.size();++i)
    //fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;

    fout.close();

}

void opeplotcontour_withintloop()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    //double s = 8.72;
    double eps = 0.09;
    double imS = 0.001;
    double sreal = 8.72;
    comp s = sreal + ii*imS;
    
    comp scomp = (comp) s;//+ ii*imS;
    comp sigmapprime = sigmab(a,m) ;//+ ii*0.02;

    double sigpRinitial = 3.6;
    double sigpRfinal = 4.0;
    double delsigpR = abs(sigpRinitial-sigpRfinal)/250.0;
    double sigpIinitial = -0.1;
    double sigpIfinal = 0.1;
    double delsigpI = abs(sigpIinitial - sigpIfinal)/250.0;
    
    //double sreal = real(s);

     

    string filename =   "OPEcontour_a_" + to_string((int)a) 
                        + "_sreal_" + to_string(sreal) 
                        + "_ims_" + to_string(imS)
                        + "_eps_" + to_string(eps)
                        + "_sigpsb_" + to_string(real(sigmab(a,m))) 
                        + ".dat";

    int size = 250;
    int isize = size;
    int jsize = size;

    fout.open(filename.c_str());

    //for(double sigmapR=sigpRinitial;sigmapR<sigpRfinal;sigmapR=sigmapR + delsigpR)
    for(int i=0;i<isize;++i)
    {
        
        //for(double sigmapI=sigpIinitial;sigmapI<sigpIfinal;sigmapI=sigmapI + delsigpI)
        for(int j=0;j<jsize;++j)
        //double someeps = 1.0e-3;
        {
            comp sigmapR = (comp)sigpRinitial + (comp)i*delsigpR;
            comp sigmapI = (comp)sigpIinitial + (comp)j*delsigpI; 
            comp sigmapcomp = (comp) sigmapR + ii*sigmapI;

            comp ope = GS_1(scomp,sigmapcomp,sigmapprime,m,eps);
            //comp ope = GS_complex_s(sreal,sigmapcomp,sigmapprime,m,eps,imS);

            fout<<setprecision(16)<<real(sigmapR)<<'\t'<<setprecision(16)<<real(sigmapI)<<'\t'<<setprecision(16)<<real(ope)<<'\t'<<setprecision(16)<<imag(ope)<<endl;
            cout<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
            //printf("%.16f\t %.16f\t %.16f\t %.16f\n",sigmapR,sigmapI,real(ope),imag(ope));
            
        }
        fout<<endl;
    }   

    fout.close();

    string filename1 = "thresh_" + filename;//OPEcontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename1.c_str());

    //fout<<real(sigpplus(s,sigmab(a,m),m))<<'\t'<<real(sigpminus(s,sigmab(a,m),m))<<'\t'<<real(sigc1(s,sigmab(a,m),m))<<'\t'<<real(sigc2(s,sigmab(a,m),m))<<'\t'<<pow(sqrt(s)-m,2.0)<<'\t'<<0<<'\t'<<0<<endl;
    fout<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<real((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\t'
        <<imag((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\t'
        <<0<<endl;
    cout<<"resigplus:"<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"imsigplus:"<<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"resigminus:"<<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"imsigminus:"<<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"resigmax:"<<real((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\n'
        <<"imsigmax:"<<imag((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\n'
        <<0<<endl;
    fout.close();

    //string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";
    string filename2 = "sigvec_" + filename;
    fout.open(filename2.c_str());

    vector<comp> sigvec;
    sigma_vector_maker_6(sigvec,s,0.0,sigmax(s,m),a,m,1.0e-8,1000,200,1.0);
    for(int i=0;i<sigvec.size();++i)
    fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;

    fout.close();

    cout<<"file1:"<<filename<<'\n'
        <<"file2:"<<filename1<<'\n'
        <<"file3:"<<filename2<<endl;

}

void opeplotcontour_momrep_withintloop()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    //double s = 8.72;
    double eps = 0.001;
    double imS = 0.0;
    double sreal = 8.72;
    comp s = sreal + ii*imS;
    
    comp scomp = (comp) s;//+ ii*imS;
    comp sigb = sigmab(a,m) ;//+ ii*0.02;
    comp q = pmom(s,sigb,m);

    double pRinitial = -0.3;
    double pRfinal = 0.3;
    double delpR = abs(pRinitial-pRfinal)/250.0;
    double pIinitial = -0.3;
    double pIfinal = 0.3;
    double delpI = abs(pIinitial - pIfinal)/250.0;
    
    //double sreal = real(s);

     

    string filename =   "OPEcontour_momrep_a_" + to_string((int)a) 
                        + "_sreal_" + to_string(sreal) 
                        + "_ims_" + to_string(imS)
                        + "_eps_" + to_string(eps)
                        + "_sigpsb_" + to_string(real(sigmab(a,m))) 
                        + ".dat";

    int size = 250;
    int isize = size;
    int jsize = size;

    fout.open(filename.c_str());

    //for(double sigmapR=sigpRinitial;sigmapR<sigpRfinal;sigmapR=sigmapR + delsigpR)
    for(int i=0;i<isize;++i)
    {
        
        //for(double sigmapI=sigpIinitial;sigmapI<sigpIfinal;sigmapI=sigmapI + delsigpI)
        for(int j=0;j<jsize;++j)
        //double someeps = 1.0e-3;
        {
            comp pR = (comp)pRinitial + (comp)i*delpR;
            comp pI = (comp)pIinitial + (comp)j*delpI; 
            comp pcomp = (comp) pR + ii*pI;

            comp ope = GS_pk(s,q,pcomp,m,eps);
            //GS_1(scomp,sigmapcomp,sigmapprime,m,eps);
            //comp ope = GS_complex_s(sreal,sigmapcomp,sigmapprime,m,eps,imS);

            fout<<setprecision(16)<<real(pR)<<'\t'<<setprecision(16)<<real(pI)<<'\t'<<setprecision(16)<<real(ope)<<'\t'<<setprecision(16)<<imag(ope)<<endl;
            cout<<"q:"<<q<<endl;
            cout<<"pR:"<<pR<<'\t'<<"pI:"<<pI<<'\t'<<"OPE:"<<real(ope)<<'\t'<<"+i "<<imag(ope)<<endl;
            //printf("%.16f\t %.16f\t %.16f\t %.16f\n",sigmapR,sigmapI,real(ope),imag(ope));
            
        }
        fout<<endl;
    }   

    fout.close();

    /*

    string filename1 = "thresh_" + filename;//OPEcontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename1.c_str());

    //fout<<real(sigpplus(s,sigmab(a,m),m))<<'\t'<<real(sigpminus(s,sigmab(a,m),m))<<'\t'<<real(sigc1(s,sigmab(a,m),m))<<'\t'<<real(sigc2(s,sigmab(a,m),m))<<'\t'<<pow(sqrt(s)-m,2.0)<<'\t'<<0<<'\t'<<0<<endl;
    fout<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<real((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\t'
        <<imag((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\t'
        <<0<<endl;
    cout<<"resigplus:"<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"imsigplus:"<<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"resigminus:"<<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"imsigminus:"<<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"resigmax:"<<real((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\n'
        <<"imsigmax:"<<imag((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\n'
        <<0<<endl;
    fout.close();

    //string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";
    string filename2 = "sigvec_" + filename;
    fout.open(filename2.c_str());

    vector<comp> sigvec;
    sigma_vector_maker_6(sigvec,s,0.0,sigmax(s,m),a,m,1.0e-8,1000,200,1.0);
    for(int i=0;i<sigvec.size();++i)
    fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;

    fout.close();

    cout<<"file1:"<<filename<<'\n'
        <<"file2:"<<filename1<<'\n'
        <<"file3:"<<filename2<<endl;
    */
}

void opeplotcontour_momrep_withintloop_withepsloop()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    //double s = 8.72;
    double eps = 0.001;
    double imS = -0.1;
    double sreal = 8.72;
    comp s = sreal + ii*imS;
    
    comp scomp = (comp) s;//+ ii*imS;
    comp sigb = sigmab(a,m) ;//+ ii*0.02;
    comp q = pmom(s,sigb,m);

    double pRinitial = -0.3;
    double pRfinal = 0.3;
    double delpR = abs(pRinitial-pRfinal)/250.0;
    double pIinitial = -0.3;
    double pIfinal = 0.3;
    double delpI = abs(pIinitial - pIfinal)/250.0;
    
    //double sreal = real(s);

    for(eps=0.01;eps<=0.1;eps=eps+0.01)
    { 

    string filename =   "OPEcontour_momrep_a_" + to_string((int)a) 
                        + "_sreal_" + to_string(sreal) 
                        + "_ims_" + to_string(imS)
                        + "_eps_" + to_string(eps)
                        + "_sigpsb_" + to_string(real(sigmab(a,m))) 
                        + ".dat";

    int size = 250;
    int isize = size;
    int jsize = size;

    fout.open(filename.c_str());

    //for(double sigmapR=sigpRinitial;sigmapR<sigpRfinal;sigmapR=sigmapR + delsigpR)
    for(int i=0;i<isize;++i)
    {
        
        //for(double sigmapI=sigpIinitial;sigmapI<sigpIfinal;sigmapI=sigmapI + delsigpI)
        for(int j=0;j<jsize;++j)
        //double someeps = 1.0e-3;
        {
            comp pR = (comp)pRinitial + (comp)i*delpR;
            comp pI = (comp)pIinitial + (comp)j*delpI; 
            comp pcomp = (comp) pR + ii*pI;

            comp ope = GS_pk(s,q,pcomp,m,eps);
            //GS_1(scomp,sigmapcomp,sigmapprime,m,eps);
            //comp ope = GS_complex_s(sreal,sigmapcomp,sigmapprime,m,eps,imS);

            fout<<setprecision(16)<<real(pR)<<'\t'<<setprecision(16)<<real(pI)<<'\t'<<setprecision(16)<<real(ope)<<'\t'<<setprecision(16)<<imag(ope)<<endl;
            cout<<"kmin+:"<<pmom(s,(sqrt(s)-m)*(sqrt(s)-m),m)<<endl;
            cout<<"kmin-:"<<pmom(s,(mysqrt(s)-m)*(mysqrt(s)-m),m)<<endl;
            cout<<"q:"<<q<<endl;
            cout<<"pR:"<<pR<<'\t'<<"pI:"<<pI<<'\t'<<"OPE:"<<real(ope)<<'\t'<<"+i "<<imag(ope)<<endl;
            //printf("%.16f\t %.16f\t %.16f\t %.16f\n",sigmapR,sigmapI,real(ope),imag(ope));
            
        }
        fout<<endl;
    }   

    fout.close();

    }
    /*

    string filename1 = "thresh_" + filename;//OPEcontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename1.c_str());

    //fout<<real(sigpplus(s,sigmab(a,m),m))<<'\t'<<real(sigpminus(s,sigmab(a,m),m))<<'\t'<<real(sigc1(s,sigmab(a,m),m))<<'\t'<<real(sigc2(s,sigmab(a,m),m))<<'\t'<<pow(sqrt(s)-m,2.0)<<'\t'<<0<<'\t'<<0<<endl;
    fout<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<real((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\t'
        <<imag((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\t'
        <<0<<endl;
    cout<<"resigplus:"<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"imsigplus:"<<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"resigminus:"<<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"imsigminus:"<<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"resigmax:"<<real((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\n'
        <<"imsigmax:"<<imag((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\n'
        <<0<<endl;
    fout.close();

    //string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";
    string filename2 = "sigvec_" + filename;
    fout.open(filename2.c_str());

    vector<comp> sigvec;
    sigma_vector_maker_6(sigvec,s,0.0,sigmax(s,m),a,m,1.0e-8,1000,200,1.0);
    for(int i=0;i<sigvec.size();++i)
    fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;

    fout.close();

    cout<<"file1:"<<filename<<'\n'
        <<"file2:"<<filename1<<'\n'
        <<"file3:"<<filename2<<endl;
    */
}

void m2kplotcontour_momrep_withintloop_withepsloop()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    //double s = 8.72;
    double eps = 0.001;
    double imS = 0.0;//-0.1;
    double sreal = 8.72;
    comp s = sreal + ii*imS;
    
    comp scomp = (comp) s;//+ ii*imS;
    comp sigb = sigmab(a,m) ;//+ ii*0.02;
    comp q = pmom(s,sigb,m);
    double N = 1000.0;

    double pRinitial = -2.0;
    double pRfinal = 2.0;
    double delpR = abs(pRinitial-pRfinal)/N;
    double pIinitial = -0.5;
    double pIfinal = 0.5;
    double delpI = abs(pIinitial - pIfinal)/N;
    
    //double sreal = real(s);

    for(eps=0.01;eps<=0.01;eps=eps+0.01)
    { 

    string filename =   "M2kcontour_momrep_a_" + to_string((int)a) 
                        + "_sreal_" + to_string(sreal) 
                        + "_ims_" + to_string(imS)
                        + "_eps_" + to_string(eps)
                        + "_sigpsb_" + to_string(real(sigmab(a,m))) 
                        + ".dat";

    int size = (int) N;
    int isize = size;
    int jsize = size;

    fout.open(filename.c_str());

    //for(double sigmapR=sigpRinitial;sigmapR<sigpRfinal;sigmapR=sigmapR + delsigpR)
    for(int i=0;i<isize;++i)
    {
        
        //for(double sigmapI=sigpIinitial;sigmapI<sigpIfinal;sigmapI=sigmapI + delsigpI)
        for(int j=0;j<jsize;++j)
        //double someeps = 1.0e-3;
        {
            comp pR = (comp)pRinitial + (comp)i*delpR;
            comp pI = (comp)pIinitial + (comp)j*delpI; 
            comp pcomp = (comp) pR + ii*pI;

            comp sigp = sigma_p(scomp,pcomp,m);

            comp ope = M2kfunc(a,sigp,m,0.0);
            //GS_pk(s,q,pcomp,m,eps);
            //GS_1(scomp,sigmapcomp,sigmapprime,m,eps);
            //comp ope = GS_complex_s(sreal,sigmapcomp,sigmapprime,m,eps,imS);

            fout<<setprecision(16)<<real(pR)<<'\t'<<setprecision(16)<<real(pI)<<'\t'<<setprecision(16)<<real(ope)<<'\t'<<setprecision(16)<<imag(ope)<<endl;
            cout<<"kmin+:"<<pmom(s,(sqrt(s)-m)*(sqrt(s)-m),m)<<endl;
            cout<<"kmin-:"<<pmom(s,(mysqrt(s)-m)*(mysqrt(s)-m),m)<<endl;
            cout<<"q:"<<q<<endl;
            cout<<"sigp:"<<sigp<<endl;
            cout<<"pR:"<<pR<<'\t'<<"pI:"<<pI<<'\t'<<"OPE:"<<real(ope)<<'\t'<<"+i "<<imag(ope)<<endl;
            //printf("%.16f\t %.16f\t %.16f\t %.16f\n",sigmapR,sigmapI,real(ope),imag(ope));
            
        }
        fout<<endl;
    }   

    fout.close();

    cout<<"bc left plus:"<<M2kbranchcut_left_momrep_plus(scomp,m)<<endl;
    cout<<"bc left minus:"<<M2kbranchcut_left_momrep_minus(scomp,m)<<endl;
    cout<<"bc right plus:"<<M2kbranchcut_right_momrep_plus(scomp,m)<<endl;
    cout<<"bc right minus:"<<M2kbranchcut_right_momrep_minus(scomp,m)<<endl;
    
    cout<<"filename:"<<filename<<endl;

    }
    /*

    string filename1 = "thresh_" + filename;//OPEcontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename1.c_str());

    //fout<<real(sigpplus(s,sigmab(a,m),m))<<'\t'<<real(sigpminus(s,sigmab(a,m),m))<<'\t'<<real(sigc1(s,sigmab(a,m),m))<<'\t'<<real(sigc2(s,sigmab(a,m),m))<<'\t'<<pow(sqrt(s)-m,2.0)<<'\t'<<0<<'\t'<<0<<endl;
    fout<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<real((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\t'
        <<imag((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\t'
        <<0<<endl;
    cout<<"resigplus:"<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"imsigplus:"<<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"resigminus:"<<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"imsigminus:"<<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"resigmax:"<<real((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\n'
        <<"imsigmax:"<<imag((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\n'
        <<0<<endl;
    fout.close();

    //string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";
    string filename2 = "sigvec_" + filename;
    fout.open(filename2.c_str());

    vector<comp> sigvec;
    sigma_vector_maker_6(sigvec,s,0.0,sigmax(s,m),a,m,1.0e-8,1000,200,1.0);
    for(int i=0;i<sigvec.size();++i)
    fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;

    fout.close();

    cout<<"file1:"<<filename<<'\n'
        <<"file2:"<<filename1<<'\n'
        <<"file3:"<<filename2<<endl;
    */
}

void omegaplotcontour_momrep_withintloop_withepsloop()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    //double s = 8.72;
    double eps = 0.001;
    double imS = 0.0;//-0.1;
    double sreal = 8.72;
    comp s = sreal + ii*imS;
    
    comp scomp = (comp) s;//+ ii*imS;
    comp sigb = sigmab(a,m) ;//+ ii*0.02;
    comp q = pmom(s,sigb,m);
    double N = 1000.0;

    double pRinitial = -2.0;
    double pRfinal = 2.0;
    double delpR = abs(pRinitial-pRfinal)/N;
    double pIinitial = -2.0;
    double pIfinal = 2.0;
    double delpI = abs(pIinitial - pIfinal)/N;
    
    //double sreal = real(s);

    for(eps=0.01;eps<=0.01;eps=eps+0.01)
    { 

    string filename =   "omegacontour_momrep_a_" + to_string((int)a) 
                        + "_sreal_" + to_string(sreal) 
                        + "_ims_" + to_string(imS)
                        + "_eps_" + to_string(eps)
                        + "_sigpsb_" + to_string(real(sigmab(a,m))) 
                        + ".dat";

    int size = (int) N;
    int isize = size;
    int jsize = size;

    fout.open(filename.c_str());

    //for(double sigmapR=sigpRinitial;sigmapR<sigpRfinal;sigmapR=sigmapR + delsigpR)
    for(int i=0;i<isize;++i)
    {
        
        //for(double sigmapI=sigpIinitial;sigmapI<sigpIfinal;sigmapI=sigmapI + delsigpI)
        for(int j=0;j<jsize;++j)
        //double someeps = 1.0e-3;
        {
            comp pR = (comp)pRinitial + (comp)i*delpR;
            comp pI = (comp)pIinitial + (comp)j*delpI; 
            comp pcomp = (comp) pR + ii*pI;

            comp sigp = sigma_p(scomp,pcomp,m);

            double pi = acos(-1.0);

            comp ope = omega_comp(pcomp,m);
            //(pcomp*pcomp/(pow(2.0*pi,2.0)omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps);
            //M2kfunc(a,sigp,m,0.0);
            //GS_pk(s,q,pcomp,m,eps);
            //GS_1(scomp,sigmapcomp,sigmapprime,m,eps);
            //comp ope = GS_complex_s(sreal,sigmapcomp,sigmapprime,m,eps,imS);

            fout<<setprecision(16)<<real(pR)<<'\t'<<setprecision(16)<<real(pI)<<'\t'<<setprecision(16)<<real(ope)<<'\t'<<setprecision(16)<<imag(ope)<<endl;
            cout<<"kmin+:"<<pmom(s,(sqrt(s)-m)*(sqrt(s)-m),m)<<endl;
            cout<<"kmin-:"<<pmom(s,(mysqrt(s)-m)*(mysqrt(s)-m),m)<<endl;
            cout<<"q:"<<q<<endl;
            cout<<"sigp:"<<sigp<<endl;
            cout<<"pR:"<<pR<<'\t'<<"pI:"<<pI<<'\t'<<"OPE:"<<real(ope)<<'\t'<<"+i "<<imag(ope)<<endl;
            //printf("%.16f\t %.16f\t %.16f\t %.16f\n",sigmapR,sigmapI,real(ope),imag(ope));
            
        }
        fout<<endl;
    }   

    fout.close();

    cout<<"bc left plus:"<<M2kbranchcut_left_momrep_plus(scomp,m)<<endl;
    cout<<"bc left minus:"<<M2kbranchcut_left_momrep_minus(scomp,m)<<endl;
    cout<<"bc right plus:"<<M2kbranchcut_right_momrep_plus(scomp,m)<<endl;
    cout<<"bc right minus:"<<M2kbranchcut_right_momrep_minus(scomp,m)<<endl;
    
    cout<<"filename:"<<filename<<endl;

    }
    /*

    string filename1 = "thresh_" + filename;//OPEcontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename1.c_str());

    //fout<<real(sigpplus(s,sigmab(a,m),m))<<'\t'<<real(sigpminus(s,sigmab(a,m),m))<<'\t'<<real(sigc1(s,sigmab(a,m),m))<<'\t'<<real(sigc2(s,sigmab(a,m),m))<<'\t'<<pow(sqrt(s)-m,2.0)<<'\t'<<0<<'\t'<<0<<endl;
    fout<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<real((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\t'
        <<imag((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\t'
        <<0<<endl;
    cout<<"resigplus:"<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"imsigplus:"<<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"resigminus:"<<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"imsigminus:"<<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"resigmax:"<<real((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\n'
        <<"imsigmax:"<<imag((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\n'
        <<0<<endl;
    fout.close();

    //string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";
    string filename2 = "sigvec_" + filename;
    fout.open(filename2.c_str());

    vector<comp> sigvec;
    sigma_vector_maker_6(sigvec,s,0.0,sigmax(s,m),a,m,1.0e-8,1000,200,1.0);
    for(int i=0;i<sigvec.size();++i)
    fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;

    fout.close();

    cout<<"file1:"<<filename<<'\n'
        <<"file2:"<<filename1<<'\n'
        <<"file3:"<<filename2<<endl;
    */
}

void kernelplotcontour_momrep_withintloop_withepsloop()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    //double s = 8.72;
    double eps = 0.001;
    double imS = 0.0;//-0.1;
    double sreal = 8.72;
    comp s = sreal + ii*imS;
    
    comp scomp = (comp) s;//+ ii*imS;
    comp sigb = sigmab(a,m) ;//+ ii*0.02;
    comp q = pmom(s,sigb,m);
    double N = 1000.0;

    double pRinitial = -0.2;
    double pRfinal = 0.2;
    double delpR = abs(pRinitial-pRfinal)/N;
    double pIinitial = -0.5;
    double pIfinal = 0.5;
    double delpI = abs(pIinitial - pIfinal)/N;
    
    //double sreal = real(s);

    for(eps=0.01;eps<=0.01;eps=eps+0.01)
    { 

    string filename =   "kernelcontour_momrep_a_" + to_string((int)a) 
                        + "_sreal_" + to_string(sreal) 
                        + "_ims_" + to_string(imS)
                        + "_eps_" + to_string(eps)
                        + "_sigpsb_" + to_string(real(sigmab(a,m))) 
                        + ".dat";

    int size = (int) N;
    int isize = size;
    int jsize = size;

    fout.open(filename.c_str());

    //for(double sigmapR=sigpRinitial;sigmapR<sigpRfinal;sigmapR=sigmapR + delsigpR)
    for(int i=0;i<isize;++i)
    {
        
        //for(double sigmapI=sigpIinitial;sigmapI<sigpIfinal;sigmapI=sigmapI + delsigpI)
        for(int j=0;j<jsize;++j)
        //double someeps = 1.0e-3;
        {
            comp pR = (comp)pRinitial + (comp)i*delpR;
            comp pI = (comp)pIinitial + (comp)j*delpI; 
            comp pcomp = (comp) pR + ii*pI;

            comp sigp = sigma_p(scomp,pcomp,m);

            double pi = acos(-1.0);

            comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps);
            //omega_comp(pcomp,m);
            //
            //M2kfunc(a,sigp,m,0.0);
            //GS_pk(s,q,pcomp,m,eps);
            //GS_1(scomp,sigmapcomp,sigmapprime,m,eps);
            //comp ope = GS_complex_s(sreal,sigmapcomp,sigmapprime,m,eps,imS);

            fout<<setprecision(16)<<real(pR)<<'\t'<<setprecision(16)<<real(pI)<<'\t'<<setprecision(16)<<real(ope)<<'\t'<<setprecision(16)<<imag(ope)<<endl;
            cout<<"kmin+:"<<pmom(s,(sqrt(s)-m)*(sqrt(s)-m),m)<<endl;
            cout<<"kmin-:"<<pmom(s,(mysqrt(s)-m)*(mysqrt(s)-m),m)<<endl;
            cout<<"q:"<<q<<endl;
            cout<<"sigp:"<<sigp<<endl;
            cout<<"pR:"<<pR<<'\t'<<"pI:"<<pI<<'\t'<<"OPE:"<<real(ope)<<'\t'<<"+i "<<imag(ope)<<endl;
            //printf("%.16f\t %.16f\t %.16f\t %.16f\n",sigmapR,sigmapI,real(ope),imag(ope));
            
        }
        fout<<endl;
    }   

    fout.close();

    cout<<"bc left plus:"<<M2kbranchcut_left_momrep_plus(scomp,m)<<endl;
    cout<<"bc left minus:"<<M2kbranchcut_left_momrep_minus(scomp,m)<<endl;
    cout<<"bc right plus:"<<M2kbranchcut_right_momrep_plus(scomp,m)<<endl;
    cout<<"bc right minus:"<<M2kbranchcut_right_momrep_minus(scomp,m)<<endl;
    cout<<"q:"<<q<<endl;
    
    cout<<"filename:"<<filename<<endl;

    comp kmax = pmom(s,0.0,m);
    //string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";
    string filename2 = "qvec_opposite_" + filename;
    fout.open(filename2.c_str());

    vector<comp> sigvec;
    mom_vector_maker_opposite_1(sigvec,s,0.0,kmax,a,m,0.01,200);
    for(int i=0;i<sigvec.size();++i)
    fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;

    fout.close();

    }

    
    /*

    string filename1 = "thresh_" + filename;//OPEcontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename1.c_str());

    //fout<<real(sigpplus(s,sigmab(a,m),m))<<'\t'<<real(sigpminus(s,sigmab(a,m),m))<<'\t'<<real(sigc1(s,sigmab(a,m),m))<<'\t'<<real(sigc2(s,sigmab(a,m),m))<<'\t'<<pow(sqrt(s)-m,2.0)<<'\t'<<0<<'\t'<<0<<endl;
    fout<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<real((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\t'
        <<imag((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\t'
        <<0<<endl;
    cout<<"resigplus:"<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"imsigplus:"<<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"resigminus:"<<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"imsigminus:"<<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"resigmax:"<<real((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\n'
        <<"imsigmax:"<<imag((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\n'
        <<0<<endl;
    fout.close();

    //string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";
    string filename2 = "sigvec_" + filename;
    fout.open(filename2.c_str());

    vector<comp> sigvec;
    sigma_vector_maker_6(sigvec,s,0.0,sigmax(s,m),a,m,1.0e-8,1000,200,1.0);
    for(int i=0;i<sigvec.size();++i)
    fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;

    fout.close();

    cout<<"file1:"<<filename<<'\n'
        <<"file2:"<<filename1<<'\n'
        <<"file3:"<<filename2<<endl;
    */
}

void kernelplotcontour_momrep_withintloop_withepsloop_drawpcut()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    //double s = 8.72;
    double eps = 0.001;
    double imS = 0.0;//-0.1;
    double sreal = 8.88;
    comp s = sreal + ii*imS;
    
    comp scomp = (comp) s;//+ ii*imS;
    comp sigb = sigmab(a,m) ;//+ ii*0.02;
    
    double N = 250.0;

    double pRinitial = -0.2;
    double pRfinal = 0.2;
    double delpR = abs(pRinitial-pRfinal)/N;
    double pIinitial = -0.5;
    double pIfinal = 0.5;
    double delpI = abs(pIinitial - pIfinal)/N;
    
    //double sreal = real(s);
    comp q;
    for(eps=0.001;eps<=0.001;eps=eps+0.001)
    { 

    q = pmom(s,sigb-ii*eps,m);
    string filename =   "kernelcontour_momrep_a_" + to_string((int)a) 
                        + "_sreal_" + to_string(sreal) 
                        + "_ims_" + to_string(imS)
                        + "_eps_" + to_string(eps)
                        + "_sigpsb_" + to_string(real(sigmab(a,m))) 
                        + ".dat";

    int size = (int) N;
    int isize = size;
    int jsize = size;

    fout.open(filename.c_str());

    //for(double sigmapR=sigpRinitial;sigmapR<sigpRfinal;sigmapR=sigmapR + delsigpR)
    for(int i=0;i<isize;++i)
    {
        
        //for(double sigmapI=sigpIinitial;sigmapI<sigpIfinal;sigmapI=sigmapI + delsigpI)
        for(int j=0;j<jsize;++j)
        //double someeps = 1.0e-3;
        {
            comp pR = (comp)pRinitial + (comp)i*delpR;
            comp pI = (comp)pIinitial + (comp)j*delpI; 
            comp pcomp = (comp) pR + ii*pI;

            comp sigp = sigma_p(scomp,pcomp,m);

            double pi = acos(-1.0);

            comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps);
            //omega_comp(pcomp,m);
            //
            //M2kfunc(a,sigp,m,0.0);
            //GS_pk(s,q,pcomp,m,eps);
            //GS_1(scomp,sigmapcomp,sigmapprime,m,eps);
            //comp ope = GS_complex_s(sreal,sigmapcomp,sigmapprime,m,eps,imS);

            fout<<setprecision(16)<<real(pR)<<'\t'<<setprecision(16)<<real(pI)<<'\t'<<setprecision(16)<<real(ope)<<'\t'<<setprecision(16)<<imag(ope)<<endl;
            cout<<"kmin+:"<<pmom(s,(sqrt(s)-m)*(sqrt(s)-m),m)<<endl;
            cout<<"kmin-:"<<pmom(s,(mysqrt(s)-m)*(mysqrt(s)-m),m)<<endl;
            cout<<"q:"<<q<<endl;
            cout<<"sigp:"<<sigp<<endl;
            cout<<"pR:"<<pR<<'\t'<<"pI:"<<pI<<'\t'<<"OPE:"<<real(ope)<<'\t'<<"+i "<<imag(ope)<<endl;
            //printf("%.16f\t %.16f\t %.16f\t %.16f\n",sigmapR,sigmapI,real(ope),imag(ope));
            
        }
        fout<<endl;
    }   

    fout.close();

    cout<<"bc left plus:"<<M2kbranchcut_left_momrep_plus(scomp,m)<<endl;
    cout<<"bc left minus:"<<M2kbranchcut_left_momrep_minus(scomp,m)<<endl;
    cout<<"bc right plus:"<<M2kbranchcut_right_momrep_plus(scomp,m)<<endl;
    cout<<"bc right minus:"<<M2kbranchcut_right_momrep_minus(scomp,m)<<endl;
    cout<<"q:"<<q<<endl;
    cout<<"qc1:"<<q_c1(s,a,m,eps)<<endl;
    cout<<"qc2:"<<q_c2(s,a,m,eps)<<endl;
    cout<<"qplus:"<<q_plus(s,a,m,eps)<<endl;
    cout<<"qminus:"<<q_minus(s,a,m,eps)<<endl;
    cout<<"filename:"<<filename<<endl;

    comp kmax = pmom(s,0.0,m);
    //string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";
    string filename2 = "pcut_" + filename;
    fout.open(filename2.c_str());

    for(double x=-1.0;x<=1.0;x=x+0.001)
    {
        fout<<real(pcut_plus(x,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(x,s,q,m,eps))<<'\t'
            <<real(pcut_minus(x,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(x,s,q,m,eps))<<'\t'
            <<0<<endl;
    }
    
    fout.close();

    string filename3 = "sigvec2_" + filename;
    fout.open(filename3.c_str());

    double rval = 0.001;
    vector<comp> qvec;
    mom_vector_maker_2(qvec,s,0.0,kmax,a,m,eps,eps,20,rval); //there are two eps here, one for ope, another for m2k
    for(int i=0;i<qvec.size();++i)
    {
        fout<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
    }
    fout.close();
    
    

    cout<<"file1:"<<filename<<'\n'
        <<"file2:"<<filename2<<'\n'
        <<"file3:"<<filename3<<endl;
    

    }

    
}


void kernelplotcontour_momrep_withintloop_withsloop()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    //double s = 8.72;
    double eps = 0.001;
    double imS = 0.0;//-0.1;
    double sreal = 8.72;
    //comp s = sreal + ii*imS;
    
    //comp scomp = (comp) s;//+ ii*imS;
    comp sigb = sigmab(a,m) ;//+ ii*0.02;
    
    double N = 250.0;

    double pRinitial = -0.2;
    double pRfinal = 0.2;
    double delpR = abs(pRinitial-pRfinal)/N;
    double pIinitial = -0.5;
    double pIfinal = 0.5;
    double delpI = abs(pIinitial - pIfinal)/N;

    double sinitial = 8.72;
    double sfinal = abs(phibthreshold(a,m));
    double dels = abs(sinitial-sfinal)/100.0;
    
    //double sreal = real(s);

    //for(eps=0.01;eps<=0.01;eps=eps+0.01)
    int count = 1;
    for(double s=sinitial;s<=sfinal;s=s+dels)
    { 
        sreal = s;
        imS = 0.0;
        comp scomp = s;

        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;

    string filename =   "kernelcontour_momrep_a_" + to_string((int)a) 
                        + "count_" + to_string(count)
                        + "_eps_" + to_string(eps)
                        + "_sigpsb_" + to_string(real(sigmab(a,m))) 
                        + ".dat";

    int size = (int) N;
    int isize = size;
    int jsize = size;

    comp q = pmom(s,sigb,m);

    fout.open(filename.c_str());

    //for(double sigmapR=sigpRinitial;sigmapR<sigpRfinal;sigmapR=sigmapR + delsigpR)
    for(int i=0;i<isize;++i)
    {
        
        //for(double sigmapI=sigpIinitial;sigmapI<sigpIfinal;sigmapI=sigmapI + delsigpI)
        for(int j=0;j<jsize;++j)
        //double someeps = 1.0e-3;
        {
            comp pR = (comp)pRinitial + (comp)i*delpR;
            comp pI = (comp)pIinitial + (comp)j*delpI; 
            comp pcomp = (comp) pR + ii*pI;

            comp sigp = sigma_p(scomp,pcomp,m);

            double pi = acos(-1.0);

            comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps);
            //omega_comp(pcomp,m);
            //
            //M2kfunc(a,sigp,m,0.0);
            //GS_pk(s,q,pcomp,m,eps);
            //GS_1(scomp,sigmapcomp,sigmapprime,m,eps);
            //comp ope = GS_complex_s(sreal,sigmapcomp,sigmapprime,m,eps,imS);

            fout<<setprecision(16)<<real(pR)<<'\t'<<setprecision(16)<<real(pI)<<'\t'<<setprecision(16)<<real(ope)<<'\t'<<setprecision(16)<<imag(ope)<<endl;
            //cout<<"kmin+:"<<pmom(s,(sqrt(s)-m)*(sqrt(s)-m),m)<<endl;
            //cout<<"kmin-:"<<pmom(s,(mysqrt(s)-m)*(mysqrt(s)-m),m)<<endl;
            //cout<<"q:"<<q<<endl;
            //cout<<"sigp:"<<sigp<<endl;
            //cout<<"pR:"<<pR<<'\t'<<"pI:"<<pI<<'\t'<<"OPE:"<<real(ope)<<'\t'<<"+i "<<imag(ope)<<endl;
            //printf("%.16f\t %.16f\t %.16f\t %.16f\n",sigmapR,sigmapI,real(ope),imag(ope));
            
        }
        fout<<endl;
    }   

    fout.close();

    cout<<"bc left plus:"<<M2kbranchcut_left_momrep_plus(scomp,m)<<endl;
    cout<<"bc left minus:"<<M2kbranchcut_left_momrep_minus(scomp,m)<<endl;
    cout<<"bc right plus:"<<M2kbranchcut_right_momrep_plus(scomp,m)<<endl;
    cout<<"bc right minus:"<<M2kbranchcut_right_momrep_minus(scomp,m)<<endl;
    cout<<"q:"<<q<<endl;

    cout<<"qc1:"<<q_c1(scomp,a,m,eps)<<endl;
    cout<<"qc2:"<<q_c2(scomp,a,m,eps)<<endl;
    cout<<"qplus:"<<q_plus(scomp,a,m,eps)<<endl;
    cout<<"qminus:"<<q_minus(scomp,a,m,eps)<<endl;


    string filename1 = "threshold_" + filename;
    fout.open(filename1.c_str());
    fout<<s<<'\t'
        <<real(M2kbranchcut_left_momrep_plus(scomp,m))<<'\t'
        <<imag(M2kbranchcut_left_momrep_plus(scomp,m))<<'\t'
        <<real(M2kbranchcut_left_momrep_minus(scomp,m))<<'\t'
        <<imag(M2kbranchcut_left_momrep_minus(scomp,m))<<'\t'
        <<real(M2kbranchcut_right_momrep_plus(scomp,m))<<'\t'
        <<imag(M2kbranchcut_right_momrep_plus(scomp,m))<<'\t'
        <<real(M2kbranchcut_right_momrep_minus(scomp,m))<<'\t'
        <<imag(M2kbranchcut_right_momrep_minus(scomp,m))<<'\t'
        <<real(q)<<'\t'
        <<imag(q)<<'\t'
        <<real(q_c1(scomp,a,m,eps))<<'\t'
        <<imag(q_c1(scomp,a,m,eps))<<'\t'
        <<real(q_c2(scomp,a,m,eps))<<'\t'
        <<imag(q_c2(scomp,a,m,eps))<<'\t'
        <<real(q_plus(scomp,a,m,eps))<<'\t'
        <<imag(q_plus(scomp,a,m,eps))<<'\t'
        <<real(q_minus(scomp,a,m,eps))<<'\t'
        <<imag(q_minus(scomp,a,m,eps))<<endl;

    fout.close();
    
    //cout<<"filename:"<<filename<<endl;

    comp kmax = pmom(s,0.0,m);
    //string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";
    string filename2 = "qvec_opposite_" + filename;
    fout.open(filename2.c_str());

    vector<comp> sigvec;
    mom_vector_maker_opposite_1(sigvec,s,0.0,kmax,a,m,eps,200);
    for(int i=0;i<sigvec.size();++i)
    fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;

    fout.close();

    cout<<"filename0:"<<filename<<endl;
    cout<<"filename1:"<<filename1<<endl;
    cout<<"filename2:"<<filename2<<endl;

    count = count + 1;

    }

    
   
}

void kernelplotcontour_momrep_withintloop_with_multiplier()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    //double s = 8.72;
    double eps = 0.001;
    double imS = 0.0;//-0.1;
    double sreal = 8.72;
    comp s = sreal + ii*imS;
    
    comp scomp = (comp) s;//+ ii*imS;
    comp sigb = sigmab(a,m) ;//+ ii*0.02;
    
    double N = 250.0;

    double pRinitial = -0.2;
    double pRfinal = 0.2;
    double delpR = abs(pRinitial-pRfinal)/N;
    double pIinitial = -0.5;
    double pIfinal = 0.5;
    double delpI = abs(pIinitial - pIfinal)/N;

    double sinitial = 8.72;
    double sfinal = abs(phibthreshold(a,m));
    double dels = abs(sinitial-sfinal)/100.0;
    
    //double sreal = real(s);

    //for(eps=0.01;eps<=0.01;eps=eps+0.01)
    int count = 1;
    //for(double s=sinitial;s<=sfinal;s=s+dels)
    for(int multiplier=-1;multiplier>=-1000;multiplier=multiplier*(10))
    { 
        //sreal = s;
        //imS = 0.0;
        //comp scomp = s;

        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;

    string filename =   "kernelcontour_momrep_a_" + to_string((int)a) 
                        + "_sreal_" + to_string(sreal)
                        + "_imS_" + to_string(imS)
                        + "_eps_" + to_string(eps)
                        + "_sigpsb_" + to_string(real(sigmab(a,m)))
                        + "_multiplier_" + to_string((int)multiplier) 
                        + ".dat";

    int size = (int) N;
    int isize = size;
    int jsize = size;

    comp q = pmom(s,sigb,m);

    fout.open(filename.c_str());

    //for(double sigmapR=sigpRinitial;sigmapR<sigpRfinal;sigmapR=sigmapR + delsigpR)
    for(int i=0;i<isize;++i)
    {
        
        //for(double sigmapI=sigpIinitial;sigmapI<sigpIfinal;sigmapI=sigmapI + delsigpI)
        for(int j=0;j<jsize;++j)
        //double someeps = 1.0e-3;
        {
            comp pR = (comp)pRinitial + (comp)i*delpR;
            comp pI = (comp)pIinitial + (comp)j*delpI; 
            comp pcomp = (comp) pR + ii*pI;

            comp sigp = sigma_p(scomp,pcomp,m);

            double pi = acos(-1.0);

            comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps);
            //omega_comp(pcomp,m);
            //
            //M2kfunc(a,sigp,m,0.0);
            //GS_pk(s,q,pcomp,m,eps);
            //GS_1(scomp,sigmapcomp,sigmapprime,m,eps);
            //comp ope = GS_complex_s(sreal,sigmapcomp,sigmapprime,m,eps,imS);

            fout<<setprecision(16)<<real(pR)<<'\t'<<setprecision(16)<<real(pI)<<'\t'<<setprecision(16)<<real(ope)<<'\t'<<setprecision(16)<<imag(ope)<<endl;
            //cout<<"kmin+:"<<pmom(s,(sqrt(s)-m)*(sqrt(s)-m),m)<<endl;
            //cout<<"kmin-:"<<pmom(s,(mysqrt(s)-m)*(mysqrt(s)-m),m)<<endl;
            //cout<<"q:"<<q<<endl;
            //cout<<"sigp:"<<sigp<<endl;
            //cout<<"pR:"<<pR<<'\t'<<"pI:"<<pI<<'\t'<<"OPE:"<<real(ope)<<'\t'<<"+i "<<imag(ope)<<endl;
            //printf("%.16f\t %.16f\t %.16f\t %.16f\n",sigmapR,sigmapI,real(ope),imag(ope));
            
        }
        fout<<endl;
    }   

    fout.close();

    cout<<"bc left plus:"<<M2kbranchcut_left_momrep_plus(scomp,m)<<endl;
    cout<<"bc left minus:"<<M2kbranchcut_left_momrep_minus(scomp,m)<<endl;
    cout<<"bc right plus:"<<M2kbranchcut_right_momrep_plus(scomp,m)<<endl;
    cout<<"bc right minus:"<<M2kbranchcut_right_momrep_minus(scomp,m)<<endl;
    cout<<"q:"<<q<<endl;

    cout<<"qc1:"<<q_c1(scomp,a,m,eps)<<endl;
    cout<<"qc2:"<<q_c2(scomp,a,m,eps)<<endl;
    cout<<"qplus:"<<q_plus(scomp,a,m,eps)<<endl;
    cout<<"qminus:"<<q_minus(scomp,a,m,eps)<<endl;


    string filename1 = "threshold_" + filename;
    fout.open(filename1.c_str());
    fout<<s<<'\t'
        <<real(M2kbranchcut_left_momrep_plus(scomp,m))<<'\t'
        <<imag(M2kbranchcut_left_momrep_plus(scomp,m))<<'\t'
        <<real(M2kbranchcut_left_momrep_minus(scomp,m))<<'\t'
        <<imag(M2kbranchcut_left_momrep_minus(scomp,m))<<'\t'
        <<real(M2kbranchcut_right_momrep_plus(scomp,m))<<'\t'
        <<imag(M2kbranchcut_right_momrep_plus(scomp,m))<<'\t'
        <<real(M2kbranchcut_right_momrep_minus(scomp,m))<<'\t'
        <<imag(M2kbranchcut_right_momrep_minus(scomp,m))<<'\t'
        <<real(q)<<'\t'
        <<imag(q)<<'\t'
        <<real(q_c1(scomp,a,m,eps))<<'\t'
        <<imag(q_c1(scomp,a,m,eps))<<'\t'
        <<real(q_c2(scomp,a,m,eps))<<'\t'
        <<imag(q_c2(scomp,a,m,eps))<<'\t'
        <<real(q_plus(scomp,a,m,eps))<<'\t'
        <<imag(q_plus(scomp,a,m,eps))<<'\t'
        <<real(q_minus(scomp,a,m,eps))<<'\t'
        <<imag(q_minus(scomp,a,m,eps))<<endl;

    fout.close();
    
    //cout<<"filename:"<<filename<<endl;

    comp kmax = pmom(s,0.0,m);
    //string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";
    string filename2 = "qvec_" + filename;
    fout.open(filename2.c_str());

    vector<comp> sigvec;
    mom_vector_maker_1(sigvec,s,0.0,kmax,a,m,eps,200,multiplier);
    for(int i=0;i<sigvec.size();++i)
    fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;

    fout.close();

    cout<<"filename0:"<<filename<<endl;
    cout<<"filename1:"<<filename1<<endl;
    cout<<"filename2:"<<filename2<<endl;

    count = count + 1;

    }

    
   
}

void kernelplotcontour_momrep_withintloop_glockletest()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    //double s = 8.72;
    double eps = 0.001;
    double imS = 0.0;//-0.1;
    double sreal = 8.72;
    comp s = sreal + ii*imS;
    
    comp scomp = (comp) s;//+ ii*imS;
    comp sigb = sigmab(a,m) ;//+ ii*0.02;
    
    double N = 250.0;

    double pRinitial = -0.2;
    double pRfinal = 0.2;
    double delpR = abs(pRinitial-pRfinal)/N;
    double pIinitial = -0.5;
    double pIfinal = 0.5;
    double delpI = abs(pIinitial - pIfinal)/N;

    double sinitial = 8.72;
    double sfinal = abs(phibthreshold(a,m));
    double dels = abs(sinitial-sfinal)/100.0;
    
    //double sreal = real(s);

    //for(eps=0.01;eps<=0.01;eps=eps+0.01)
    int count = 1;


    comp kmax = pmom(s,0.0,m);
    //string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";
    string filename2 = "qvec_glockletest_a_" + to_string((int)a) 
                        + "_sreal_" + to_string(sreal)
                        + "_imS_" + to_string(imS)
                        + "_eps_" + to_string(eps)
                        + "_sigb_" + to_string(real(sigmab(a,m)))
                        + ".dat";
    fout.open(filename2.c_str());

    vector<comp> sigvec;
    mom_vector_maker_1(sigvec,s,0.0,kmax,a,m,eps,25,1.0);
    for(int i=0;i<sigvec.size();++i)
    fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;

    fout.close();
    //for(double s=sinitial;s<=sfinal;s=s+dels)
    for(int k=0;k<sigvec.size();++k)
    { 
        //sreal = s;
        //imS = 0.0;
        //comp scomp = s;
    comp q = pmom(s,sigb-ii*eps,m);
    comp qmom = sigvec[k];
    cout<<"s:"<<s<<'\t'<<"qmom:"<<qmom<<'\t'<<"q:"<<q<<endl;
    double reqmom = real(qmom);
    double imqmom = imag(qmom);
    
    string filename =   "kernelcontour_momrep_glockletest_a_" + to_string((int)a) 
                        + "_sreal_" + to_string(sreal)
                        + "_imS_" + to_string(imS)
                        + "_eps_" + to_string(eps)
                        + "_sigpsb_" + to_string(real(sigmab(a,m)))
                        + "_qvec_" + to_string(k)
                        + ".dat";

    

    

    int size = (int) N;
    int isize = size;
    int jsize = size;

    
    

    fout.open(filename.c_str());

    //for(double sigmapR=sigpRinitial;sigmapR<sigpRfinal;sigmapR=sigmapR + delsigpR)
    for(int i=0;i<isize;++i)
    {
        
        //for(double sigmapI=sigpIinitial;sigmapI<sigpIfinal;sigmapI=sigmapI + delsigpI)
        for(int j=0;j<jsize;++j)
        //double someeps = 1.0e-3;
        {
            comp pR = (comp)pRinitial + (comp)i*delpR;
            comp pI = (comp)pIinitial + (comp)j*delpI; 
            comp pcomp = (comp) pR + ii*pI;

            comp sigp = sigma_p(scomp,pcomp,m);

            double pi = acos(-1.0);

            comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,qmom,pcomp,m,eps)*M2kfunc(a,sigp,m,eps);
            //omega_comp(pcomp,m);
            //
            //M2kfunc(a,sigp,m,0.0);
            //GS_pk(s,q,pcomp,m,eps);
            //GS_1(scomp,sigmapcomp,sigmapprime,m,eps);
            //comp ope = GS_complex_s(sreal,sigmapcomp,sigmapprime,m,eps,imS);

            fout<<setprecision(16)<<real(pR)<<'\t'<<setprecision(16)<<real(pI)<<'\t'<<setprecision(16)<<real(ope)<<'\t'<<setprecision(16)<<imag(ope)<<endl;
            //cout<<"kmin+:"<<pmom(s,(sqrt(s)-m)*(sqrt(s)-m),m)<<endl;
            //cout<<"kmin-:"<<pmom(s,(mysqrt(s)-m)*(mysqrt(s)-m),m)<<endl;
            //cout<<"q:"<<q<<endl;
            //cout<<"sigp:"<<sigp<<endl;
            //cout<<"pR:"<<pR<<'\t'<<"pI:"<<pI<<'\t'<<"OPE:"<<real(ope)<<'\t'<<"+i "<<imag(ope)<<endl;
            //printf("%.16f\t %.16f\t %.16f\t %.16f\n",sigmapR,sigmapI,real(ope),imag(ope));
            
        }
        fout<<endl;
    }   

    fout.close();

    cout<<"bc left plus:"<<M2kbranchcut_left_momrep_plus(scomp,m)<<endl;
    cout<<"bc left minus:"<<M2kbranchcut_left_momrep_minus(scomp,m)<<endl;
    cout<<"bc right plus:"<<M2kbranchcut_right_momrep_plus(scomp,m)<<endl;
    cout<<"bc right minus:"<<M2kbranchcut_right_momrep_minus(scomp,m)<<endl;
    cout<<"q:"<<q<<endl;

    cout<<"qc1:"<<q_c1(scomp,a,m,eps)<<endl;
    cout<<"qc2:"<<q_c2(scomp,a,m,eps)<<endl;
    cout<<"qplus:"<<q_plus(scomp,a,m,eps)<<endl;
    cout<<"qminus:"<<q_minus(scomp,a,m,eps)<<endl;


    string filename1 = "threshold_" + filename;
    fout.open(filename1.c_str());
    fout<<real(s)<<'\t'
        <<imag(s)<<'\t'
        <<real(M2kbranchcut_left_momrep_plus(scomp,m))<<'\t'
        <<imag(M2kbranchcut_left_momrep_plus(scomp,m))<<'\t'
        <<real(M2kbranchcut_left_momrep_minus(scomp,m))<<'\t'
        <<imag(M2kbranchcut_left_momrep_minus(scomp,m))<<'\t'
        <<real(M2kbranchcut_right_momrep_plus(scomp,m))<<'\t'
        <<imag(M2kbranchcut_right_momrep_plus(scomp,m))<<'\t'
        <<real(M2kbranchcut_right_momrep_minus(scomp,m))<<'\t'
        <<imag(M2kbranchcut_right_momrep_minus(scomp,m))<<'\t'
        <<real(q)<<'\t'
        <<imag(q)<<'\t'
        <<real(q_c1(scomp,a,m,eps))<<'\t'
        <<imag(q_c1(scomp,a,m,eps))<<'\t'
        <<real(q_c2(scomp,a,m,eps))<<'\t'
        <<imag(q_c2(scomp,a,m,eps))<<'\t'
        <<real(q_plus_k(scomp,qmom,m,eps))<<'\t'
        <<imag(q_plus_k(scomp,qmom,m,eps))<<'\t'
        <<real(q_minus_k(scomp,qmom,m,eps))<<'\t'
        <<imag(q_minus_k(scomp,qmom,m,eps))<<endl;

    fout.close();
    
    //cout<<"filename:"<<filename<<endl;

    

    cout<<"filename0:"<<filename<<endl;
    cout<<"filename1:"<<filename1<<endl;
    cout<<"filename2:"<<filename2<<endl;

    count = count + 1;

    }

    
   
}

void branchcuts_and_thresholds_sloop()
{
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    double eps = 0.00001;


    double sinitial = 8.72;
    double sfinal = abs(phibthreshold(a,m));
    cout<<"threshold:"<<sfinal<<endl;
    double dels = abs(sinitial-sfinal)/1000.0;
    comp sigb = sigmab(a,m);

    ofstream fout;

    string filename =   "bcandthresholds_a_" + to_string(a)
                        + "_eps_" + to_string(eps) 
                        + "_sloop.dat";

    fout.open(filename.c_str());
    for(double s=sinitial;s<sfinal;s=s+dels)
    {
        comp qplus = q_plus(s,a,m,eps);
        comp qminus = q_minus(s,a,m,eps);
        comp q = pmom(s,sigb-ii*eps,m);


        fout<<s<<'\t'
            <<real(qplus)<<'\t'
            <<imag(qplus)<<'\t'
            <<real(qminus)<<'\t'
            <<imag(qminus)<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<endl;
        
        /*cout<<s<<'\t'
            <<real(qplus)<<'\t'
            <<imag(qplus)<<'\t'
            <<real(qminus)<<'\t'
            <<imag(qminus)<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<endl;*/

    }

    fout.close();
    
}

void branchcuts_and_thresholds_qvecloop()
{
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    double eps = 0.001;


    double sinitial = 8.72;
    double sfinal = abs(phibthreshold(a,m));
    cout<<"threshold:"<<sfinal<<endl;
    double dels = abs(sinitial-sfinal)/1000.0;
    double s = sinitial;
    comp sigb = sigmab(a,m);

    comp kmax = pmom(s,0.0,m);
    //string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";
    
    

    vector<comp> sigvec;
    mom_vector_maker_1(sigvec,s,0.0,kmax,a,m,eps,200,1.0);
    
    

    ofstream fout;

    string filename =   "bcandthresholds_a_" + to_string(a)
                        + "_eps_" + to_string(eps) 
                        + "_qvecloop.dat";

    fout.open(filename.c_str());
    for(int i=0;i<sigvec.size();++i)
    {

        comp qmom = sigvec[i];
        comp sigq = sigma_p(s,qmom,m);
        comp qplus = q_plus_k(s,sigq,m,eps);
        comp qminus = q_minus_k(s,sigq,m,eps);
        comp q = pmom(s,sigb-ii*eps,m);


        fout<<i<<'\t'
            <<real(qplus)<<'\t'
            <<imag(qplus)<<'\t'
            <<real(qminus)<<'\t'
            <<imag(qminus)<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<endl;
        
        /*cout<<s<<'\t'
            <<real(qplus)<<'\t'
            <<imag(qplus)<<'\t'
            <<real(qminus)<<'\t'
            <<imag(qminus)<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<endl;*/

    }

    fout.close();
    
}

void glockletest()
{
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    
    double sreal = 8.72;
    double imS = 0.0;
    comp s = sreal + ii*imS;

    double qinitial = 0.0;
    double qfinal = 1.0;
    double delq = abs(qfinal - qinitial)/1000.0;
    double eps = 0.1;
    double x = 1.0;

    comp kmax = pmom(s,0.0,m);
    
    //vector<comp> sigvec;
    //mom_vector_maker_1(sigvec,s,0.0,kmax,a,m,eps,250,1.0);
    
    ofstream fout;
    string filename =   "glockletest_qreal_sreal_" + to_string(sreal) 
                            + "_imS_" + to_string(imS) 
                            + "_eps_" + to_string(eps) + ".dat";
    fout.open(filename.c_str());
    

    for(double x=-1;x<=1;x=x+0.01)
    {
    for(double q=qinitial;q<=qfinal;q=q+delq)
    //for(int i=0;i<sigvec.size();++i)
    {
        //comp q = sigvec[i];
        fout<<real(q)<<'\t'
            <<imag(q)<<'\t'
            <<real(pcut_plus(x,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(x,s,q,m,eps))<<'\t'
            <<real(pcut_minus(x,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(x,s,q,m,eps))<<endl;
    }

    
    }
    fout.close();

    /*cout<<pcut_minus(1,s,1,m,0.0)<<endl;
    cout<<pcut_minus(-1,s,1,m,0.0)<<endl;
    cout<<q_plus_k(s,1,m,0.0)<<endl;
    cout<<-q_minus_k(s,1,m,0.0)<<endl;*/
}


void glockletest_C1()
{
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    
    double sreal = 8.72;
    double imS = 0.0;
    comp s = sreal + ii*imS;

    double qinitial = 0.0;
    double qfinal = 1.0;
    double delq = abs(qfinal - qinitial)/1000.0;
    double eps = 0.001;
    double x = 1.0;

    comp kmax = pmom(s,0.0,m);
    double rval = 0.001;
    vector<comp> sigvec;
    mom_vector_maker_2(sigvec,s,0.0,kmax,a,m,eps,eps,250,rval);//there are two eps here, one for ope another for m2k
    
    

    for(int x=-1;x<=1;x=x+1)
    {
        ofstream fout;
        string filename =   "glockletest_qC1_sreal_" + to_string(sreal) 
                            + "_imS_" + to_string(imS) 
                            + "_x_" + to_string(x) + "_eps_" + to_string(eps) + "_multiplier_-20.dat";
        fout.open(filename.c_str());
    //for(double q=qinitial;q<=qfinal;q=q+delq)
    for(int i=0;i<sigvec.size();++i)
    {
        comp q = sigvec[i];
        fout<<real(q)<<'\t'
            <<imag(q)<<'\t'
            <<real(pcut_plus(x,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(x,s,q,m,eps))<<'\t'
            <<real(pcut_minus(x,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(x,s,q,m,eps))<<endl;
    }

    fout.close();
    }

    ofstream fout;
    string filename =   "glockletest_qC1_sreal_" + to_string(sreal) 
                            + "_imS_" + to_string(imS) 
                            + "_eps_" + to_string(eps) + "_r_" + to_string(rval) + ".dat";
    fout.open(filename.c_str());

    for(double x=-1;x<=1;x=x+0.01)
    {
        
    //for(double q=qinitial;q<=qfinal;q=q+delq)
    for(int i=0;i<sigvec.size();++i)
    {
        comp q = sigvec[i];
        fout<<real(q)<<'\t'
            <<imag(q)<<'\t'
            <<real(pcut_plus(x,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(x,s,q,m,eps))<<'\t'
            <<real(pcut_minus(x,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(x,s,q,m,eps))<<endl;
    }

    
    }
    fout.close();

    comp sigb = sigmab(a,m);
    comp qb = pmom(s,sigb-ii*eps,m);
    cout<<qb<<endl;
}

void opeplotcontour_withintloop_withloopeps()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    //double s = 8.72;
    //double eps = 0.09;
    double imS = -0.001;
    double sreal = 8.72;
    comp s = sreal + ii*imS;
    
    comp scomp = (comp) s;//+ ii*imS;
    comp sigmapprime = sigmab(a,m) ;//+ ii*0.02;

    double sigpRinitial = 3.6;
    double sigpRfinal = 4.0;
    double delsigpR = abs(sigpRinitial-sigpRfinal)/250.0;
    double sigpIinitial = -0.1;
    double sigpIfinal = 0.1;
    double delsigpI = abs(sigpIinitial - sigpIfinal)/250.0;
    
    //double sreal = real(s);

    for(double eps=0.01;eps<=0.1;eps=eps+0.01)
    { 

    string filename =   "OPEcontour_a_" + to_string((int)a) 
                        + "_sreal_" + to_string(sreal) 
                        + "_ims_" + to_string(imS)
                        + "_eps_" + to_string(eps)
                        + "_sigpsb_" + to_string(real(sigmab(a,m))) 
                        + ".dat";

    int size = 250;
    int isize = size;
    int jsize = size;

    fout.open(filename.c_str());

    //for(double sigmapR=sigpRinitial;sigmapR<sigpRfinal;sigmapR=sigmapR + delsigpR)
    for(int i=0;i<isize;++i)
    {
        
        //for(double sigmapI=sigpIinitial;sigmapI<sigpIfinal;sigmapI=sigmapI + delsigpI)
        for(int j=0;j<jsize;++j)
        //double someeps = 1.0e-3;
        {
            comp sigmapR = (comp)sigpRinitial + (comp)i*delsigpR;
            comp sigmapI = (comp)sigpIinitial + (comp)j*delsigpI; 
            comp sigmapcomp = (comp) sigmapR + ii*sigmapI;

            comp ope = GS_1(scomp,sigmapcomp,sigmapprime,m,eps);
            //comp ope = GS_complex_s(sreal,sigmapcomp,sigmapprime,m,eps,imS);

            fout<<setprecision(16)<<real(sigmapR)<<'\t'<<setprecision(16)<<real(sigmapI)<<'\t'<<setprecision(16)<<real(ope)<<'\t'<<setprecision(16)<<imag(ope)<<endl;
            cout<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
            //printf("%.16f\t %.16f\t %.16f\t %.16f\n",sigmapR,sigmapI,real(ope),imag(ope));
            
        }
        fout<<endl;
    }   

    fout.close();

    string filename1 = "thresh_" + filename;//OPEcontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename1.c_str());

    //fout<<real(sigpplus(s,sigmab(a,m),m))<<'\t'<<real(sigpminus(s,sigmab(a,m),m))<<'\t'<<real(sigc1(s,sigmab(a,m),m))<<'\t'<<real(sigc2(s,sigmab(a,m),m))<<'\t'<<pow(sqrt(s)-m,2.0)<<'\t'<<0<<'\t'<<0<<endl;
    fout<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<real((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\t'
        <<imag((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\t'
        <<0<<endl;
    cout<<"resigplus:"<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"imsigplus:"<<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"resigminus:"<<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"imsigminus:"<<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\n'
        <<"resigmax:"<<real((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\n'
        <<"imsigmax:"<<imag((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\n'
        <<0<<endl;
    fout.close();

    //string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";
    string filename2 = "sigvec_" + filename;
    fout.open(filename2.c_str());

    vector<comp> sigvec;
    sigma_vector_maker_6(sigvec,s,0.0,sigmax(s,m),a,m,1.0e-8,1000,200,1.0);
    for(int i=0;i<sigvec.size();++i)
    fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;

    fout.close();

    cout<<"file1:"<<filename<<'\n'
        <<"file2:"<<filename1<<'\n'
        <<"file3:"<<filename2<<endl;

    }
}

void zpk_energyfactor_contour_sigp_vs_sigk_withintloop()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    //double s = 8.72;
    double eps = 0.0;
    double imS = 0.0;
    double sreal = 8.72;
    comp s = sreal + ii*imS;
    
    comp scomp = (comp) s;//+ ii*imS;
    //comp sigmapprime = sigmab(a,m) ;//+ ii*0.02;

    double sigpinitial = 0.0;
    double sigpfinal = 4.0;
    double delsigp = abs(sigpinitial-sigpfinal)/250.0;
    double sigkinitial = 0.0;
    double sigkfinal = 4.0;
    double delsigk = abs(sigkinitial - sigkfinal)/250.0;
    
    //double sreal = real(s);

     

    string filename =   "zpk_energyfactor_contour_a_" + to_string((int)a) 
                        + "_sreal_" + to_string(sreal) 
                        + "_eps_" + to_string(eps)
                        + ".dat";

    int size = 250;
    int isize = size;
    int jsize = size;

    fout.open(filename.c_str());

    //for(double sigmapR=sigpRinitial;sigmapR<sigpRfinal;sigmapR=sigmapR + delsigpR)
    for(int i=0;i<isize;++i)
    {
        
        //for(double sigmapI=sigpIinitial;sigmapI<sigpIfinal;sigmapI=sigmapI + delsigpI)
        for(int j=0;j<jsize;++j)
        //double someeps = 1.0e-3;
        {
            comp sigmap = (comp)sigpinitial + (comp)i*delsigp;
            comp sigmak = (comp)sigkinitial + (comp)j*delsigk; 
            //comp sigmapcomp = (comp) sigmapR + ii*sigmapI;

            comp result = zpk_energyfactor(scomp,sigmap,sigmak,m,eps);

            fout<<setprecision(16)<<real(sigmap)<<'\t'
                <<setprecision(16)<<real(sigmak)<<'\t'
                <<setprecision(16)<<real(result)<<'\t'
                <<setprecision(16)<<imag(result)<<endl;
            cout<<sigmap<<'\t'<<sigmak<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
            //printf("%.16f\t %.16f\t %.16f\t %.16f\n",sigmapR,sigmapI,real(ope),imag(ope));
            
        }
        fout<<endl;
    }   

    fout.close();

    

}

void zpk_energyfactor_contour_sigp_vs_s_withintloop()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    //double s = 8.72;
    double eps = 0.0;
    double imS = 0.0;
    //double sreal = 8.72;
    //comp s = sreal + ii*imS;
    
    //comp scomp = (comp) s;//+ ii*imS;
    //comp sigmapprime = sigmab(a,m) ;//+ ii*0.02;

    double sigpinitial = 0.0;
    double sigpfinal = 4.0;
    double delsigp = abs(sigpinitial-sigpfinal)/250.0;
    double sigkinitial = 0.0;
    double sigkfinal = 4.0;
    double delsigk = abs(sigkinitial - sigkfinal)/250.0;
    double sinitial = 8.72;
    double sfinal = 9.0;
    double dels = abs(sinitial-sfinal)/250.0;
    
    //double sreal = real(s);
    comp sigb = sigmab(a,m);

     

    string filename =   "zpk_energyfactor_contour_a_" + to_string((int)a) 
                        + "_sigk_" + to_string(real(sigb)) 
                        + "_eps_" + to_string(eps)
                        + ".dat";

    int size = 250;
    int isize = size;
    int jsize = size;

    fout.open(filename.c_str());

    //for(double sigmapR=sigpRinitial;sigmapR<sigpRfinal;sigmapR=sigmapR + delsigpR)
    for(int i=0;i<isize;++i)
    {
        
        //for(double sigmapI=sigpIinitial;sigmapI<sigpIfinal;sigmapI=sigmapI + delsigpI)
        for(int j=0;j<jsize;++j)
        //double someeps = 1.0e-3;
        {
            comp scomp = sinitial + (comp)i*dels;
            comp sigmap = (comp)sigpinitial + (comp)j*delsigp;
            //comp sigmak = (comp)sigkinitial + (comp)j*delsigk; 
            comp sigmak = sigmab(a,m);
            double sreal = real(scomp);
             
            //comp sigmapcomp = (comp) sigmapR + ii*sigmapI;

            comp result = zpk_energyfactor(scomp,sigmap,sigmak,m,eps);

            fout<<setprecision(16)<<real(scomp)<<'\t'
                <<setprecision(16)<<real(sigmap)<<'\t'
                <<setprecision(16)<<real(result)<<'\t'
                <<setprecision(16)<<imag(result)<<endl;
            cout<<sreal<<'\t'<<sigmap<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
            //printf("%.16f\t %.16f\t %.16f\t %.16f\n",sigmapR,sigmapI,real(ope),imag(ope));
            
        }
        fout<<endl;
    }   

    fout.close();

    

}



void kernelplotcontour_withintloop()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    //double s = 8.72;
    comp s = 8.72 - ii*0.01;
    double eps = 0.0;//1.0e-5;
    double imS = 0.07;
    comp scomp = (comp) s;//+ ii*imS;
    comp sigmapprime = sigmab(a,m) ;//+ ii*0.02;

    double sigpRinitial = 3.6;
    double sigpRfinal = 4.0;
    double delsigpR = abs(sigpRinitial-sigpRfinal)/250.0;
    double sigpIinitial = -0.1;
    double sigpIfinal = 0.1;
    double delsigpI = abs(sigpIinitial - sigpIfinal)/250.0;
    
    double sreal = real(s);

     

    string filename = "kernelcontour_a_" + to_string((int)a) + "_s-i_" + to_string(sreal) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    int size = 250;
    int isize = size;
    int jsize = size;

    fout.open(filename.c_str());

    //for(double sigmapR=sigpRinitial;sigmapR<sigpRfinal;sigmapR=sigmapR + delsigpR)
    for(int i=0;i<isize;++i)
    {
        
        //for(double sigmapI=sigpIinitial;sigmapI<sigpIfinal;sigmapI=sigmapI + delsigpI)
        for(int j=0;j<jsize;++j)
        //double someeps = 1.0e-3;
        {
            comp sigmapR = (comp)sigpRinitial + (comp)i*delsigpR;
            comp sigmapI = (comp)sigpIinitial + (comp)j*delsigpI; 
            comp sigmapcomp = (comp) sigmapR + ii*sigmapI;

            comp ope = GS(scomp,sigmapprime,sigmapcomp,m,eps)*tau1(scomp,sigmapcomp,m)*M2kfunc(a,sigmapcomp,m,eps);

            fout<<setprecision(16)<<real(sigmapR)<<'\t'<<setprecision(16)<<real(sigmapI)<<'\t'<<setprecision(16)<<real(ope)<<'\t'<<setprecision(16)<<imag(ope)<<endl;
            cout<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
            //printf("%.16f\t %.16f\t %.16f\t %.16f\n",sigmapR,sigmapI,real(ope),imag(ope));
            
        }
        fout<<endl;
    }   

    fout.close();

    string filename1 = "thresh_" + filename;//OPEcontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename1.c_str());

    //fout<<real(sigpplus(s,sigmab(a,m),m))<<'\t'<<real(sigpminus(s,sigmab(a,m),m))<<'\t'<<real(sigc1(s,sigmab(a,m),m))<<'\t'<<real(sigc2(s,sigmab(a,m),m))<<'\t'<<pow(sqrt(s)-m,2.0)<<'\t'<<0<<'\t'<<0<<endl;
    fout<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<real((sqrt(s)-m)*(sqrt(s)-m))<<'\t'
        <<imag((sqrt(s)-m)*(sqrt(s)-m))<<'\t'
        <<0<<endl;
    fout.close();

    //string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";
    string filename2 = "sigvec_" + filename;
    fout.open(filename2.c_str());

    vector<comp> sigvec;
    sigma_vector_maker_8(sigvec,s,0.0,sigmax(s,m),a,m,1.0e-8,1000,200,0.0);
    for(int i=0;i<sigvec.size();++i)
    fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;

    fout.close();
    cout<<"files: 1)"<<filename<<'\n'
        <<"       2)"<<filename1<<'\n'
        <<"       3)"<<filename2<<endl;

}

void kernelplotcontour_withintloop_vs_s3()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    //double s = 8.72;
    //comp s = 8.72;// - ii*0.01;
    double eps = 1.0e-5;
    double imS = 0.07;
    //comp scomp = (comp) s;//+ ii*imS;
    comp sigmapprime = sigmab(a,m) ;//+ ii*0.02;

    double sigpRinitial = 3.6;
    double sigpRfinal = 4.0;
    double delsigpR = abs(sigpRinitial-sigpRfinal)/250.0;
    double sigpIinitial = -0.1;
    double sigpIfinal = 0.1;
    double delsigpI = abs(sigpIinitial - sigpIfinal)/250.0;
    
    //double sreal = real(s);
    double sinitial = 8.77;
    double sfinal = 8.79;
    double dels = abs(sfinal-sinitial)/50.0;
    int count = 0;

    for(double s=sinitial;s<=sfinal;s=s+dels)
    { 
        comp scomp = (comp) s;

        //string filename = "kernelcontour_a_" + to_string((int)a) + "_s_" + to_string(sreal) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";
        string filename =   "kernelcontour_a_" + to_string((int)a) 
                            + "_sigpsb_" + to_string(real(sigmab(a,m))) 
                            + "_sinitial_" + to_string(sinitial) 
                            + "_dels_" + to_string(dels)
                            + "_" + to_string(count) 
                            + ".dat";


        int size = 250;
        int isize = size;
        int jsize = size;

        fout.open(filename.c_str());

        //for(double sigmapR=sigpRinitial;sigmapR<sigpRfinal;sigmapR=sigmapR + delsigpR)
        for(int i=0;i<isize;++i)
        {
        
            //for(double sigmapI=sigpIinitial;sigmapI<sigpIfinal;sigmapI=sigmapI + delsigpI)
            for(int j=0;j<jsize;++j)
            //double someeps = 1.0e-3;
            {
                comp sigmapR = (comp)sigpRinitial + (comp)i*delsigpR;
                comp sigmapI = (comp)sigpIinitial + (comp)j*delsigpI; 
                comp sigmapcomp = (comp) sigmapR + ii*sigmapI;

                comp ope = GS(scomp,sigmapprime,sigmapcomp,m,eps)*tau1(scomp,sigmapcomp,m)*M2kfunc(a,sigmapcomp,m,eps);

                fout<<setprecision(16)<<real(sigmapR)<<'\t'<<setprecision(16)<<real(sigmapI)<<'\t'<<setprecision(16)<<real(ope)<<'\t'<<setprecision(16)<<imag(ope)<<endl;
                //cout<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
                //printf("%.16f\t %.16f\t %.16f\t %.16f\n",sigmapR,sigmapI,real(ope),imag(ope));
            
            }
            fout<<endl;
        }   

        fout.close();

        string filename1 = "thresh_" + filename;//OPEcontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

        fout.open(filename1.c_str());

        //fout<<real(sigpplus(s,sigmab(a,m),m))<<'\t'<<real(sigpminus(s,sigmab(a,m),m))<<'\t'<<real(sigc1(s,sigmab(a,m),m))<<'\t'<<real(sigc2(s,sigmab(a,m),m))<<'\t'<<pow(sqrt(s)-m,2.0)<<'\t'<<0<<'\t'<<0<<endl;
        fout<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
            <<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
            <<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
            <<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
            <<real((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\t'
            <<imag((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\t'
            <<0<<endl;
        fout.close();

        //string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";
        string filename2 = "sigvec_" + filename;
        fout.open(filename2.c_str());

        vector<comp> sigvec;
        sigma_vector_maker_6(sigvec,scomp,0.0,sigmax(scomp,m),a,m,1.0e-8,1000,200,1.0);
        for(int i=0;i<sigvec.size();++i)
        fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;

        fout.close();
        cout<<"s3: "<<s<<" dels: "<<dels<<endl; 
        cout<<"files: 1)"<<filename<<'\n'
            <<"       2)"<<filename1<<'\n'
            <<"       3)"<<filename2<<endl;
        cout<<"---------------------------------------------"<<endl;
        count = count + 1;
    }
}


void kernelplot_withsigvec()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    //double s = 8.72;
    double imS = 0.00001;
    comp s = 8.72-ii*imS ;
    double eps = 0.0;//1.0e-3;
    
    comp scomp = (comp) s;//+ ii*imS;
    comp sigmapprime = sigmab(a,m) ;//+ ii*0.02;

    //double sigpRinitial = 3.6;
    //double sigpRfinal = 4.0;
    //double delsigpR = abs(sigpRinitial-sigpRfinal)/250.0;
    //double sigpIinitial = -0.1;
    //double sigpIfinal = 0.1;
    //double delsigpI = abs(sigpIinitial - sigpIfinal)/250.0;
    
    double sreal = real(s);
    comp sigb = sigmab(a,m);

    
    for(int box=1;box<=1;box=box*10) 
    {
        double N1 = 5000;
        double N2 = 500;
        vector<comp> sigvec;
        double eps1 = 0.0;//abs((1.0e-2)*eps*eps_energy_factor_minus(s,sigb,m)); 
        //sigma_vector_maker_5(sigvec,s,0.0,sigmax(s,m),a,m,eps1,N1,N2,box);
        sigma_vector_maker_8(sigvec,s,0.0,sigmax(s,m),a,m,eps1,N1,N2,0);

        /*string filename =   "kernel_sigvec_a_" + to_string((int)a) 
                            + "_s_" + to_string(sreal) 
                            + "_sigpsb_" + to_string(real(sigmab(a,m))) 
                            + "_epspow_" + to_string((int) abs(log10(eps)))
                            + "_N1_" + to_string((int)N1)
                            + "_N2_" + to_string((int)N2)
                            + "_box_" + to_string(box)
                            + "_pow2eneps.dat";
        */

        string filename =   "kernel_sigvec_a_" + to_string((int)a) 
                            + "_s_" + to_string(sreal) 
                            + "_sigpsb_" + to_string(real(sigmab(a,m))) 
                            + "_eps2_" + to_string(imS)
                            + "_N1_" + to_string((int)N1)
                            + "_N2_" + to_string((int)N2)
                            + "_box_" + to_string(box)
                            + "_pow2eneps.dat";

        int size = 250;
        int isize = size;
        int jsize = size;

        fout.open(filename.c_str());

        for(int i=0;i<sigvec.size();++i)
        {
            comp ker = GS(s,sigb,sigvec[i],m,eps)*tau1(s,sigvec[i],m)*M2kfunc(a,sigvec[i],m,eps);
            fout<<i<<'\t'<<real(ker)<<'\t'<<imag(ker)<<endl;
        }

        fout.close();


        //fout<<real(sigpplus(s,sigmab(a,m),m))<<'\t'<<real(sigpminus(s,sigmab(a,m),m))<<'\t'<<real(sigc1(s,sigmab(a,m),m))<<'\t'<<real(sigc2(s,sigmab(a,m),m))<<'\t'<<pow(sqrt(s)-m,2.0)<<'\t'<<0<<'\t'<<0<<endl;
        cout<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
            <<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
            <<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
            <<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
            <<(sqrt(s)-m)*(sqrt(s)-m)<<'\t'<<0<<endl;
    
    }

}

void ope_tau_m2_plot_withsigvec()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    //double s = 8.72;
    double imS = 0.0;
    comp s = 8.72-ii*imS ;
    double eps = 1.0e-5;
    
    comp scomp = (comp) s;//+ ii*imS;
    comp sigmapprime = sigmab(a,m) ;//+ ii*0.02;

    //double sigpRinitial = 3.6;
    //double sigpRfinal = 4.0;
    //double delsigpR = abs(sigpRinitial-sigpRfinal)/250.0;
    //double sigpIinitial = -0.1;
    //double sigpIfinal = 0.1;
    //double delsigpI = abs(sigpIinitial - sigpIfinal)/250.0;
    
    double sreal = real(s);
    comp sigb = sigmab(a,m);

    
    for(int box=1;box<=1;box=box*10) 
    {
        double N1 = 5000;
        double N2 = 500;
        vector<comp> sigvec;
        double eps1 = abs((1.0e-2)*eps*eps_energy_factor_minus(s,sigb,m)); 
        //sigma_vector_maker_5(sigvec,s,0.0,sigmax(s,m),a,m,eps1,N1,N2,box);
        sigma_vector_maker_6(sigvec,s,0.0,sigmax(s,m),a,m,eps1,N1,N2,1);

        /*string filename =   "kernel_sigvec_a_" + to_string((int)a) 
                            + "_s_" + to_string(sreal) 
                            + "_sigpsb_" + to_string(real(sigmab(a,m))) 
                            + "_epspow_" + to_string((int) abs(log10(eps)))
                            + "_N1_" + to_string((int)N1)
                            + "_N2_" + to_string((int)N2)
                            + "_box_" + to_string(box)
                            + "_pow2eneps.dat";
        */

        string filename =   "deltaGS_sigvec_a_" + to_string((int)a) 
                            + "_s_" + to_string(sreal) 
                            + "_sigpsb_" + to_string(real(sigmab(a,m))) 
                            + "_eps2_" + to_string(imS)
                            + "_N1_" + to_string((int)N1)
                            + "_N2_" + to_string((int)N2)
                            + "_box_" + to_string(box)
                            + "_pow2eneps.dat";

        int size = 250;
        int isize = size;
        int jsize = size;

        fout.open(filename.c_str());

        for(int i=0;i<sigvec.size();++i)
        {
            comp sigpoint = sigpplus_witheps(s,sigb,m,eps);
            //cout<<"sigpoint:"<<sigpoint<<endl;
            comp sigpoint_final = real(sigpoint) + ii*eps1;
            comp res1 = deltaGS_1(s,sigb,sigvec[i],sigpoint_final,m,eps); 
            //GS(s,sigb,sigvec[i],m,eps);
            comp res2 = GS(s,sigb,sigvec[i],m,eps);//tau1(s,sigvec[i],m);
            comp res3 = M2kfunc(a,sigvec[i],m,eps);
            fout<<i<<'\t'<<real(res1)<<'\t'<<imag(res1)
                   <<'\t'<<real(res2)<<'\t'<<imag(res2)
                   <<'\t'<<real(res3)<<'\t'<<imag(res3)<<endl;
        }

        fout.close();


        //fout<<real(sigpplus(s,sigmab(a,m),m))<<'\t'<<real(sigpminus(s,sigmab(a,m),m))<<'\t'<<real(sigc1(s,sigmab(a,m),m))<<'\t'<<real(sigc2(s,sigmab(a,m),m))<<'\t'<<pow(sqrt(s)-m,2.0)<<'\t'<<0<<'\t'<<0<<endl;
        cout<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
            <<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
            <<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
            <<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
            <<(sqrt(s)-m)*(sqrt(s)-m)<<'\t'<<0<<endl;
    
    }

}


void opeplotcontour_withintloop_multisigvec()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    double s = 8.72;
    double eps = 1.0e-3;
    double imS = 0.07;
    comp scomp = (comp) s;//+ ii*imS;
    comp sigmapprime = sigmab(a,m) ;//+ ii*0.02;

    double sigpRinitial = 3.6;
    double sigpRfinal = 4.0;
    double delsigpR = abs(sigpRinitial-sigpRfinal)/250.0;
    double sigpIinitial = -0.1;
    double sigpIfinal = 0.1;
    double delsigpI = abs(sigpIinitial - sigpIfinal)/250.0;
    

     

    string filename =   "kernelcontour_a_" + to_string((int)a) 
                        + "_s_" + to_string(s) 
                        + "_sigpsb_" + to_string(real(sigmab(a,m))) 
                        + "_epspow_" + to_string((int)abs(log10(eps)))
                        + "_pow2eneps.dat";

    int size = 250;
    int isize = size;
    int jsize = size;

    fout.open(filename.c_str());

    //for(double sigmapR=sigpRinitial;sigmapR<sigpRfinal;sigmapR=sigmapR + delsigpR)
    for(int i=0;i<isize;++i)
    {
        
        //for(double sigmapI=sigpIinitial;sigmapI<sigpIfinal;sigmapI=sigmapI + delsigpI)
        for(int j=0;j<jsize;++j)
        //double someeps = 1.0e-3;
        {
            comp sigmapR = (comp)sigpRinitial + (comp)i*delsigpR;
            comp sigmapI = (comp)sigpIinitial + (comp)j*delsigpI; 
            comp sigmapcomp = (comp) sigmapR + ii*sigmapI;

            comp ope = GS(scomp,sigmapcomp,sigmapprime,m,eps)*tau1(scomp,sigmapcomp,m)*M2kfunc(a,sigmapcomp,m,eps);
            //GS(s,sigb,sigvec[i],m,eps)*tau1(s,sigvec[i],m)*M2kfunc(a,sigvec[i],m,eps);

            fout<<setprecision(16)<<real(sigmapR)<<'\t'<<setprecision(16)<<real(sigmapI)<<'\t'<<setprecision(16)<<real(ope)<<'\t'<<setprecision(16)<<imag(ope)<<endl;
            cout<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
            //printf("%.16f\t %.16f\t %.16f\t %.16f\n",sigmapR,sigmapI,real(ope),imag(ope));
            
        }
        fout<<endl;
    }   

    fout.close();

    string filename1 = "thresh_" + filename;//OPEcontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename1.c_str());

    //fout<<real(sigpplus(s,sigmab(a,m),m))<<'\t'<<real(sigpminus(s,sigmab(a,m),m))<<'\t'<<real(sigc1(s,sigmab(a,m),m))<<'\t'<<real(sigc2(s,sigmab(a,m),m))<<'\t'<<pow(sqrt(s)-m,2.0)<<'\t'<<0<<'\t'<<0<<endl;
    fout<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<(sqrt(s)-m)*(sqrt(s)-m)<<'\t'<<0<<endl;
    fout.close();

    double eps1 = abs((1.0e-2)*eps*eps_energy_factor_minus(s,sigmab(a,m),m));
    //string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";
    
    for(int box=1;box<=100000;box=box*10)
    {
        string filename2 = "sigvec_box_" + to_string(box) + "_"  + filename;
        fout.open(filename2.c_str());

        vector<comp> sigvec;
        sigma_vector_maker_5(sigvec,s,0.0,sigmax(s,m),a,m,eps1,5000,500,box);
        for(int i=0;i<sigvec.size();++i)
        {
            fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;
            //fout<<endl;
        }

        fout.close();
    }

}


void kernelplotcontour_withintloop_vs_sigvec_delta()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    double s = 8.72;
    double eps = 1.0e-8;
    double imS = 0.07;
    comp scomp = (comp) s;//+ ii*imS;
    comp sigb = sigmab(a,m) ;//+ ii*0.02;

    /*
    double sigpRinitial = 3.6;
    double sigpRfinal = 4.0;
    double delsigpR = abs(sigpRinitial-sigpRfinal)/250.0;
    double sigpIinitial = -0.1;
    double sigpIfinal = 0.1;
    double delsigpI = abs(sigpIinitial - sigpIfinal)/250.0;
    */

     

    string filename =   "kernel_delta_sigvec_contour_a_" + to_string((int)a) 
                        + "_s_" + to_string(s) 
                        + "_sigpsb_" + to_string(real(sigmab(a,m))) 
                        + "_epspow_" + to_string((int)abs(log10(eps)))
                        + ".dat";

    int size = 500;
    int isize = size;
    int jsize = 4*size;

    double initialdelta = 1.0e-9;
    double finaldelta = 1.0;
    double deldelta = abs(finaldelta-initialdelta)/((double) size);

    //double initial
    


    fout.open(filename.c_str());

    //for(double sigmapR=sigpRinitial;sigmapR<sigpRfinal;sigmapR=sigmapR + delsigpR)
    for(int i=0;i<isize;++i)
    {
        double eps1 = initialdelta + ((double) i)*deldelta;
        double delta = abs((eps1)*eps*eps_energy_factor_minus(s,sigmab(a,m),m));
        
        vector<comp> sigvec;
        sigma_vector_maker_6(sigvec,s,0.0,sigmax(s,m),a,m,delta,500,250,10);
        
        //for(double sigmapI=sigpIinitial;sigmapI<sigpIfinal;sigmapI=sigmapI + delsigpI)
        for(int j=0;j<jsize;++j)
        //double someeps = 1.0e-3;
        {
            /*comp sigmapR = (comp)sigpRinitial + (comp)i*delsigpR;
            comp sigmapI = (comp)sigpIinitial + (comp)j*delsigpI; 
            comp sigmapcomp = (comp) sigmapR + ii*sigmapI;
            */

            comp ope = GS(scomp,sigb,sigvec[j],m,eps)*tau1(scomp,sigvec[j],m)*M2kfunc(a,sigvec[j],m,eps);
            //GS(s,sigb,sigvec[i],m,eps)*tau1(s,sigvec[i],m)*M2kfunc(a,sigvec[i],m,eps);

            //fout<<setprecision(16)<<real(sigmapR)<<'\t'<<setprecision(16)<<real(sigmapI)<<'\t'<<setprecision(16)<<real(ope)<<'\t'<<setprecision(16)<<imag(ope)<<endl;

            //cout<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
            //printf("%.16f\t %.16f\t %.16f\t %.16f\n",sigmapR,sigmapI,real(ope),imag(ope));

            fout<<setprecision(16)<<eps1<<'\t'
                <<setprecision(16)<<delta<<'\t'
                <<setprecision(16)<<j<<'\t'
                <<setprecision(16)<<real(ope)<<'\t'
                <<setprecision(16)<<imag(ope)<<endl;

            
            
            cout<<"delta:"<<delta<<'\t'
                <<"sigvec:"<<sigvec[j]<<'\t'
                <<"ker:"<<ope<<endl;
            
        }
        fout<<endl;
    }   

    fout.close();

    

}

void tauplotcontour_withintloop()
{
    ofstream fout;

    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    double s = 8.6;
    double eps = 0.0;//1.0e-7;
    double imS = 0.03;
    comp scomp = (comp) s+ ii*imS;
    comp sigmapprime = sigmab(a,m) + ii*0.02;

    double sigpRinitial = 0;
    double sigpRfinal = 20.0;
    double delsigpR = abs(sigpRinitial-sigpRfinal)/250.0;
    double sigpIinitial = -0.1;
    double sigpIfinal = 0.1;
    double delsigpI = abs(sigpIinitial - sigpIfinal)/250.0;
    

     

    string filename = "TAUcontour_a_" + to_string((int)a) + "_s+i_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    int size = 250;
    int isize = size;
    int jsize = size;

    fout.open(filename.c_str());

    //for(double sigmapR=sigpRinitial;sigmapR<sigpRfinal;sigmapR=sigmapR + delsigpR)
    for(int i=0;i<isize;++i)
    {
        
        //for(double sigmapI=sigpIinitial;sigmapI<sigpIfinal;sigmapI=sigmapI + delsigpI)
        for(int j=0;j<jsize;++j)
        //double someeps = 1.0e-3;
        {
            comp sigmapR = (comp)sigpRinitial + (comp)i*delsigpR;
            comp sigmapI = (comp)sigpIinitial + (comp)j*delsigpI; 
            comp sigmapcomp = (comp) sigmapR + ii*sigmapI;

            comp ope = tau2(scomp,sigmapcomp,m);//GS(scomp,sigmapcomp,sigmapprime,m,eps);

            fout<<setprecision(16)<<real(sigmapR)<<'\t'<<setprecision(16)<<real(sigmapI)<<'\t'<<setprecision(16)<<real(ope)<<'\t'<<setprecision(16)<<imag(ope)<<endl;
            cout<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
            //printf("%.16f\t %.16f\t %.16f\t %.16f\n",sigmapR,sigmapI,real(ope),imag(ope));
            
        }
        fout<<endl;
    }   

    fout.close();
    cout<<(sqrt(scomp)-m)*(sqrt(scomp)-m)<<'\t'<<(sqrt(scomp)+m)*(sqrt(scomp)+m)<<endl;

    string filename1 = "thresh" + filename;//OPEcontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename1.c_str());
    fout<<real((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\t'<<imag((sqrt(scomp)-m)*(sqrt(scomp)-m))<<'\t'<<real((sqrt(scomp)+m)*(sqrt(scomp)+m))<<'\t'<<imag((sqrt(scomp)+m)*(sqrt(scomp)+m))<<endl;

    //fout<<real(sigpplus(s,sigmab(a,m),m))<<'\t'<<real(sigpminus(s,sigmab(a,m),m))<<'\t'<<real(sigc1(s,sigmab(a,m),m))<<'\t'<<real(sigc2(s,sigmab(a,m),m))<<'\t'<<pow(sqrt(s)-m,2.0)<<'\t'<<0<<'\t'<<0<<endl;
    /*fout<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<endl;
    */
    fout.close();

    //string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    //fout.open(filename2.c_str());

    //vector<comp> sigvec;
    //sigma_vector_maker(sigvec,s,0.0,sigmax(s,m),a,m,1000.0,1.0e-2);
    //for(int i=0;i<sigvec.size();++i)
    //fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;

    //fout.close();

}


void kernelplotcontour()
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
    double delsigpR = abs(sigpRinitial-sigpRfinal)/500.0;
    double sigpIinitial = -0.1;
    double sigpIfinal = 0.1;
    double delsigpI = abs(sigpIinitial - sigpIfinal)/500.0;
    

     

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

    string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename2.c_str());

    vector<comp> sigvec;
    //sigma_vector_maker(sigvec,s,0.0,sigmax(s,m),a,m,1000.0,1.0e-5);
    for(int i=0;i<sigvec.size();++i)
    fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;

    fout.close();
}

void kernelplotcontour_vs_s()
{
    ofstream fout;

    double a = 16.0;
    double m = 1.0;
    double s = 8.65;
    double eps = 1.0e-5;
    comp scomp = (comp) s;
    comp sigmapprime = sigmab(a,m);
    double sigpRinitial = 3.55;
    double sigpRfinal = 5.0;
    double delsigpR = abs(sigpRinitial-sigpRfinal)/500.0;
    double sigpIinitial = -0.1;
    double sigpIfinal = 0.1;
    double delsigpI = abs(sigpIinitial - sigpIfinal)/500.0;
    

     

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

    string filename2 = "sigveccontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename2.c_str());

    vector<comp> sigvec;
    sigma_vector_maker(sigvec,s,0.0,sigmax(s,m),a,m,1000.0,1.0e-2);
    for(int i=0;i<sigvec.size();++i)
    fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;

    fout.close();
}

void mommaker()
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

double sbfuncconst(	double a,
					double r,
					double m		)
{
	double k2k = 1.0/a + r/(2.0*a*a);

	return 4.0*(-k2k*k2k + m*m);
}

comp KallenTriangleComplex(	comp x,
							comp y,
							comp z	)
{

	comp A = x*x + y*y + z*z - 2.0*(x*y + y*z + z*x);
	
	return A;
}

//This gives out the Boundstate Momentum
comp qfunc(     double E,
                double a,
                double r,
                double m     )
{
    comp Ecomp = (comp) E;
    comp sb = (comp) sbfuncconst(a,r,m);
    comp mcomp = (comp) m;
	comp sqt = sqrt(KallenTriangleComplex(Ecomp*Ecomp,sb,mcomp*mcomp));
    comp q1 = sqt/(2.0*Ecomp);

    return q1;
}


comp OPE_before_PW( comp s,
                    comp p,
                    comp k,
                    comp x,
                    double m,
                    double eps  )
{
    comp ii = {0.0,1.0};
    comp sigp = sigma_p(s,p,m);
    comp sigk = sigma_p(s,k,m);
    comp A = pow(sqrt(s)-omega_comp(p,m)-omega_comp(k,m),2.0);
    return -Hfunc(sigp,sigk,m)/(A - p*p - k*k - 2.0*p*k*x - m*m + ii*eps);
}

void only_pcut_x()
{
    double a = 16.0;
    double m = 1.0;
    double sreal = 8.72;
    double simag = 0.1;
    comp ii = {0.0,1.0};
    //comp s = sreal + ii*simag;
    double eps = 0.0;

    comp sigb = sigmab(a,m);
    //comp q = pmom(s,sigb-ii*eps,m);
    //comp kmax = pmom(s,0.0,m);


    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;

    double delxreal = abs(xrealinitial-xrealfinal)/(double)points;
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double sinitial = 8.72;
    double sfinal = (double)real(phibthreshold(a,m));
    double dels = abs(sinitial-sfinal)/10.0;

    int count = 1;
    for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    {
        comp s = s3 + ii*simag;
        comp q = pmom(s,sigb-ii*eps,m);
    
        ofstream fout;
        string filename1 = "pcutonly_test_a_" + to_string(a) 
                          + "_eps_" + to_string(eps)
                          + "_count_" + to_string(count)
                          + ".dat";

        fout.open(filename1.c_str());                   
        for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        {
            comp pcutp = pcut_plus_comp(x,s,q,m,eps);
            comp pcutm = pcut_minus_comp(x,s,q,m,eps);
            fout<<x<<'\t'<<real(pcutp)<<'\t'<<imag(pcutp)<<'\t'<<real(pcutm)<<'\t'<<imag(pcutm)<<'\t'<<real(q)<<'\t'<<imag(q)<<endl;
        }

        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;
        cout<<"file:"<<filename1<<endl;
        count = count + 1;
    }
}

void only_pcut_x_comp()
{
    double a = 16.0;
    double m = 1.0;
    double sreal = 8.72;
    double simag = -0.1;
    comp ii = {0.0,1.0};
    //comp s = sreal + ii*simag;
    double eps = 0.0;

    comp sigb = sigmab(a,m);
    //comp q = pmom(s,sigb-ii*eps,m);
    //comp kmax = pmom(s,0.0,m);


    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 2000;

    double delxreal = abs(xrealinitial-xrealfinal)/(double)points;
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double sinitial = 8.72;
    double sfinal = (double)real(phibthreshold(a,m));
    double dels = abs(sinitial-sfinal)/10.0;

    int count = 1;
    for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    {
        comp s = s3 + ii*simag;
        comp q = pmom(s,sigb-ii*eps,m);
        double jmultiplier = 1.0;
    
        ofstream fout;
        string filename1 = "pcutonly_test_a_" + to_string(a) 
                          + "_eps_" + to_string(eps)
                          + "_count_" + to_string(count)
                          + ".dat";

        fout.open(filename1.c_str());                   
        for(double j=0.0;j<=jmultiplier;j=j+0.001)
        {
            comp x = -1.0 - ii*j;
            comp pcutp = pcut_plus_comp(x,s,q,m,eps);
            comp pcutm = pcut_minus_comp(x,s,q,m,eps);
            fout<<real(x)<<'\t'<<imag(x)<<'\t'<<real(pcutp)<<'\t'<<imag(pcutp)<<'\t'<<real(pcutm)<<'\t'<<imag(pcutm)<<'\t'<<real(q)<<'\t'<<imag(q)<<endl;
        }
        for(double xval=xrealinitial;xval<=xrealfinal;xval=xval+delxreal)
        {
            comp x = - ii*jmultiplier + xval;
            comp pcutp = pcut_plus_comp(x,s,q,m,eps);
            comp pcutm = pcut_minus_comp(x,s,q,m,eps);
            fout<<real(x)<<'\t'<<imag(x)<<'\t'<<real(pcutp)<<'\t'<<imag(pcutp)<<'\t'<<real(pcutm)<<'\t'<<imag(pcutm)<<'\t'<<real(q)<<'\t'<<imag(q)<<endl;
        }
        for(double j=0.0;j<=jmultiplier;j=j+0.001)
        {
            comp x = 1.0 - ii*jmultiplier + ii*j;
            
            comp pcutp = pcut_plus_comp(x,s,q,m,eps);
            comp pcutm = pcut_minus_comp(x,s,q,m,eps);
            fout<<real(x)<<'\t'<<imag(x)<<'\t'<<real(pcutp)<<'\t'<<imag(pcutp)<<'\t'<<real(pcutm)<<'\t'<<imag(pcutm)<<'\t'<<real(q)<<'\t'<<imag(q)<<endl;
        }

        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;
        cout<<"file:"<<filename1<<endl;
        count = count + 1;
    }
}



/*  

    This is the massive data generator to perform Glockle-like 
    tests. For a given s3, we plot the contour plot of the kernel
    with one of the momentum set at the bound state pole. We then
    make a contour vector q_vec, that we also plot in the contour. 
    For each momentum in the q_vec, we generate the pcut data. We 
    repeat this for consecutive s3's as well. All of this will be 
    plotted using the contour_maker.ipynb notebook.

*/
void kernel_pcut_x_glockletest_vs_s3_qvec()
{
    ofstream fout;
    double a = 16.0;
    double m = 1.0;
    double sreal = 8.72;
    double simag = -0.01;
    comp ii = {0.0,1.0};
    double eps = 0.01;
    double eps1 = 0.01;

    comp sigb = sigmab(a,m);
    


    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int delspoints = 10;
    int qvecpoints = 50;
    double qvec_r = 0.01;

    double delxreal = abs(xrealinitial-xrealfinal)/(double)points;
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double sinitial = 8.72;
    double sfinal = (double)real(phibthreshold(a,m));
    double dels = abs(sinitial-sfinal)/(double)delspoints;

    int scount = 0;
    for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    {
        int pcount = 0;
        comp s = s3 + ii*simag;
        comp q = pmom(s,sigb-ii*eps1,m);
        comp kmax = pmom(s,0.0,m);

    
        cout<<"s = "<<s<<'\t'<<"dels = "<<dels<<endl;

        vector<comp> qvec;
        mom_vector_maker_3(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        


        //This part writes the qvec in to a file //
        string qvecfile =   "qvec_a_" + to_string(a)
                          + "_eps_" + to_string(eps)
                          + "_scount_" + to_string(scount)
                          + ".dat";


        fout.open(qvecfile.c_str());
        for(int i=0;i<qvec.size();++i)
        {
            fout<<i<<'\t'<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
        }
        fout.close();
        cout<<"qvecfile : "<<qvecfile<<endl;

        //This part writes the kernel with k=q
        double pRinitial = -0.2;
        double pRfinal = 0.2;
        double delpR = abs(pRinitial-pRfinal)/points;
        double pIinitial = -0.5;
        double pIfinal = 0.5;
        double delpI = abs(pIinitial - pIfinal)/points;
        int isize = (int) points;
        int jsize = (int) points;

        string kernelfile =   "kernel_a_" + to_string(a)
                            + "_eps_" + to_string(eps)
                            + "_scount_" + to_string(scount)
                            + ".dat";

        fout.open(kernelfile.c_str());
        for(int i=0;i<isize;++i)
        {
            for(int j=0;j<jsize;++j)
            {
                comp pR = (comp)pRinitial + (comp)i*delpR;
                comp pI = (comp)pIinitial + (comp)j*delpI; 
                comp pcomp = (comp) pR + ii*pI;

                comp sigp = sigma_p(s,pcomp,m);

                double pi = acos(-1.0);

                comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            
                fout<<setprecision(16)<<real(pR)<<'\t'
                    <<setprecision(16)<<real(pI)<<'\t'
                    <<setprecision(16)<<real(ope)<<'\t'
                    <<setprecision(16)<<imag(ope)<<endl;
            
            }
        
        }
        fout.close();
        cout<<"kernelfile : "<<kernelfile<<endl;

        //Here we print the relevant thresholds for 
        //for the kernel and trace the pcut for the OPE

        string thresholdfile = "threshold_a_" + to_string(a)
                            + "_eps_" + to_string(eps)
                            + "_scount_" + to_string(scount)
                            + ".dat";

        
        fout.open(thresholdfile.c_str());
        fout<<real(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<'\t'
            <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<dels<<endl;

        fout.close();
        cout<<"thresholdfile : "<<thresholdfile<<endl;

        string opecutfile = "opecuttracer_a_" + to_string(a) 
                          + "_eps_" + to_string(eps)
                          + "_scount_" + to_string(scount)
                          + ".dat";

        fout.open(opecutfile.c_str());                   
        for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        {
            comp pcutp = pcut_plus_comp(x,s,q,m,eps);
            comp pcutm = pcut_minus_comp(x,s,q,m,eps);
            fout<<x<<'\t'
                <<real(pcutp)<<'\t'
                <<imag(pcutp)<<'\t'
                <<real(pcutm)<<'\t'
                <<imag(pcutm)<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<abs(-q-pcutp)<<'\t'
                <<abs(q-pcutm)<<endl;
        }
        fout.close();

        cout<<"opecutfile : "<<opecutfile<<endl;

        for(int k=0;k<qvec.size();++k)
        {
            string pcutfile = "pcut_a_" + to_string(a) 
                          + "_eps_" + to_string(eps)
                          + "_scount_" + to_string(scount)
                          + "_pcount_" + to_string(pcount)
                          + ".dat";

            fout.open(pcutfile.c_str());                   
            for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
            {
                comp pcutp = pcut_plus_comp(x,s,qvec[k],m,eps);
                comp pcutm = pcut_minus_comp(x,s,qvec[k],m,eps);
                fout<<x<<'\t'
                    <<real(pcutp)<<'\t'
                    <<imag(pcutp)<<'\t'
                    <<real(pcutm)<<'\t'
                    <<imag(pcutm)<<'\t'
                    <<real(q)<<'\t'
                    <<imag(q)<<'\t'
                    <<abs(-q-pcutp)<<'\t'
                    <<abs(q-pcutm)<<endl;
            }
            fout.close();
            cout<<"qvec["<<k<<"] = "<<qvec[k]<<'\t'<<"pcount = "<<pcount<<endl;
            cout<<"pcutfile : "<<pcutfile<<endl;

            //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;
            //cout<<"file:"<<filename1<<endl;
            pcount = pcount + 1;
        }
        scount = scount + 1;
    }
}


//This momentum vector makers are defined according
//to sebastians definitions found in the paper
void mom_vector_maker_seba_imsneg(  vector<comp> &qvec,
                                    vector<comp> &weights, 
                                    comp s,
                                    comp kmin,
                                    comp kmax,
                                    double a,
                                    double m,
                                    double eps,
                                    double eps_for_m2k,
                                    double points   )
{
    //initial setup//
    comp ii = {0.0,1.0};
    comp q = pmom(s,sigmab(a,m),m);
    comp qc2 = q_c2(s,a,m,eps);
    comp pplus = pcut_plus(1.0,s,q,m,eps);
    comp pminus = pcut_plus(-1.0,s,q,m,eps);
    comp qc1 = q_c1(s,a,m,eps);
    double rad = real(qc1);
    //cout<<"radius = "<<rad<<endl;
    
    if(rad>0.0)
    {
        //cout<<"contour chosen"<<endl;
        //firstnode - secondnode
        comp firstnode = kmin;
        comp secondnode = -2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2);
        line_maker_with_weights(qvec,weights,firstnode,secondnode,2.0*points/10.0);

        //secondnode - thirdnode
        comp thirdnode = 0.5*real(pplus - pminus) - abs(qc2)*ii;
        line_maker_with_weights(qvec,weights,secondnode,thirdnode,2.0*points/10.0);

        //thirdnode - fourthnode 
        comp fourthnode = (2.0/3.0)*(1.0 - 2.0*ii)*abs(qc2);
        line_maker_with_weights(qvec,weights,thirdnode,fourthnode,2.0*points/10.0);

        //fourthnode - fifthnode
        comp fifthnode = (3.0/2.0)*abs(qc2);
        line_maker_with_weights(qvec,weights,fourthnode,fifthnode,2.0*points/10.0);

        //fifthnode - sixthnode 
        comp sixthnode = kmax;
        line_maker_with_weights(qvec,weights,fifthnode,sixthnode,2.0*points/10.0);
    }
    else
    {
        //cout<<"straight line chosen"<<endl;
        line_maker_with_weights(qvec,weights,kmin,kmax,points);     
    }

}

//This momentum vector makers are defined according
//to sebastians definitions found in the paper
void mom_vector_maker_seba_imspos(  vector<comp> &qvec,
                                    vector<comp> &weights, 
                                    comp s,
                                    comp kmin,
                                    comp kmax,
                                    double a,
                                    double m,
                                    double eps,
                                    double eps_for_m2k,
                                    double points         )
{
    int tag1 =0;
                                    int tag2 = 0;
                                    int switch_for_gvec_fixer = 0;
    //initial setup//
    comp ii = {0.0,1.0};
    comp q = pmom(s,sigmab(a,m),m);
    comp qc2 = q_c2(s,a,m,eps);
    comp pplus = pcut_plus(1.0,s,q,m,eps);
    comp pminus = pcut_plus(-1.0,s,q,m,eps);
    comp qc1 = q_c1(real(s),a,m,eps);
    comp sigb = sigmab(a,m);
    double rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    
    if(rad>0.0)
    {
        cout<<"contour chosen"<<endl;
        //firstnode - secondnode
        comp firstnode = kmin;
        comp secondnode = -2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2);
        line_maker_with_weights(qvec,weights,firstnode,secondnode,4.0*points/20.0);

        //secondnode - thirdnode
        comp sigmac2 = sigc2(s,sigb,m);
        comp p = pmom(s,sigmac2,m);
        comp k = q;
        comp alphapk = ((sqrt(s)-omega_comp(p,m)-omega_comp(k,m))*(sqrt(s)-omega_comp(p,m)-omega_comp(k,m)) - p*p - k*k - m*m)/(2.0*p*k);
        //comp zpk = 2.0*s*(sigb + ii*eps) - (s+sigb - m*m)*(s+m*m-sigmac2);
        //comp zpk = 2.0*s*(sigmac2 + ii*eps) - (s+sigmac2 - m*m)*(s+m*m-sigb);
        
        double x0 = -1.0*real(alphapk);
        //double x0 = real(zpk);
        cout<<"sigma c2 = "<<sigmac2<<endl;
        
        cout<<"x0 = "<<x0<<endl;
        comp thirdnode = -pcut_minus(-x0,s,q,m,eps);
        line_maker_with_weights(qvec,weights,secondnode,thirdnode,4.0*points/20.0);

        tag1 = qvec.size() - 1;
        //thirdnode - fourthnode 
        comp fourthnode = pcut_minus(x0,s,q,m,eps);
        line_maker_with_weights(qvec,weights,thirdnode,fourthnode,1.0*points/20.0);

        tag2 = qvec.size() - 1;
        //fourthnode - fifthnode
        comp fifthnode = (2.0/3.0)*(1.0-2.0*ii)*abs(qc2);
        line_maker_with_weights(qvec,weights,fourthnode,fifthnode,3.0*points/20.0);

        //fifthnode - sixthnode 
        comp sixthnode = (3.0/2.0)*abs(qc2);
        line_maker_with_weights(qvec,weights,fifthnode,sixthnode,4.0*points/20.0);

        //sixthnode - seventhnode
        comp seventhnode = kmax;
        line_maker_with_weights(qvec,weights,sixthnode,seventhnode,4.0*points/20.0);

        switch_for_gvec_fixer = 0;

    }
    else
    {
        cout<<"straight line chosen"<<endl;
        line_maker_with_weights(qvec,weights,kmin,kmax,points);   
        switch_for_gvec_fixer = 1;  
    }

}


//instead of going through the lowest point between
//the two cuts, it goes through the cut from little bit
//above

void mom_vector_maker_seba_imspos_1(  vector<comp> &qvec,
                                    vector<comp> &weights, 
                                    comp s,
                                    comp kmin,
                                    comp kmax,
                                    double a,
                                    double m,
                                    double eps,
                                    double eps_for_m2k,
                                    double points   )
{
    //initial setup//
    comp ii = {0.0,1.0};
    comp q = pmom(s,sigmab(a,m),m);
    comp qc2 = q_c2(s,a,m,eps);
    comp pplus = pcut_plus(1.0,s,q,m,eps);
    comp pminus = pcut_plus(-1.0,s,q,m,eps);
    comp qc1 = q_c1(s,a,m,eps);
    comp sigb = sigmab(a,m);
    double rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    
    if(rad>0.0)
    {
        cout<<"contour chosen"<<endl;
        //firstnode - secondnode
        comp firstnode = kmin;
        comp secondnode = -2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2);
        
        //secondnode - thirdnode
        comp sigmac2 = sigc2(s,sigb,m);
        comp p = pmom(s,sigmac2,m);
        comp k = q;
        comp alphapk = ((sqrt(s)-omega_comp(p,m)-omega_comp(k,m))*(sqrt(s)-omega_comp(p,m)-omega_comp(k,m)) - p*p - k*k - m*m)/(2.0*p*k);
        //comp zpk = 2.0*s*(sigb + ii*eps) - (s+sigb - m*m)*(s+m*m-sigmac2);
        //comp zpk = 2.0*s*(sigmac2 + ii*eps) - (s+sigmac2 - m*m)*(s+m*m-sigb);
        
        double x0 = -1.0*real(alphapk);
        //double x0 = real(zpk);
        cout<<"sigma c2 = "<<sigmac2<<endl;
        
        cout<<"x0 = "<<x0<<endl;
        comp thirdnode = -pcut_minus(-x0,s,q,m,eps);//-pcut_minus((-x0-1.0)/2.0,s,q,m,eps);
        
        //tag1 = qvec.size() - 1;
        //thirdnode - fourthnode 
        comp fourthnode = pcut_minus(x0,s,q,m,eps);
        
        comp tempsecondnode;

        int tempcount1, tempcount2, tempcount3;
        tempcount1 = 1;
        tempcount2 = 1;
        tempcount3 = 1;
        //0 means active, 1 means off 

        if(imag(thirdnode)>imag(secondnode))
        {
            tempsecondnode = real(secondnode) + ii*imag(-secondnode);
            if(real(thirdnode)>real(secondnode))
            {
                secondnode = -secondnode;
                tempcount1 = 0;
                tempcount3 = 0; 
            }
            else 
            {
                secondnode = tempsecondnode;
                tempcount1 = 1;
                tempcount3 = 0;
            }
        }

        if(tempcount3==0)
        {
            if(tempcount1==0)
            {
                line_maker_with_weights(qvec,weights,firstnode,secondnode,1.0*points/10.0);
                comp intermediatenode = -secondnode;
                
                line_maker_with_weights(qvec,weights,secondnode,thirdnode,1.0*points/10.0);
                
                line_maker_with_weights(qvec,weights,thirdnode,intermediatenode,2.0*points/10.0);
                
                line_maker_with_weights(qvec,weights,intermediatenode,fourthnode,2.0*points/10.0);

            }
            else
            {
                line_maker_with_weights(qvec,weights,firstnode,secondnode,1.0*points/10.0);
                comp intermediatenode = real(tempsecondnode) + ii*imag(-tempsecondnode);
                line_maker_with_weights(qvec,weights,secondnode,thirdnode,1.0*points/10.0);
                
                line_maker_with_weights(qvec,weights,thirdnode,intermediatenode,2.0*points/10.0);

                line_maker_with_weights(qvec,weights,intermediatenode,fourthnode,2.0*points/10.0);

        
                
            }
            
        }
        else
        {
            line_maker_with_weights(qvec,weights,firstnode,secondnode,2.0*points/10.0);

            line_maker_with_weights(qvec,weights,secondnode,thirdnode,2.0*points/10.0);

            line_maker_with_weights(qvec,weights,thirdnode,fourthnode,2.0*points/10.0);

        
        }
        
        
        //tag2 = qvec.size() - 1;
        //fourthnode - fifthnode
        comp fifthnode = (2.0/3.0)*(1.0-2.0*ii)*abs(qc2);//real(qc1) + ii*imag(fourthnode);
        line_maker_with_weights(qvec,weights,fourthnode,fifthnode,2.0*points/10.0);

        //fifthnode - sixthnode 
        comp sixthnode = (3.0/2.0)*abs(qc2);
        line_maker_with_weights(qvec,weights,fifthnode,sixthnode,1.0*points/10.0);

        //sixthnode - seventhnode
        comp seventhnode = kmax;
        line_maker_with_weights(qvec,weights,sixthnode,seventhnode,1.0*points/10.0);

    }
    else
    {
        cout<<"straight line chosen"<<endl;
        line_maker_with_weights(qvec,weights,kmin,kmax,points);     
    }

}

void mom_vector_maker_seba_imspos_2(  vector<comp> &qvec,
                                    vector<comp> &weights, 
                                    comp s,
                                    comp kmin,
                                    comp kmax,
                                    double a,
                                    double m,
                                    double eps,
                                    double eps_for_m2k,
                                    double points,
                                    int &tag1,
                                    int &tag2,
                                    int &switch_for_gvec_fixer   )
{
    //initial setup//
    comp ii = {0.0,1.0};
    comp q = pmom(s,sigmab(a,m),m);
    comp qc2 = q_c2(s,a,m,eps);
    comp pplus = pcut_plus(1.0,s,q,m,eps);
    comp pminus = pcut_plus(-1.0,s,q,m,eps);
    comp qc1 = q_c1(real(s),a,m,eps);
    comp sigb = sigmab(a,m);
    double rad = abs(real(qc1));
    cout<<"radius = "<<rad<<endl;
    
    if(rad>0.0)
    {
        switch_for_gvec_fixer = 0;
        cout<<"contour chosen"<<endl;
        //firstnode - secondnode
        comp firstnode = kmin;
        comp secondnode = -rad/3.0;//-2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2);
        //1
        line_maker_with_weights(qvec,weights,firstnode,secondnode,1.0*points/4.0);

        //secondnode - thirdnode
        comp sigmac2 = sigc2(s,sigb,m);
        comp p = (q_c1(s,a,m,eps) + qc2)/2.0;//pmom(s,sigmac2,m);
        comp k = q;
        comp alphapk = ((sqrt(s)-omega_comp(p,m)-omega_comp(k,m))*(sqrt(s)-omega_comp(p,m)-omega_comp(k,m)) - p*p - k*k - m*m)/(2.0*p*k);
        //comp zpk = 2.0*s*(sigb + ii*eps) - (s+sigb - m*m)*(s+m*m-sigmac2);
        //comp zpk = 2.0*s*(sigmac2 + ii*eps) - (s+sigmac2 - m*m)*(s+m*m-sigb);
        
        double x0 = -1.0*real(alphapk);
        //double x0 = real(zpk);
        cout<<"sigma c2 = "<<sigmac2<<endl;
        
        cout<<"x0 = "<<x0<<endl;
        //comp thirdnode = -pcut_minus((-x0-0.77)/2.0,s,q,m,eps);
        
        
        //tag1 = qvec.size() - 1;
        //thirdnode - fourthnode 
        comp f1 = pcut_minus((x0+1.0)/2.0,s,q,m,eps);
        
        comp f2 = 0.5*(real(pplus - pminus)) - ii*abs(qc2);

        comp temp_f3_1 = pcut_minus(-1,s,q,m,eps);
        comp temp_f3_2 = pcut_minus(1,s,q,m,eps);
        comp temp_f3_3 = pcut_plus(-1,s,q,m,eps);
        comp temp_f3_4 = pcut_plus(1,s,q,m,eps);
        comp f3 = 0.5*(real(pplus - pminus)) - ii*5.0/4.0*abs(qc2);;
        /*if(imag(temp_f3_1)<imag(temp_f3_2) && imag(temp_f3_1)<imag(temp_f3_3) && imag(temp_f3_1)<imag(temp_f3_4))
        {
            f3 = pcut_minus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_2)<imag(temp_f3_1) && imag(temp_f3_2)<imag(temp_f3_3) && imag(temp_f3_2)<imag(temp_f3_4))
        {
            f3 = pcut_minus(1,s,q,m,eps);
        }
        else if(imag(temp_f3_3)<imag(temp_f3_1) && imag(temp_f3_3)<imag(temp_f3_2) && imag(temp_f3_3)<imag(temp_f3_4))
        {
            f3 = pcut_plus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_4)<imag(temp_f3_1) && imag(temp_f3_4)<imag(temp_f3_2) && imag(temp_f3_4)<imag(temp_f3_3))
        {
            f3 = pcut_plus(1,s,q,m,eps);
        }*/

        comp thirdnode = real(secondnode) + ii*imag((f2 + f3)/2.0);//-pcut_minus(-x0,s,q,m,eps);
        //2
        line_maker_with_weights(qvec,weights,secondnode,thirdnode,1.0*points/4.0);
        tag1 = qvec.size() - 1;
        
        //comp f3 = pcut_minus(1,s,q,m,eps);
        //comp fourthnode = real(thirdnode) + ii*imag((f2 + f3)/2.0);//real(thirdnode) + ii*imag(f1)*0.75;//pcut_minus(2.0*x0/3.0,s,q,m,eps);
        //3
        //line_maker_with_weights(qvec,weights,thirdnode,fourthnode,2.0*points/10.0);

        //tag2 = qvec.size() - 1;
        //fourthnode - fifthnode
        //comp fifthnode = (f2 + f3)/2.0;//real(f1) + ii*imag(f1)*0.75;//real(qc1) + ii*imag(fourthnode);//(2.0/3.0)*(1.0-2.0*ii)*abs(qc2);
        //4
        //line_maker_with_weights(qvec,weights,fourthnode,fifthnode,2.0*points/10.0);

        //fifthnode - sixthnode 
        comp sixthnode = real(q_c1(s,a,m,eps)) + (ii*imag((f2+f3)/2.0));//(3.0/2.0)*abs(qc2);
        //5
        //line_maker_with_weights(qvec,weights,fifthnode,sixthnode,1.0*points/4.0);
        line_maker_with_weights(qvec,weights,thirdnode,sixthnode,1.0*points/4.0);
        tag2 = qvec.size() - 1;
        //sixthnode - seventhnode
        comp seventhnode = kmax;
        //6
        line_maker_with_weights(qvec,weights,sixthnode,seventhnode,1.0*points/4.0);

    }
    else
    {
        switch_for_gvec_fixer = 1;
        cout<<"straight line chosen"<<endl;
        line_maker_with_weights(qvec,weights,kmin,kmax,points);     
    }

}


void mom_vector_maker_seba_imspos_2_with_contour47(  vector<comp> &qvec,
                                    vector<comp> &weights, 
                                    comp s,
                                    comp kmin,
                                    comp kmax,
                                    double a,
                                    double m,
                                    double eps,
                                    double eps_for_m2k,
                                    double points,
                                    int &tag1,
                                    int &tag2,
                                    int &switch_for_gvec_fixer   )
{
    //initial setup//
    comp ii = {0.0,1.0};
    comp q = pmom(s,sigmab(a,m),m);
    comp qc2 = q_c2(s,a,m,eps);
    comp pplus = pcut_plus(1.0,s,q,m,eps);
    comp pminus = pcut_plus(-1.0,s,q,m,eps);
    comp qc1 = q_c1(real(s),a,m,eps);
    comp sigb = sigmab(a,m);
    double rad = abs(real(qc1));
    cout<<"radius = "<<rad<<endl;

    //this section is only for looking into if the 
    //pcuts have crossed the real axis or not

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    if(axis_flag==0)
    //if(rad>0.0)
    {
        switch_for_gvec_fixer = 0;
        cout<<"contour chosen"<<endl;
        //firstnode - secondnode
        comp firstnode = kmin;
        comp secondnode = -rad/1.50;//-2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2);
        //1
        line_maker_with_weights(qvec,weights,firstnode,secondnode,1.0*points/4.0);

        //secondnode - thirdnode
        comp sigmac2 = sigc2(s,sigb,m);
        comp p = (q_c1(s,a,m,eps) + qc2)/2.0;//pmom(s,sigmac2,m);
        comp k = q;
        comp alphapk = ((sqrt(s)-omega_comp(p,m)-omega_comp(k,m))*(sqrt(s)-omega_comp(p,m)-omega_comp(k,m)) - p*p - k*k - m*m)/(2.0*p*k);
        //comp zpk = 2.0*s*(sigb + ii*eps) - (s+sigb - m*m)*(s+m*m-sigmac2);
        //comp zpk = 2.0*s*(sigmac2 + ii*eps) - (s+sigmac2 - m*m)*(s+m*m-sigb);
        
        double x0 = -1.0*real(alphapk);
        //double x0 = real(zpk);
        cout<<"sigma c2 = "<<sigmac2<<endl;
        
        cout<<"x0 = "<<x0<<endl;
        //comp thirdnode = -pcut_minus((-x0-0.77)/2.0,s,q,m,eps);
        
        
        //tag1 = qvec.size() - 1;
        //thirdnode - fourthnode 
        comp f1 = pcut_minus((x0+1.0)/2.0,s,q,m,eps);
        
        comp f2 = 0.5*(real(pplus - pminus)) - ii*abs(qc2);

        comp temp_f3_1 = pcut_minus(-1,s,q,m,eps);
        comp temp_f3_2 = pcut_minus(1,s,q,m,eps);
        comp temp_f3_3 = pcut_plus(-1,s,q,m,eps);
        comp temp_f3_4 = pcut_plus(1,s,q,m,eps);
        comp f3 = 0.5*(real(pplus - pminus)) - ii*5.0/4.0*abs(qc2);;
        /*if(imag(temp_f3_1)<imag(temp_f3_2) && imag(temp_f3_1)<imag(temp_f3_3) && imag(temp_f3_1)<imag(temp_f3_4))
        {
            f3 = pcut_minus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_2)<imag(temp_f3_1) && imag(temp_f3_2)<imag(temp_f3_3) && imag(temp_f3_2)<imag(temp_f3_4))
        {
            f3 = pcut_minus(1,s,q,m,eps);
        }
        else if(imag(temp_f3_3)<imag(temp_f3_1) && imag(temp_f3_3)<imag(temp_f3_2) && imag(temp_f3_3)<imag(temp_f3_4))
        {
            f3 = pcut_plus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_4)<imag(temp_f3_1) && imag(temp_f3_4)<imag(temp_f3_2) && imag(temp_f3_4)<imag(temp_f3_3))
        {
            f3 = pcut_plus(1,s,q,m,eps);
        }*/

        comp thirdnode = real(secondnode) + ii*imag((f2 + f3)/2.0);//-pcut_minus(-x0,s,q,m,eps);
        //2
        line_maker_with_weights(qvec,weights,secondnode,thirdnode,1.0*points/4.0);
        tag1 = qvec.size() - 1;
        
        //comp f3 = pcut_minus(1,s,q,m,eps);
        //comp fourthnode = real(thirdnode) + ii*imag((f2 + f3)/2.0);//real(thirdnode) + ii*imag(f1)*0.75;//pcut_minus(2.0*x0/3.0,s,q,m,eps);
        //3
        //line_maker_with_weights(qvec,weights,thirdnode,fourthnode,2.0*points/10.0);

        //tag2 = qvec.size() - 1;
        //fourthnode - fifthnode
        //comp fifthnode = (f2 + f3)/2.0;//real(f1) + ii*imag(f1)*0.75;//real(qc1) + ii*imag(fourthnode);//(2.0/3.0)*(1.0-2.0*ii)*abs(qc2);
        //4
        //line_maker_with_weights(qvec,weights,fourthnode,fifthnode,2.0*points/10.0);

        //fifthnode - sixthnode 
        comp sixthnode = real(q_c1(s,a,m,eps)) + (ii*imag((f2+f3)/2.0));//(3.0/2.0)*abs(qc2);
        //5
        //line_maker_with_weights(qvec,weights,fifthnode,sixthnode,1.0*points/4.0);
        line_maker_with_weights(qvec,weights,thirdnode,sixthnode,1.0*points/4.0);
        tag2 = qvec.size() - 1;
        //sixthnode - seventhnode
        comp seventhnode = kmax;
        //6
        line_maker_with_weights(qvec,weights,sixthnode,seventhnode,1.0*points/4.0);

    }
    else
    {
        switch_for_gvec_fixer = 1;
        cout<<"straight line chosen"<<endl;
        line_maker_with_weights(qvec,weights,kmin,kmax,points);     
    }

}

void rightmost_left_pcut(   vector<comp> &pcutpvec,
                            vector<comp> &pcutmvec,
                            int flag_pcut,
                            comp &value    )
{
    comp ii = {0.0,1.0};
    double real_part_pcut = 0.0;
    double imag_part_pcut = 0.0;
    double temp_real_pcut = 0.0;
    double temp_imag_pcut = 0.0;
    double actual_real_pcut = 0.0;
    double actual_imag_pcut = 0.0;
    //if flag_pcut==0, then pcutm is to the left
    //if flag_pcut==1, then pcutp is to the left
    int somecounter = 0;
    if(flag_pcut==0)
    {
        for(int i=0;i<pcutmvec.size();++i)
        {
            real_part_pcut = real(pcutmvec[i]);
            imag_part_pcut = imag(pcutmvec[i]);

            if(imag_part_pcut<=0.0)
            {
                if(somecounter==0)
                {
                    actual_real_pcut = real(pcutmvec[i]);
                    actual_imag_pcut = imag(pcutmvec[i]);
                    somecounter=1;
                    continue;
                }

                if(real_part_pcut>actual_real_pcut)
                {
                    actual_real_pcut = real_part_pcut;
                    actual_imag_pcut = imag_part_pcut;
                }
                
            }
        }
    }
    else
    {
        for(int i=0;i<pcutpvec.size();++i)
        {
            real_part_pcut = real(pcutpvec[i]);
            imag_part_pcut = imag(pcutpvec[i]);

            if(imag_part_pcut<=0.0)
            {
                if(somecounter==0)
                {
                    actual_real_pcut = real(pcutpvec[i]);
                    actual_imag_pcut = imag(pcutpvec[i]);
                    somecounter=1;
                    continue;
                }

                if(real_part_pcut>actual_real_pcut)
                {
                    actual_real_pcut = real_part_pcut;
                    actual_imag_pcut = imag_part_pcut;
                }

            }
        }
    }

    value = actual_real_pcut + ii*actual_imag_pcut ;
}


void leftmost_left_pcut(   vector<comp> &pcutpvec,
                            vector<comp> &pcutmvec,
                            int flag_pcut,
                            comp &value    )
{
    comp ii = {0.0,1.0};
    double real_part_pcut = 0.0;
    double imag_part_pcut = 0.0;
    double temp_real_pcut = 0.0;
    double temp_imag_pcut = 0.0;
    double actual_real_pcut = 0.0;
    double actual_imag_pcut = 0.0;
    //if flag_pcut==0, then pcutm is to the left
    //if flag_pcut==1, then pcutp is to the left
    int somecounter = 0;
    if(flag_pcut==0)
    {
        for(int i=0;i<pcutmvec.size();++i)
        {
            real_part_pcut = real(pcutmvec[i]);
            imag_part_pcut = imag(pcutmvec[i]);


            if(somecounter==0)
            {
                actual_real_pcut = real_part_pcut;
                actual_imag_pcut = imag_part_pcut;
                somecounter=1;
                continue;
            }

            if(real_part_pcut<actual_real_pcut)
            {
                actual_real_pcut = real_part_pcut;
                actual_imag_pcut = imag_part_pcut;
            }
            
            
        }
    }
    else
    {
        for(int i=0;i<pcutpvec.size();++i)
        {
            real_part_pcut = real(pcutpvec[i]);
            imag_part_pcut = imag(pcutpvec[i]);


            if(somecounter==0)
            {
                actual_real_pcut = real_part_pcut;
                actual_imag_pcut = imag_part_pcut;
                somecounter=1;
                continue;
            }

            if(real_part_pcut<actual_real_pcut)
            {
                actual_real_pcut = real_part_pcut;
                actual_imag_pcut = imag_part_pcut;
            }
        }
    }

    value = actual_real_pcut + ii*actual_imag_pcut ;
}

void lowest_left_pcut(   vector<comp> &pcutpvec,
                            vector<comp> &pcutmvec,
                            int flag_pcut,
                            comp &value    )
{
    comp ii = {0.0,1.0};
    double real_part_pcut = 0.0;
    double imag_part_pcut = 0.0;
    double temp_real_pcut = 0.0;
    double temp_imag_pcut = 0.0;
    double actual_real_pcut = 0.0;
    double actual_imag_pcut = 0.0;
    //if flag_pcut==0, then pcutm is to the left
    //if flag_pcut==1, then pcutp is to the left
    int somecounter = 0;
    if(flag_pcut==0)
    {
        for(int i=0;i<pcutmvec.size();++i)
        {
            real_part_pcut = real(pcutmvec[i]);
            imag_part_pcut = imag(pcutmvec[i]);


            if(somecounter==0)
            {
                actual_real_pcut = real_part_pcut;
                actual_imag_pcut = imag_part_pcut;
                somecounter=1;
                continue;
            }

            if(imag_part_pcut<actual_imag_pcut)
            {
                actual_real_pcut = real_part_pcut;
                actual_imag_pcut = imag_part_pcut;
            }
            
            
        }
    }
    else
    {
        for(int i=0;i<pcutpvec.size();++i)
        {
            real_part_pcut = real(pcutpvec[i]);
            imag_part_pcut = imag(pcutpvec[i]);


            if(somecounter==0)
            {
                actual_real_pcut = real_part_pcut;
                actual_imag_pcut = imag_part_pcut;
                somecounter=1;
                continue;
            }

            if(imag_part_pcut<actual_imag_pcut)
            {
                actual_real_pcut = real_part_pcut;
                actual_imag_pcut = imag_part_pcut;
            }
        }
    }

    value = actual_real_pcut + ii*actual_imag_pcut ;
}


//the new vector maker to solve for mphib 3d for ims>0
void mom_vector_maker_seba_imspos_3(  vector<comp> &qvec,
                                    vector<comp> &weights, 
                                    comp s,
                                    comp kmin,
                                    comp kmax,
                                    double a,
                                    double m,
                                    double eps,
                                    double eps_for_m2k,
                                    double points,
                                    int &tag1,
                                    int &tag2,
                                    int &switch_for_gvec_fixer   )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    comp qc2 = q_c2(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    
    cout<<"axis_flag = "<<axis_flag<<endl;
    if(axis_flag==0)
    cout<<"cut to the right crosses real line"<<endl;
    else
    cout<<"cut to the right doesnot cross real line"<<endl;

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;

    if(axis_flag==0)
    {
        cout<<"deforming the contour backward first 1 "<<endl;
        comp q = pmom(s,sigmab(a,m),m);
        comp qc2 = q_c2(s,a,m,eps);
        comp pplus = pcut_plus(-1.0,s,q,m,eps);
        comp pminus = pcut_plus(+1.0,s,q,m,eps);
        cout<<"pplus = "<<pplus<<'\t'<<"pminus = "<<pminus<<endl;
        cout<<"pplus-pminus = "<<0.5*real(pplus - pminus)<<endl;
        comp qc1 = q_c1(s,a,m,eps);
        double rad = real(qc1);

        comp firstnode = kmin;
        comp secondnode;
        if(flag_pcut==0)
        {
            secondnode = pcutmvec[pcutmvec.size()/2];
        }
        else
        {
            secondnode = pcutpvec[pcutpvec.size()/2];
        }

        line_maker_with_weights(qvec,weights,firstnode,secondnode/2.0,points);

        comp thirdnode;
        if(flag_pcut==0)
        {
            thirdnode = pcutmvec[9*pcutmvec.size()/10];
        }
        else
        {
            thirdnode = pcutpvec[9*pcutpvec.size()/10];
        }

        line_maker_with_weights(qvec,weights,qvec[qvec.size()-1],thirdnode,points);

        tag1 = qvec.size() - 1;

        comp tempforthnode1;
        if(flag_pcut==0)
        {
            tempforthnode1 = pcutmvec[lowest_m_index];
        }
        else
        {
            tempforthnode1 = pcutpvec[lowest_p_index];
        }

        comp tempforthnode2;
        if(flag_pcut==0)
        {
            tempforthnode2 = pcutpvec[lowest_p_index];
        }
        else
        {
            tempforthnode2 = pcutmvec[lowest_m_index];
        }

        //comp forthnode = real(thirdnode) + ii*0.45*imag((tempforthnode1 + tempforthnode2)/2.0);
        comp forthnode = real(thirdnode) + ii*imag((tempforthnode1));// + tempforthnode2)/2.0);

        line_maker_with_weights(qvec,weights,thirdnode,forthnode,points);

        comp fifthnode;

        comp pth = phibthreshold(a,m);
        double multiplier = (1.0 - real(s)/real(pth));

        if(flag_pcut==0)
        {
            fifthnode = pcutpvec[multiplier*9.5*pcutpvec.size()/10];
        }
        else
        {
            fifthnode = pcutmvec[multiplier*9.5*pcutmvec.size()/10];
        }

        line_maker_with_weights(qvec,weights,forthnode,fifthnode,points);

        tag2 = qvec.size() - 1;

        comp tempsixthnode;

        if(flag_pcut==0)
        {
            tempsixthnode = pcutpvec[highest_p_index];
        }
        else
        {
            tempsixthnode = pcutmvec[highest_m_index];
        }

        comp sixthnode = real(tempsixthnode) + ii*imag(fifthnode);
        line_maker_with_weights(qvec,weights,fifthnode,sixthnode,points);

        comp seventhnode = kmax;

        line_maker_with_weights(qvec,weights,sixthnode,seventhnode,points);


         




    }
}

void mom_vector_maker_seba_imspos_4(  vector<comp> &qvec,
                                    vector<comp> &weights, 
                                    comp s,
                                    comp kmin,
                                    comp kmax,
                                    double a,
                                    double m,
                                    double eps,
                                    double eps_for_m2k,
                                    double points,
                                    int &tag1,
                                    int &tag2,
                                    int &switch_for_gvec_fixer   )
{
    //initial setup//
    comp ii = {0.0,1.0};
    comp q = pmom(s,sigmab(a,m),m);
    comp qc2 = q_c2(s,a,m,eps);
    comp pplus = pcut_plus(1.0,s,q,m,eps);
    comp pminus = pcut_plus(-1.0,s,q,m,eps);
    comp qc1 = q_c1(real(s),a,m,eps);
    comp sigb = sigmab(a,m);
    double rad = abs(real(qc1));
    cout<<"radius = "<<rad<<endl;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);

    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    if(axis_flag==0)
    {
        switch_for_gvec_fixer = 0;
        cout<<"contour chosen"<<endl;
        
        //firstnode - secondnode
        comp firstnode = kmin;

        comp radend;
        comp radstart;
        rightmost_left_pcut(pcutpvec,pcutmvec,flag_pcut,radend);
        leftmost_left_pcut(pcutpvec,pcutmvec,flag_pcut,radstart);

        double newrad = (real(radend) + real(radstart))/2.0;

        comp secondnode = newrad;//-rad/3.0;//-2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2);
        //1
        line_maker_with_weights(qvec,weights,firstnode,secondnode,1.0*points/4.0);

        //secondnode - thirdnode
        comp sigmac2 = sigc2(s,sigb,m);
        comp p = (q_c1(s,a,m,eps) + qc2)/2.0;//pmom(s,sigmac2,m);
        comp k = q;
        comp alphapk = ((sqrt(s)-omega_comp(p,m)-omega_comp(k,m))*(sqrt(s)-omega_comp(p,m)-omega_comp(k,m)) - p*p - k*k - m*m)/(2.0*p*k);
        //comp zpk = 2.0*s*(sigb + ii*eps) - (s+sigb - m*m)*(s+m*m-sigmac2);
        //comp zpk = 2.0*s*(sigmac2 + ii*eps) - (s+sigmac2 - m*m)*(s+m*m-sigb);
        
        double x0 = -1.0*real(alphapk);
        //double x0 = real(zpk);
        cout<<"sigma c2 = "<<sigmac2<<endl;
        
        cout<<"x0 = "<<x0<<endl;
        //comp thirdnode = -pcut_minus((-x0-0.77)/2.0,s,q,m,eps);
        
        
        //tag1 = qvec.size() - 1;
        //thirdnode - fourthnode 
        comp f1 = pcut_minus((x0+1.0)/2.0,s,q,m,eps);
        
        comp f2 = 0.5*(real(pplus - pminus)) - ii*abs(qc2);

        comp temp_f3_1 = pcut_minus(-1,s,q,m,eps);
        comp temp_f3_2 = pcut_minus(1,s,q,m,eps);
        comp temp_f3_3 = pcut_plus(-1,s,q,m,eps);
        comp temp_f3_4 = pcut_plus(1,s,q,m,eps);
        comp f3 = 0.5*(real(pplus - pminus)) - ii*5.0/4.0*abs(qc2);;
        /*if(imag(temp_f3_1)<imag(temp_f3_2) && imag(temp_f3_1)<imag(temp_f3_3) && imag(temp_f3_1)<imag(temp_f3_4))
        {
            f3 = pcut_minus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_2)<imag(temp_f3_1) && imag(temp_f3_2)<imag(temp_f3_3) && imag(temp_f3_2)<imag(temp_f3_4))
        {
            f3 = pcut_minus(1,s,q,m,eps);
        }
        else if(imag(temp_f3_3)<imag(temp_f3_1) && imag(temp_f3_3)<imag(temp_f3_2) && imag(temp_f3_3)<imag(temp_f3_4))
        {
            f3 = pcut_plus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_4)<imag(temp_f3_1) && imag(temp_f3_4)<imag(temp_f3_2) && imag(temp_f3_4)<imag(temp_f3_3))
        {
            f3 = pcut_plus(1,s,q,m,eps);
        }*/

        comp thirdnode = real(secondnode) + ii*0.85*imag((f2 + f3)/2.0);//-pcut_minus(-x0,s,q,m,eps);
        //2
        line_maker_with_weights(qvec,weights,secondnode,thirdnode,1.0*points/4.0);
        tag1 = qvec.size() - 1;
        
        //comp f3 = pcut_minus(1,s,q,m,eps);
        //comp fourthnode = real(thirdnode) + ii*imag((f2 + f3)/2.0);//real(thirdnode) + ii*imag(f1)*0.75;//pcut_minus(2.0*x0/3.0,s,q,m,eps);
        //3
        //line_maker_with_weights(qvec,weights,thirdnode,fourthnode,2.0*points/10.0);

        //tag2 = qvec.size() - 1;
        //fourthnode - fifthnode
        //comp fifthnode = (f2 + f3)/2.0;//real(f1) + ii*imag(f1)*0.75;//real(qc1) + ii*imag(fourthnode);//(2.0/3.0)*(1.0-2.0*ii)*abs(qc2);
        //4
        //line_maker_with_weights(qvec,weights,fourthnode,fifthnode,2.0*points/10.0);

        //fifthnode - sixthnode 
        comp sixthnode = real(q_c1(s,a,m,eps)) + (ii*imag((f2+f3)/2.0));//(3.0/2.0)*abs(qc2);
        //5
        //line_maker_with_weights(qvec,weights,fifthnode,sixthnode,1.0*points/4.0);
        line_maker_with_weights(qvec,weights,thirdnode,sixthnode,1.0*points/4.0);
        tag2 = qvec.size() - 1;
        //sixthnode - seventhnode
        comp seventhnode = kmax;
        //6
        line_maker_with_weights(qvec,weights,sixthnode,seventhnode,1.0*points/4.0);

    }
    else
    {
        switch_for_gvec_fixer = 1;
        cout<<"straight line chosen"<<endl;
        line_maker_with_weights(qvec,weights,kmin,kmax,points);     
    }

}

void mom_vector_maker_seba_imspos_5(  vector<comp> &qvec,
                                    vector<comp> &weights, 
                                    comp s,
                                    comp kmin,
                                    comp kmax,
                                    double a,
                                    double m,
                                    double eps,
                                    double eps_for_m2k,
                                    double points,
                                    int &tag1,
                                    int &tag2,
                                    int &switch_for_gvec_fixer   )
{
    //initial setup//
    comp ii = {0.0,1.0};
    comp q = pmom(s,sigmab(a,m),m);
    comp qc2 = q_c2(s,a,m,eps);
    comp pplus = pcut_plus(1.0,s,q,m,eps);
    comp pminus = pcut_plus(-1.0,s,q,m,eps);
    comp qc1 = q_c1(real(s),a,m,eps);
    comp sigb = sigmab(a,m);
    double rad = abs(real(qc1));
    //cout<<"radius = "<<rad<<endl;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);

    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    if(axis_flag==0)
    {
        switch_for_gvec_fixer = 0;
        //cout<<"contour chosen"<<endl;
        
        //firstnode - secondnode
        comp firstnode = kmin;

        comp radend;
        comp radstart;
        rightmost_left_pcut(pcutpvec,pcutmvec,flag_pcut,radend);
        leftmost_left_pcut(pcutpvec,pcutmvec,flag_pcut,radstart);

        double newrad = (real(radend) + real(radstart))/2.0;

        comp secondnode = newrad;//-rad/3.0;//-2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2);
        //1
        line_maker_with_weights(qvec,weights,firstnode,secondnode,1.0*points/4.0);

        //secondnode - thirdnode
        comp sigmac2 = sigc2(s,sigb,m);
        comp p = (q_c1(s,a,m,eps) + qc2)/2.0;//pmom(s,sigmac2,m);
        comp k = q;
        comp alphapk = ((sqrt(s)-omega_comp(p,m)-omega_comp(k,m))*(sqrt(s)-omega_comp(p,m)-omega_comp(k,m)) - p*p - k*k - m*m)/(2.0*p*k);
        //comp zpk = 2.0*s*(sigb + ii*eps) - (s+sigb - m*m)*(s+m*m-sigmac2);
        //comp zpk = 2.0*s*(sigmac2 + ii*eps) - (s+sigmac2 - m*m)*(s+m*m-sigb);
        
        double x0 = -1.0*real(alphapk);
        //double x0 = real(zpk);
        //cout<<"sigma c2 = "<<sigmac2<<endl;
        
        //cout<<"x0 = "<<x0<<endl;
        //comp thirdnode = -pcut_minus((-x0-0.77)/2.0,s,q,m,eps);
        
        
        //tag1 = qvec.size() - 1;
        //thirdnode - fourthnode 
        comp f1 = pcut_minus((x0+1.0)/2.0,s,q,m,eps);
        
        comp f2 = 0.5*(real(pplus - pminus)) - ii*abs(qc2);

        comp temp_f3_1 = pcut_minus(-1,s,q,m,eps);
        comp temp_f3_2 = pcut_minus(1,s,q,m,eps);
        comp temp_f3_3 = pcut_plus(-1,s,q,m,eps);
        comp temp_f3_4 = pcut_plus(1,s,q,m,eps);
        comp f3 = 0.5*(real(pplus - pminus)) - ii*5.0/4.0*abs(qc2);;
        /*if(imag(temp_f3_1)<imag(temp_f3_2) && imag(temp_f3_1)<imag(temp_f3_3) && imag(temp_f3_1)<imag(temp_f3_4))
        {
            f3 = pcut_minus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_2)<imag(temp_f3_1) && imag(temp_f3_2)<imag(temp_f3_3) && imag(temp_f3_2)<imag(temp_f3_4))
        {
            f3 = pcut_minus(1,s,q,m,eps);
        }
        else if(imag(temp_f3_3)<imag(temp_f3_1) && imag(temp_f3_3)<imag(temp_f3_2) && imag(temp_f3_3)<imag(temp_f3_4))
        {
            f3 = pcut_plus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_4)<imag(temp_f3_1) && imag(temp_f3_4)<imag(temp_f3_2) && imag(temp_f3_4)<imag(temp_f3_3))
        {
            f3 = pcut_plus(1,s,q,m,eps);
        }*/

        comp lowest_left;
        lowest_left_pcut(pcutpvec,pcutmvec,flag_pcut,lowest_left);

        //comp thirdnode = real(secondnode) + ii*0.85*imag((f2 + f3)/2.0);//-pcut_minus(-x0,s,q,m,eps);
        comp thirdnode = real(secondnode) + ii*imag(lowest_left) - ii*0.01;//-pcut_minus(-x0,s,q,m,eps);
        
        
        //2
        line_maker_with_weights(qvec,weights,secondnode,thirdnode,1.0*points/4.0);
        tag1 = qvec.size() - 1;
        
        //comp f3 = pcut_minus(1,s,q,m,eps);
        //comp fourthnode = real(thirdnode) + ii*imag((f2 + f3)/2.0);//real(thirdnode) + ii*imag(f1)*0.75;//pcut_minus(2.0*x0/3.0,s,q,m,eps);
        //3
        //line_maker_with_weights(qvec,weights,thirdnode,fourthnode,2.0*points/10.0);

        //tag2 = qvec.size() - 1;
        //fourthnode - fifthnode
        //comp fifthnode = (f2 + f3)/2.0;//real(f1) + ii*imag(f1)*0.75;//real(qc1) + ii*imag(fourthnode);//(2.0/3.0)*(1.0-2.0*ii)*abs(qc2);
        //4
        //line_maker_with_weights(qvec,weights,fourthnode,fifthnode,2.0*points/10.0);

        //fifthnode - sixthnode 
        //comp sixthnode = real(q_c1(s,a,m,eps)) + (ii*imag((f2+f3)/2.0));//(3.0/2.0)*abs(qc2);
        comp sixthnode = real(q_c1(s,a,m,eps)) + (ii*imag(thirdnode));//(3.0/2.0)*abs(qc2);
        
        
        //5
        //line_maker_with_weights(qvec,weights,fifthnode,sixthnode,1.0*points/4.0);
        line_maker_with_weights(qvec,weights,thirdnode,sixthnode,1.0*points/4.0);
        tag2 = qvec.size() - 1;
        //sixthnode - seventhnode
        comp seventhnode = kmax;
        //6
        line_maker_with_weights(qvec,weights,sixthnode,seventhnode,1.0*points/4.0);

    }
    else
    {
        switch_for_gvec_fixer = 1;
        //cout<<"straight line chosen"<<endl;
        line_maker_with_weights(qvec,weights,kmin,kmax,points);     
    }

}

void mom_vector_maker_seba_imspos_5_ext1(  vector<comp> &qvec,
                                    vector<comp> &weights, 
                                    comp s,
                                    comp kmin,
                                    comp kmax,
                                    double a,
                                    double m,
                                    double eps,
                                    double eps_for_m2k,
                                    double points,
                                    int &tag1,
                                    int &tag2,
                                    int &switch_for_gvec_fixer   )
{
    //initial setup//
    comp ii = {0.0,1.0};
    comp q = pmom(s,sigmab(a,m),m);
    comp qc2 = q_c2(s,a,m,eps);
    comp pplus = pcut_plus(1.0,s,q,m,eps);
    comp pminus = pcut_plus(-1.0,s,q,m,eps);
    comp qc1 = q_c1(real(s),a,m,eps);
    comp sigb = sigmab(a,m);
    double rad = abs(real(qc1));
    //cout<<"radius = "<<rad<<endl;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);

    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    if(axis_flag==0)
    {
        switch_for_gvec_fixer = 0;
        //cout<<"contour chosen"<<endl;
        
        //firstnode - secondnode
        comp firstnode = kmin;

        comp radend;
        comp radstart;
        rightmost_left_pcut(pcutpvec,pcutmvec,flag_pcut,radend);
        leftmost_left_pcut(pcutpvec,pcutmvec,flag_pcut,radstart);

        double newrad = (real(radend) + real(radstart))/2.0;

        comp secondnode = newrad;//-rad/3.0;//-2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2);
        //1
        line_maker_with_weights(qvec,weights,firstnode,secondnode,1.0*points/4.0);

        //secondnode - thirdnode
        comp sigmac2 = sigc2(s,sigb,m);
        comp p = (q_c1(s,a,m,eps) + qc2)/2.0;//pmom(s,sigmac2,m);
        comp k = q;
        comp alphapk = ((sqrt(s)-omega_comp(p,m)-omega_comp(k,m))*(sqrt(s)-omega_comp(p,m)-omega_comp(k,m)) - p*p - k*k - m*m)/(2.0*p*k);
        //comp zpk = 2.0*s*(sigb + ii*eps) - (s+sigb - m*m)*(s+m*m-sigmac2);
        //comp zpk = 2.0*s*(sigmac2 + ii*eps) - (s+sigmac2 - m*m)*(s+m*m-sigb);
        
        double x0 = -1.0*real(alphapk);
        //double x0 = real(zpk);
        //cout<<"sigma c2 = "<<sigmac2<<endl;
        
        //cout<<"x0 = "<<x0<<endl;
        //comp thirdnode = -pcut_minus((-x0-0.77)/2.0,s,q,m,eps);
        
        
        //tag1 = qvec.size() - 1;
        //thirdnode - fourthnode 
        comp f1 = pcut_minus((x0+1.0)/2.0,s,q,m,eps);
        
        comp f2 = 0.5*(real(pplus - pminus)) - ii*abs(qc2);

        comp temp_f3_1 = pcut_minus(-1,s,q,m,eps);
        comp temp_f3_2 = pcut_minus(1,s,q,m,eps);
        comp low_f3_1;
        if(imag(temp_f3_1)<imag(temp_f3_2))
        {
            low_f3_1 = temp_f3_1;
        }
        else 
        {
            low_f3_1 = temp_f3_2;
        }
        comp temp_f3_3 = pcut_plus(-1,s,q,m,eps);
        comp temp_f3_4 = pcut_plus(1,s,q,m,eps);
        comp low_f3_2;
        if(imag(temp_f3_3)<imag(temp_f3_4))
        {
            low_f3_2 = temp_f3_3;
        }
        else 
        {
            low_f3_2 = temp_f3_4;
        }
        comp low_f3;
        if(imag(low_f3_1)<imag(low_f3_2))
        {
            low_f3 = low_f3_1;
        }
        else 
        {
            low_f3 = low_f3_2;
        }

        //cout<<"branchcut endpoints 1 = "<<temp_f3_1<<endl;
        //cout<<"branchcut endpoints 2 = "<<temp_f3_2<<endl;
        //cout<<"branchcut endpoints 3 = "<<temp_f3_3<<endl;
        //cout<<"branchcut endpoints 4 = "<<temp_f3_4<<endl;
        //cout<<"chosen = "<<low_f3<<endl;

        comp f3 = 0.5*(real(pplus - pminus)) - ii*5.0/4.0*abs(qc2);;
        /*if(imag(temp_f3_1)<imag(temp_f3_2) && imag(temp_f3_1)<imag(temp_f3_3) && imag(temp_f3_1)<imag(temp_f3_4))
        {
            f3 = pcut_minus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_2)<imag(temp_f3_1) && imag(temp_f3_2)<imag(temp_f3_3) && imag(temp_f3_2)<imag(temp_f3_4))
        {
            f3 = pcut_minus(1,s,q,m,eps);
        }
        else if(imag(temp_f3_3)<imag(temp_f3_1) && imag(temp_f3_3)<imag(temp_f3_2) && imag(temp_f3_3)<imag(temp_f3_4))
        {
            f3 = pcut_plus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_4)<imag(temp_f3_1) && imag(temp_f3_4)<imag(temp_f3_2) && imag(temp_f3_4)<imag(temp_f3_3))
        {
            f3 = pcut_plus(1,s,q,m,eps);
        }*/

        comp lowest_left;
        lowest_left_pcut(pcutpvec,pcutmvec,flag_pcut,lowest_left);

        comp low_qc1;
        if(imag(qc1)<imag(-qc1))
        {
            low_qc1 = qc1;
        }
        else 
        {
            low_qc1 = -qc1;
        }
        
        comp low_qc2;
        if(imag(qc2)<imag(-qc2))
        {
            low_qc2 = qc2;
        }
        else 
        {
            low_qc2 = -qc2;
        }

        //cout<<"qc 1 = "<<qc1<<endl;
        //cout<<"low qc1 = "<<low_qc1<<endl;

        //cout<<"qc 2 = "<<qc2<<endl;
        //cout<<"low qc2 = "<<low_qc2<<endl;
        //comp thirdnode = real(secondnode) + ii*0.85*imag((f2 + f3)/2.0);//-pcut_minus(-x0,s,q,m,eps);
        comp thirdnode = real(secondnode) + ii*imag((low_f3 + low_qc2)/2.0);//real(secondnode) + ii*imag(lowest_left) - ii*0.01;//-pcut_minus(-x0,s,q,m,eps);
        
        
        //2
        line_maker_with_weights(qvec,weights,secondnode,thirdnode,1.0*points/4.0);
        tag1 = qvec.size() - 1;
        
        //comp f3 = pcut_minus(1,s,q,m,eps);
        //comp fourthnode = real(thirdnode) + ii*imag((f2 + f3)/2.0);//real(thirdnode) + ii*imag(f1)*0.75;//pcut_minus(2.0*x0/3.0,s,q,m,eps);
        //3
        //line_maker_with_weights(qvec,weights,thirdnode,fourthnode,2.0*points/10.0);

        //tag2 = qvec.size() - 1;
        //fourthnode - fifthnode
        //comp fifthnode = (f2 + f3)/2.0;//real(f1) + ii*imag(f1)*0.75;//real(qc1) + ii*imag(fourthnode);//(2.0/3.0)*(1.0-2.0*ii)*abs(qc2);
        //4
        //line_maker_with_weights(qvec,weights,fourthnode,fifthnode,2.0*points/10.0);

        //fifthnode - sixthnode 
        //comp sixthnode = real(q_c1(s,a,m,eps)) + (ii*imag((f2+f3)/2.0));//(3.0/2.0)*abs(qc2);
        comp sixthnode = real(q_c1(s,a,m,eps)) + (ii*imag(thirdnode));//(3.0/2.0)*abs(qc2);
        
        
        //5
        //line_maker_with_weights(qvec,weights,fifthnode,sixthnode,1.0*points/4.0);
        line_maker_with_weights(qvec,weights,thirdnode,sixthnode,1.0*points/4.0);
        tag2 = qvec.size() - 1;
        //sixthnode - seventhnode
        comp seventhnode = kmax;
        //6
        line_maker_with_weights(qvec,weights,sixthnode,seventhnode,1.0*points/4.0);

    }
    else
    {
        switch_for_gvec_fixer = 1;
        //cout<<"straight line chosen"<<endl;
        line_maker_with_weights(qvec,weights,kmin,kmax,points);     
    }

}


void contour_for_resonance( vector<comp> &qvec,
                            vector<comp> &weights,
                            comp kmin, 
                            comp kmax,
                            double shift,
                            double qvecpoints   )
{
    
    comp ii = {0.0,1.0};
    comp midpoint = (kmax + kmin)/2.0;
    midpoint = midpoint + ii*shift;
    
    
    line_maker_with_weights(qvec,weights,kmin,midpoint,qvecpoints/2.0);
    line_maker_with_weights(qvec,weights,midpoint,kmax,qvecpoints/2.0);
}

void contour_for_resonance_1( vector<comp> &qvec,
                            vector<comp> &weights,
                            comp q,
                            comp kmin, 
                            comp kmax,
                            double shift,
                            double qvecpoints   )
{
    
    comp ii = {0.0,1.0};
    comp midpoint = (kmax + kmin)/2.0;
    midpoint = midpoint + ii*shift;
    
    comp firstpoint = real(q) - abs(shift);
    //cout<<kmin<<'\t'<<firstpoint<<endl;
    comp secondpoint = real(firstpoint) + ii*imag(q) + ii*shift;
    comp thirdpoint = real(secondpoint) + 2.0*abs(shift) + ii*imag(secondpoint);
    comp forthpoint = real(thirdpoint);

    if(imag(q)<=0.0)
    {
        line_maker_with_weights(qvec,weights,kmin,firstpoint,qvecpoints/5.0);
        line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints/5.0);
        line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints/5.0);
        line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints/5.0);
        line_maker_with_weights(qvec,weights,forthpoint,kmax,qvecpoints/5.0);
    }
    else
    {
        line_maker_with_weights(qvec,weights,kmin,kmax,qvecpoints);
    }
}

void contour_for_resonance_2( vector<comp> &qvec,
                            vector<comp> &weights,
                            comp kmin, 
                            double reshift,
                            double imshift,
                            comp kmax,
                            double qvecpoints   )
{
    
    comp ii = {0.0,1.0};

    comp shift = reshift + ii*imshift;
    
    
    line_maker_with_weights(qvec,weights,kmin,shift,qvecpoints/2.0);
    line_maker_with_weights(qvec,weights,shift,kmax,qvecpoints/2.0);
}

comp M2kfunc_1(   double a,
                comp sigk,
                double m,
                double epsilon    )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    //cout<<"sigk:"<<sigk<<endl;
    comp num = 16.0*pi*sqrt(sigk);
    comp denom = -1.0/a - ii*sqrt((sigk+ii*epsilon)/4.0 - m*m);

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

    comp rho = mysqrt(sigk - 4.0*m*m)/(32.0*pi*sqrt(sigk));

    return m2k/(1.0 + 2.0*ii*rho*m2k);
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

void contour_for_resonance_3( vector<comp> &qvec,
                            vector<comp> &weights,
                            comp s,
                            double m,
                            comp kmin, 
                            comp kmax,
                            double eps,
                            double shift,
                            double qvecpoints   )
{
    
    comp ii = {0.0,1.0};


    comp m2k_cut_1 = M2kbranchcut_right_momrep_plus_eps(s,m,eps);
    comp m2k_cut_2 = M2kbranchcut_right_momrep_minus_eps(s,m,eps);
    comp selected_m2k_cut;
    if(real(m2k_cut_1)>real(m2k_cut_2))
    {
        selected_m2k_cut = m2k_cut_1;
    }
    else
    {
        selected_m2k_cut = m2k_cut_2;
    }

    
    
    if(imag(selected_m2k_cut)<0.0)
    {
        comp firstpoint = kmin;
        comp secondpoint = ii*imag(selected_m2k_cut) + ii*shift;

        line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints/4.0);

        comp thirdpoint = ii*imag(secondpoint) + real(selected_m2k_cut) + real(selected_m2k_cut)/2.0;

        line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints/2.0);

        comp forthpoint = kmax;

        line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints/4.0);
    }
    else 
    {
        comp firstpoint = kmin;
        comp secondpoint = ii*shift; 
        
        line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints/4.0);
        
        comp thirdpoint = ii*imag(secondpoint) + real(selected_m2k_cut) + real(selected_m2k_cut)/2.0;

        line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints/2.0);
        
        comp forthpoint = kmax;
        
        line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints/4.0);

    }

}

//this contour is for re(s)<9 and im(s)<0 
void contour_for_resonance_4( vector<comp> &qvec,
                            vector<comp> &weights,
                            comp s,
                            double m,
                            comp kmin, 
                            comp kmax,
                            double eps,
                            double qvecpoints,
                            int &tag_for_m2k   )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);

    vector<comp> m2kcut;
    for(int i=0;i<500;++i)
    {
        double delsigk = abs(1.0)/500;
        double delsig = i*delsigk;
        comp sigk = 4.0 + (comp)delsig;

        comp m2kcut_in_p = pmom(s,sigk,m);
        m2kcut.push_back(m2kcut_in_p);
    }

    comp firstpoint = kmin;
    comp secondpoint = -0.15 -0.35*ii;
    line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints/4.0);

    comp thirdpoint = m2kcut[2*m2kcut.size()/3];
    line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints/4.0);

    tag_for_m2k = qvec.size() - 1;

    comp forthpoint = kmax; 
    line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints/2.0);


}

//this contour is for sheet -1 of the three body amplitude
//this is to look for mirror poles 
void contour_for_resonance_for_sheet_minus1_imspos( vector<comp> &qvec,
                            vector<comp> &weights,
                            comp s,
                            double m,
                            comp kmin, 
                            comp kmax,
                            double eps,
                            double qvecpoints,
                            int &tag_for_m2k   )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);

    vector<comp> m2kcut;
    for(int i=0;i<500;++i)
    {
        double delsigk = abs(1.0)/500;
        double delsig = i*delsigk;
        comp sigk = 4.0 + (comp)delsig;

        comp m2kcut_in_p = pmom(real(s) - ii*imag(s),sigk,m);
        m2kcut.push_back(m2kcut_in_p);
    }

    comp firstpoint = kmin;
    comp secondpoint = -0.15 -0.35*ii;
    line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints/4.0);

    comp thirdpoint = m2kcut[2*m2kcut.size()/3];
    line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints/4.0);

    tag_for_m2k = qvec.size() - 1;

    comp forthpoint = kmax; 
    line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints/2.0);


}


void contour_for_resonance_5( vector<comp> &qvec,
                            vector<comp> &weights,
                            comp s,
                            double m,
                            comp kmin, 
                            comp kmax,
                            double eps,
                            double qvecpoints,
                            int &tag_for_m2k   )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);

    vector<comp> m2kcut;
    for(int i=0;i<500;++i)
    {
        double delsigk = abs(1.0)/500;
        double delsig = i*delsigk;
        comp sigk = 4.0 + (comp)delsig;

        comp m2kcut_in_p = pmom(s,sigk,m);
        if(imag(m2kcut_in_p)>0.0)
        {
            comp temp = -real(m2kcut_in_p) - ii*imag(m2kcut_in_p);
            m2kcut_in_p = temp;
        }
        m2kcut.push_back(m2kcut_in_p);
    }

    comp firstpoint = kmin;
    comp secondpoint = -0.05 -0.05*ii;
    line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints/4.0);

    comp thirdpoint = secondpoint;//real(secondpoint) + ii*imag(m2kcut[m2kcut.size()/2]);
    //line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints/1.0);
    
    comp forthpoint = m2kcut[m2kcut.size()/2];
    line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints/4.0);
    
    tag_for_m2k = qvec.size() - 1;

    comp fifthpoint = ii*imag(qvec[qvec.size()-1]) + real(kmax);
    line_maker_with_weights(qvec,weights,forthpoint,fifthpoint,qvecpoints/4.0);

    comp sixthpoint = kmax; 
    line_maker_with_weights(qvec,weights,fifthpoint,sixthpoint,qvecpoints/4.0);


}

//this contour is for second sheet of dS the top part is s plane
//this will work with simag>0.0 part only
void contour_for_resonance_6( vector<comp> &qvec,
                            vector<comp> &weights,
                            double a,
                            comp s,
                            double m,
                            comp kmin,
                            comp kmax,
                            double eps,
                            double qvecpoints,
                            int &tag_for_m2k   )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);

    vector<comp> m2kcut;
    for(int i=0;i<500;++i)
    {
        double delsigk = abs(1.0)/500;
        double delsig = i*delsigk;
        comp sigk = 4.0 + (comp)delsig;

        comp m2kcut_in_p = -pmom(s,sigk,m);
        m2kcut.push_back(m2kcut_in_p);
    }

    comp sigb = sigmab(a,m);
    comp qpole = -pmom(s,sigb,m);
    comp m2kbc = -M2kbranchcut_right_momrep_plus_eps(s,m,eps);
    
    double shiftreal = 0.05;
    //if(real(s)>9.0)
    int set = 0;

    if(set==0)
    {
        comp firstpoint = kmin;
        comp secondpoint = real(m2kbc) + 0.005*ii;//-0.005 + 0.005*ii;
        //line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints*2.0/10.0);

        //comp thirdpoint = real(qpole) - shiftreal + ii*imag(qpole);
        
        comp thirdpoint = (qpole + m2kbc)/2.0;
        
        
        //line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints*2.0/10.0);
        line_maker_with_weights(qvec,weights,firstpoint,thirdpoint,qvecpoints*2.0/10.0);

        //comp forthpoint = real(thirdpoint) + ii*imag(m2kcut[2*m2kcut.size()/3])/2.0;
        comp forthpoint = real(thirdpoint) + ii*imag(m2kcut[2*m2kcut.size()/3]);
        line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints*2.0/10.0);

        comp fifthpoint = m2kcut[2*m2kcut.size()/3];
        line_maker_with_weights(qvec,weights,forthpoint,fifthpoint,qvecpoints*2.0/10.0);

        tag_for_m2k = qvec.size() - 1;

        comp sixthpoint = real(kmax);
        line_maker_with_weights(qvec,weights,fifthpoint,sixthpoint,qvecpoints*2.0/10.0);

        comp seventhpoint = kmax;
        line_maker_with_weights(qvec,weights,sixthpoint,seventhpoint,qvecpoints*2.0/10.0);


    }
    else 
    {
        
        comp firstpoint = kmin;
        comp secondpoint = real(m2kbc);// + 0.005*ii;//-0.005 + 0.005*ii;
        line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints*2.0/10.0);

        //comp thirdpoint = real(qpole) - shiftreal + ii*imag(qpole);
        
        comp thirdpoint = (qpole + m2kbc)/2.0;
        
        
        line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints*1.0/10.0);
        //line_maker_with_weights(qvec,weights,firstpoint,thirdpoint,qvecpoints*2.0/10.0);

        //comp forthpoint = real(thirdpoint) + ii*imag(m2kcut[2*m2kcut.size()/3])/2.0;
        comp forthpoint = real(thirdpoint) + ii*imag(m2kcut[2*m2kcut.size()/3]);
        line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints*1.0/10.0);

        comp fifthpoint = m2kcut[2*m2kcut.size()/3];
        line_maker_with_weights(qvec,weights,forthpoint,fifthpoint,qvecpoints*2.0/10.0);

        tag_for_m2k = qvec.size() - 1;

        comp sixthpoint = real(kmax);
        line_maker_with_weights(qvec,weights,fifthpoint,sixthpoint,qvecpoints*2.0/10.0);

        comp seventhpoint = kmax;
        line_maker_with_weights(qvec,weights,sixthpoint,seventhpoint,qvecpoints*2.0/10.0);


    }
    
    
    

}

//this contour is for sheet -1 to look for mirror poles 
void contour_for_resonance_for_sheet_minus1_imsneg( vector<comp> &qvec,
                            vector<comp> &weights,
                            double a,
                            comp s,
                            double m,
                            comp kmin,
                            comp kmax,
                            double eps,
                            double qvecpoints,
                            int &tag_for_m2k   )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);

    vector<comp> m2kcut;
    for(int i=0;i<500;++i)
    {
        double delsigk = abs(1.0)/500;
        double delsig = i*delsigk;
        comp sigk = 4.0 + (comp)delsig;

        comp m2kcut_in_p = -pmom(real(s)-ii*imag(s),sigk,m);
        m2kcut.push_back(m2kcut_in_p);
    }

    comp sigb = sigmab(a,m);
    comp qpole = -pmom(s,sigb,m);
    comp m2kbc = -M2kbranchcut_right_momrep_plus_eps(s,m,eps);
    
    double shiftreal = 0.05;
    //if(real(s)>9.0)
    int set = 0;

    if(set==0)
    {
        comp firstpoint = kmin;
        comp secondpoint = real(m2kbc) + 0.005*ii;//-0.005 + 0.005*ii;
        //line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints*2.0/10.0);

        //comp thirdpoint = real(qpole) - shiftreal + ii*imag(qpole);
        
        comp thirdpoint = (qpole + m2kbc)/2.0;
        
        
        //line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints*2.0/10.0);
        line_maker_with_weights(qvec,weights,firstpoint,thirdpoint,qvecpoints*2.0/10.0);

        //comp forthpoint = real(thirdpoint) + ii*imag(m2kcut[2*m2kcut.size()/3])/2.0;
        comp forthpoint = real(thirdpoint) + ii*imag(m2kcut[2*m2kcut.size()/3]);
        line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints*2.0/10.0);

        comp fifthpoint = m2kcut[2*m2kcut.size()/3];
        line_maker_with_weights(qvec,weights,forthpoint,fifthpoint,qvecpoints*2.0/10.0);

        tag_for_m2k = qvec.size() - 1;

        comp sixthpoint = real(kmax);
        line_maker_with_weights(qvec,weights,fifthpoint,sixthpoint,qvecpoints*2.0/10.0);

        comp seventhpoint = kmax;
        line_maker_with_weights(qvec,weights,sixthpoint,seventhpoint,qvecpoints*2.0/10.0);


    }
    else 
    {
        
        comp firstpoint = kmin;
        comp secondpoint = real(m2kbc);// + 0.005*ii;//-0.005 + 0.005*ii;
        line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints*2.0/10.0);

        //comp thirdpoint = real(qpole) - shiftreal + ii*imag(qpole);
        
        comp thirdpoint = (qpole + m2kbc)/2.0;
        
        
        line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints*1.0/10.0);
        //line_maker_with_weights(qvec,weights,firstpoint,thirdpoint,qvecpoints*2.0/10.0);

        //comp forthpoint = real(thirdpoint) + ii*imag(m2kcut[2*m2kcut.size()/3])/2.0;
        comp forthpoint = real(thirdpoint) + ii*imag(m2kcut[2*m2kcut.size()/3]);
        line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints*1.0/10.0);

        comp fifthpoint = m2kcut[2*m2kcut.size()/3];
        line_maker_with_weights(qvec,weights,forthpoint,fifthpoint,qvecpoints*2.0/10.0);

        tag_for_m2k = qvec.size() - 1;

        comp sixthpoint = real(kmax);
        line_maker_with_weights(qvec,weights,fifthpoint,sixthpoint,qvecpoints*2.0/10.0);

        comp seventhpoint = kmax;
        line_maker_with_weights(qvec,weights,sixthpoint,seventhpoint,qvecpoints*2.0/10.0);


    }
    
    
    

}


//this contour is made for the third sheet of dS
//the bottom simag panel 
void contour_for_resonance_7( vector<comp> &qvec,
                            vector<comp> &weights,
                            double a,
                            comp s,
                            double m,
                            comp kmin,
                            comp kmax,
                            double eps,
                            double qvecpoints,
                            int &tag_for_m2k_1,
                            int &tag_for_m2k_2   )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    double shiftreal = 0.05;
    double shiftcontour_from_m2kbc = 0.015;

    vector<comp> m2kcut;
    for(int i=0;i<500;++i)
    {
        double delsigk = abs(1.0)/500;
        double delsig = i*delsigk;
        comp sigk = 4.0 + (comp)delsig;

        comp m2kcut_in_p;
        if(imag(s)<0.0)
        m2kcut_in_p = pmom(s,sigk,m);
        else 
        m2kcut_in_p = -pmom(s,sigk,m);

        m2kcut.push_back(m2kcut_in_p);
    }

    comp sigb = sigmab(a,m);
    comp qpole;
    comp m2kbc;

    if(imag(s)<0.0)
    {
        qpole = pmom(s,sigb,m);
        m2kbc = M2kbranchcut_right_momrep_plus_eps(s,m,eps);
    }
    else
    {
        qpole = -pmom(s,sigb,m);
        m2kbc = -M2kbranchcut_right_momrep_plus_eps(s,m,eps);
    }
    //if(real(s)>9.0)
    {
        comp firstpoint = kmin;
        comp m2k_point = m2kcut[5];//m2kcut[1*m2kcut.size()/50];
        comp secondpoint = -real(m2kbc) + ii*imag(kmin + m2k_point)/2.0;//-0.05 - 0.05*ii;
        
        line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints*1.0/10.0);

        
        comp thirdpoint = m2k_point;//m2kcut[1*m2kcut.size()/20];//(qpole + m2kbc)/2.0;
        
        line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints*2.0/10.0);
        tag_for_m2k_1 = qvec.size() - 1;

        double shifts = 4.0;

        comp forthpoint = 10.0*real((qpole + m2kbc)/shifts) + ii*imag(thirdpoint);
        
        line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints*1.0/10.0);

        comp fifthpoint = real(forthpoint) - ii*imag((qpole + m2kbc)/shifts);//(qpole + m2kbc)/shifts;
        
        line_maker_with_weights(qvec,weights,forthpoint,fifthpoint,qvecpoints*1.0/10.0);

        comp sixthpoint = real(m2kbc) + ii*shiftcontour_from_m2kbc;
        //line_maker_with_weights(qvec,weights,fifthpoint,sixthpoint,qvecpoints*1.0/10.0);

        comp seventhpoint = (-fifthpoint + m2kbc)/2.0;//(kmin + fifthpoint)/2.0;
        
        //line_maker_with_weights(qvec,weights,sixthpoint,seventhpoint,qvecpoints*1.0/10.0);

        comp eighthpoint = thirdpoint;
        //line_maker_with_weights(qvec,weights,seventhpoint,eighthpoint,qvecpoints*1.0/10.0);
        
        tag_for_m2k_2 = qvec.size() - 1;

        comp ninthpoint = real(kmax) + ii*imag(seventhpoint);
        //line_maker_with_weights(qvec,weights,eighthpoint,ninthpoint,qvecpoints*1.0/10.0);
        
        comp tenthpoint = kmax;
        //line_maker_with_weights(qvec,weights,ninthpoint,tenthpoint,qvecpoints*1.0/10.0);
        

    }

}

/*  This is the same massive data generator but it puts all the values
    of individual qvec points in to one file, this helps us to see if 
    for a given energy, the qvec contour we formed, for each point of 
    the qvec, the pcut does not crosses the BS pole in momentum repre-
    sentation. 

*/

void kernel_pcut_x_glockletest_vs_s3_allqvec()
{
    //cout<<"here"<<endl;
    ofstream fout;
    double a = -10.1;//was 16
    double m = 1.0;
    double sreal = -3.0;//3.0;
    double simag = -0.01;//-0.05;
    comp ii = {0.0,1.0};
    double eps = 0.0;//1.0e-3;//1.0e-4;
    double eps1 = eps;//0.0;

    //comp sigb = sigmab(a,m);
    comp sigb = 2.0*m*m;//sigmab(a,m);
    
    

    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int delspoints = 30;
    int qvecpoints = 500;
    double qvec_r = 0.01;

    double delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double rad = 0.16;
    double pRinitial = -1.221;//-abs(abs(rad)+0.015);
    double pRfinal = 1.22;//abs(abs(rad)+0.015);
    double delpR = abs(pRinitial-pRfinal)/points;
    double pIinitial = -1.221;//-abs(abs(rad)+0.015);
    double pIfinal = 1.22;//abs(abs(rad)+0.015);
    double delpI = abs(pIinitial - pIfinal)/points;

    double sinitial = 8.98;//8.7;//(double)real(phibthreshold(a,m));//8.72;
    double sfinal = 9.01;//(double)real(phibthreshold(a,m));//9.0;//(double)real(phibthreshold(a,m));//8.0;
    double dels = abs(sinitial-sfinal)/(double)delspoints;

    int scount = 0;

    for(int i=0;i<delspoints+1;++i)
    //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    //int i=0;
    {
        
        //if(i!=0) break;
        cout<<"am = "<<a<<endl;
        double s3 = sinitial + i*dels;
        //double s3 = 8.7;//sfinal - i*dels;
        //if(scount!=0) break;
        int pcount = 0;
        
        //s3 = 8.99;
        double delsimag = abs(-0.01 - 0.01)/delspoints;
        //simag = -0.01 + i*delsimag;
        simag = -0.001;
        comp s = s3 + ii*simag;

        //This section only focuses on the third node for the contour
        comp q1 = pmom(s,sigmab(a,m),m);
        //comp q1 = pmom(s,2.0,m);
        
        comp qc2 = q_c2(s,a,m,eps);
        comp pplus = pcut_plus(1.0,s,q1,m,eps);
        comp pminus = pcut_plus(-1.0,s,q1,m,eps);
        comp thirdnode = 0.5*real(pplus - pminus) - abs(qc2)*ii;
        comp qc1 = q_c1(s,a,m,eps);
        double rad = abs(qc1);
        comp fifthnode = (2.0/3.0)*(1.0-2.0*ii)*abs(qc2);
        cout<<"radius = "<<rad<<endl;
        
        /*if(rad>0.0)
        {
            pRinitial = -abs(abs(rad) + abs(rad)/3.0);
            pRfinal = abs(abs(rad) + abs(rad + 0.01057943)/3.0);
            delpR = abs(pRinitial-pRfinal)/points;
            pIinitial = -abs(abs(rad) + abs(rad)/3.0);
            pIfinal = abs(abs(rad) + abs(rad + 0.0103793)/3.0);
            delpI = abs(pIinitial - pIfinal)/points;
        }
        else
        {
            pRinitial = -1.017;
            pRfinal = 1.0;
            delpR = abs(pRinitial-pRfinal)/points;
            pIinitial = -1.019;
            pIfinal = 1.0;
            delpI = abs(pIinitial - pIfinal)/points;
        }*/
        
        comp q = pmom(s,sigb-ii*eps1,m);
        comp kmax = pmom(s,0.0,m);
        cout<<"run = "<<scount+1<<endl;
    
        cout<<"s = "<<s<<'\t'<<"dels = "<<dels<<endl;
        cout<<"eps = "<<eps<<'\t'<<"eps_for_m2k = "<<eps1<<endl;

        vector<comp> qvec;
        vector<comp> weights;
        ifstream fin;
        string input_contour = "contur_res.dat";
        fin.open(input_contour.c_str());
        //int tag1,tag2; 
        if(simag<0.0)
        {
            //mom_vector_maker_4(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            comp qc2 = q_c2(s,a,m,eps);
            cout<<"qc2 = "<<qc2<<endl;
            //mom_vector_maker_43_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,qvec_r);
            //mom_vector_maker_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints);
            qvecpoints = 500;
            //mom_vector_maker_43_with_weights_with_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,qvec_r);
            double reshift = 0.1;
            double imshift = -0.1;
            //contour_for_resonance_2(qvec,weights,0.0,reshift,imshift,kmax,qvecpoints);
            //contour_for_resonance_3(qvec,weights,s,m,0.0,kmax,eps,-0.1,qvecpoints);
            int tag_for_m2k;
            int tag_for_m2k_1;
            int tag_for_m2k_2;
            //contour_for_resonance_7(qvec,weights,a,s,m,0.0,kmax,eps,qvecpoints,tag_for_m2k_1,tag_for_m2k_2);
            
            //contour_for_resonance_5(qvec,weights,s,m,0.0,kmax,eps,qvecpoints,tag_for_m2k);
            
            /*double some_a, some_b;
            while(fin>>some_a>>some_b)
            {
                comp temp_contour = some_a + ii*some_b;
                cout<<temp_contour<<endl;
                qvec.push_back(temp_contour);
            }*/

            /*for(int qind=0;qind<qvec.size();++qind)
            {
                comp reqvec = -real(qvec[qind]);
                comp imqvec = -imag(qvec[qind]);
                qvec[qind] = reqvec + ii*imqvec;
            }*/


            //contour_for_resonance_6(qvec,weights,a,s,m,0.0,kmax,eps,qvecpoints,tag_for_m2k);
            contour_for_resonance_for_sheet_minus1_imsneg(qvec,weights,a,s,m,0.0,kmax,eps,qvecpoints,tag_for_m2k);
            //this is the last one i was using
            //contour_for_resonance_4(qvec,weights,s,m,0.0,kmax,eps,qvecpoints,tag_for_m2k);
            
            //line_maker_with_weights(qvec,weights,0.0,kmax,qvecpoints);
                
            //mom_vector_maker_41(qvec,#)
            for(int qind=0;qind<qvec.size();++qind)
            {
                comp reqvec = real(qvec[qind]);
                comp imqvec = -imag(qvec[qind]);
                qvec[qind] = reqvec + ii*imqvec;
            }
        }
        else if(simag>=0.0)
        {
            //mom_vector_maker_4(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            //mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            //mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            int gvec_fixer_tag = 0;
            int tag1 = 0;
            int tag2 = 0;
            //mom_vector_maker_seba_imspos_2(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,tag1,tag2,gvec_fixer_tag);
            //mom_vector_maker_seba_imspos_4(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,tag1,tag2,gvec_fixer_tag);
            //mom_vector_maker_seba_imspos_5(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,tag1,tag2,gvec_fixer_tag);
            //mom_vector_maker_seba_imspos_5_ext1(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,tag1,tag2,gvec_fixer_tag);
            
            //line_maker_with_weights(qvec,weights,0.0,kmax,qvecpoints);
            //contour_for_resonance(qvec,weights,0.0,kmax,0.1,qvecpoints);
            cout<<M2kbranchcut_left_momrep_plus(s,m)<<endl;
            cout<<M2kbranchcut_left_momrep_minus(s,m)<<endl;
            cout<<M2kbranchcut_right_momrep_plus_eps(s,m,eps1)<<endl;
            cout<<M2kbranchcut_right_momrep_minus_eps(s,m,eps1)<<endl;

            double reshift = 0.1;
            double imshift = -0.1;
            //contour_for_resonance_2(qvec,weights,0.0,reshift,imshift,kmax,qvecpoints);
            int tag_for_m2k;
            int tag_for_m2k_1;
            int tag_for_m2k_2;
            
            //contour_for_resonance_4(qvec,weights,s,m,0.0,kmax,eps,qvecpoints,tag_for_m2k);
            contour_for_resonance_for_sheet_minus1_imspos(qvec,weights,s,m,0.0,kmax,eps,qvecpoints,tag_for_m2k);
            //this is the last one we were using
            //contour_for_resonance_6(qvec,weights,a,s,m,0.0,kmax,eps,qvecpoints,tag_for_m2k);
            
            //in test
            //contour_for_resonance_7(qvec,weights,a,s,m,0.0,kmax,eps,qvecpoints,tag_for_m2k_1,tag_for_m2k_2);
            
            //contour_for_resonance_1(qvec,weights,M2kbranchcut_right_momrep_plus_eps(s,m,eps1),0.0,kmax,-abs(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))/10,qvecpoints);
            //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints);
            //mom_vector_maker_seba_imspos_2_with_contour47(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,tag1,tag2,gvec_fixer_tag);
            //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints);
            
            for(int qind=0;qind<qvec.size();++qind)
            {
                comp reqvec = real(qvec[qind]);
                comp imqvec = -imag(qvec[qind]);
                qvec[qind] = reqvec + ii*imqvec;
            }
        }
        fin.close();
        //mom_vector_maker_61(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r,tag1,tag2);
        
        cout<<"qvec size = "<<qvec.size()<<endl;
        //cout<<"tag1 = "<<tag1<<'\t'<<"tag2 = "<<tag2<<endl;
        //for(int i=0;i<qvec.size();++i) cout<<i<<'\t'<<qvec[i]<<endl;

        //This part writes the qvec in to a file //
        string qvecfile =   //"qvec_a_" + to_string(a)
                          //+ "_eps_" + to_string(eps)
                           "qvec_scount_" + to_string(scount)
                          + ".dat";


        fout.open(qvecfile.c_str());
        for(int i=0;i<qvec.size();++i)
        {
            fout<<i<<'\t'<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
        }
        fout.close();
        cout<<"qvecfile : "<<qvecfile<<endl;

        //This part writes the kernel with k=q
        
        int isize = (int) points;
        int jsize = (int) points;

        string kernelfile =  // "kernel_" //a_" + to_string(a)
                            //+ "_eps_" + to_string(eps)
                             "kernel_scount_" + to_string(scount)
                            + ".dat";

        fout.open(kernelfile.c_str());
        for(int i=0;i<isize;++i)
        {
            for(int j=0;j<jsize;++j)
            {
                comp pR = (comp)pRinitial + (comp)i*delpR;
                comp pI = (comp)pIinitial + (comp)j*delpI; 
                comp pcomp = (comp) pR + ii*pI;

                comp sigp = sigma_p(s,pcomp,m);

                double pi = acos(-1.0);

                //comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
                //comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*M2kfunc(a,sigp,m,eps1);
                double theta = 0.0;//-90*pi/180.0;
                //comp ope = M2kfunc_rotatedcut(a,sigp,m,eps1,theta);//M2kfunc_secondsheet(a,sigp,m,eps1);//M2kfunc_1(a,sigp,m,eps1);
                //comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc_rotatedcut(a,sigp,m,eps1,theta);//M2kfunc_1(a,sigp,m,eps1);
                
                comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
                
                //comp ope = GS_pk(s,q,pcomp,m,eps);
                //comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)))*GS_pk(s,q,pcomp,m,eps);
                
                fout<<setprecision(16)<<real(pR)<<'\t'
                    <<setprecision(16)<<real(pI)<<'\t'
                    <<setprecision(16)<<real(ope)<<'\t'
                    <<setprecision(16)<<imag(ope)<<endl;
            
            }
        
        }
        fout.close();
        cout<<"kernelfile : "<<kernelfile<<endl;

        //Here we print the relevant thresholds for 
        //for the kernel and trace the pcut for the OPE

        string thresholdfile = //"threshold_a_" + to_string(a)
                            //+ "_eps_" + to_string(eps)
                            + "threshold_scount_" + to_string(scount)
                            + ".dat";

        comp qpole = pmom(s,sigmab(a,m),m);

        fout.open(thresholdfile.c_str());
        fout<<real(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<real(qpole)<<'\t'
            <<imag(qpole)<<'\t'
            <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<eps<<'\t'
            <<a<<endl;

        fout.close();
        cout<<"thresholdfile : "<<thresholdfile<<endl;

        delxreal = abs(xrealinitial-xrealfinal)/(20*points);

        string opecutfile = //"opecuttracer_a_" + to_string(a) 
                          //+ "_eps_" + to_string(eps)
                           "opecuttracer_scount_" + to_string(scount)
                          + ".dat";

        fout.open(opecutfile.c_str());  

        vector<double> xvec;
        vector<comp> pcutpvec;
        vector<comp> pcutmvec;

        for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        {
            xvec.push_back(x);
            comp pcutp = pcut_plus(x,s,q,m,eps);
            comp pcutm = pcut_minus(x,s,q,m,eps);
            pcutpvec.push_back(pcutp);
            pcutmvec.push_back(pcutm);

        }
        pcut_fixer(pcutpvec,pcutmvec);

        //for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        for(int i=0;i<pcutpvec.size();++i)
        {
            double x = xvec[i];
            comp pcutp = pcutpvec[i];//pcut_plus(x,s,q,m,eps);
            comp pcutm = pcutmvec[i];//pcut_minus(x,s,q,m,eps);

            fout<<x<<'\t'
                <<real(pcutp)<<'\t'
                <<imag(pcutp)<<'\t'
                <<real(pcutm)<<'\t'
                <<imag(pcutm)<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<abs(-q-pcutp)<<'\t'
                <<abs(q-pcutm)<<endl;
        }
        fout.close();

        cout<<"opecutfile : "<<opecutfile<<endl;

        string m2kcutfile =     "m2kcuttracer_scount_" + to_string(scount)
                            +   ".dat";
        fout.open(m2kcutfile.c_str());

        for(int i=0;i<5000;++i)
        {
            comp ii = {0.0,1.0};
            double pi = acos(-1.0);
            double theta = 0.0;//-90*pi/180;
            //since the cut is rotated -90 degree
            //we can just write it straightforward
            double delsigk = abs(3.0)/5000;
            //double delsig = -i*delsigk;
            
            double delsig = i*delsigk;
            
            //comp sigk = 4.0 + (comp)ii*delsig;
            comp sigk = 4.0 + (comp)delsig;
            
            
            //comp m2kcut = m2k_cut_tracer(sigk, m, theta);

            comp m2kcut_in_p = pmom(s,sigk,m);
            if(imag(m2kcut_in_p)>0.0)
            {
                comp temp = -real(m2kcut_in_p) - ii*imag(m2kcut_in_p);
                m2kcut_in_p = temp;
            }
            //cout<<sigk<<'\t'<<m2kcut_in_p<<endl;
            fout<<real(m2kcut_in_p)<<'\t'<<imag(m2kcut_in_p)<<endl;
        }
        //cout<<sigma_p(s,ii*m,m)<<endl;
        fout.close();

        cout<<"m2k_cutfile : "<<m2kcutfile<<endl;

        string pcutfile = //"pcut_a_" + to_string(a) 
                          //+ "_eps_" + to_string(eps)
                           "pcut_scount_" + to_string(scount)
                          + ".dat";

        fout.open(pcutfile.c_str());  

        pcount = 0;

        delxreal = abs(xrealinitial-xrealfinal)/((double)points);
        for(int k=0;k<qvec.size();++k)
        {
                             
            for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
            {
                comp pcutp = pcut_plus(x,s,qvec[k],m,eps);
                comp pcutm = pcut_minus(x,s,qvec[k],m,eps);
                fout<<x<<'\t'
                    <<real(pcutp)<<'\t'
                    <<imag(pcutp)<<'\t'
                    <<real(pcutm)<<'\t'
                    <<imag(pcutm)<<'\t'
                    <<real(q)<<'\t'
                    <<imag(q)<<'\t'
                    <<abs(-q-pcutp)<<'\t'
                    <<abs(q-pcutm)<<endl;
            }
            
            //cout<<"qvec["<<k<<"] = "<<qvec[k]<<'\t'<<"pcount = "<<pcount<<endl;
            //cout<<"pcutfile : "<<pcutfile<<endl;

            //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;
            //cout<<"file:"<<filename1<<endl;
            pcount = pcount + 1;
            //cout<<"loop completed = "<<pcount<<endl;
        }
        fout.close();
        cout<<"pcutfile : "<<pcutfile<<endl;
        cout<<endl;
        scount = scount + 1;
    }
}




void kernel_pcut_x_glockletest_vs_s3_allqvec_sigma_space()
{
    //cout<<"here"<<endl;
    ofstream fout;
    double a = -2;
    double m = 1.0;
    double sreal = 3.0;//3.0;
    double simag = 0.0;//-0.05;
    comp ii = {0.0,1.0};
    double eps = 0.0;//1.0e-3;//1.0e-4;
    double eps1 = eps;//0.0;

    //comp sigb = sigmab(a,m);
    comp sigb = 2.0*m*m;//sigmab(a,m);
    
    

    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int delspoints = 500;
    int qvecpoints = 500;
    double qvec_r = 0.01;

    double delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double rad = 0.16;
    double pRinitial = -1.6101;//-0.221;//-abs(abs(rad)+0.015);
    double pRfinal = 1.61;//0.22;//abs(abs(rad)+0.015);
    double delpR = abs(pRinitial-pRfinal)/points;
    double pIinitial = -1.51;//-0.221;//-abs(abs(rad)+0.015);
    double pIfinal = 1.5;//0.22;//abs(abs(rad)+0.015);
    double delpI = abs(pIinitial - pIfinal)/points;

    double sigRinitial = 3.901;
    double sigRfinal = 4.10;
    double delsigR = abs(sigRinitial-sigRfinal)/points;
    double sigIinitial = -0.1;
    double sigIfinal = 0.101;
    double delsigI = abs(sigIinitial-sigIfinal)/points;

    double sinitial = 8.60;//(double)real(phibthreshold(a,m));//8.72;
    double sfinal = (double)real(phibthreshold(a,m));//9.0;//(double)real(phibthreshold(a,m));//8.0;
    double dels = abs(sinitial-sfinal)/(double)delspoints;

    int scount = 0;

    //for(int i=0;i<delspoints+1;++i)
    //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    //int i=0;
    {
        
        //if(i!=0) break;
        cout<<"am = "<<a<<endl;
        //double s3 = sinitial + i*dels;
        double s3 = 0.0;//sfinal - i*dels;
        //if(scount!=0) break;
        int pcount = 0;
        
        s3 = 8.85;
        double delsimag = abs(-0.01 - 0.01)/delspoints;
        //simag = -0.01 + i*delsimag;
        simag = +0.00;
        comp s = s3 + ii*simag;

        //This section only focuses on the third node for the contour
        //comp q1 = pmom(s,sigmab(a,m),m);
        comp q1 = pmom(s,2.0,m);
        
        comp qc2 = q_c2(s,a,m,eps);
        comp pplus = pcut_plus(1.0,s,q1,m,eps);
        comp pminus = pcut_plus(-1.0,s,q1,m,eps);
        comp thirdnode = 0.5*real(pplus - pminus) - abs(qc2)*ii;
        comp qc1 = q_c1(s,a,m,eps);
        double rad = abs(qc1);
        comp fifthnode = (2.0/3.0)*(1.0-2.0*ii)*abs(qc2);
        cout<<"radius = "<<rad<<endl;
        /*
        if(rad>0.0)
        {
            pRinitial = -abs(abs(rad) + abs(rad)/3.0);
            pRfinal = abs(abs(rad) + abs(rad)/3.0);
            delpR = abs(pRinitial-pRfinal)/points;
            pIinitial = -abs(abs(rad) + abs(rad)/3.0);
            pIfinal = abs(abs(rad) + abs(rad)/3.0);
            delpI = abs(pIinitial - pIfinal)/points;
        }
        else
        {
            pRinitial = -1.01;
            pRfinal = 1.0;
            delpR = abs(pRinitial-pRfinal)/points;
            pIinitial = -1.01;
            pIfinal = 1.0;
            delpI = abs(pIinitial - pIfinal)/points;
        }
        */
        comp q = pmom(s,sigb-ii*eps1,m);
        comp kmax = pmom(s,0.0,m);
        cout<<"run = "<<scount+1<<endl;
    
        cout<<"s = "<<s<<'\t'<<"dels = "<<dels<<endl;
        cout<<"eps = "<<eps<<'\t'<<"eps_for_m2k = "<<eps1<<endl;

        vector<comp> qvec;
        vector<comp> weights;
        //int tag1,tag2;
        if(simag<0.0)
        {
            //mom_vector_maker_4(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            comp qc2 = q_c2(s,a,m,eps);
            cout<<"qc2 = "<<qc2<<endl;
            //mom_vector_maker_43_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,qvec_r);
            //mom_vector_maker_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints);
            qvecpoints = 500;
            //mom_vector_maker_43_with_weights_with_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,qvec_r);
            double reshift = 0.1;
            double imshift = -0.1;
            contour_for_resonance_2(qvec,weights,0.0,reshift,imshift,kmax,qvecpoints);
            //line_maker_with_weights(qvec,weights,0.0,kmax,qvecpoints);
                
            //mom_vector_maker_41(qvec,#)
        }
        else if(simag>=0.0)
        {
            //mom_vector_maker_4(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            //mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            //mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            int gvec_fixer_tag = 0;
            int tag1 = 0;
            int tag2 = 0;
            //mom_vector_maker_seba_imspos_2(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,tag1,tag2,gvec_fixer_tag);
            //mom_vector_maker_seba_imspos_4(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,tag1,tag2,gvec_fixer_tag);
            //mom_vector_maker_seba_imspos_5(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,tag1,tag2,gvec_fixer_tag);
            //line_maker_with_weights(qvec,weights,0.0,kmax,qvecpoints);
            //contour_for_resonance(qvec,weights,0.0,kmax,0.1,qvecpoints);
            cout<<M2kbranchcut_left_momrep_plus(s,m)<<endl;
            cout<<M2kbranchcut_left_momrep_minus(s,m)<<endl;
            cout<<M2kbranchcut_right_momrep_plus_eps(s,m,eps1)<<endl;
            cout<<M2kbranchcut_right_momrep_minus_eps(s,m,eps1)<<endl;

            double reshift = 0.1;
            double imshift = -0.1;
            contour_for_resonance_2(qvec,weights,0.0,reshift,imshift,kmax,qvecpoints);
            //contour_for_resonance_1(qvec,weights,M2kbranchcut_right_momrep_plus_eps(s,m,eps1),0.0,kmax,-abs(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))/10,qvecpoints);
            //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints);
            //mom_vector_maker_seba_imspos_2_with_contour47(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,tag1,tag2,gvec_fixer_tag);
            //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints);
            
            /*for(int qind=0;qind<qvec.size();++qind)
            {
                comp reqvec = real(qvec[qind]);
                comp imqvec = -imag(qvec[qind]);
                qvec[qind] = reqvec + ii*imqvec;
            }*/
        }
        //mom_vector_maker_61(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r,tag1,tag2);
        
        cout<<"qvec size = "<<qvec.size()<<endl;
        //cout<<"tag1 = "<<tag1<<'\t'<<"tag2 = "<<tag2<<endl;
        //for(int i=0;i<qvec.size();++i) cout<<i<<'\t'<<qvec[i]<<endl;

        //This part writes the qvec in to a file //
        string qvecfile =   //"qvec_a_" + to_string(a)
                          //+ "_eps_" + to_string(eps)
                           "qvec_scount_" + to_string(scount)
                          + ".dat";


        fout.open(qvecfile.c_str());
        for(int i=0;i<qvec.size();++i)
        {
            fout<<i<<'\t'<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
        }
        fout.close();
        cout<<"qvecfile : "<<qvecfile<<endl;

        //This part writes the kernel with k=q
        
        int isize = (int) points;
        int jsize = (int) points;

        string kernelfile =  // "kernel_" //a_" + to_string(a)
                            //+ "_eps_" + to_string(eps)
                             "kernel_scount_" + to_string(scount)
                            + ".dat";

        fout.open(kernelfile.c_str());
        for(int i=0;i<isize;++i)
        {
            for(int j=0;j<jsize;++j)
            {
                comp pR = (comp)pRinitial + (comp)i*delpR;
                comp pI = (comp)pIinitial + (comp)j*delpI; 
                comp pcomp = (comp) pR + ii*pI;

                comp sigR = (comp)sigRinitial + (comp)i*delsigR;
                comp sigI = (comp)sigIinitial + (comp)j*delsigI;

                comp sigp = sigR + ii*sigI;

                //comp sigp = sigma_p(s,pcomp,m);

                double pi = acos(-1.0);

                //comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
                //comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*M2kfunc(a,sigp,m,eps1);

                double theta = -90.0*pi/180.0;

                comp ope = M2kfunc_rotatedcut(a,sigp,m,eps1,theta);
                
                //comp ope = M2kfunc_secondsheet(a,sigp,m,eps1);//M2kfunc_1(a,sigp,m,eps1);
                
                //comp ope = GS_pk(s,q,pcomp,m,eps);
                //comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)))*GS_pk(s,q,pcomp,m,eps);
                
                fout<<setprecision(16)<<real(sigR)<<'\t'
                    <<setprecision(16)<<real(sigI)<<'\t'
                    <<setprecision(16)<<real(ope)<<'\t'
                    <<setprecision(16)<<imag(ope)<<endl;
            
            }
        
        }
        fout.close();
        cout<<"kernelfile : "<<kernelfile<<endl;

        //Here we print the relevant thresholds for 
        //for the kernel and trace the pcut for the OPE

        string thresholdfile = //"threshold_a_" + to_string(a)
                            //+ "_eps_" + to_string(eps)
                            + "threshold_scount_" + to_string(scount)
                            + ".dat";

        
        fout.open(thresholdfile.c_str());
        fout<<real(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<'\t'
            <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<eps<<'\t'
            <<a<<endl;

        fout.close();
        cout<<"thresholdfile : "<<thresholdfile<<endl;

        delxreal = abs(xrealinitial-xrealfinal)/(20*points);

        string opecutfile = //"opecuttracer_a_" + to_string(a) 
                          //+ "_eps_" + to_string(eps)
                           "opecuttracer_scount_" + to_string(scount)
                          + ".dat";

        fout.open(opecutfile.c_str());  

        vector<double> xvec;
        vector<comp> pcutpvec;
        vector<comp> pcutmvec;

        for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        {
            xvec.push_back(x);
            comp pcutp = pcut_plus(x,s,q,m,eps);
            comp pcutm = pcut_minus(x,s,q,m,eps);
            pcutpvec.push_back(pcutp);
            pcutmvec.push_back(pcutm);

        }
        pcut_fixer(pcutpvec,pcutmvec);

        //for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        for(int i=0;i<pcutpvec.size();++i)
        {
            double x = xvec[i];
            comp pcutp = pcutpvec[i];//pcut_plus(x,s,q,m,eps);
            comp pcutm = pcutmvec[i];//pcut_minus(x,s,q,m,eps);

            fout<<x<<'\t'
                <<real(pcutp)<<'\t'
                <<imag(pcutp)<<'\t'
                <<real(pcutm)<<'\t'
                <<imag(pcutm)<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<abs(-q-pcutp)<<'\t'
                <<abs(q-pcutm)<<endl;
        }
        fout.close();

        cout<<"opecutfile : "<<opecutfile<<endl;

        string pcutfile = //"pcut_a_" + to_string(a) 
                          //+ "_eps_" + to_string(eps)
                           "pcut_scount_" + to_string(scount)
                          + ".dat";

        fout.open(pcutfile.c_str());  

        pcount = 0;

        delxreal = abs(xrealinitial-xrealfinal)/((double)points);
        for(int k=0;k<qvec.size();++k)
        {
                             
            for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
            {
                comp pcutp = pcut_plus(x,s,qvec[k],m,eps);
                comp pcutm = pcut_minus(x,s,qvec[k],m,eps);
                fout<<x<<'\t'
                    <<real(pcutp)<<'\t'
                    <<imag(pcutp)<<'\t'
                    <<real(pcutm)<<'\t'
                    <<imag(pcutm)<<'\t'
                    <<real(q)<<'\t'
                    <<imag(q)<<'\t'
                    <<abs(-q-pcutp)<<'\t'
                    <<abs(q-pcutm)<<endl;
            }
            
            //cout<<"qvec["<<k<<"] = "<<qvec[k]<<'\t'<<"pcount = "<<pcount<<endl;
            //cout<<"pcutfile : "<<pcutfile<<endl;

            //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;
            //cout<<"file:"<<filename1<<endl;
            pcount = pcount + 1;
            //cout<<"loop completed = "<<pcount<<endl;
        }
        fout.close();
        cout<<"pcutfile : "<<pcutfile<<endl;
        cout<<endl;
        scount = scount + 1;
    }
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

void ope_qq_plot_vs_real_s()
{
    //cout<<"here"<<endl;
    ofstream fout;
    double a = 2.0;
    double m = 1.0;
    double sreal = 3.0;//3.0;
    double simag = +1.0e-5;//-0.05;
    comp ii = {0.0,1.0};
    double eps = 0.0;//1.0e-3;//1.0e-4;
    double eps1 = eps;//0.0;

    comp sigb = sigmab(a,m);
    


    int points = 250;
    int delspoints = 100;
    int qvecpoints = 1000;
    double qvec_r = 0.01;

    
    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50

    double sinitial = 3.9;//(double)real(phibthreshold(a,m));//8.72;
    double sfinal = 7.1;//(double)real(phibthreshold(a,m));//9.0;//(double)real(phibthreshold(a,m));//8.0;
    double dels = abs(sinitial-sfinal)/(double)delspoints;

    int scount = 0;

    string kernelfile =  // "kernel_" //a_" + to_string(a)
                            //+ "_eps_" + to_string(eps)
                             "ope_a_" + to_string((int)a)
                            + ".dat";

    fout.open(kernelfile.c_str());

    for(int i=0;i<delspoints+1;++i)
    //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    //int i=0;
    {
        
        //if(i!=0) break;
        cout<<"am = "<<a<<endl;
        double s3 = sinitial + i*dels;
        //if(scount!=0) break;
        int pcount = 0;
        comp s = s3 + ii*simag;

        //This section only focuses on the third node for the contour
        //comp q1 = pmom(s,sigmab(a,m),m);
        comp q1 = pmom(s,2.0,m);
        
        comp qc2 = q_c2(s,a,m,eps);
        comp pplus = pcut_plus(1.0,s,q1,m,eps);
        comp pminus = pcut_plus(-1.0,s,q1,m,eps);
        comp thirdnode = 0.5*real(pplus - pminus) - abs(qc2)*ii;
        comp qc1 = q_c1(real(s),a,m,eps);
        double rad = real(qc1);
        comp fifthnode = (2.0/3.0)*(1.0-2.0*ii)*abs(qc2);
        cout<<"radius = "<<rad<<endl;
        

        comp q = pmom(s,sigb-ii*eps1,m);
        comp kmax = pmom(s,0.0,m);
        cout<<"run = "<<scount+1<<endl;
    
        cout<<"s = "<<s<<'\t'<<"dels = "<<dels<<endl;
        cout<<"eps = "<<eps<<'\t'<<"eps_for_m2k = "<<eps1<<endl;

        comp sigp = sigma_p(s,q,m);
        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;

        comp ope = -gsq*GS_pk(s,q,q,m,eps);

        fout<<setprecision(16)<<real(s)<<'\t'
            <<setprecision(16)<<real(ope)<<'\t'
            <<setprecision(16)<<imag(ope)<<endl;
        
        
    }
    cout<<"kernelfile : "<<kernelfile<<endl;
    fout.close();
        
}



void mom_vector_maker_43_2( double somepoint,   
                            vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - somepoint*rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);
            
            cout<<"turningpoint 0 = "<<qvec[0]<<endl;
            startingpoint = qvec[qvec.size()-1];
            cout<<"turningpoint 1 = "<<startingpoint<<endl;
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            cout<<"turningpoint 2 = "<<startingpoint<<endl;
            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[lowest_p_index] - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[lowest_m_index] - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            if(lowest_index_yval<=imag(startingpoint))
            {
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 3 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0); 

                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 4 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = max;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
                cout<<"end point 5 = "<<qvec[qvec.size() - 1]<<endl;
            }
            else
            {
                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 3 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
                
                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 4 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = max;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
                cout<<"end point 5 = "<<qvec[qvec.size() - 1]<<endl;
            }
        }
        else
        {
            cout<<"deforming the contour downward first"<<endl;
            //comp startingpoint = min;
            //comp delq = abs(startingpoint - rad/2.0)/(points);
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,1);

            //startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            //comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/5.0);


            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            lowest_index_yval = imag(pcutpvec[lowest_p_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
            highest_index_xval = real(pcutpvec[highest_p_index]);
        }
        else
        {
            lowest_index_yval = imag(pcutmvec[lowest_m_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
            highest_index_xval = real(pcutpvec[highest_m_index]);
        }


        if(abs(lowest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour downward first"<<endl;
            comp startingpoint = min;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);


        }
        else 
        {
            cout<<"taking the straight line"<<endl;
            //cout<<"prob here with axis_flag="<<axis_flag<<endl;
            comp startingpoint = min;
            double totpoints = points;// + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);

        }
    }



    
}


//same as the previous function but for single 
//three body energy, s, value 
void kernel_pcut_x_glockletest_vs_single_s3_allqvec()
{
    
    ofstream fout;
    double a = 16.0;
    double m = 1.0;
    double sreal =  8.8;//(double) real(Grightbranchpoint(a,m)) + 1.0e-9;//real(phibthreshold(a,m)) - 1.0e-5; //8.94;//// 8.9999;//3.0;
    double simag = 1.0e-3;//0.5;//-0.05;
    comp ii = {0.0,1.0};
    double eps = 0.0;//1.0e-4;
    double eps1 = 0.0;

    comp sigb = sigmab(a,m);
    

    //initial setup

    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int delspoints = 300;
    int qvecpoints = 100;
    double qvec_r = 0.0;

    double delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    //contour x and y axis limits
    double pRinitial = -0.25101;
    double pRfinal = 0.251;
    double delpR = abs(pRinitial-pRfinal)/points;
    double pIinitial = -0.25101;
    double pIfinal = 0.251;
    double delpI = abs(pIinitial - pIfinal)/points;

    

    int scount = 0;

    
    int pcount = 0;
    comp s = sreal + ii*simag;
    comp q = pmom(s,sigb-ii*eps1,m);
    //comp q = pmom(s,2.0*m*m,m);
    
    comp kmax = pmom(s,0.0,m);
    cout<<"run = "<<scount+1<<endl;
    
    cout<<"s = "<<s<<endl;

    comp qq1 = q_c1(s,a,m,eps);
    comp qq2 = q_c2(s,a,m,eps);
    cout<<"qq1 = "<<qq1<<'\t'<<" qq2 = "<<qq2<<endl;

    vector<comp> qvec;
    vector<comp> weights;
    //int tag1,tag2;

    //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
    if(simag<=0.0)
    {
        //mom_vector_maker_4(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        mom_vector_maker_43_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,qvec_r);
        //mom_vector_maker_41(qvec,#)
        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        comp startingpoint = kmin;
        comp endingpoint = kmax;
        //line_maker(qvec,startingpoint,endingpoint,qvecpoints);
    }
    else if(simag>0.0)
    {
        //mom_vector_maker_4(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //int tag1 = 0;
        //int tag2 = 0;
        //mom_vector_maker_47_2_3_1(0.1,qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r,tag1,tag2);
        //mom_vector_maker_seba_imspos(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints);
        int gvec_fixer_tag = 0;
        int tag1 = 0;
        int tag2 = 0;
        //mom_vector_maker_seba_imspos_2(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,tag1,tag2,gvec_fixer_tag);
        //mom_vector_maker_seba_imsneg(qvec,weights,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps,(double)qvecpoints);
        mom_vector_maker_43_with_weights_with_seba_imsneg(qvec,weights,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps,(double)qvecpoints,qvec_r);
        //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        
        for(int qind=0;qind<qvec.size();++qind)
        {
            comp reqvec = real(qvec[qind]);
            comp imqvec = -imag(qvec[qind]);
            qvec[qind] = reqvec + ii*imqvec;
        }
    }
    //mom_vector_maker_61(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r,tag1,tag2);

    
        
    cout<<"qvec size = "<<qvec.size()<<endl;
    //cout<<"tag1 = "<<tag1<<'\t'<<"tag2 = "<<tag2<<endl;
    //for(int i=0;i<qvec.size();++i) cout<<i<<'\t'<<qvec[i]<<endl;

    //This part writes the qvec in to a file //
    string qvecfile =   //"qvec_" //+ to_string(a)
                        //+ "_eps_" + to_string(eps)
                        + "qvec_scount_" + to_string(scount)
                        + ".dat";


    fout.open(qvecfile.c_str());
    for(int i=0;i<qvec.size();++i)
    {
        fout<<i<<'\t'<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
    }
    fout.close();
    cout<<"qvecfile : "<<qvecfile<<endl;

    //This part writes the kernel with k=q
      
    int isize = (int) points;
    int jsize = (int) points;

    string kernelfile =   //"kernel_"// + to_string(a)
                        //+ "_eps_" + to_string(eps)
                         "kernel_scount_" + to_string(scount)
                        + ".dat";

    fout.open(kernelfile.c_str());
    for(int i=0;i<isize;++i)
    {
        for(int j=0;j<jsize;++j)
        {
            comp pR = (comp)pRinitial + (comp)i*delpR;
            comp pI = (comp)pIinitial + (comp)j*delpI; 
            comp pcomp = (comp) pR + ii*pI;
            comp sigp = sigma_p(s,pcomp,m);

            double pi = acos(-1.0);
            comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            //comp ope = GS_pk(s,q,pcomp,m,eps);//(pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            
            fout<<setprecision(16)<<real(pR)<<'\t'
                <<setprecision(16)<<real(pI)<<'\t'
                <<setprecision(16)<<real(ope)<<'\t'
                <<setprecision(16)<<imag(ope)<<endl;
            
        }
        
    }
    fout.close();
    cout<<"kernelfile : "<<kernelfile<<endl;

    //Here we print the relevant thresholds for 
    //for the kernel and trace the pcut for the OPE

    string thresholdfile = //"threshold_" //+ to_string(a)
                        //+ "_eps_" + to_string(eps)
                            "threshold_scount_" + to_string(scount)
                            + ".dat";

        
    fout.open(thresholdfile.c_str());
    fout<<real(pcut_plus(1.0,s,q,m,eps))<<'\t'
        <<imag(pcut_plus(1.0,s,q,m,eps))<<'\t'
        <<real(pcut_plus(-1.0,s,q,m,eps))<<'\t'
        <<imag(pcut_plus(-1.0,s,q,m,eps))<<'\t'
        <<real(pcut_minus(1.0,s,q,m,eps))<<'\t'
        <<imag(pcut_minus(1.0,s,q,m,eps))<<'\t'
        <<real(pcut_minus(-1.0,s,q,m,eps))<<'\t'
        <<imag(pcut_minus(-1.0,s,q,m,eps))<<'\t'
        <<real(q)<<'\t'
        <<imag(q)<<'\t'
        <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
        <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
        <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
        <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
        <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
        <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
        <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
        <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
        <<real(s)<<'\t'
        <<imag(s)<<'\t'
        <<eps<<'\t'
        <<a<<endl;

    fout.close();
    cout<<"thresholdfile : "<<thresholdfile<<endl;

    delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);

    string opecutfile = //"opecuttracer_"// + to_string(a) 
                          //+ "_eps_" + to_string(eps)
                          + "opecuttracer_scount_" + to_string(scount)
                          + ".dat";

    fout.open(opecutfile.c_str());  

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
    {
        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);
        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);

    }
    pcut_fixer(pcutpvec,pcutmvec);

    //for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
    for(int i=0;i<pcutpvec.size();++i)
    {
        double x = xvec[i];
        comp pcutp = pcutpvec[i];//pcut_plus(x,s,q,m,eps);
        comp pcutm = pcutmvec[i];//pcut_minus(x,s,q,m,eps);

        fout<<x<<'\t'
            <<real(pcutp)<<'\t'
            <<imag(pcutp)<<'\t'
            <<real(pcutm)<<'\t'
            <<imag(pcutm)<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<'\t'
            <<abs(-q-pcutp)<<'\t'
            <<abs(q-pcutm)<<endl;
    }
    fout.close();

    cout<<"opecutfile : "<<opecutfile<<endl;

    string pcutfile = //"pcut_"// + to_string(a) 
                    //+ "_eps_" + to_string(eps)
                    + "pcut_scount_" + to_string(scount)
                    + ".dat";

    fout.open(pcutfile.c_str());  

    pcount = 0;

    delxreal = abs(xrealinitial-xrealfinal)/((double)points);
    for(int k=0;k<qvec.size();++k)
    {
                             
        for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        {
            comp pcutp = pcut_plus(x,s,qvec[k],m,eps);
            comp pcutm = pcut_minus(x,s,qvec[k],m,eps);
            fout<<x<<'\t'
                <<real(pcutp)<<'\t'
                <<imag(pcutp)<<'\t'
                <<real(pcutm)<<'\t'
                <<imag(pcutm)<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<abs(-q-pcutp)<<'\t'
                <<abs(q-pcutm)<<endl;
        }
            
            //cout<<"qvec["<<k<<"] = "<<qvec[k]<<'\t'<<"pcount = "<<pcount<<endl;
            //cout<<"pcutfile : "<<pcutfile<<endl;

            //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;
            //cout<<"file:"<<filename1<<endl;
            pcount = pcount + 1;
            //cout<<"loop completed = "<<pcount<<endl;
    }
    fout.close();
    cout<<"pcutfile : "<<pcutfile<<endl;
    cout<<endl;
    scount = scount + 1;
    
}

void kernel_pcut_x_glockletest_vs_single_s3_singleqvec()
{
    
    ofstream fout;
    double a = 16.0;
    double m = 1.0;
    double sreal = 8.8;//3.0;
    double simag = 1.0e-5;//0.1;//0.05;//-0.05;
    comp ii = {0.0,1.0};
    double eps = 0.0;//1.0e-4;
    double eps1 = 0.0;

    comp sigb = sigmab(a,m);
    

    //initial setup

    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int delspoints = 300;
    int qvecpoints = 500;
    double qvec_r = 0.01;

    double delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    //contour x and y axis limits
    double pRinitial = -0.2051;
    double pRfinal = 0.205;
    double delpR = abs(pRinitial-pRfinal)/points;
    double pIinitial = -0.2051;
    double pIfinal = 0.205;
    double delpI = abs(pIinitial - pIfinal)/points;

    

    int scount = 0;

    
    int pcount = 0;
    comp s = sreal + ii*simag;
    comp q = pmom(s,sigb-ii*eps1,m);
    //comp q = pmom(s,2.0*m*m,m);
    
    comp kmax = pmom(s,0.0,m);
    cout<<"run = "<<scount+1<<endl;
    
    cout<<"s = "<<s<<endl;

    comp qq1 = q_c1(s,a,m,eps);
    comp qq2 = q_c2(s,a,m,eps);
    cout<<"qq1 = "<<qq1<<'\t'<<" qq2 = "<<qq2<<endl;

    vector<comp> qvec;
    vector<comp> weights; //added the weights vector last time
    //int tag1,tag2;
    int gvec_fixer_tag = 0;
    int tag1 = 0;
    int tag2 = 0;

    //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
    if(simag<=0.0)
    {
        //mom_vector_maker_4(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_41(qvec,#)
        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        comp startingpoint = kmin;
        comp endingpoint = kmax;
        line_maker(qvec,startingpoint,endingpoint,qvecpoints);
    }
    else if(simag>0.0)
    {
        //mom_vector_maker_4(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //int tag1 = 0;
        //int tag2 = 0;
        //mom_vector_maker_47_2_2(0.065,qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r,tag1,tag2);
        //mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_47_2_A1(1000,qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r,tag1,tag2);
        //mom_vector_maker_seba_imspos(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints);
        //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints);
        
        mom_vector_maker_seba_imspos_2(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,tag1,tag2,gvec_fixer_tag);
        //mom_vector_maker_seba_imspos_2_with_contour47(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,tag1,tag2,gvec_fixer_tag);
        
        //mom_vector_maker_43_with_weights_with_seba_imsneg(qvec,weights,real(s)-ii*abs(imag(s)),0.0,kmax,a,m,eps,eps,(double)qvecpoints,qvec_r);

        //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints);
        /*for(int qind=0;qind<qvec.size();++qind)
        {
            comp reqvec = real(qvec[qind]);
            comp imqvec = -imag(qvec[qind]);
            qvec[qind] = reqvec + ii*imqvec;
        }*/
    }
    //mom_vector_maker_61(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r,tag1,tag2);
    
    /*int size = qvec.size();
    Eigen::VectorXcd Gvec(size);
    Eigen::VectorXcd dsol(size);
    comp qval = q;

            
    if(simag<0.0)
    {
        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
    }
    else if(simag>=0.0)
    {
        Gvec_maker_momrep_withtags_1(Gvec,s,qvec,qval,m,eps,tag1,tag2);
    }

    int switch_for_gvec_fixer=0;
    if(switch_for_gvec_fixer==0)
    {
        //Gvec_fixer_1(Gvec,Gvec,qvec,s,qval,m,eps,tag1,tag2); 
        //Gvec_fixer_2(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
        Gvec_fixer_3(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
        tag_fixer_for_gvec_fixer3(Gvec,qvec,s,qval,m,eps,tag1);
               
    }*/

    cout<<"tag fixer = "<<tag1<<'\t'<<tag2<<endl;
        
    cout<<"qvec size = "<<qvec.size()<<endl;
    cout<<"tags = "<<tag1<<'\t'<<tag2<<endl;
    cout<<"tagged qvec = "<<qvec[tag1]<<'\t'<<qvec[tag2]<<endl;
    //cout<<"tag1 = "<<tag1<<'\t'<<"tag2 = "<<tag2<<endl;
    //for(int i=0;i<qvec.size();++i) cout<<i<<'\t'<<qvec[i]<<endl;

    //This part writes the qvec in to a file //
    string qvecfile =   //"qvec_" //+ to_string(a)
                        //+ "_eps_" + to_string(eps)
                        + "qvec_scount_" + to_string(scount)
                        + ".dat";


    fout.open(qvecfile.c_str());
    for(int i=0;i<qvec.size();++i)
    {
        fout<<i<<'\t'<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
    }
    fout.close();
    cout<<"qvecfile : "<<qvecfile<<endl;

    //This part writes the kernel with k=q
      
    int isize = (int) points;
    int jsize = (int) points;

    string kernelfile =   //"kernel_"// + to_string(a)
                        //+ "_eps_" + to_string(eps)
                         "kernel_scount_" + to_string(scount)
                        + ".dat";

    fout.open(kernelfile.c_str());
    for(int i=0;i<isize;++i)
    {
        for(int j=0;j<jsize;++j)
        {
            comp pR = (comp)pRinitial + (comp)i*delpR;
            comp pI = (comp)pIinitial + (comp)j*delpI; 
            comp pcomp = (comp) pR + ii*pI;
            comp sigp = sigma_p(s,pcomp,m);

            double pi = acos(-1.0);
            comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            //comp ope = GS_pk(s,q,pcomp,m,eps);//(pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            
            fout<<setprecision(16)<<real(pR)<<'\t'
                <<setprecision(16)<<real(pI)<<'\t'
                <<setprecision(16)<<real(ope)<<'\t'
                <<setprecision(16)<<imag(ope)<<endl;
            
        }
        
    }
    fout.close();
    cout<<"kernelfile : "<<kernelfile<<endl;

    //Here we print the relevant thresholds for 
    //for the kernel and trace the pcut for the OPE

    string thresholdfile = //"threshold_" //+ to_string(a)
                        //+ "_eps_" + to_string(eps)
                            "threshold_scount_" + to_string(scount)
                            + ".dat";

        
    fout.open(thresholdfile.c_str());
    fout<<real(pcut_plus(1.0,s,q,m,eps))<<'\t'
        <<imag(pcut_plus(1.0,s,q,m,eps))<<'\t'
        <<real(pcut_plus(-1.0,s,q,m,eps))<<'\t'
        <<imag(pcut_plus(-1.0,s,q,m,eps))<<'\t'
        <<real(pcut_minus(1.0,s,q,m,eps))<<'\t'
        <<imag(pcut_minus(1.0,s,q,m,eps))<<'\t'
        <<real(pcut_minus(-1.0,s,q,m,eps))<<'\t'
        <<imag(pcut_minus(-1.0,s,q,m,eps))<<'\t'
        <<real(q)<<'\t'
        <<imag(q)<<'\t'
        <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
        <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
        <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
        <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
        <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
        <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
        <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
        <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
        <<real(s)<<'\t'
        <<imag(s)<<'\t'
        <<eps<<'\t'
        <<a<<endl;

    fout.close();
    cout<<"thresholdfile : "<<thresholdfile<<endl;

    delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);

    string opecutfile = //"opecuttracer_"// + to_string(a) 
                          //+ "_eps_" + to_string(eps)
                          + "opecuttracer_scount_" + to_string(scount)
                          + ".dat";

    fout.open(opecutfile.c_str());  

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
    {
        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);
        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);

    }
    pcut_fixer(pcutpvec,pcutmvec);

    //for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
    for(int i=0;i<pcutpvec.size();++i)
    {
        double x = xvec[i];
        comp pcutp = pcutpvec[i];//pcut_plus(x,s,q,m,eps);
        comp pcutm = pcutmvec[i];//pcut_minus(x,s,q,m,eps);

        fout<<x<<'\t'
            <<real(pcutp)<<'\t'
            <<imag(pcutp)<<'\t'
            <<real(pcutm)<<'\t'
            <<imag(pcutm)<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<'\t'
            <<abs(-q-pcutp)<<'\t'
            <<abs(q-pcutm)<<endl;
    }
    fout.close();

    cout<<"opecutfile : "<<opecutfile<<endl;

    

    pcount = 0;

    delxreal = abs(xrealinitial-xrealfinal)/((double)points);
    for(int k=0;k<qvec.size();++k)
    {
        string pcutfile = //"pcut_"// + to_string(a) 
                    //+ "_eps_" + to_string(eps)
                    + "pcut_scount_" + to_string(scount)
                    + "_pcount_" + to_string(pcount)
                    + ".dat";

        fout.open(pcutfile.c_str());                          
        for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        {
            comp pcutp = pcut_plus(x,s,qvec[k],m,eps);
            comp pcutm = pcut_minus(x,s,qvec[k],m,eps);
            fout<<x<<'\t'
                <<real(pcutp)<<'\t'
                <<imag(pcutp)<<'\t'
                <<real(pcutm)<<'\t'
                <<imag(pcutm)<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<abs(-q-pcutp)<<'\t'
                <<abs(q-pcutm)<<endl;
        }
            
        //cout<<"qvec["<<k<<"] = "<<qvec[k]<<'\t'<<"pcount = "<<pcount<<endl;
        //cout<<"pcutfile : "<<pcutfile<<endl;

        //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;
        //cout<<"file:"<<filename1<<endl;
        pcount = pcount + 1;
        //cout<<"loop completed = "<<pcount<<endl;
        fout.close();
        cout<<"pcutfile : "<<pcutfile<<endl;
    }
    
    cout<<endl;
    scount = scount + 1;
    
}


void kernel_pcut_x_glockletest_vs_single_s3_singleqvec_nonanalyticity_singleplot()
{
    
    ofstream fout;
    double a = 16.0;
    double m = 1.0;
    double sreal = 8.5;//3.0;
    double simag = 0.1;//0.1;//0.05;//-0.05;
    comp ii = {0.0,1.0};
    double eps = 0.0;//1.0e-4;
    double eps1 = 0.0;

    comp sigb = sigmab(a,m);
    

    //initial setup

    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int delspoints = 300;
    int qvecpoints = 1500;
    double qvec_r = 0.01;

    double delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    //contour x and y axis limits
    double pRinitial = -0.2051;
    double pRfinal = 0.205;
    double delpR = abs(pRinitial-pRfinal)/points;
    double pIinitial = -0.2051;
    double pIfinal = 0.205;
    double delpI = abs(pIinitial - pIfinal)/points;

    

    int scount = 0;

    
    int pcount = 0;
    comp s = sreal + ii*simag;
    comp q = pmom(s,sigb-ii*eps1,m);
    //comp q = pmom(s,2.0*m*m,m);
    
    comp kmax = pmom(s,0.0,m);
    cout<<"run = "<<scount+1<<endl;
    
    cout<<"s = "<<s<<endl;

    comp qq1 = q_c1(s,a,m,eps);
    comp qq2 = q_c2(s,a,m,eps);
    cout<<"qq1 = "<<qq1<<'\t'<<" qq2 = "<<qq2<<endl;

    vector<comp> qvec;
    vector<comp> weights; //added the weights vector last time
    int tag1,tag2;

    //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
    if(simag<=0.0)
    {
        //mom_vector_maker_4(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_41(qvec,#)
        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        comp startingpoint = kmin;
        comp endingpoint = kmax;
        line_maker(qvec,startingpoint,endingpoint,qvecpoints);
    }
    else if(simag>0.0)
    {
        //mom_vector_maker_4(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //int tag1 = 0;
        //int tag2 = 0;
        //mom_vector_maker_47_2_2(0.065,qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r,tag1,tag2);
        //mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_47_2_A1(1000,qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r,tag1,tag2);
        //mom_vector_maker_seba_imspos(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints);
        //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints);
        int gvec_fixer_tag = 0;
        int tag1 = 0;
        int tag2 = 0;
        //mom_vector_maker_seba_imspos_2(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,tag1,tag2,gvec_fixer_tag);
        //mom_vector_maker_seba_imspos_2_with_contour47(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,tag1,tag2,gvec_fixer_tag);
        //mom_vector_maker_seba_imspos_2_with_contour47(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,tag1,tag2,gvec_fixer_tag);
        //mom_vector_maker_seba_imspos_4(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,tag1,tag2,gvec_fixer_tag);
        
        mom_vector_maker_seba_imspos_5(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints,tag1,tag2,gvec_fixer_tag);
        //mom_vector_maker_43_with_weights_with_seba_imsneg(qvec,weights,real(s)-ii*abs(imag(s)),0.0,kmax,a,m,eps,eps,(double)qvecpoints,qvec_r);

        //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)qvecpoints);
        /*for(int qind=0;qind<qvec.size();++qind)
        {
            comp reqvec = real(qvec[qind]);
            comp imqvec = -imag(qvec[qind]);
            qvec[qind] = reqvec + ii*imqvec;
        }*/
    }
    //mom_vector_maker_61(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r,tag1,tag2);


    //------------------------------------------
    //this portion is to add contours by hand
    /*ifstream fin;
    string contour_file = "for_digonto_sr8.8_si0.001.txt";
    fin.open(contour_file.c_str());
    double qvecx = 0.0;
    double qvecy = 0.0;
    vector<comp> contour;
    while(fin>>qvecx>>qvecy)
    {
        comp contour_from_file = qvecx + ii*qvecy;
        contour.push_back(contour_from_file); 
    }
    fin.close();
    qvec = contour;
    */
    //--------------------------------------------
        
    cout<<"qvec size = "<<qvec.size()<<endl;
    //cout<<"tag1 = "<<tag1<<'\t'<<"tag2 = "<<tag2<<endl;
    //for(int i=0;i<qvec.size();++i) cout<<i<<'\t'<<qvec[i]<<endl;

    //This part writes the qvec in to a file //
    string qvecfile =   //"qvec_" //+ to_string(a)
                        //+ "_eps_" + to_string(eps)
                        + "qvec_scount_" + to_string(scount)
                        + ".dat";


    fout.open(qvecfile.c_str());
    for(int i=0;i<qvec.size();++i)
    {
        fout<<i<<'\t'<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
    }
    fout.close();
    cout<<"qvecfile : "<<qvecfile<<endl;

    //This part writes the kernel with k=q
      
    int isize = (int) points;
    int jsize = (int) points;

    string kernelfile =   //"kernel_"// + to_string(a)
                        //+ "_eps_" + to_string(eps)
                         "kernel_scount_" + to_string(scount)
                        + ".dat";

    fout.open(kernelfile.c_str());
    for(int i=0;i<isize;++i)
    {
        for(int j=0;j<jsize;++j)
        {
            comp pR = (comp)pRinitial + (comp)i*delpR;
            comp pI = (comp)pIinitial + (comp)j*delpI; 
            comp pcomp = (comp) pR + ii*pI;
            comp sigp = sigma_p(s,pcomp,m);

            double pi = acos(-1.0);
            comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            //comp ope = GS_pk(s,q,pcomp,m,eps);//(pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            
            fout<<setprecision(16)<<real(pR)<<'\t'
                <<setprecision(16)<<real(pI)<<'\t'
                <<setprecision(16)<<real(ope)<<'\t'
                <<setprecision(16)<<imag(ope)<<endl;
            
        }
        
    }
    fout.close();
    cout<<"kernelfile : "<<kernelfile<<endl;

    //Here we print the relevant thresholds for 
    //for the kernel and trace the pcut for the OPE

    string thresholdfile = //"threshold_" //+ to_string(a)
                        //+ "_eps_" + to_string(eps)
                            "threshold_scount_" + to_string(scount)
                            + ".dat";

        
    fout.open(thresholdfile.c_str());
    fout<<real(pcut_plus(1.0,s,q,m,eps))<<'\t'
        <<imag(pcut_plus(1.0,s,q,m,eps))<<'\t'
        <<real(pcut_plus(-1.0,s,q,m,eps))<<'\t'
        <<imag(pcut_plus(-1.0,s,q,m,eps))<<'\t'
        <<real(pcut_minus(1.0,s,q,m,eps))<<'\t'
        <<imag(pcut_minus(1.0,s,q,m,eps))<<'\t'
        <<real(pcut_minus(-1.0,s,q,m,eps))<<'\t'
        <<imag(pcut_minus(-1.0,s,q,m,eps))<<'\t'
        <<real(q)<<'\t'
        <<imag(q)<<'\t'
        <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
        <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
        <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
        <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
        <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
        <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
        <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
        <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
        <<real(s)<<'\t'
        <<imag(s)<<'\t'
        <<eps<<'\t'
        <<a<<endl;

    fout.close();
    cout<<"thresholdfile : "<<thresholdfile<<endl;

    delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);

    string opecutfile = //"opecuttracer_"// + to_string(a) 
                          //+ "_eps_" + to_string(eps)
                          + "opecuttracer_scount_" + to_string(scount)
                          + ".dat";

    fout.open(opecutfile.c_str());  

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
    {
        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);
        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);

    }
    pcut_fixer(pcutpvec,pcutmvec);

    //for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
    for(int i=0;i<pcutpvec.size();++i)
    {
        double x = xvec[i];
        comp pcutp = pcutpvec[i];//pcut_plus(x,s,q,m,eps);
        comp pcutm = pcutmvec[i];//pcut_minus(x,s,q,m,eps);

        fout<<x<<'\t'
            <<real(pcutp)<<'\t'
            <<imag(pcutp)<<'\t'
            <<real(pcutm)<<'\t'
            <<imag(pcutm)<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<'\t'
            <<abs(-q-pcutp)<<'\t'
            <<abs(q-pcutm)<<endl;
    }
    fout.close();

    cout<<"opecutfile : "<<opecutfile<<endl;

    

    pcount = 0;

    string pcutfile = //"pcut_"// + to_string(a) 
                    //+ "_eps_" + to_string(eps)
                    + "pcut_scount_" + to_string(scount)
                    //+ "_pcount_" + to_string(pcount)
                    + ".dat";

    fout.open(pcutfile.c_str());     

    delxreal = abs(xrealinitial-xrealfinal)/((double)points);
    for(int k=0;k<qvec.size();++k)
    {
        for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        {
            comp pcutp = pcut_plus(x,s,qvec[k],m,eps);
            comp pcutm = pcut_minus(x,s,qvec[k],m,eps);
            fout<<x<<'\t'
                <<real(pcutp)<<'\t'
                <<imag(pcutp)<<'\t'
                <<real(pcutm)<<'\t'
                <<imag(pcutm)<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<abs(-q-pcutp)<<'\t'
                <<abs(q-pcutm)<<endl;
        }
            
        //cout<<"qvec["<<k<<"] = "<<qvec[k]<<'\t'<<"pcount = "<<pcount<<endl;
        //cout<<"pcutfile : "<<pcutfile<<endl;

        //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;
        //cout<<"file:"<<filename1<<endl;
        pcount = pcount + 1;
        //cout<<"loop completed = "<<pcount<<endl;
        
        cout<<"pcutfile : "<<pcutfile<<endl;
    }
    fout.close();

    cout<<endl;
    scount = scount + 1;
    
}



void kernel_pcut_x_glockletest_vs_a_allqvec()
{
    //cout<<"here"<<endl;
    ofstream fout;
    //double a = 1.42;
    double m = 1.0;
    double sreal = 3.0;//3.0;
    double simag = 0.0;//-0.05;
    comp ii = {0.0,1.0};
    double eps = 1.0e-4;
    double eps1 = 0.0;//eps;//0.0;

    
    


    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int delspoints = 300;
    int qvecpoints = 50;
    double qvec_r = 0.01;

    double delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double pRinitial = -1.0;
    double pRfinal = 1.0;
    double delpR = abs(pRinitial-pRfinal)/points;
    double pIinitial = -1.0;
    double pIfinal = 1.0;
    double delpI = abs(pIinitial - pIfinal)/points;

    //double sinitial = 3.0;//8.90667;//8.72;//(double)real(phibthreshold(a,m));//8.72;
    //double sfinal = (double)real(phibthreshold(a,m));//9.0;//(double)real(phibthreshold(a,m));//8.0;
    //double dels = abs(sinitial-sfinal)/(double)delspoints;

    double delapoints = 50.0;
    double initiala = 1.4;
    double finala = 2.0;
    double dela = abs(initiala - finala)/delapoints;

    int scount = 0;

    for(int i=0;i<delapoints+1;++i)
    //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    {

        double a = initiala + i*dela;
        cout<<"am = "<<a<<endl;
        
        comp sigb = sigmab(a,m);
        //if(i!=0) break;
        //double s3 = sinitial + i*dels;
        //if(scount!=0) break;
        int pcount = 0;
        comp s = sreal + ii*simag;
        comp q = pmom(s,sigb-ii*eps1,m);
        comp kmax = pmom(s,0.0,m);
        cout<<"run = "<<scount+1<<endl;
    
        cout<<"s = "<<s<<endl;//'\t'<<"dels = "<<dels<<endl;
        cout<<"eps = "<<eps<<'\t'<<"eps_for_m2k = "<<eps1<<endl;

        vector<comp> qvec;
        int tag1,tag2;
        if(simag<=0.0)
        //mom_vector_maker_4(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_41(qvec,#)
        else if(simag>0.0)
        {
            //mom_vector_maker_4(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            for(int qind=0;qind<qvec.size();++qind)
            {
                comp reqvec = real(qvec[qind]);
                comp imqvec = -imag(qvec[qind]);
                qvec[qind] = reqvec + ii*imqvec;
            }
        }
        //mom_vector_maker_61(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r,tag1,tag2);
        
        cout<<"qvec size = "<<qvec.size()<<endl;
        //cout<<"tag1 = "<<tag1<<'\t'<<"tag2 = "<<tag2<<endl;
        //for(int i=0;i<qvec.size();++i) cout<<i<<'\t'<<qvec[i]<<endl;

        //This part writes the qvec in to a file //
        string qvecfile =   //"qvec_a_" + to_string(a)
                          //+ "_eps_" + to_string(eps)
                           "qvec_scount_" + to_string(scount)
                          + ".dat";


        fout.open(qvecfile.c_str());
        for(int i=0;i<qvec.size();++i)
        {
            fout<<i<<'\t'<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
        }
        fout.close();
        cout<<"qvecfile : "<<qvecfile<<endl;

        //This part writes the kernel with k=q
        
        int isize = (int) points;
        int jsize = (int) points;

        string kernelfile =  // "kernel_" //a_" + to_string(a)
                            //+ "_eps_" + to_string(eps)
                             "kernel_scount_" + to_string(scount)
                            + ".dat";

        fout.open(kernelfile.c_str());
        for(int i=0;i<isize;++i)
        {
            for(int j=0;j<jsize;++j)
            {
                comp pR = (comp)pRinitial + (comp)i*delpR;
                comp pI = (comp)pIinitial + (comp)j*delpI; 
                comp pcomp = (comp) pR + ii*pI;

                comp sigp = sigma_p(s,pcomp,m);

                double pi = acos(-1.0);

                comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            
                fout<<setprecision(16)<<real(pR)<<'\t'
                    <<setprecision(16)<<real(pI)<<'\t'
                    <<setprecision(16)<<real(ope)<<'\t'
                    <<setprecision(16)<<imag(ope)<<endl;
            
            }
        
        }
        fout.close();
        cout<<"kernelfile : "<<kernelfile<<endl;

        //Here we print the relevant thresholds for 
        //for the kernel and trace the pcut for the OPE

        string thresholdfile = //"threshold_a_" + to_string(a)
                            //+ "_eps_" + to_string(eps)
                            + "threshold_scount_" + to_string(scount)
                            + ".dat";

        
        fout.open(thresholdfile.c_str());
        fout<<real(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<'\t'
            <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<eps<<'\t'
            <<a<<endl;

        fout.close();
        cout<<"thresholdfile : "<<thresholdfile<<endl;

        delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);

        string opecutfile = //"opecuttracer_a_" + to_string(a) 
                          //+ "_eps_" + to_string(eps)
                           "opecuttracer_scount_" + to_string(scount)
                          + ".dat";

        fout.open(opecutfile.c_str());  

        vector<double> xvec;
        vector<comp> pcutpvec;
        vector<comp> pcutmvec;

        for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        {
            xvec.push_back(x);
            comp pcutp = pcut_plus(x,s,q,m,eps);
            comp pcutm = pcut_minus(x,s,q,m,eps);
            pcutpvec.push_back(pcutp);
            pcutmvec.push_back(pcutm);

        }
        pcut_fixer(pcutpvec,pcutmvec);

        //for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        for(int i=0;i<pcutpvec.size();++i)
        {
            double x = xvec[i];
            comp pcutp = pcutpvec[i];//pcut_plus(x,s,q,m,eps);
            comp pcutm = pcutmvec[i];//pcut_minus(x,s,q,m,eps);

            fout<<x<<'\t'
                <<real(pcutp)<<'\t'
                <<imag(pcutp)<<'\t'
                <<real(pcutm)<<'\t'
                <<imag(pcutm)<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<abs(-q-pcutp)<<'\t'
                <<abs(q-pcutm)<<endl;
        }
        fout.close();

        cout<<"opecutfile : "<<opecutfile<<endl;

        string pcutfile = //"pcut_a_" + to_string(a) 
                          //+ "_eps_" + to_string(eps)
                           "pcut_scount_" + to_string(scount)
                          + ".dat";

        fout.open(pcutfile.c_str());  

        pcount = 0;

        delxreal = abs(xrealinitial-xrealfinal)/((double)points);
        for(int k=0;k<qvec.size();++k)
        {
                             
            for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
            {
                comp pcutp = pcut_plus(x,s,qvec[k],m,eps);
                comp pcutm = pcut_minus(x,s,qvec[k],m,eps);
                fout<<x<<'\t'
                    <<real(pcutp)<<'\t'
                    <<imag(pcutp)<<'\t'
                    <<real(pcutm)<<'\t'
                    <<imag(pcutm)<<'\t'
                    <<real(q)<<'\t'
                    <<imag(q)<<'\t'
                    <<abs(-q-pcutp)<<'\t'
                    <<abs(q-pcutm)<<endl;
            }
            
            //cout<<"qvec["<<k<<"] = "<<qvec[k]<<'\t'<<"pcount = "<<pcount<<endl;
            //cout<<"pcutfile : "<<pcutfile<<endl;

            //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;
            //cout<<"file:"<<filename1<<endl;
            pcount = pcount + 1;
            //cout<<"loop completed = "<<pcount<<endl;
        }
        fout.close();
        cout<<"pcutfile : "<<pcutfile<<endl;
        cout<<endl;
        scount = scount + 1;
    }
}


void kernel_pcut_x_glockletest_vs_single_a_allqvec()
{
    //cout<<"here"<<endl;
    ofstream fout;
    //double a = 1.42;
    double m = 1.0;
    double sreal = 3.0;//3.0;
    double simag = 0.0;//-0.05;
    comp ii = {0.0,1.0};
    double eps = 1.0e-3;//1.0e-4;
    double eps1 = 0.0;//eps;//0.0;

    
    


    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int delspoints = 300;
    int qvecpoints = 50;
    double qvec_r = 0.01;

    double delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double pRinitial = -1.0;
    double pRfinal = 1.0;
    double delpR = abs(pRinitial-pRfinal)/points;
    double pIinitial = -1.0;
    double pIfinal = 1.0;
    double delpI = abs(pIinitial - pIfinal)/points;

    //double sinitial = 3.0;//8.90667;//8.72;//(double)real(phibthreshold(a,m));//8.72;
    //double sfinal = (double)real(phibthreshold(a,m));//9.0;//(double)real(phibthreshold(a,m));//8.0;
    //double dels = abs(sinitial-sfinal)/(double)delspoints;

    double delapoints = 50.0;
    double initiala = 1.4;
    double finala = 2.0;
    double dela = abs(initiala - finala)/delapoints;

    int scount = 0;

    //for(int i=0;i<delapoints+1;++i)
    //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    //{

        //double a = initiala + i*dela;

        double a = 2.0;
        cout<<"am = "<<a<<endl;
        
        comp sigb = sigmab(a,m);
        //if(i!=0) break;
        //double s3 = sinitial + i*dels;
        //if(scount!=0) break;
        int pcount = 0;
        comp s = sreal + ii*simag;
        comp q = pmom(s,sigb-ii*eps1,m);
        comp kmax = pmom(s,0.0,m);
        cout<<"run = "<<scount+1<<endl;
    
        cout<<"s = "<<s<<endl;//'\t'<<"dels = "<<dels<<endl;
        cout<<"eps = "<<eps<<'\t'<<"eps_for_m2k = "<<eps1<<endl;

        vector<comp> qvec;
        int tag1,tag2;
        if(simag<=0.0)
        //mom_vector_maker_4(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_41(qvec,#)
        else if(simag>0.0)
        {
            //mom_vector_maker_4(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
            //for(int qind=0;qind<qvec.size();++qind)
            //{
            //    comp reqvec = real(qvec[qind]);
            //    comp imqvec = -imag(qvec[qind]);
            //    qvec[qind] = reqvec + ii*imqvec;
            //}

        }
        //mom_vector_maker_61(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r,tag1,tag2);
        
        cout<<"qvec size = "<<qvec.size()<<endl;
        //cout<<"tag1 = "<<tag1<<'\t'<<"tag2 = "<<tag2<<endl;
        //for(int i=0;i<qvec.size();++i) cout<<i<<'\t'<<qvec[i]<<endl;

        //This part writes the qvec in to a file //
        string qvecfile =   //"qvec_a_" + to_string(a)
                          //+ "_eps_" + to_string(eps)
                           "qvec_scount_" + to_string(scount)
                          + ".dat";


        fout.open(qvecfile.c_str());
        for(int i=0;i<qvec.size();++i)
        {
            fout<<i<<'\t'<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
        }
        fout.close();
        cout<<"qvecfile : "<<qvecfile<<endl;

        //This part writes the kernel with k=q
        
        int isize = (int) points;
        int jsize = (int) points;

        string kernelfile =  // "kernel_" //a_" + to_string(a)
                            //+ "_eps_" + to_string(eps)
                             "kernel_scount_" + to_string(scount)
                            + ".dat";

        fout.open(kernelfile.c_str());
        for(int i=0;i<isize;++i)
        {
            for(int j=0;j<jsize;++j)
            {
                comp pR = (comp)pRinitial + (comp)i*delpR;
                comp pI = (comp)pIinitial + (comp)j*delpI; 
                comp pcomp = (comp) pR + ii*pI;

                comp sigp = sigma_p(s,pcomp,m);

                double pi = acos(-1.0);

                comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            
                fout<<setprecision(16)<<real(pR)<<'\t'
                    <<setprecision(16)<<real(pI)<<'\t'
                    <<setprecision(16)<<real(ope)<<'\t'
                    <<setprecision(16)<<imag(ope)<<endl;
            
            }
        
        }
        fout.close();
        cout<<"kernelfile : "<<kernelfile<<endl;

        //Here we print the relevant thresholds for 
        //for the kernel and trace the pcut for the OPE

        string thresholdfile = //"threshold_a_" + to_string(a)
                            //+ "_eps_" + to_string(eps)
                            + "threshold_scount_" + to_string(scount)
                            + ".dat";

        
        fout.open(thresholdfile.c_str());
        fout<<real(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<'\t'
            <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<eps<<'\t'
            <<a<<endl;

        fout.close();
        cout<<"thresholdfile : "<<thresholdfile<<endl;

        delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);

        string opecutfile = //"opecuttracer_a_" + to_string(a) 
                          //+ "_eps_" + to_string(eps)
                           "opecuttracer_scount_" + to_string(scount)
                          + ".dat";

        fout.open(opecutfile.c_str());  

        vector<double> xvec;
        vector<comp> pcutpvec;
        vector<comp> pcutmvec;

        for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        {
            xvec.push_back(x);
            comp pcutp = pcut_plus(x,s,q,m,eps);
            comp pcutm = pcut_minus(x,s,q,m,eps);
            pcutpvec.push_back(pcutp);
            pcutmvec.push_back(pcutm);

        }
        pcut_fixer(pcutpvec,pcutmvec);

        //for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        for(int i=0;i<pcutpvec.size();++i)
        {
            double x = xvec[i];
            comp pcutp = pcutpvec[i];//pcut_plus(x,s,q,m,eps);
            comp pcutm = pcutmvec[i];//pcut_minus(x,s,q,m,eps);

            fout<<x<<'\t'
                <<real(pcutp)<<'\t'
                <<imag(pcutp)<<'\t'
                <<real(pcutm)<<'\t'
                <<imag(pcutm)<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<abs(-q-pcutp)<<'\t'
                <<abs(q-pcutm)<<endl;
        }
        fout.close();

        cout<<"opecutfile : "<<opecutfile<<endl;

        string pcutfile = //"pcut_a_" + to_string(a) 
                          //+ "_eps_" + to_string(eps)
                           "pcut_scount_" + to_string(scount)
                          + ".dat";

        fout.open(pcutfile.c_str());  

        pcount = 0;

        delxreal = abs(xrealinitial-xrealfinal)/((double)points);
        for(int k=0;k<qvec.size();++k)
        {
                             
            for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
            {
                comp pcutp = pcut_plus(x,s,qvec[k],m,eps);
                comp pcutm = pcut_minus(x,s,qvec[k],m,eps);
                fout<<x<<'\t'
                    <<real(pcutp)<<'\t'
                    <<imag(pcutp)<<'\t'
                    <<real(pcutm)<<'\t'
                    <<imag(pcutm)<<'\t'
                    <<real(q)<<'\t'
                    <<imag(q)<<'\t'
                    <<abs(-q-pcutp)<<'\t'
                    <<abs(q-pcutm)<<endl;
            }
            
            //cout<<"qvec["<<k<<"] = "<<qvec[k]<<'\t'<<"pcount = "<<pcount<<endl;
            //cout<<"pcutfile : "<<pcutfile<<endl;

            //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;
            //cout<<"file:"<<filename1<<endl;
            pcount = pcount + 1;
            //cout<<"loop completed = "<<pcount<<endl;
        }
        fout.close();
        cout<<"pcutfile : "<<pcutfile<<endl;
        cout<<endl;
        scount = scount + 1;
    //}
}



void kernel_pcut_x_glockletest_vs_s3_allqvec_omp()
{
    //cout<<"here"<<endl;
    
    double a = 16.0;
    double m = 1.0;
    double sreal = 3.0;//3.0;
    double simag = -0.05;
    comp ii = {0.0,1.0};
    double eps = 0.0;//1.0e-5;
    double eps1 = 0.0;

    comp sigb = sigmab(a,m);
    


    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int delspoints = 300;
    int qvecpoints = 50;
    double qvec_r = 0.01;

    double delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double sinitial = 8.72;//(double)real(phibthreshold(a,m));//8.72;
    double sfinal = 9.0;//(double)real(phibthreshold(a,m));//9.0;//(double)real(phibthreshold(a,m));//8.0;
    double dels = abs(sinitial-sfinal)/(double)delspoints;

    int scount = 0;
    int icount = 0;

    //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)

    omp_set_num_threads(6);
    ofstream fout[(int)delspoints+1];

    #pragma omp parallel for 
    for(icount=0;icount<(int)delspoints+1;++icount)
    {
        
        //if(scount!=0) break;
        double s3 = sinitial + icount*dels;
        scount = icount;
        int pcount = 0;
        comp s = s3 + ii*simag;
        comp q = pmom(s,sigb-ii*eps1,m);
        comp kmax = pmom(s,0.0,m);


        

        vector<comp> qvec;
        int tag1,tag2;
        //mom_vector_maker_4(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        mom_vector_maker_6(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r,tag1,tag2);
        
        
        //for(int i=0;i<qvec.size();++i) cout<<i<<'\t'<<qvec[i]<<endl;

        //This part writes the qvec in to a file //
        string qvecfile =   "qvec_a_" + to_string(a)
                          //+ "_eps_" + to_string(eps)
                          + "_scount_" + to_string(scount)
                          + ".dat";


        fout[icount].open(qvecfile.c_str());
        for(int i=0;i<qvec.size();++i)
        {
            fout[icount]<<i<<'\t'<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
        }
        fout[icount].close();
        

        //This part writes the kernel with k=q
        double pRinitial = -0.2;
        double pRfinal = 0.2;
        double delpR = abs(pRinitial-pRfinal)/points;
        double pIinitial = -0.5;
        double pIfinal = 0.5;
        double delpI = abs(pIinitial - pIfinal)/points;
        int isize = (int) points;
        int jsize = (int) points;

        string kernelfile =   "kernel_a_" + to_string(a)
                            //+ "_eps_" + to_string(eps)
                            + "_scount_" + to_string(scount)
                            + ".dat";

        fout[icount].open(kernelfile.c_str());
        for(int i=0;i<isize;++i)
        {
            for(int j=0;j<jsize;++j)
            {
                comp pR = (comp)pRinitial + (comp)i*delpR;
                comp pI = (comp)pIinitial + (comp)j*delpI; 
                comp pcomp = (comp) pR + ii*pI;

                comp sigp = sigma_p(s,pcomp,m);

                double pi = acos(-1.0);

                comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            
                fout[icount]<<setprecision(16)<<real(pR)<<'\t'
                    <<setprecision(16)<<real(pI)<<'\t'
                    <<setprecision(16)<<real(ope)<<'\t'
                    <<setprecision(16)<<imag(ope)<<endl;
            
            }
        
        }
        fout[icount].close();
        

        //Here we print the relevant thresholds for 
        //for the kernel and trace the pcut for the OPE

        string thresholdfile = "threshold_a_" + to_string(a)
                            //+ "_eps_" + to_string(eps)
                            + "_scount_" + to_string(scount)
                            + ".dat";

        
        fout[icount].open(thresholdfile.c_str());
        fout[icount]<<real(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<'\t'
            <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<eps<<endl;

        fout[icount].close();
        

        delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);

        string opecutfile = "opecuttracer_a_" + to_string(a) 
                          //+ "_eps_" + to_string(eps)
                          + "_scount_" + to_string(scount)
                          + ".dat";

        fout[icount].open(opecutfile.c_str());  

        vector<double> xvec;
        vector<comp> pcutpvec;
        vector<comp> pcutmvec;

        for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        {
            xvec.push_back(x);
            comp pcutp = pcut_plus(x,s,q,m,eps);
            comp pcutm = pcut_minus(x,s,q,m,eps);
            pcutpvec.push_back(pcutp);
            pcutmvec.push_back(pcutm);

        }
        pcut_fixer(pcutpvec,pcutmvec);

        //for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        for(int i=0;i<pcutpvec.size();++i)
        {
            double x = xvec[i];
            comp pcutp = pcutpvec[i];//pcut_plus(x,s,q,m,eps);
            comp pcutm = pcutmvec[i];//pcut_minus(x,s,q,m,eps);

            fout[icount]<<x<<'\t'
                <<real(pcutp)<<'\t'
                <<imag(pcutp)<<'\t'
                <<real(pcutm)<<'\t'
                <<imag(pcutm)<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<abs(-q-pcutp)<<'\t'
                <<abs(q-pcutm)<<endl;
        }
        fout[icount].close();

        

        string pcutfile = "pcut_a_" + to_string(a) 
                          //+ "_eps_" + to_string(eps)
                          + "_scount_" + to_string(scount)
                          + ".dat";

        fout[icount].open(pcutfile.c_str());  

        pcount = 0;

        delxreal = abs(xrealinitial-xrealfinal)/((double)points);
        for(int k=0;k<qvec.size();++k)
        {
                             
            for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
            {
                comp pcutp = pcut_plus(x,s,qvec[k],m,eps);
                comp pcutm = pcut_minus(x,s,qvec[k],m,eps);
                fout[icount]<<x<<'\t'
                    <<real(pcutp)<<'\t'
                    <<imag(pcutp)<<'\t'
                    <<real(pcutm)<<'\t'
                    <<imag(pcutm)<<'\t'
                    <<real(q)<<'\t'
                    <<imag(q)<<'\t'
                    <<abs(-q-pcutp)<<'\t'
                    <<abs(q-pcutm)<<endl;
            }
            
            //cout<<"qvec["<<k<<"] = "<<qvec[k]<<'\t'<<"pcount = "<<pcount<<endl;
            //cout<<"pcutfile : "<<pcutfile<<endl;

            //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;
            //cout<<"file:"<<filename1<<endl;
            //pcount = pcount + 1;
            //cout<<"loop completed = "<<pcount<<endl;
        }
        fout[icount].close();

        if(omp_get_thread_num()==0)
        {
            std::cout<<"run = "<<scount+1<<endl;
    
            std::cout<<"s = "<<s<<'\t'<<"dels = "<<dels<<endl;
            std::cout<<"qvec size = "<<qvec.size()<<endl;
            //std::cout<<"tag1 = "<<tag1<<'\t'<<"tag2 = "<<tag2<<endl;
            std::cout<<"kernelfile : "<<kernelfile<<endl;
            std::cout<<"qvecfile : "<<qvecfile<<endl;
            std::cout<<"opecutfile : "<<opecutfile<<endl;
            std::cout<<"thresholdfile : "<<thresholdfile<<endl;
            std::cout<<"pcutfile : "<<pcutfile<<endl;
            std::cout<<endl;
            //scount = scount + 1;
        }
    }
}

void kernel_pcut_x_glockletest_vs_s3imag_allqvec()
{
    //cout<<"here"<<endl;
    ofstream fout;
    double a = 16.0;
    double m = 1.0;
    double sreal = 8.72;//3.0;
    double simag = 0.0;
    comp ii = {0.0,1.0};
    double eps = 0.0;//1.0e-5;
    double eps1 = 0.0;

    comp sigb = sigmab(a,m);
    


    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int delspoints = 100;
    int qvecpoints = 10;
    double qvec_r = 0.01;

    double delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double sinitial = -0.1;//8.72;//(double)real(phibthreshold(a,m));//8.72;
    double sfinal = 0.1;//(double)real(phibthreshold(a,m));//9.0;//(double)real(phibthreshold(a,m));//8.0;
    double dels = abs(sinitial-sfinal)/(double)delspoints;

    int scount = 0;
    for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    {
        int pcount = 0;
        comp s = sreal + ii*s3;// + ii*simag;
        comp q = pmom(s,sigb-ii*eps1,m);
        comp kmax = pmom(s,0.0,m);
        cout<<"run = "<<scount+1<<endl;
    
        cout<<"s = "<<s<<'\t'<<"dels = "<<dels<<endl;

        vector<comp> qvec;
        
        //mom_vector_maker_4(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        int tag1,tag2;
        //mom_vector_maker_4(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        mom_vector_maker_6(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r,tag1,tag2);
        
        cout<<"qvec size = "<<qvec.size()<<endl;

        //This part writes the qvec in to a file //
        string qvecfile =   "qvec_a_" + to_string(a)
                          //+ "_eps_" + to_string(eps)
                          + "_scount_" + to_string(scount)
                          + ".dat";


        fout.open(qvecfile.c_str());
        for(int i=0;i<qvec.size();++i)
        {
            fout<<i<<'\t'<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
        }
        fout.close();
        cout<<"qvecfile : "<<qvecfile<<endl;

        //This part writes the kernel with k=q
        double pRinitial = -0.2;
        double pRfinal = 0.2;
        double delpR = abs(pRinitial-pRfinal)/points;
        double pIinitial = -0.5;
        double pIfinal = 0.5;
        double delpI = abs(pIinitial - pIfinal)/points;
        int isize = (int) points;
        int jsize = (int) points;

        string kernelfile =   "kernel_a_" + to_string(a)
                            //+ "_eps_" + to_string(eps)
                            + "_scount_" + to_string(scount)
                            + ".dat";

        fout.open(kernelfile.c_str());
        for(int i=0;i<isize;++i)
        {
            for(int j=0;j<jsize;++j)
            {
                comp pR = (comp)pRinitial + (comp)i*delpR;
                comp pI = (comp)pIinitial + (comp)j*delpI; 
                comp pcomp = (comp) pR + ii*pI;

                comp sigp = sigma_p(s,pcomp,m);

                double pi = acos(-1.0);

                comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            
                fout<<setprecision(16)<<real(pR)<<'\t'
                    <<setprecision(16)<<real(pI)<<'\t'
                    <<setprecision(16)<<real(ope)<<'\t'
                    <<setprecision(16)<<imag(ope)<<endl;
            
            }
        
        }
        fout.close();
        cout<<"kernelfile : "<<kernelfile<<endl;

        //Here we print the relevant thresholds for 
        //for the kernel and trace the pcut for the OPE

        string thresholdfile = "threshold_a_" + to_string(a)
                            //+ "_eps_" + to_string(eps)
                            + "_scount_" + to_string(scount)
                            + ".dat";

        
        fout.open(thresholdfile.c_str());
        fout<<real(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<'\t'
            <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<eps<<endl;

        fout.close();
        cout<<"thresholdfile : "<<thresholdfile<<endl;

        delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);

        string opecutfile = "opecuttracer_a_" + to_string(a) 
                          //+ "_eps_" + to_string(eps)
                          + "_scount_" + to_string(scount)
                          + ".dat";

        fout.open(opecutfile.c_str());  

        vector<double> xvec;
        vector<comp> pcutpvec;
        vector<comp> pcutmvec;

        for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        {
            xvec.push_back(x);
            comp pcutp = pcut_plus(x,s,q,m,eps);
            comp pcutm = pcut_minus(x,s,q,m,eps);
            pcutpvec.push_back(pcutp);
            pcutmvec.push_back(pcutm);

        }
        pcut_fixer(pcutpvec,pcutmvec);

        //for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        for(int i=0;i<pcutpvec.size();++i)
        {
            double x = xvec[i];
            comp pcutp = pcutpvec[i];//pcut_plus(x,s,q,m,eps);
            comp pcutm = pcutmvec[i];//pcut_minus(x,s,q,m,eps);

            fout<<x<<'\t'
                <<real(pcutp)<<'\t'
                <<imag(pcutp)<<'\t'
                <<real(pcutm)<<'\t'
                <<imag(pcutm)<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<abs(-q-pcutp)<<'\t'
                <<abs(q-pcutm)<<endl;
        }
        fout.close();

        cout<<"opecutfile : "<<opecutfile<<endl;

        string pcutfile = "pcut_a_" + to_string(a) 
                          //+ "_eps_" + to_string(eps)
                          + "_scount_" + to_string(scount)
                          + ".dat";

        fout.open(pcutfile.c_str());  

        pcount = 0;

        delxreal = abs(xrealinitial-xrealfinal)/((double)points);
        for(int k=0;k<qvec.size();++k)
        {
                             
            for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
            {
                comp pcutp = pcut_plus(x,s,qvec[k],m,eps);
                comp pcutm = pcut_minus(x,s,qvec[k],m,eps);
                fout<<x<<'\t'
                    <<real(pcutp)<<'\t'
                    <<imag(pcutp)<<'\t'
                    <<real(pcutm)<<'\t'
                    <<imag(pcutm)<<'\t'
                    <<real(q)<<'\t'
                    <<imag(q)<<'\t'
                    <<abs(-q-pcutp)<<'\t'
                    <<abs(q-pcutm)<<endl;
            }
            
            //cout<<"qvec["<<k<<"] = "<<qvec[k]<<'\t'<<"pcount = "<<pcount<<endl;
            //cout<<"pcutfile : "<<pcutfile<<endl;

            //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;
            //cout<<"file:"<<filename1<<endl;
            pcount = pcount + 1;
            //cout<<"loop completed = "<<pcount<<endl;
        }
        fout.close();
        cout<<"pcutfile : "<<pcutfile<<endl;
        cout<<endl;
        scount = scount + 1;
    }
}


void kernel_pcut_x_glockletest_vs_eps_allqvec_temp()
{
    ofstream fout;
    double a = 16.0;
    double m = 1.0;
    double sreal = 3.0;//3.0;
    double simag = 0.0;
    comp ii = {0.0,1.0};
    double eps = 0.0;
    double eps1 = 0.0;

    comp sigb = sigmab(a,m);
    


    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int delspoints = 100;
    int qvecpoints = 4000;
    double qvec_r = 0.01;

    double delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double sinitial = 8.72;
    double sfinal = (double)real(phibthreshold(a,m));//8.0;
    double dels = abs(sinitial-sfinal)/(double)delspoints;

    int epscount = 0;
    double s3 = 8.72;


    //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    double epsinitial = 0.001;
    double epsfinal = 0.0;
    double epspoints = 300.0;
    double deleps = abs(epsinitial - epsfinal)/epspoints;

    for(double eps=epsinitial;eps>=epsfinal;eps=eps-deleps)
    {
        int pcount = 0;
        comp s = s3 + ii*simag;
        comp q = pmom(s,sigb-ii*eps1,m);
        comp kmax = pmom(s,0.0,m);

    
        cout<<"s = "<<s<<'\t'<<"dels = "<<dels<<endl;

        vector<comp> qvec;
        mom_vector_maker_4(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        
        cout<<"qvec size = "<<qvec.size()<<endl;

        //This part writes the qvec in to a file //
        string qvecfile =   "qvec_a_" + to_string(a)
                          + "_epscount_" + to_string(epscount)
                          + ".dat";


        fout.open(qvecfile.c_str());
        for(int i=0;i<qvec.size();++i)
        {
            fout<<i<<'\t'<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
        }
        fout.close();
        cout<<"qvecfile : "<<qvecfile<<endl;

        //This part writes the kernel with k=q
        double pRinitial = -0.2;
        double pRfinal = 0.2;
        double delpR = abs(pRinitial-pRfinal)/points;
        double pIinitial = -0.5;
        double pIfinal = 0.5;
        double delpI = abs(pIinitial - pIfinal)/points;
        int isize = (int) points;
        int jsize = (int) points;

        string kernelfile =   "kernel_a_" + to_string(a)
                            + "_epscount_" + to_string(epscount)
                            + ".dat";

        fout.open(kernelfile.c_str());
        for(int i=0;i<isize;++i)
        {
            for(int j=0;j<jsize;++j)
            {
                comp pR = (comp)pRinitial + (comp)i*delpR;
                comp pI = (comp)pIinitial + (comp)j*delpI; 
                comp pcomp = (comp) pR + ii*pI;

                comp sigp = sigma_p(s,pcomp,m);

                double pi = acos(-1.0);

                comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            
                fout<<setprecision(16)<<real(pR)<<'\t'
                    <<setprecision(16)<<real(pI)<<'\t'
                    <<setprecision(16)<<real(ope)<<'\t'
                    <<setprecision(16)<<imag(ope)<<endl;
            
            }
        
        }
        fout.close();
        cout<<"kernelfile : "<<kernelfile<<endl;

        //Here we print the relevant thresholds for 
        //for the kernel and trace the pcut for the OPE

        string thresholdfile = "threshold_a_" + to_string(a)
                            + "_epscount_" + to_string(epscount)
                            + ".dat";

        
        fout.open(thresholdfile.c_str());
        fout<<real(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<'\t'
            <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<eps<<endl;

        fout.close();
        cout<<"thresholdfile : "<<thresholdfile<<endl;

        delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);

        string opecutfile = "opecuttracer_a_" + to_string(a) 
                          + "_epscount_" + to_string(epscount)
                          + ".dat";

        fout.open(opecutfile.c_str());                   
        for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        {
            comp pcutp = pcut_plus(x,s,q,m,eps);
            comp pcutm = pcut_minus(x,s,q,m,eps);
            fout<<x<<'\t'
                <<real(pcutp)<<'\t'
                <<imag(pcutp)<<'\t'
                <<real(pcutm)<<'\t'
                <<imag(pcutm)<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<abs(-q-pcutp)<<'\t'
                <<abs(q-pcutm)<<endl;
        }
        fout.close();

        cout<<"opecutfile : "<<opecutfile<<endl;

        string pcutfile = "pcut_a_" + to_string(a) 
                          + "_epscount_" + to_string(epscount)
                          + ".dat";

        fout.open(pcutfile.c_str());  

        pcount = 0;

        delxreal = abs(xrealinitial-xrealfinal)/((double)points);
        for(int k=0;k<qvec.size();++k)
        {
                             
            for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
            {
                comp pcutp = pcut_plus(x,s,qvec[k],m,eps);
                comp pcutm = pcut_minus(x,s,qvec[k],m,eps);
                fout<<x<<'\t'
                    <<real(pcutp)<<'\t'
                    <<imag(pcutp)<<'\t'
                    <<real(pcutm)<<'\t'
                    <<imag(pcutm)<<'\t'
                    <<real(q)<<'\t'
                    <<imag(q)<<'\t'
                    <<abs(-q-pcutp)<<'\t'
                    <<abs(q-pcutm)<<endl;
            }
            
            //cout<<"qvec["<<k<<"] = "<<qvec[k]<<'\t'<<"pcount = "<<pcount<<endl;
            //cout<<"pcutfile : "<<pcutfile<<endl;

            //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;
            //cout<<"file:"<<filename1<<endl;
            pcount = pcount + 1;
            //cout<<"loop completed = "<<pcount<<endl;
        }
        fout.close();
        cout<<"pcutfile : "<<pcutfile<<endl;
        cout<<endl;
        epscount = epscount + 1;
    }
}


void kernel_pcut_x_glockletest_vs_eps_allqvec()
{
    //cout<<"here"<<endl;
    ofstream fout;
    double a = 16.0;
    double m = 1.0;
    double sreal = 3.0;//3.0;
    double s3 = 8.72;
    double simag = 0.01;//-0.05;
    comp ii = {0.0,1.0};
    double eps = 0.0;//1.0e-5;
    double eps1 = 0.0;

    comp sigb = sigmab(a,m);
    


    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int delspoints = 300;
    int qvecpoints = 50;
    double qvec_r = 0.01;

    double delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double sinitial = 8.72;//(double)real(phibthreshold(a,m));//8.72;
    double sfinal = (double)real(phibthreshold(a,m));//9.0;//(double)real(phibthreshold(a,m));//8.0;
    double dels = abs(sinitial-sfinal)/(double)delspoints;

    int delepspoints = 300;

    double epsinitial = 0.1;
    double epsfinal = 0.0;
    double deleps = abs(epsinitial - epsfinal)/delepspoints;

    int scount = 0;

    

    for(int i=0;i<delepspoints;++i)
    //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    {
        eps = epsinitial - i*deleps;
        //double s3 = sfinal - i*dels;
        //if(scount!=0) break;
        int pcount = 0;
        comp s = s3 + ii*simag;
        comp q = pmom(s,sigb-ii*eps1,m);
        comp kmax = pmom(s,0.0,m);
        cout<<"run = "<<scount+1<<endl;
    
        cout<<"s = "<<s<<'\t'<<"dels = "<<dels<<endl;
        cout<<"eps = "<<eps<<endl;

        vector<comp> qvec;
        int tag1,tag2;
        mom_vector_maker_4(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        //mom_vector_maker_61(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r,tag1,tag2);
        
        cout<<"qvec size = "<<qvec.size()<<endl;
        //cout<<"tag1 = "<<tag1<<'\t'<<"tag2 = "<<tag2<<endl;
        //for(int i=0;i<qvec.size();++i) cout<<i<<'\t'<<qvec[i]<<endl;

        //This part writes the qvec in to a file //
        string qvecfile =   "qvec_a_" + to_string(a)
                          //+ "_eps_" + to_string(eps)
                          + "_scount_" + to_string(scount)
                          + ".dat";


        fout.open(qvecfile.c_str());
        for(int i=0;i<qvec.size();++i)
        {
            fout<<i<<'\t'<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
        }
        fout.close();
        cout<<"qvecfile : "<<qvecfile<<endl;

        //This part writes the kernel with k=q
        double pRinitial = -0.2;
        double pRfinal = 0.2;
        double delpR = abs(pRinitial-pRfinal)/points;
        double pIinitial = -0.5;
        double pIfinal = 0.5;
        double delpI = abs(pIinitial - pIfinal)/points;
        int isize = (int) points;
        int jsize = (int) points;

        string kernelfile =   "kernel_a_" + to_string(a)
                            //+ "_eps_" + to_string(eps)
                            + "_scount_" + to_string(scount)
                            + ".dat";

        fout.open(kernelfile.c_str());
        for(int i=0;i<isize;++i)
        {
            for(int j=0;j<jsize;++j)
            {
                comp pR = (comp)pRinitial + (comp)i*delpR;
                comp pI = (comp)pIinitial + (comp)j*delpI; 
                comp pcomp = (comp) pR + ii*pI;

                comp sigp = sigma_p(s,pcomp,m);

                double pi = acos(-1.0);

                comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            
                fout<<setprecision(16)<<real(pR)<<'\t'
                    <<setprecision(16)<<real(pI)<<'\t'
                    <<setprecision(16)<<real(ope)<<'\t'
                    <<setprecision(16)<<imag(ope)<<endl;
            
            }
        
        }
        fout.close();
        cout<<"kernelfile : "<<kernelfile<<endl;

        //Here we print the relevant thresholds for 
        //for the kernel and trace the pcut for the OPE

        string thresholdfile = "threshold_a_" + to_string(a)
                            //+ "_eps_" + to_string(eps)
                            + "_scount_" + to_string(scount)
                            + ".dat";

        
        fout.open(thresholdfile.c_str());
        fout<<real(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<'\t'
            <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<eps<<endl;

        fout.close();
        cout<<"thresholdfile : "<<thresholdfile<<endl;

        delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);

        string opecutfile = "opecuttracer_a_" + to_string(a) 
                          //+ "_eps_" + to_string(eps)
                          + "_scount_" + to_string(scount)
                          + ".dat";

        fout.open(opecutfile.c_str());  

        vector<double> xvec;
        vector<comp> pcutpvec;
        vector<comp> pcutmvec;

        for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        {
            xvec.push_back(x);
            comp pcutp = pcut_plus(x,s,q,m,eps);
            comp pcutm = pcut_minus(x,s,q,m,eps);
            pcutpvec.push_back(pcutp);
            pcutmvec.push_back(pcutm);

        }
        pcut_fixer(pcutpvec,pcutmvec);

        //for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        for(int i=0;i<pcutpvec.size();++i)
        {
            double x = xvec[i];
            comp pcutp = pcutpvec[i];//pcut_plus(x,s,q,m,eps);
            comp pcutm = pcutmvec[i];//pcut_minus(x,s,q,m,eps);

            fout<<x<<'\t'
                <<real(pcutp)<<'\t'
                <<imag(pcutp)<<'\t'
                <<real(pcutm)<<'\t'
                <<imag(pcutm)<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<abs(-q-pcutp)<<'\t'
                <<abs(q-pcutm)<<endl;
        }
        fout.close();

        cout<<"opecutfile : "<<opecutfile<<endl;

        string pcutfile = "pcut_a_" + to_string(a) 
                          //+ "_eps_" + to_string(eps)
                          + "_scount_" + to_string(scount)
                          + ".dat";

        fout.open(pcutfile.c_str());  

        pcount = 0;

        delxreal = abs(xrealinitial-xrealfinal)/((double)points);
        for(int k=0;k<qvec.size();++k)
        {
                             
            for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
            {
                comp pcutp = pcut_plus(x,s,qvec[k],m,eps);
                comp pcutm = pcut_minus(x,s,qvec[k],m,eps);
                fout<<x<<'\t'
                    <<real(pcutp)<<'\t'
                    <<imag(pcutp)<<'\t'
                    <<real(pcutm)<<'\t'
                    <<imag(pcutm)<<'\t'
                    <<real(q)<<'\t'
                    <<imag(q)<<'\t'
                    <<abs(-q-pcutp)<<'\t'
                    <<abs(q-pcutm)<<endl;
            }
            
            //cout<<"qvec["<<k<<"] = "<<qvec[k]<<'\t'<<"pcount = "<<pcount<<endl;
            //cout<<"pcutfile : "<<pcutfile<<endl;

            //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;
            //cout<<"file:"<<filename1<<endl;
            pcount = pcount + 1;
            //cout<<"loop completed = "<<pcount<<endl;
        }
        fout.close();
        cout<<"pcutfile : "<<pcutfile<<endl;
        cout<<endl;
        scount = scount + 1;
    }
}


void kernel_pcut_x_glockletest_vs_s3_eps()
{
    ofstream fout;
    double a = 16.0;
    double m = 1.0;
    double sreal = 8.72;
    double simag = 0.10;
    comp ii = {0.0,1.0};
    double eps = 0.0;
    double eps1 = 0.0;

    comp sigb = sigmab(a,m);
    


    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int delspoints = 10;
    int qvecpoints = 100;
    int epspoints = 50;
    double qvec_r = 0.01;

    double delxreal = abs(xrealinitial-xrealfinal)/(double)points;
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double sinitial = 8.72;
    double sfinal = (double)real(phibthreshold(a,m));
    double simaginitial = 0.01;
    double simagfinal = 0.1;
    double dels = abs(sinitial-sfinal)/(double)delspoints;
    double delsimag = 0.01;//abs(simaginitial-simagfinal)/(double)delspoints;

    int scount = 0;
    for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    {
        int pcount = 0;
        comp s = s3 + ii*simag;
        
        comp kmax = pmom(s,0.0,m);

    
        cout<<"s = "<<s<<'\t'<<"dels = "<<dels<<endl;

        

        //This part writes the kernel with k=q
        double pRinitial = -0.2;
        double pRfinal = 0.2;
        double delpR = abs(pRinitial-pRfinal)/points;
        double pIinitial = -0.5;
        double pIfinal = 0.5;
        double delpI = abs(pIinitial - pIfinal)/points;
        int isize = (int) points;
        int jsize = (int) points;

        int epscount = 0;
        double epsinitial = 0.0;
        double epsfinal =0.1;
        double deleps = abs(epsinitial - epsfinal)/epspoints;

        for(epscount=0;epscount<epspoints;++epscount)
        {
            eps = epsinitial + epscount*deleps;
            double tempeps = eps - simag*eps*eps;
            eps = tempeps;
            //eps = eps;
            eps1 = eps;
            comp q = pmom(s,sigb-ii*eps1,m);
            string kernelfile =   "kernel_a_" + to_string(a)
                                + "_epscount_" + to_string(epscount)
                                + "_scount_" + to_string(scount)
                                + ".dat";

            fout.open(kernelfile.c_str());
            for(int i=0;i<isize;++i)
            {
                for(int j=0;j<jsize;++j)
                {
                    comp pR = (comp)pRinitial + (comp)i*delpR;
                    comp pI = (comp)pIinitial + (comp)j*delpI; 
                    comp pcomp = (comp) pR + ii*pI;

                    comp sigp = sigma_p(s,pcomp,m);

                    double pi = acos(-1.0);

                    comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            
                    fout<<setprecision(16)<<real(pR)<<'\t'
                        <<setprecision(16)<<real(pI)<<'\t'
                        <<setprecision(16)<<real(ope)<<'\t'
                        <<setprecision(16)<<imag(ope)<<endl;
            
                }
        
            }
            fout.close();
            cout<<"kernelfile : "<<kernelfile<<endl;

            //Here we print the relevant thresholds for 
            //for the kernel and trace the pcut for the OPE

            string thresholdfile = "threshold_a_" + to_string(a)
                                    + "_epscount_" + to_string(epscount)
                                    + "_scount_" + to_string(scount)
                                    + ".dat";

        
            fout.open(thresholdfile.c_str());
            fout<<real(pcut_plus_comp(1.0,s,q,m,eps))<<'\t'
                <<imag(pcut_plus_comp(1.0,s,q,m,eps))<<'\t'
                <<real(pcut_plus_comp(-1.0,s,q,m,eps))<<'\t'
                <<imag(pcut_plus_comp(-1.0,s,q,m,eps))<<'\t'
                <<real(pcut_minus_comp(1.0,s,q,m,eps))<<'\t'
                <<imag(pcut_minus_comp(1.0,s,q,m,eps))<<'\t'
                <<real(pcut_minus_comp(-1.0,s,q,m,eps))<<'\t'
                <<imag(pcut_minus_comp(-1.0,s,q,m,eps))<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
                <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
                <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
                <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
                <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
                <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
                <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
                <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
                <<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<eps<<endl;

            fout.close();
            cout<<"thresholdfile : "<<thresholdfile<<endl;

            string opecutfile = "opecuttracer_a_" + to_string(a) 
                                + "_epscount_" + to_string(epscount)
                                + "_scount_" + to_string(scount)
                                + ".dat";

            fout.open(opecutfile.c_str());                   
            for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
            {
                comp pcutp = pcut_plus_comp(x,s,q,m,eps);
                comp pcutm = pcut_minus_comp(x,s,q,m,eps);
                fout<<x<<'\t'
                    <<real(pcutp)<<'\t'
                    <<imag(pcutp)<<'\t'
                    <<real(pcutm)<<'\t'
                    <<imag(pcutm)<<'\t'
                    <<real(q)<<'\t'
                    <<imag(q)<<'\t'
                    <<abs(-q-pcutp)<<'\t'
                    <<abs(q-pcutm)<<endl;
            }
            fout.close();

            cout<<"opecutfile : "<<opecutfile<<endl;

        
            cout<<endl;
        }
        scount = scount + 1;
    }
}

void kernel_pcut_x_glockletest_vs_s3imag_eps_criticalepstest()
{
    ofstream fout;
    double a = 16.0;
    double m = 1.0;
    double sreal = 8.72;
    double simag = 0.10;
    comp ii = {0.0,1.0};
    double eps = 0.0;
    double eps1 = 0.0;

    comp sigb = sigmab(a,m);
    


    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int delspoints = 10;
    int qvecpoints = 100;
    int epspoints = 50;
    double qvec_r = 0.01;

    double delxreal = abs(xrealinitial-xrealfinal)/(double)points;
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double sinitial = 8.72;
    double sfinal = (double)real(phibthreshold(a,m));
    double simaginitial = 0.01;
    double simagfinal = 0.1;
    double dels = abs(sinitial-sfinal)/(double)delspoints;
    double delsimag = 0.01;//abs(simaginitial-simagfinal)/(double)delspoints;

    int scount = 0;
    for(double s3=simaginitial;s3<=simagfinal;s3=s3+delsimag)
    {
        int pcount = 0;
        comp s = sreal + ii*s3;
        
        comp kmax = pmom(s,0.0,m);

    
        cout<<"s = "<<s<<'\t'<<"dels = "<<dels<<endl;

        

        //This part writes the kernel with k=q
        double pRinitial = -0.2;
        double pRfinal = 0.2;
        double delpR = abs(pRinitial-pRfinal)/points;
        double pIinitial = -0.5;
        double pIfinal = 0.5;
        double delpI = abs(pIinitial - pIfinal)/points;
        int isize = (int) points;
        int jsize = (int) points;

        int epscount = 0;
        double epsinitial = 0.0;
        double epsfinal =0.002;
        double deleps = abs(epsinitial - epsfinal)/epspoints;

        for(epscount=0;epscount<epspoints;++epscount)
        {
            //if(epscount!=0) break;
            eps = epsinitial + epscount*deleps;
            //double tempeps = eps - simag*eps*eps;
            //eps = tempeps;
            //eps = eps;
            eps1 = 0.0;//eps;
            comp q = pmom(s,sigb-ii*eps1,m);
            string kernelfile =   "kernel_a_" + to_string(a)
                                + "_epscount_" + to_string(epscount)
                                + "_scount_" + to_string(scount)
                                + ".dat";

            fout.open(kernelfile.c_str());
            for(int i=0;i<isize;++i)
            {
                for(int j=0;j<jsize;++j)
                {
                    comp pR = (comp)pRinitial + (comp)i*delpR;
                    comp pI = (comp)pIinitial + (comp)j*delpI; 
                    comp pcomp = (comp) pR + ii*pI;

                    comp sigp = sigma_p(s,pcomp,m);

                    double pi = acos(-1.0);

                    comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            
                    fout<<setprecision(16)<<real(pR)<<'\t'
                        <<setprecision(16)<<real(pI)<<'\t'
                        <<setprecision(16)<<real(ope)<<'\t'
                        <<setprecision(16)<<imag(ope)<<endl;
            
                }
        
            }
            fout.close();
            cout<<"kernelfile : "<<kernelfile<<endl;

            //Here we print the relevant thresholds for 
            //for the kernel and trace the pcut for the OPE

            string thresholdfile = "threshold_a_" + to_string(a)
                                    + "_epscount_" + to_string(epscount)
                                    + "_scount_" + to_string(scount)
                                    + ".dat";

        
            fout.open(thresholdfile.c_str());
            fout<<real(pcut_plus_comp(1.0,s,q,m,eps))<<'\t'
                <<imag(pcut_plus_comp(1.0,s,q,m,eps))<<'\t'
                <<real(pcut_plus_comp(-1.0,s,q,m,eps))<<'\t'
                <<imag(pcut_plus_comp(-1.0,s,q,m,eps))<<'\t'
                <<real(pcut_minus_comp(1.0,s,q,m,eps))<<'\t'
                <<imag(pcut_minus_comp(1.0,s,q,m,eps))<<'\t'
                <<real(pcut_minus_comp(-1.0,s,q,m,eps))<<'\t'
                <<imag(pcut_minus_comp(-1.0,s,q,m,eps))<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
                <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
                <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
                <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
                <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
                <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
                <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
                <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
                <<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<eps<<endl;

            fout.close();
            cout<<"thresholdfile : "<<thresholdfile<<endl;

            string opecutfile = "opecuttracer_a_" + to_string(a) 
                                + "_epscount_" + to_string(epscount)
                                + "_scount_" + to_string(scount)
                                + ".dat";

            fout.open(opecutfile.c_str());                   
            for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
            {
                comp pcutp = pcut_plus_comp(x,s,q,m,eps);
                comp pcutm = pcut_minus_comp(x,s,q,m,eps);
                fout<<x<<'\t'
                    <<real(pcutp)<<'\t'
                    <<imag(pcutp)<<'\t'
                    <<real(pcutm)<<'\t'
                    <<imag(pcutm)<<'\t'
                    <<real(q)<<'\t'
                    <<imag(q)<<'\t'
                    <<abs(-q-pcutp)<<'\t'
                    <<abs(q-pcutm)<<endl;
            }
            fout.close();

            cout<<"opecutfile : "<<opecutfile<<endl;

        
            cout<<endl;
        }
        scount = scount + 1;
    }
}

void kernel_pcut_x_glockletest_vs_s3_eps_forsigk2m()
{
    ofstream fout;
    double a = 16.0;
    double m = 1.0;
    double sreal = 8.72;
    double simag = 0.0;
    comp ii = {0.0,1.0};
    double eps = 0.0;
    double eps1 = 0.0;

    comp sigb = 2.0*m*m;//sigmab(a,m);
    


    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int delspoints = 50;
    int qvecpoints = 100;
    int epspoints = 50;
    double qvec_r = 0.01;

    double delxreal = abs(xrealinitial-xrealfinal)/(double)points;
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double sinitial = 8.72;
    double sfinal = 9.0;//(double)real(phibthreshold(a,m));
    double dels = abs(sinitial-sfinal)/(double)delspoints;

    int scount = 0;
    for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    {
        int pcount = 0;
        comp s = s3 + ii*simag;
        
        //s = {8.78287,0.0};
        comp kmax = pmom(s,0.0,m);

    
        cout<<"s = "<<s<<'\t'<<"dels = "<<dels<<endl;

        

        //This part writes the kernel with k=q
        double pRinitial = -1;
        double pRfinal = 1;
        double delpR = abs(pRinitial-pRfinal)/points;
        double pIinitial = -1;
        double pIfinal = 1;
        double delpI = abs(pIinitial - pIfinal)/points;
        int isize = (int) points;
        int jsize = (int) points;

        int epscount = 0;
        double epsinitial = 0.0;
        double epsfinal =0.1;
        double deleps = abs(epsinitial - epsfinal)/epspoints;

        for(epscount=0;epscount<epspoints;++epscount)
        {
            if(epscount!=0) break;
            eps = epsinitial + epscount*deleps;
            //double tempeps = eps - simag*eps*eps;
            //eps = tempeps;
            //eps = eps;
            eps1 = 0.0;//eps;
            comp q = pmom(s,sigb-ii*eps1,m);
            string kernelfile =   "kernel_a_" + to_string(a)
                                + "_epscount_" + to_string(epscount)
                                + "_scount_" + to_string(scount)
                                + ".dat";

            fout.open(kernelfile.c_str());
            for(int i=0;i<isize;++i)
            {
                for(int j=0;j<jsize;++j)
                {
                    comp pR = (comp)pRinitial + (comp)i*delpR;
                    comp pI = (comp)pIinitial + (comp)j*delpI; 
                    comp pcomp = (comp) pR + ii*pI;

                    comp sigp = sigma_p(s,pcomp,m);

                    double pi = acos(-1.0);

                    comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            
                    fout<<setprecision(16)<<real(pR)<<'\t'
                        <<setprecision(16)<<real(pI)<<'\t'
                        <<setprecision(16)<<real(ope)<<'\t'
                        <<setprecision(16)<<imag(ope)<<endl;
            
                }
        
            }
            fout.close();
            cout<<"kernelfile : "<<kernelfile<<endl;

            //Here we print the relevant thresholds for 
            //for the kernel and trace the pcut for the OPE

            string thresholdfile = "threshold_a_" + to_string(a)
                                    + "_epscount_" + to_string(epscount)
                                    + "_scount_" + to_string(scount)
                                    + ".dat";

        
            fout.open(thresholdfile.c_str());
            fout<<real(pcut_plus_comp(1.0,s,q,m,eps))<<'\t'
                <<imag(pcut_plus_comp(1.0,s,q,m,eps))<<'\t'
                <<real(pcut_plus_comp(-1.0,s,q,m,eps))<<'\t'
                <<imag(pcut_plus_comp(-1.0,s,q,m,eps))<<'\t'
                <<real(pcut_minus_comp(1.0,s,q,m,eps))<<'\t'
                <<imag(pcut_minus_comp(1.0,s,q,m,eps))<<'\t'
                <<real(pcut_minus_comp(-1.0,s,q,m,eps))<<'\t'
                <<imag(pcut_minus_comp(-1.0,s,q,m,eps))<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
                <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
                <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
                <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
                <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
                <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
                <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
                <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
                <<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(phibthreshold(a,m))<<'\t'
                <<eps<<endl;

            fout.close();
            cout<<"thresholdfile : "<<thresholdfile<<endl;


            vector<double> xvec;
            vector<comp> pcutpvec;
            vector<comp> pcutmvec;

            for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
            {
                xvec.push_back(x);
                comp pcutp = pcut_plus(x,s,q,m,eps);
                comp pcutm = pcut_minus(x,s,q,m,eps);
                pcutpvec.push_back(pcutp);
                pcutmvec.push_back(pcutm);

            }
            pcut_fixer(pcutpvec,pcutmvec);

            string opecutfile = "opecuttracer_a_" + to_string(a) 
                                + "_epscount_" + to_string(epscount)
                                + "_scount_" + to_string(scount)
                                + ".dat";

            fout.open(opecutfile.c_str());                   
            //for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
            for(int somei=0;somei<pcutpvec.size();++somei)
            {
                double x = xvec[somei];
                comp pcutp = pcutpvec[somei];//pcut_plus_comp(x,s,q,m,eps);
                comp pcutm = pcutmvec[somei];//pcut_minus_comp(x,s,q,m,eps);
                fout<<x<<'\t'
                    <<real(pcutp)<<'\t'
                    <<imag(pcutp)<<'\t'
                    <<real(pcutm)<<'\t'
                    <<imag(pcutm)<<'\t'
                    <<real(q)<<'\t'
                    <<imag(q)<<'\t'
                    <<abs(-q-pcutp)<<'\t'
                    <<abs(q-pcutm)<<endl;
            }
            fout.close();

            cout<<"opecutfile : "<<opecutfile<<endl;
            cout<<"phib threshold = "<<(double)real(phibthreshold(a,m))<<endl;
        
            cout<<endl;
        }
        scount = scount + 1;
    }
}




void only_pcut_x_comp_qvec()
{
    double a = 16.0;
    double m = 1.0;
    double sreal = 8.72;
    double simag = -0.1;
    comp ii = {0.0,1.0};
    //comp s = sreal + ii*simag;
    double eps = 0.0;

    comp sigb = sigmab(a,m);
    //comp q = pmom(s,sigb-ii*eps,m);
    //comp kmax = pmom(s,0.0,m);


    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 2000;

    double delxreal = abs(xrealinitial-xrealfinal)/(double)points;
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double sinitial = 8.72;
    double sfinal = (double)real(phibthreshold(a,m));
    double dels = abs(sinitial-sfinal)/10.0;

    int count = 1;
    for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    {
        comp s = s3 + ii*simag;
        comp q = pmom(s,sigb-ii*eps,m);
        double jmultiplier = 1.0;
    
        ofstream fout;
        string filename1 = "pcutonly_test_a_" + to_string(a) 
                          + "_eps_" + to_string(eps)
                          + "_count_" + to_string(count)
                          + ".dat";

        fout.open(filename1.c_str());                   
        for(double j=0.0;j<=jmultiplier;j=j+0.001)
        {
            comp x = -1.0 - ii*j;
            comp pcutp = pcut_plus_comp(x,s,q,m,eps);
            comp pcutm = pcut_minus_comp(x,s,q,m,eps);
            fout<<real(x)<<'\t'<<imag(x)<<'\t'<<real(pcutp)<<'\t'<<imag(pcutp)<<'\t'<<real(pcutm)<<'\t'<<imag(pcutm)<<'\t'<<real(q)<<'\t'<<imag(q)<<endl;
        }
        for(double xval=xrealinitial;xval<=xrealfinal;xval=xval+delxreal)
        {
            comp x = - ii*jmultiplier + xval;
            comp pcutp = pcut_plus_comp(x,s,q,m,eps);
            comp pcutm = pcut_minus_comp(x,s,q,m,eps);
            fout<<real(x)<<'\t'<<imag(x)<<'\t'<<real(pcutp)<<'\t'<<imag(pcutp)<<'\t'<<real(pcutm)<<'\t'<<imag(pcutm)<<'\t'<<real(q)<<'\t'<<imag(q)<<endl;
        }
        for(double j=0.0;j<=jmultiplier;j=j+0.001)
        {
            comp x = 1.0 - ii*jmultiplier + ii*j;
            
            comp pcutp = pcut_plus_comp(x,s,q,m,eps);
            comp pcutm = pcut_minus_comp(x,s,q,m,eps);
            fout<<real(x)<<'\t'<<imag(x)<<'\t'<<real(pcutp)<<'\t'<<imag(pcutp)<<'\t'<<real(pcutm)<<'\t'<<imag(pcutm)<<'\t'<<real(q)<<'\t'<<imag(q)<<endl;
        }

        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;
        cout<<"file:"<<filename1<<endl;
        count = count + 1;
    }
}


void OPE_PW_poletesting_x()
{
    double a = 16.0;
    double m = 1.0;
    double sreal = 8.72;
    double simag = 1.0;
    comp ii = {0.0,1.0};
    comp s = sreal + ii*simag;
    double eps = 0.0;

    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps,m);
    comp kmax = pmom(s,0.0,m);


    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;

    double delxreal = abs(xrealinitial-xrealfinal)/(double)points;
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    double pfinal = real(kmax);
    double delp = abs(pinitial-pfinal)/50.0;

    int count = 1;
    for(double p=pinitial;p<=pfinal;p=p+delp)
    {
        ofstream fout;
        string filename1 = "OPE_PW_test_a_" + to_string(a) 
                          + "_sreal_" + to_string(sreal)
                          + "_simag_" + to_string(simag)
                          + "_eps_" + to_string(eps)
                          + "_count_" + to_string(count)
                          + ".dat";

        fout.open(filename1.c_str());                   
        for(int i=0;i<points;++i)
        {
            for(int j=0;j<points;++j)
            {
                double xvalreal = xrealinitial + i*delxreal;
                double xvalimag = ximaginitial + j*delximag;
                comp x = xvalreal + ii*xvalimag;

                comp ope = OPE_before_PW(s,p,q,x,m,eps);

                fout<<xvalreal<<'\t'<<xvalimag<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
            }
        }

        cout<<"p:"<<p<<'\t'<<"delp:"<<delp<<endl;
        cout<<"file:"<<filename1<<endl;
        count = count + 1;
    }
}

void test_lowest_point()
{
    double s = 8.72;
    double a = 16.0;
    double m = 1.0;
    double eps = 1.0e-8;
    double delta = abs((1.0e-2)*eps*eps_energy_factor_minus(s,sigmab(a,m),m));
        
    comp sigb = sigmab(a,m) ;

    vector<comp> sigvec;
    sigma_vector_maker_6(sigvec,s,0.0,sigmax(s,m),a,m,delta,5000,500,10);

    double temp = 0.0;
    comp tempsig = {0.0,0.0};

    for(int i=0;i<sigvec.size();++i)
    {
        comp ope = GS(s,sigb,sigvec[i],m,eps)*tau1(s,sigvec[i],m)*M2kfunc(a,sigvec[i],m,eps);

        if(imag(ope)<temp)
        { 
            temp = imag(ope);
            tempsig = sigvec[i];
        }

    }

    cout<<tempsig<<'\t'<<temp<<endl;
    cout<<sigpplus_witheps(s,sigb,m,eps)<<endl;
}


void test_sigma_plus_reflection()
{
    comp ii = {0.0,1.0};
    double s = 8.72;
    double a = 16.0;
    double m = 1.0;
    double eps = 1.0e-5;
    double delta = abs((1.0e-2)*eps*eps_energy_factor_minus(s,sigmab(a,m),m));
        
    comp sigb = sigmab(a,m) ;

    vector<comp> sigvec;
    sigma_vector_maker_6(sigvec,s,0.0,sigmax(s,m),a,m,delta,5000,500,10);

    double temp = 0.0;
    comp tempsig = {0.0,0.0};

    int points = 10000;

    double Resigmap = real(sigpplus_witheps(s,sigb,m,eps)) - 0.001;
    double branchpoint = imag(sigpplus_witheps(s,sigb,m,eps));
    //cout<<sigpplus_witheps(s,sigb,m,eps)<<endl;

    double delsigI = branchpoint/((double)points);

    for(int i=0;i<=2*points;++i)
    {
        comp sigp = Resigmap + ii*((double)i)*delsigI;
        comp ope = GS(s,sigb,sigp,m,eps)*tau1(s,sigp,m)*M2kfunc(a,sigp,m,eps);

        cout<<real(sigp)<<'\t'<<imag(sigp)<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;

    }

    //cout<<tempsig<<'\t'<<temp<<endl;
   
}


void testing_sigp_qfunc()
{
    comp ii = {0.0,1.0};

    double a = 16.0;
    double m = 1.0;
    comp s = 8.72 + ii*0.0;
    double eps = 0.01;
    comp sigb = sigmab(a,m);

    comp q = pmom(s,sigb,m);

    comp sigq = sigma_p(s,q,m);

    cout<<sigq<<endl;

    comp q1 = sqrt(sigb/4.0 - m*m);
    comp check = -1.0/a - ii*q1;

    cout<<"phibthreshold:"<<pow(sqrt(sigmab(16,1)) + 1.0,2.0)<<endl;

    cout<<q<<endl;
    cout<<q1<<endl;
    cout<<check<<endl;

    //for(double sreal=8.72;sreal<=real(pow(sqrt(sigmab(16,1)) + 1.0,2.0));sreal=sreal+0.01)
    //{
    //    cout<<sreal<<'\t'<<pmom(sreal,sigb,m)<<endl;
    //}

    /*for(double sigk=0.0;sigk<=4.0;sigk=sigk+0.01)
    {
        cout<<sigk<<'\t'<<pmom(s,sigk,m)<<endl;
    }

    cout<<(sqrt(s)-m)*(sqrt(s)-m)<<endl;
    */

    cout<<q_plus(s,a,m,eps)<<endl;
    cout<<q_minus(s,a,m,eps)<<endl;
    cout<<q_c1(s,a,m,eps)<<endl;
    cout<<q_c2(s,a,m,eps)<<endl;
    
}


void testing_momvec_maker_2()
{
    double a = 16.0;
    double m = 1.0;
    double phib = abs(real(phibthreshold(a,m)));
    double eps = 0.01;


    for(double s=8.72;s<=phib;s=s+0.01)
    {
        vector<comp> qvec;
        comp kmax = pmom(s,0.0,m);
        mom_vector_maker_2(qvec,s,0.0,kmax,a,m,eps,eps,1000,0.01); // there are two eps here, one for ope another for m2k
        cout<<s<<'\t'<<qvec.size()<<endl;
    }

}

void testing_pcuts()
{
    double a = 16.0;
    double m = 1.0;
    comp s = {8.905,0.1};
    comp p = {0.075,0.147};
    double eps = 0.01;

    for(double x=-1.0;x<=1.0;x=x+0.1)
    {
        cout<<"x:"<<x<<endl;
        comp pcutp = pcut_plus_comp(x,s,p,m,eps);
        comp pcutm = pcut_minus_comp(x,s,p,m,eps);
        cout<<x<<'\t'
            <<real(pcutp)<<'\t'
            <<imag(pcutp)<<'\t'
            <<real(pcutm)<<'\t'
            <<imag(pcutm)<<endl;
    }
}

void testing_pcuts_and_kernel( comp s )
{
    double a = 16.0;
    double m = 1.0;
    //comp s = {8.72,0.1};
    comp p = pmom(s,sigmab(a,m),m);
    double eps = 0.01;

    for(double x=-1.0;x<=1.0;x=x+0.1)
    {
        cout<<"x:"<<x<<endl;
        comp pcutp = pcut_plus_withprinter(x,s,p,m,eps,0);
        comp pcutm = pcut_minus_withprinter(x,s,p,m,eps,0);
        cout<<x<<'\t'
            <<real(pcutp)<<'\t'
            <<imag(pcutp)<<'\t'
            <<real(pcutm)<<'\t'
            <<imag(pcutm)<<endl;
    }

    ofstream fout;
    //double a = 16.0;
    //double m = 1.0;
    double sreal = 8.72;
    double simag = 0.1;
    comp ii = {0.0,1.0};
    //double eps = 0.01;
    double eps1 = 0.0;

    comp sigb = sigmab(a,m);
    


    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int delspoints = 100;
    int qvecpoints = 100;
    double qvec_r = 0.01;

    double delxreal = abs(xrealinitial-xrealfinal)/(double)points;
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double sinitial = 8.72;
    double sfinal = (double)real(phibthreshold(a,m));
    double dels = abs(sinitial-sfinal)/(double)delspoints;

    int scount = 0;
    
        int pcount = 0;
        //comp s = s3 + ii*simag;
        comp q = pmom(s,sigb-ii*eps1,m);
        comp kmax = pmom(s,0.0,m);

        
    
        cout<<"s = "<<s<<'\t'<<"dels = "<<dels<<endl;

        vector<comp> qvec;
        mom_vector_maker_3(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        


        //This part writes the qvec in to a file //
        string qvecfile =   "qvec_a_" + to_string(a)
                          + "_eps_" + to_string(eps)
                          + "_scount_" + to_string(scount)
                          + ".dat";


        fout.open(qvecfile.c_str());
        for(int i=0;i<qvec.size();++i)
        {
            fout<<i<<'\t'<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
        }
        fout.close();
        cout<<"qvecfile : "<<qvecfile<<endl;

        //This part writes the kernel with k=q
        double pRinitial = -0.2;
        double pRfinal = 0.2;
        double delpR = abs(pRinitial-pRfinal)/points;
        double pIinitial = -0.5;
        double pIfinal = 0.5;
        double delpI = abs(pIinitial - pIfinal)/points;
        int isize = (int) points;
        int jsize = (int) points;

        string kernelfile =   "kernel_a_" + to_string(a)
                            + "_eps_" + to_string(eps)
                            + "_scount_" + to_string(scount)
                            + ".dat";

        fout.open(kernelfile.c_str());
        for(int i=0;i<isize;++i)
        {
            for(int j=0;j<jsize;++j)
            {
                comp pR = (comp)pRinitial + (comp)i*delpR;
                comp pI = (comp)pIinitial + (comp)j*delpI; 
                comp pcomp = (comp) pR + ii*pI;

                comp sigp = sigma_p(s,pcomp,m);

                double pi = acos(-1.0);

                comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            
                fout<<setprecision(16)<<real(pR)<<'\t'
                    <<setprecision(16)<<real(pI)<<'\t'
                    <<setprecision(16)<<real(ope)<<'\t'
                    <<setprecision(16)<<imag(ope)<<endl;
            
            }
        
        }
        fout.close();
        cout<<"kernelfile : "<<kernelfile<<endl;

        //Here we print the relevant thresholds for 
        //for the kernel and trace the pcut for the OPE

        string thresholdfile = "threshold_a_" + to_string(a)
                            + "_eps_" + to_string(eps)
                            + "_scount_" + to_string(scount)
                            + ".dat";

        
        fout.open(thresholdfile.c_str());
        fout<<real(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<'\t'
            <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<dels<<endl;

        fout.close();
        cout<<"thresholdfile : "<<thresholdfile<<endl;

        string opecutfile = "opecuttracer_a_" + to_string(a) 
                          + "_eps_" + to_string(eps)
                          + "_scount_" + to_string(scount)
                          + ".dat";

        fout.open(opecutfile.c_str());                   
        for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        {
            comp pcutp = pcut_plus_withprinter(x,s,q,m,eps,1);
            comp pcutm = pcut_minus_withprinter(x,s,q,m,eps,1);
            fout<<x<<'\t'
                <<real(pcutp)<<'\t'
                <<imag(pcutp)<<'\t'
                <<real(pcutm)<<'\t'
                <<imag(pcutm)<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<abs(-q-pcutp)<<'\t'
                <<abs(q-pcutm)<<endl;
        }
        fout.close();

        cout<<"opecutfile : "<<opecutfile<<endl;

        string pcutfile = "pcut_a_" + to_string(a) 
                          + "_eps_" + to_string(eps)
                          + "_scount_" + to_string(scount)
                          + ".dat";

        fout.open(pcutfile.c_str());  

        for(int k=0;k<qvec.size();++k)
        {
                             
            for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
            {
                comp pcutp = pcut_plus_withprinter(x,s,qvec[k],m,eps,1);
                comp pcutm = pcut_minus_withprinter(x,s,qvec[k],m,eps,1);
                fout<<x<<'\t'
                    <<real(pcutp)<<'\t'
                    <<imag(pcutp)<<'\t'
                    <<real(pcutm)<<'\t'
                    <<imag(pcutm)<<'\t'
                    <<real(q)<<'\t'
                    <<imag(q)<<'\t'
                    <<abs(-q-pcutp)<<'\t'
                    <<abs(q-pcutm)<<endl;
            }
            
            //cout<<"qvec["<<k<<"] = "<<qvec[k]<<'\t'<<"pcount = "<<pcount<<endl;
            //cout<<"pcutfile : "<<pcutfile<<endl;

            //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;
            //cout<<"file:"<<filename1<<endl;
            pcount = pcount + 1;
        }
        fout.close();
        cout<<"pcutfile : "<<pcutfile<<endl;
        cout<<endl;
        scount = scount + 1;
    
}

double critical_epsilon(    comp s,
                            comp p,
                            comp k,
                            double m    )
{
    double sR = (double) real(s);
    double sI = (double) imag(s);

    double pR = (double) real(p);
    double pI = (double) imag(p);

    double kR = (double) real(k);
    double kI = (double) imag(k);

    double r_s = sqrt(abs(sR*sR + sI*sI));
    double R_s = sqrt(r_s);
    double theta_s = atan(sI/sR);

    double A_s = cos(theta_s/2.0);
    double B_s = sin(theta_s/2.0);

    double r_omega_p = sqrt(abs(pow(pR*pR - pI*pI + m*m,2.0) + 4.0*pR*pR*pI*pI));
    double R_omega_p = sqrt(r_omega_p);
    double theta_omega_p = atan((2.0*pR*pI)/(pR*pR - pI*pI + m*m));

    double A_omega_p = cos(theta_omega_p/2.0);
    double B_omega_p = sin(theta_omega_p/2.0);

    double r_omega_k = sqrt(abs(pow(kR*kR - kI*kI + m*m,2.0) + 4.0*kR*kR*kI*kI));
    double R_omega_k = sqrt(r_omega_k);
    double theta_omega_k = atan((2.0*kR*kI)/(kR*kR - kI*kI + m*m));

    double A_omega_k = cos(theta_omega_k/2.0);
    double B_omega_k = sin(theta_omega_k/2.0);

    double res =  2.0*R_s*R_omega_p*(A_s*B_omega_p + A_omega_p*B_s)
                + 2.0*R_s*R_omega_k*(A_s*B_omega_k + A_omega_k*B_s)
                - 2.0*R_omega_p*R_omega_k*(A_omega_p*B_omega_k + A_omega_k*B_omega_p)
                + 2.0*(pR*kI + pI*kR) - sI;
    
    return res;
    //cout<<res<<endl;
}

void test_critical_epsilon_vs_pcontour()
{
    ofstream fout;
    string filename = "crit_eps_test_pcontour_1.dat";
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    comp s = {8.902,0.1};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb,m);

    double points = 250;
    double pRinitial = -0.2;
    double pRfinal = 0.2;
    double delpR = abs(pRinitial-pRfinal)/points;
    double pIinitial = -0.5;
    double pIfinal = 0.5;
    double delpI = abs(pIinitial - pIfinal)/points;
    int isize = (int) points;
    int jsize = (int) points;

    cout<<"q:"<<q<<endl;
    cout<<"s:"<<s<<endl;
    fout.open(filename.c_str());
    for(int i=0;i<isize;++i)
    {
        for(int j=0;j<jsize;++j)
        {
            comp pR = (comp)pRinitial + (comp)i*delpR;
            comp pI = (comp)pIinitial + (comp)j*delpI; 
            comp pcomp = (comp) pR + ii*pI;

            comp sigp = sigma_p(s,pcomp,m);

            double pi = acos(-1.0);

            double res = critical_epsilon(s,pcomp,q,m);

            fout<<setprecision(16)<<real(pR)<<'\t'
                <<setprecision(16)<<real(pI)<<'\t'
                <<setprecision(16)<<res<<endl;
                
            
        }
        
    }
    fout.close();
}

void test_critical_epsilon_vs_pcontour_with0points()
{
    ofstream fout,fout1;
    string filename = "crit_eps_test_pcontour_1.dat";
    string filename2 = "eps_zeroing.dat";
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    comp s = {8.902,0.1};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb,m);

    double points = 250;
    double pRinitial = -0.2;
    double pRfinal = 0.2;
    double delpR = abs(pRinitial-pRfinal)/points;
    double pIinitial = -0.5;
    double pIfinal = 0.5;
    double delpI = abs(pIinitial - pIfinal)/points;
    int isize = (int) points;
    int jsize = (int) points;
    double eps = 0.01;

    cout<<"q:"<<q<<endl;
    cout<<"s:"<<s<<endl;
    fout.open(filename.c_str());
    fout1.open(filename2.c_str());
    for(int i=0;i<isize;++i)
    {
        for(int j=0;j<jsize;++j)
        {
            comp pR = (comp)pRinitial + (comp)i*delpR;
            comp pI = (comp)pIinitial + (comp)j*delpI; 
            comp pcomp = (comp) pR + ii*pI;

            comp sigp = sigma_p(s,pcomp,m);

            double pi = acos(-1.0);

            double res = critical_epsilon(s,pcomp,q,m);

            fout<<setprecision(16)<<real(pR)<<'\t'
                <<setprecision(16)<<real(pI)<<'\t'
                <<setprecision(16)<<res<<endl;

            if(abs(eps-res)<1.0e-5)
            {
                fout1<<setprecision(16)<<real(pR)<<'\t'
                     <<setprecision(16)<<real(pI)<<'\t'
                     <<abs(eps-res)<<endl;
            }

            
        }
        
    }
    fout.close();
    fout1.close();
}

void test_critical_epsilon_vs_s3()
{
    ofstream fout;
    string filename = "crit_eps_test_simag_0.1_p=0.dat";
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    //comp s = {8.902,0.1};
    

    double points = 250;
    double pRinitial = -0.2;
    double pRfinal = 0.2;
    double delpR = abs(pRinitial-pRfinal)/points;
    double pIinitial = -0.5;
    double pIfinal = 0.5;
    double delpI = abs(pIinitial - pIfinal)/points;
    int isize = (int) points;
    int jsize = (int) points;


    double sreal_initial = 8.72;
    double sreal_final = (double) real(phibthreshold(a,m));
    double simag = 0.1;

    double dels = abs(sreal_initial-sreal_final)/points;



    //cout<<"q:"<<q<<endl;
    //cout<<"s:"<<s<<endl;

    fout.open(filename.c_str());
    for(int i=0;i<points;++i)
    {
        comp s = sreal_initial + i*dels + ii*simag;
    
        comp sigb = sigmab(a,m);
        comp q = pmom(s,sigb,m*m);
        comp kmax = pmom(s,0.0,m*m);
        comp kmin = {0.0,0.0};
        //comp pR = (comp)pRinitial + (comp)i*delpR;
        //comp pI = (comp)pIinitial + (comp)j*delpI; 
        comp pcomp = kmin;//kmin;//(comp) pR + ii*pI;

        comp sigp = sigma_p(s,pcomp,m);

        double pi = acos(-1.0);

        double res = critical_epsilon(s,pcomp,q,m);

        fout<<setprecision(6)<<real(s)<<'\t'
            <<setprecision(6)<<imag(s)<<'\t'
            <<setprecision(6)<<res<<'\t'
            <<setprecision(6)<<real(pcomp)<<'\t'
            <<setprecision(6)<<imag(pcomp)<<endl;
                
    }        
    fout.close();
}

void test_critical_epsilon_pcut_numbers()
{
    double a = 16.0;
    double m = 1.0;
    double eps = 0.0056;
    comp s = {8.72,0.10};

    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb,m);
    comp qc1 = q_c1(s,a,m,eps);
    comp qc2 = q_c2(s,a,m,eps);

    comp ope = GS_pk(s,qc2,q,m,eps);
    double tmp1 = 1000.0;
    double tmp2 = 1000.0;
    double lowestx1 = 0.0;
    double lowestx2 = 0.0;
    double lowestepsval = -2.0;
    double tmp3 = 1000.0;

    for(double eps=0.0;eps<=0.01;eps=eps+0.0001)
    {
    for(double x=-1.0;x<=1.0;x=x+0.0001)
    {
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);
        //cout<<"pcutp:"<<pcut_plus(x,s,q,m,eps)<<'\t'<<"pcutm:"<<pcut_minus(x,s,q,m,eps)<<endl;
        comp opepcutp = GS_pk(s,pcutp,q,m,eps);
        comp opepcutm = GS_pk(s,pcutm,q,m,eps);
        double pcutpminusqc2 = abs(pcutp-qc2);
        double pcutmminusqc2 = abs(pcutm-qc2);
        //cout<<"opepp:"<<opepcutp<<'\t'<<"opepm:"<<opepcutm<<endl;
        //cout<<pcutpminusqc2<<'\t'<<pcutmminusqc2<<endl;
        if(pcutpminusqc2<tmp1)
        {
            tmp1 = pcutpminusqc2;
            lowestx1 = x;
            lowestepsval = eps;
        }
        if(pcutmminusqc2<tmp2)
        {
            tmp2 = pcutmminusqc2;
            lowestx2 = x;
            lowestepsval = eps;
        }
        
    }


    }
    cout<<"lowest val "<<tmp1<<'\t'<<tmp2<<endl;
    cout<<"lowest x "<<lowestx1<<'\t'<<lowestx2<<endl;
    cout<<"lowest eps "<<lowestepsval<<endl;
    double epscalc = critical_epsilon(s,qc1,q,m);

    //cout<<epscalc<<endl;
    //cout<<eps - epscalc<<endl;
    //cout<<ope<<endl;
    //cout<<qc1<<'\t'<<qc2<<endl;
    //cout<<pcut_plus(1,s,q,m,eps)<<endl;

    
}

void test_distance_between_pcuts()
{
    comp ii = {0.0,1.0};
    double sreal = 8.72;
    double simag = 0.0;
    comp s = sreal + ii*simag;
    double a = 16.0;
    double m = 1.0;
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb,m);
    double eps = 1.0e-7;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 50000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }

    double tempdiff = 1.0e+16;
    double tempxp = 0.0;
    double tempxm = 0.0;

    for(int i=0;i<pcutpvec.size();++i)
    {
        for(int j=0;j<pcutmvec.size();++j)
        {
            double diff = abs(pcutpvec[i] - pcutmvec[j]);

            //cout<<"xp = "<<xvec[i]<<'\t'<<"xm = "<<xvec[j]<<endl;

            if(diff<=tempdiff)
            {
                //cout<<"diff minimum then previous"<<endl;
                //cout<<"diff = "<<diff<<endl;
                //if(imag(pcutpvec[i])<0.0 && imag(pcutmvec[j])<0.0)
                {
                    tempdiff = diff;
                    tempxp = xvec[i];
                    tempxm = xvec[j];
                }
            }
            
        }
    }

    cout<<"diff is lowest for pcutp = "<<pcut_plus(tempxp,s,q,m,eps)<<" with x = "<<tempxp<<endl;
    cout<<" and pcutm = "<<pcut_minus(tempxm,s,q,m,eps)<<" with x = "<<tempxm<<endl;
    cout<<" the diff is = "<<abs(pcut_plus(tempxp,s,q,m,eps)-pcut_minus(tempxm,s,q,m,eps));
}

void test_with_Mphib()
{
    double a = 16.0;
    double m = 1.0;
    double s = 8.72;
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb,m);
    double eps = 1.0e-7;
    double eps_for_m2k = 0.0;

    comp kmin = {0.0,0.0};
    comp kmax = pmom(s,0.0,m);
    cout<<"kmax = "<<kmax<<endl;

    double qvec_r = 0.01;
    double points = 1000.0;

    vector<comp> qvec;
    mom_vector_maker_4(qvec,s,kmin,kmax,a,m,eps,eps_for_m2k,points,qvec_r);

    cout<<"OPE"<<endl;
    for(int i=0;i<qvec.size();++i)
    {
        cout<<"qvec = "<<qvec[i]<<'\t'<<"GS = "<<GS_pk(s,qvec[i],q,m,eps)<<endl;
    }


}

void kernel_pcut_criticalepsilontest_vs_s3(double eps)
{
    ofstream fout;
    double a = 16.0;
    double m = 1.0;
    double sreal = 8.70;//3.0;
    double simag = 0.1;
    comp ii = {0.0,1.0};
    //double eps = 1.0e-1;
    double eps1 = 0.0;

    comp sigb = sigmab(a,m);
    


    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int xpoints = 10000;
    int delspoints = 100;
    int qvecpoints = 5000;
    double qvec_r = 0.01;

    double delxreal = abs(xrealinitial-xrealfinal)/(xpoints);
    double delximag = abs(ximaginitial-ximagfinal)/(5.0*(double)points);

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double sinitial = 8.70;
    double sfinal = (double)real(phibthreshold(a,m));//8.0;
    double dels = abs(sinitial-sfinal)/(double)delspoints;

    int scount = 0;
    for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    {
        int pcount = 0;
        comp s = s3 + ii*simag;
        comp q = pmom(s,sigb-ii*eps1,m);
        comp kmax = pmom(s,0.0,m);

    
        cout<<"s = "<<s<<'\t'<<"dels = "<<dels<<endl;
        cout<<"eps = "<<eps<<endl;

        vector<comp> qvec;
        mom_vector_maker_4(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        
        cout<<"qvec size = "<<qvec.size()<<endl;

        //This part writes the qvec in to a file //
        string qvecfile =   "qvec_a_" + to_string(a)
                          + "_scount_" + to_string(scount)
                          + ".dat";


        fout.open(qvecfile.c_str());
        for(int i=0;i<qvec.size();++i)
        {
            fout<<i<<'\t'<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
        }
        fout.close();
        cout<<"qvecfile : "<<qvecfile<<endl;

        //This part writes the kernel with k=q
        double pRinitial = -0.2;
        double pRfinal = 0.2;
        double delpR = abs(pRinitial-pRfinal)/points;
        double pIinitial = -0.5;
        double pIfinal = 0.5;
        double delpI = abs(pIinitial - pIfinal)/points;
        int isize = (int) points;
        int jsize = (int) points;

        string kernelfile =   "kernel_a_" + to_string(a)
                            + "_scount_" + to_string(scount)
                            + ".dat";

        fout.open(kernelfile.c_str());
        for(int i=0;i<isize;++i)
        {
            for(int j=0;j<jsize;++j)
            {
                comp pR = (comp)pRinitial + (comp)i*delpR;
                comp pI = (comp)pIinitial + (comp)j*delpI; 
                comp pcomp = (comp) pR + ii*pI;

                comp sigp = sigma_p(s,pcomp,m);

                double pi = acos(-1.0);

                comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            
                fout<<setprecision(16)<<real(pR)<<'\t'
                    <<setprecision(16)<<real(pI)<<'\t'
                    <<setprecision(16)<<real(ope)<<'\t'
                    <<setprecision(16)<<imag(ope)<<endl;
            
            }
        
        }
        fout.close();
        cout<<"kernelfile : "<<kernelfile<<endl;

        //Here we print the relevant thresholds for 
        //for the kernel and trace the pcut for the OPE

        string thresholdfile = "threshold_a_" + to_string(a)
                            + "_scount_" + to_string(scount)
                            + ".dat";

        
        comp pcutp_close;
        comp pcutm_close;
        double x_pcutp;
        double x_pcutm;

        vector<double> xvec;
        vector<comp> pcutpvec;
        vector<comp> pcutmvec;

        for(int i=0;i<xpoints+1;++i)
        {
            double x = xrealinitial + i*delxreal;

            xvec.push_back(x);
            
            comp pcutp = pcut_plus(x,s,q,m,eps);
            comp pcutm = pcut_minus(x,s,q,m,eps);

            pcutpvec.push_back(pcutp);
            pcutmvec.push_back(pcutm);

        }

        pcut_fixer(pcutpvec,pcutmvec);
        //closest_pcuts_bottom(s,q,m,eps,pcutp_close,pcutm_close,x_pcutp,x_pcutm);
        closest_pcuts_bottom_withpcutvec(pcutpvec,pcutmvec,xvec,s,q,m,eps,pcutp_close,pcutm_close,x_pcutp,x_pcutm);
        
        cout<<"pcutp_close = "<<pcutp_close<<endl;
        cout<<"pcutm_close = "<<pcutm_close<<endl;
        cout<<"x_pcutp = "<<x_pcutp<<endl;
        cout<<"x_pcutm = "<<x_pcutm<<endl;
        

        fout.open(thresholdfile.c_str());
        fout<<real(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<'\t'
            <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<eps<<'\t'
            <<real(pcutp_close)<<'\t'
            <<imag(pcutp_close)<<'\t'
            <<real(pcutm_close)<<'\t'
            <<imag(pcutm_close)<<endl;

        fout.close();
        cout<<"thresholdfile : "<<thresholdfile<<endl;

        string opecutfile = "opecuttracer_a_" + to_string(a) 
                          + "_scount_" + to_string(scount)
                          + ".dat";

        fout.open(opecutfile.c_str()); 
        double somexcount = 1;                  
        //for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        for(int i=0;i<pcutpvec.size();++i)
        {
            double x = xvec[i];
            comp pcutp = pcutpvec[i];//pcut_plus(x,s,q,m,eps);
            comp pcutm = pcutmvec[i];//pcut_minus(x,s,q,m,eps);
            fout<<x<<'\t'
                <<real(pcutp)<<'\t'
                <<imag(pcutp)<<'\t'
                <<real(pcutm)<<'\t'
                <<imag(pcutm)<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<abs(-q-pcutp)<<'\t'
                <<abs(q-pcutm)<<endl;
            somexcount = somexcount + 1;
        }
        fout.close();

        cout<<"number of x summed = "<<somexcount<<endl;

        cout<<"opecutfile : "<<opecutfile<<endl;

        
        
        cout<<endl;
        scount = scount + 1;
    }
}

void kernel_pcut_criticalepsilontest_vs_single_s3()
{
    ofstream fout;
    double a = 16.0;
    double m = 1.0;
    double sreal = 8.72;//3.0;
    double simag = 0.1;
    comp ii = {0.0,1.0};
    double eps = 0.01;
    double eps1 = 0.0;

    comp sigb = sigmab(a,m);
    


    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int delspoints = 100;
    int qvecpoints = 5000;
    double qvec_r = 0.01;

    double delxreal = abs(xrealinitial-xrealfinal)/500000;//(5.0*(double)points);
    double delximag = abs(ximaginitial-ximagfinal)/(5.0*(double)points);

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double sinitial = 8.70;
    double sfinal = (double)real(phibthreshold(a,m));//8.0;
    double dels = abs(sinitial-sfinal)/(double)delspoints;

    int scount = 0;
    for(double s3=sinitial;s3<sinitial+dels;s3=s3+dels)
    {
        int pcount = 0;
        comp s = s3 + ii*simag;
        comp q = pmom(s,sigb-ii*eps1,m);
        comp kmax = pmom(s,0.0,m);

    
        cout<<"s = "<<s<<'\t'<<"dels = "<<dels<<endl;

        vector<comp> qvec;
        mom_vector_maker_4(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        
        cout<<"qvec size = "<<qvec.size()<<endl;

        //This part writes the qvec in to a file //
        string qvecfile =   "qvec_a_" + to_string(a)
                          + "_eps_" + to_string(eps)
                          + "_scount_" + to_string(scount)
                          + ".dat";


        fout.open(qvecfile.c_str());
        for(int i=0;i<qvec.size();++i)
        {
            fout<<i<<'\t'<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
        }
        fout.close();
        cout<<"qvecfile : "<<qvecfile<<endl;

        //This part writes the kernel with k=q
        double pRinitial = -0.2;
        double pRfinal = 0.2;
        double delpR = abs(pRinitial-pRfinal)/points;
        double pIinitial = -0.5;
        double pIfinal = 0.5;
        double delpI = abs(pIinitial - pIfinal)/points;
        int isize = (int) points;
        int jsize = (int) points;

        string kernelfile =   "kernel_a_" + to_string(a)
                            + "_eps_" + to_string(eps)
                            + "_scount_" + to_string(scount)
                            + ".dat";

        fout.open(kernelfile.c_str());
        for(int i=0;i<isize;++i)
        {
            for(int j=0;j<jsize;++j)
            {
                comp pR = (comp)pRinitial + (comp)i*delpR;
                comp pI = (comp)pIinitial + (comp)j*delpI; 
                comp pcomp = (comp) pR + ii*pI;

                comp sigp = sigma_p(s,pcomp,m);

                double pi = acos(-1.0);

                comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            
                fout<<setprecision(16)<<real(pR)<<'\t'
                    <<setprecision(16)<<real(pI)<<'\t'
                    <<setprecision(16)<<real(ope)<<'\t'
                    <<setprecision(16)<<imag(ope)<<endl;
            
            }
        
        }
        fout.close();
        cout<<"kernelfile : "<<kernelfile<<endl;

        //Here we print the relevant thresholds for 
        //for the kernel and trace the pcut for the OPE

        string thresholdfile = "threshold_a_" + to_string(a)
                            + "_eps_" + to_string(eps)
                            + "_scount_" + to_string(scount)
                            + ".dat";

        
        comp pcutp_close;
        comp pcutm_close;
        double x_pcutp;
        double x_pcutm;

        closest_pcuts_bottom(s,q,m,eps,pcutp_close,pcutm_close,x_pcutp,x_pcutm);
        cout<<"pcutp_close = "<<pcutp_close<<endl;
        cout<<"pcutm_close = "<<pcutm_close<<endl;
        cout<<"x_pcutp = "<<x_pcutp<<endl;
        cout<<"x_pcutm = "<<x_pcutm<<endl;
        

        fout.open(thresholdfile.c_str());
        fout<<real(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<'\t'
            <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<eps<<'\t'
            <<real(pcutp_close)<<'\t'
            <<imag(pcutp_close)<<'\t'
            <<real(pcutm_close)<<'\t'
            <<imag(pcutm_close)<<endl;

        fout.close();
        cout<<"thresholdfile : "<<thresholdfile<<endl;

        string opecutfile = "opecuttracer_a_" + to_string(a) 
                          + "_eps_" + to_string(eps)
                          + "_scount_" + to_string(scount)
                          + ".dat";

        fout.open(opecutfile.c_str()); 
        double somexcount = 1;                  
        for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        {
            comp pcutp = pcut_plus(x,s,q,m,eps);
            comp pcutm = pcut_minus(x,s,q,m,eps);
            fout<<x<<'\t'
                <<real(pcutp)<<'\t'
                <<imag(pcutp)<<'\t'
                <<real(pcutm)<<'\t'
                <<imag(pcutm)<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<abs(-q-pcutp)<<'\t'
                <<abs(q-pcutm)<<endl;
            somexcount = somexcount + 1;
        }
        fout.close();

        cout<<"number of x summed = "<<somexcount<<endl;

        cout<<"opecutfile : "<<opecutfile<<endl;

        
        
        cout<<endl;
        scount = scount + 1;
    }
}

void test_minimum_distance_vs_s3_eps()
{
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    double eps1 = 0.0066;
    double eps2 = 0.0;


    double sinitial = 8.72;
    double sfinal = real(phibthreshold(a,m));
    double simag = 0.1;

    double dels = abs(sinitial - sfinal)/100.0;

    for(int i=0;i<100;++i)
    {
        double sreal = sinitial + i*dels;

        comp s = sreal + ii*simag;

        comp    pcutp_close1, pcutp_close2,
                pcutm_close1, pcutm_close2;
        double  x_pcutp1, x_pcutm1,
                x_pcutp2, x_pcutm2;
        
        comp sigb = sigmab(a,m);
        comp q = pmom(s,sigb,m);

        closest_pcuts_bottom(s,q,m,eps1,pcutp_close1,pcutm_close1,x_pcutp1,x_pcutm1);
        closest_pcuts_bottom(s,q,m,eps2,pcutp_close2,pcutm_close2,x_pcutp2,x_pcutm2);

        cout<<"s = "<<s<<endl;
        cout<<"eps1:"<<eps1<<'\t'
            <<"pcp1:"<<pcutp_close1<<'\t'
            <<"pcm1:"<<pcutm_close1<<endl;
        cout<<"eps2:"<<eps2<<'\t'
            <<"pcp2:"<<pcutp_close2<<'\t'
            <<"pcm2:"<<pcutm_close2<<endl;
        cout<<"xcp1:"<<x_pcutp1<<'\t'
            <<"xcm1:"<<x_pcutm1<<endl;
        cout<<"xcp2:"<<x_pcutp2<<'\t'
            <<"xcm2:"<<x_pcutm2<<endl;

        
    }

}


void kernel_pcut_x_glockletest_vs_sigk_allqvec()
{
    ofstream fout;
    double a = 16.0;
    double m = 1.0;
    double sreal = 3.0;//3.0;
    double simag = 0.0;
    comp ii = {0.0,1.0};
    double eps = 1.0e-5;
    double eps1 = 0.0;

    comp sigb = sigmab(a,m);
    


    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int delspoints = 100;
    int qvecpoints = 4000;
    double qvec_r = 0.01;

    double delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double sinitial = 8.72;
    double sfinal = (double)real(phibthreshold(a,m));//8.0;
    double dels = abs(sinitial-sfinal)/(double)delspoints;

    double sigkinitial = 2.0*m*m;
    double sigkfinal = 4.0*m*m;
    double delsigk = abs(sigkinitial - sigkfinal)/200.0;

    int scount = 0;
    //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    for(double sigk=sigkinitial;sigk<=sigkfinal;sigk=sigk+delsigk)
    {
        int pcount = 0;
        comp s = sinitial;//s3 + ii*simag;
        //comp q = pmom(s,sigb-ii*eps1,m);
        comp actualq = pmom(s,sigb-ii*eps1,m);
        comp q = pmom(s,sigk,m);
        comp kmax = pmom(s,0.0,m);

    
        cout<<"s = "<<s<<'\t'<<"dels = "<<dels<<endl;

        vector<comp> qvec;
        mom_vector_maker_4(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        
        cout<<"qvec size = "<<qvec.size()<<endl;

        //This part writes the qvec in to a file //
        string qvecfile =   "qvec_a_" + to_string(a)
                          + "_scount_" + to_string(scount)
                          + ".dat";


        fout.open(qvecfile.c_str());
        for(int i=0;i<qvec.size();++i)
        {
            fout<<i<<'\t'<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
        }
        fout.close();
        cout<<"qvecfile : "<<qvecfile<<endl;

        //This part writes the kernel with k=q
        double pRinitial = -1.0;
        double pRfinal = 1.0;
        double delpR = abs(pRinitial-pRfinal)/points;
        double pIinitial = -1.0;
        double pIfinal = 1.0;
        double delpI = abs(pIinitial - pIfinal)/points;
        int isize = (int) points;
        int jsize = (int) points;

        string kernelfile =   "kernel_a_" + to_string(a)
                            + "_scount_" + to_string(scount)
                            + ".dat";

        fout.open(kernelfile.c_str());
        for(int i=0;i<isize;++i)
        {
            for(int j=0;j<jsize;++j)
            {
                comp pR = (comp)pRinitial + (comp)i*delpR;
                comp pI = (comp)pIinitial + (comp)j*delpI; 
                comp pcomp = (comp) pR + ii*pI;

                comp sigp = sigma_p(s,pcomp,m);

                double pi = acos(-1.0);

                comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            
                fout<<setprecision(16)<<real(pR)<<'\t'
                    <<setprecision(16)<<real(pI)<<'\t'
                    <<setprecision(16)<<real(ope)<<'\t'
                    <<setprecision(16)<<imag(ope)<<endl;
            
            }
        
        }
        fout.close();
        cout<<"kernelfile : "<<kernelfile<<endl;

        //Here we print the relevant thresholds for 
        //for the kernel and trace the pcut for the OPE

        string thresholdfile = "threshold_a_" + to_string(a)
                            + "_scount_" + to_string(scount)
                            + ".dat";

        
        fout.open(thresholdfile.c_str());
        fout<<real(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<'\t'
            <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<eps<<'\t'
            <<sigk<<'\t'
            <<real(actualq)<<'\t'
            <<imag(actualq)<<endl;

        fout.close();
        cout<<"thresholdfile : "<<thresholdfile<<endl;

        delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);

        string opecutfile = "opecuttracer_a_" + to_string(a) 
                          + "_scount_" + to_string(scount)
                          + ".dat";

        fout.open(opecutfile.c_str());  

        vector<double> xvec;
        vector<comp> pcutpvec;
        vector<comp> pcutmvec;

        for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        {
            xvec.push_back(x);
            comp pcutp = pcut_plus(x,s,q,m,eps);
            comp pcutm = pcut_minus(x,s,q,m,eps);
            pcutpvec.push_back(pcutp);
            pcutmvec.push_back(pcutm);

        }
        pcut_fixer(pcutpvec,pcutmvec);

        //for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        for(int i=0;i<pcutpvec.size();++i)
        {
            double x = xvec[i];
            comp pcutp = pcutpvec[i];//pcut_plus(x,s,q,m,eps);
            comp pcutm = pcutmvec[i];//pcut_minus(x,s,q,m,eps);

            fout<<x<<'\t'
                <<real(pcutp)<<'\t'
                <<imag(pcutp)<<'\t'
                <<real(pcutm)<<'\t'
                <<imag(pcutm)<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<abs(-q-pcutp)<<'\t'
                <<abs(q-pcutm)<<endl;
        }
        fout.close();

        cout<<"opecutfile : "<<opecutfile<<endl;

        string pcutfile = "pcut_a_" + to_string(a) 
                          + "_scount_" + to_string(scount)
                          + ".dat";

        fout.open(pcutfile.c_str());  

        pcount = 0;

        delxreal = abs(xrealinitial-xrealfinal)/((double)points);
        for(int k=0;k<qvec.size();++k)
        {
                             
            for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
            {
                comp pcutp = pcut_plus(x,s,qvec[k],m,eps);
                comp pcutm = pcut_minus(x,s,qvec[k],m,eps);
                fout<<x<<'\t'
                    <<real(pcutp)<<'\t'
                    <<imag(pcutp)<<'\t'
                    <<real(pcutm)<<'\t'
                    <<imag(pcutm)<<'\t'
                    <<real(q)<<'\t'
                    <<imag(q)<<'\t'
                    <<abs(-q-pcutp)<<'\t'
                    <<abs(q-pcutm)<<endl;
            }
            
            //cout<<"qvec["<<k<<"] = "<<qvec[k]<<'\t'<<"pcount = "<<pcount<<endl;
            //cout<<"pcutfile : "<<pcutfile<<endl;

            //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;
            //cout<<"file:"<<filename1<<endl;
            pcount = pcount + 1;
            //cout<<"loop completed = "<<pcount<<endl;
        }
        fout.close();
        cout<<"pcutfile : "<<pcutfile<<endl;
        cout<<endl;
        scount = scount + 1;
    }
}

void kernel_pcut_x_glockletest_vs_k_allqvec()
{
    ofstream fout;
    double a = 16.0;
    double m = 1.0;
    double sreal = 3.0;//3.0;
    double simag = 0.0;
    comp ii = {0.0,1.0};
    double eps = 0.0;//1.0e-5;
    double eps1 = 0.0;

    comp sigb = sigmab(a,m);
    


    double xrealinitial = -1.0;
    double xrealfinal = 1.0;
    double ximaginitial = 2.0*xrealinitial;
    double ximagfinal = 2.0*xrealfinal;

    int points = 250;
    int delspoints = 100;
    int qvecpoints = 500;
    double qvec_r = 0.01;

    double delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);
    double delximag = abs(ximaginitial-ximagfinal)/(double)points;

    double pinitial = 0.0;
    //double pfinal = real(kmax);
    //double delp = abs(pinitial-pfinal)/50.0;

    double sinitial = 8.72;
    double sfinal = (double)real(phibthreshold(a,m));//8.0;
    double dels = abs(sinitial-sfinal)/(double)delspoints;

    double sigkinitial = 2.0*m*m;
    double sigkfinal = 4.0*m*m;
    double delsigk = abs(sigkinitial - sigkfinal)/200.0;

    double spole = 8.782848835105487;
    comp s = spole;
    double kinitial = 0.0;
    double kfinal = real(pmom(s,0.0,m));
    double kpoints = 200;
    double delk = abs(kinitial - kfinal)/kpoints;

    int scount = 0;
    //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    //for(double sigk=sigkinitial;sigk<=sigkfinal;sigk=sigk+delsigk)
    for(int i=0;i<kpoints;++i) 
    {
        double k = kinitial + i*delk;
        comp sigk = sigma_p(s,k,m);
        int pcount = 0;
        //comp s = sinitial;//s3 + ii*simag;
        //comp q = pmom(s,sigb-ii*eps1,m);
        comp actualq = pmom(s,sigb-ii*eps1,m);
        double kmin = 0.0;
        comp q = k;//pmom(s,sigk,m);
        comp kmax = pmom(s,0.0,m);

    
        cout<<"s = "<<s<<'\t'<<"dels = "<<dels<<endl;

        vector<comp> qvec;
        //mom_vector_maker_4(qvec,s,0.0,kmax,a,m,eps,eps1,(double)qvecpoints,qvec_r);
        line_maker(qvec,kmin,kmax,qvecpoints);
        
        cout<<"qvec size = "<<qvec.size()<<endl;

        //This part writes the qvec in to a file //
        string qvecfile =   "qvec_a_" + to_string(a)
                          + "_scount_" + to_string(scount)
                          + ".dat";


        fout.open(qvecfile.c_str());
        for(int i=0;i<qvec.size();++i)
        {
            fout<<i<<'\t'<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
        }
        fout.close();
        cout<<"qvecfile : "<<qvecfile<<endl;

        //This part writes the kernel with k=q
        double pRinitial = -1.0;
        double pRfinal = kfinal;
        double delpR = abs(pRinitial-pRfinal)/points;
        double pIinitial = -1.0;
        double pIfinal = 1.0;
        double delpI = abs(pIinitial - pIfinal)/points;
        int isize = (int) points;
        int jsize = (int) points;

        string kernelfile =   "kernel_a_" + to_string(a)
                            + "_scount_" + to_string(scount)
                            + ".dat";

        fout.open(kernelfile.c_str());
        for(int i=0;i<isize;++i)
        {
            for(int j=0;j<jsize;++j)
            {
                comp pR = (comp)pRinitial + (comp)i*delpR;
                comp pI = (comp)pIinitial + (comp)j*delpI; 
                comp pcomp = (comp) pR + ii*pI;

                comp sigp = sigma_p(s,pcomp,m);

                double pi = acos(-1.0);

                comp ope = (pcomp*pcomp/(pow(2.0*pi,2.0)*omega_comp(pcomp,m)))*GS_pk(s,q,pcomp,m,eps)*M2kfunc(a,sigp,m,eps1);
            
                fout<<setprecision(16)<<real(pR)<<'\t'
                    <<setprecision(16)<<real(pI)<<'\t'
                    <<setprecision(16)<<real(ope)<<'\t'
                    <<setprecision(16)<<imag(ope)<<endl;
            
            }
        
        }
        fout.close();
        cout<<"kernelfile : "<<kernelfile<<endl;

        //Here we print the relevant thresholds for 
        //for the kernel and trace the pcut for the OPE

        string thresholdfile = "threshold_a_" + to_string(a)
                            + "_scount_" + to_string(scount)
                            + ".dat";

        
        fout.open(thresholdfile.c_str());
        fout<<real(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_plus(-1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(1.0,s,q,m,eps))<<'\t'
            <<real(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<imag(pcut_minus(-1.0,s,q,m,eps))<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<'\t'
            <<real(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_plus(s,m))<<'\t'
            <<real(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<imag(M2kbranchcut_left_momrep_minus(s,m))<<'\t'
            <<real(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_plus_eps(s,m,eps1))<<'\t'
            <<real(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<imag(M2kbranchcut_right_momrep_minus_eps(s,m,eps1))<<'\t'
            <<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<eps<<'\t'
            <<real(sigk)<<'\t'
            <<imag(sigk)<<'\t'
            <<real(actualq)<<'\t'
            <<imag(actualq)<<endl;

        fout.close();
        cout<<"thresholdfile : "<<thresholdfile<<endl;

        delxreal = abs(xrealinitial-xrealfinal)/(20.0*points);

        string opecutfile = "opecuttracer_a_" + to_string(a) 
                          + "_scount_" + to_string(scount)
                          + ".dat";

        fout.open(opecutfile.c_str());  

        vector<double> xvec;
        vector<comp> pcutpvec;
        vector<comp> pcutmvec;

        for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        {
            xvec.push_back(x);
            comp pcutp = pcut_plus(x,s,q,m,eps);
            comp pcutm = pcut_minus(x,s,q,m,eps);
            pcutpvec.push_back(pcutp);
            pcutmvec.push_back(pcutm);

        }
        pcut_fixer(pcutpvec,pcutmvec);

        //for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
        for(int i=0;i<pcutpvec.size();++i)
        {
            double x = xvec[i];
            comp pcutp = pcutpvec[i];//pcut_plus(x,s,q,m,eps);
            comp pcutm = pcutmvec[i];//pcut_minus(x,s,q,m,eps);

            fout<<x<<'\t'
                <<real(pcutp)<<'\t'
                <<imag(pcutp)<<'\t'
                <<real(pcutm)<<'\t'
                <<imag(pcutm)<<'\t'
                <<real(q)<<'\t'
                <<imag(q)<<'\t'
                <<abs(-q-pcutp)<<'\t'
                <<abs(q-pcutm)<<endl;
        }
        fout.close();

        cout<<"opecutfile : "<<opecutfile<<endl;

        string pcutfile = "pcut_a_" + to_string(a) 
                          + "_scount_" + to_string(scount)
                          + ".dat";

        fout.open(pcutfile.c_str());  

        pcount = 0;

        delxreal = abs(xrealinitial-xrealfinal)/((double)points);
        for(int k=0;k<qvec.size();++k)
        {
                             
            for(double x=xrealinitial;x<=xrealfinal;x=x+delxreal)
            {
                comp pcutp = pcut_plus(x,s,qvec[k],m,eps);
                comp pcutm = pcut_minus(x,s,qvec[k],m,eps);
                fout<<x<<'\t'
                    <<real(pcutp)<<'\t'
                    <<imag(pcutp)<<'\t'
                    <<real(pcutm)<<'\t'
                    <<imag(pcutm)<<'\t'
                    <<real(q)<<'\t'
                    <<imag(q)<<'\t'
                    <<abs(-q-pcutp)<<'\t'
                    <<abs(q-pcutm)<<endl;
            }
            
            //cout<<"qvec["<<k<<"] = "<<qvec[k]<<'\t'<<"pcount = "<<pcount<<endl;
            //cout<<"pcutfile : "<<pcutfile<<endl;

            //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<endl;
            //cout<<"file:"<<filename1<<endl;
            pcount = pcount + 1;
            //cout<<"loop completed = "<<pcount<<endl;
        }
        fout.close();
        cout<<"pcutfile : "<<pcutfile<<endl;
        cout<<endl;
        scount = scount + 1;
    }
}


void testq()
{
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;

    double sreal = 8.93;
    double simaginitial = -0.1;
    double simagfinal = 0.1;
    double spoints = 100.0;

    double delsimag = abs(simaginitial - simagfinal)/spoints;

    for(int i=0;i<spoints+1;++i)
    {
        double simag = simaginitial + i*delsimag;
        comp s = sreal + ii*simag;
        comp sigb = sigmab(a,m);
        comp q = pmom(s,sigb,m);

        cout<<s<<'\t'<<q<<endl;

    }
}

void check_qvec_for_ims_pos()
{
    comp s = {8.6,0.05};
    double a = 16.0;
    double m = 1.0;
    comp min = {0.0,0.0};
    comp max = pmom(s,0.0,m);
    double eps = 0.0;
    double eps_for_m2k = eps;

    vector<comp> qvec;
    mom_vector_maker_47(qvec,s, min,max,a,m,eps,eps_for_m2k,2000,0.01);

    for(int i=0;i<qvec.size();++i){
        cout<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
    }
}

int main(int argc, char *argv[])
{
    //check_qvec_for_ims_pos();
    //double eps = atof(argv[1]);
    //cout<<"eps selected = "<<eps<<endl;
    //comp ii = {0.0,1.0};
    //comp s = 8.93 + ii*0.1;
    //double a = 16.0;
    //double m = 1.0;
    //comp qc1 = q_c1(s,a,m,eps);
    //comp rad;
    //cout<<"qc1 = "<<qc1<<endl;
    //rad = real(qc1);
    //cout<<"radius = "<<rad<<endl;
    
    //comp ii = {0.0,1.0};
    //double sigp = 3.5;
    //comp sigb = sigmab(16,1);
    //comp sigk = sigb;
    //comp s = 8.72 + ii*1.0e-5;
    //double eps1 = 1.0e-5;
    /*cout<<zppprime(s,sigp,sigk,1.0,eps1)<<'\t'
        <<zppprime_1(s,sigp,sigk,1.0,eps1)<<'\t'
        <<zppprime_complex_s(real(s),sigp,sigk,1.0,eps1,imag(s))<<endl;
    */

    //zpk_energyfactor_contour_sigp_vs_sigk_withintloop();
    //cout<<"phibthreshold:"<<pow(sqrt(sigmab(16,1)) + 1.0,2.0)<<endl;
    //zpk_energyfactor_contour_sigp_vs_s_withintloop();
    //opeplotcontour();
    //opeplotcontour_withintloop();
    //opeplotcontour_withintloop_withloopeps();//
    //m2kplotcontour_momrep_withintloop_withepsloop();
    //omegaplotcontour_momrep_withintloop_withepsloop();
    //kernelplotcontour_momrep_withintloop_withsloop();
    //kernelplotcontour_momrep_withintloop_with_multiplier();
    //kernelplotcontour_momrep_withintloop_glockletest();

    //branchcuts_and_thresholds_sloop();
    //branchcuts_and_thresholds_qvecloop();
    //glockletest();

    //kernel_pcut_x_glockletest_vs_s3_qvec();
    //testing_pcuts();

    kernel_pcut_x_glockletest_vs_s3_allqvec();//this was the last last thing we were checking

    //M2k in sigma space:
    //kernel_pcut_x_glockletest_vs_s3_allqvec_sigma_space();
    //

    //ope_qq_plot_vs_real_s(); //checking the sign of the ope branch points are correct or not for Mphib

    //kernel_pcut_x_glockletest_vs_a_allqvec();

    //kernel_pcut_x_glockletest_vs_single_a_allqvec();
    
    //kernel_pcut_x_glockletest_vs_single_s3_allqvec(); //this was the last one that we were checking

    //kernel_pcut_x_glockletest_vs_single_s3_singleqvec();// this was the last last last thing we were checking
    //kernel_pcut_x_glockletest_vs_single_s3_singleqvec_nonanalyticity_singleplot();
    //kernel_pcut_x_glockletest_vs_s3_allqvec_omp();
    //kernel_pcut_x_glockletest_vs_s3imag_allqvec();

    //testq();
   
    //kernel_pcut_x_glockletest_vs_eps_allqvec();
    //kernel_pcut_x_glockletest_vs_sigk_allqvec();
    //kernel_pcut_x_glockletest_vs_k_allqvec();


    //kernel_pcut_criticalepsilontest_vs_s3(eps);
    //kernel_pcut_criticalepsilontest_vs_single_s3();
    //test_minimum_distance_vs_s3_eps();
    //test_with_Mphib();
    //test_distance_between_pcuts();
    //test_critical_epsilon_pcut_numbers();
    //test_critical_epsilon_vs_pcontour_with0points();
    //test_critical_epsilon_vs_s3();
    
    //kernel_pcut_x_glockletest_vs_s3imag_eps_criticalepstest();
    //comp s = {8.92,0.1};
    //testing_pcuts_and_kernel(s);
    //kernel_pcut_x_glockletest_vs_s3_eps_forsigk2m();
    //kernel_pcut_x_glockletest_vs_s3_eps();
    //testing_momvec_maker_2();
    //kernelplotcontour_momrep_withintloop_withepsloop_drawpcut();
    //OPE_PW_poletesting_x();
    //only_pcut_x_comp();
    //glockletest_C1();
    //ope_tau_m2_plot_withsigvec();
    //opeplotcontour_momrep_withintloop();
    //opeplotcontour_momrep_withintloop_withepsloop();
    //kernelplotcontour_withintloop();
    //kernelplotcontour_withintloop_vs_s3();
    //kernelplot_withsigvec();
    //opeplotcontour_withintloop_multisigvec();
    //kernelplotcontour_withintloop_vs_sigvec_delta();
    //test_lowest_point();
    //test_sigma_plus_reflection();
    //tauplotcontour_withintloop();
    //opeplot();
    //kernelplotcontour();
    //mommaker();
    //testing_sigp_qfunc();

    //comp phibthres = phibthreshold()
    /*double s = 8.85;
    double a = 2.0;
    double m = 1.0;

    comp sb = sigmab(a,m);

    comp q = pmom(s,sb,m);

    cout<<q<<endl;

    cout<<qfunc(sqrt(s),a,0,m)<<endl; 
    
    //list_of_points();*/
    /*comp ii = {0.0,1.0};
    comp scomp = 8.6 + ii*0.02;
    comp sigmapprime = 3.98 + ii*0.02;
    double m = 1.0;
    double eps = 0;
    cout<<sigpplus_witheps(scomp,sigmapprime,m,eps)<<'\n'
        <<sigpminus_witheps(scomp,sigmapprime,m,eps)<<'\n'
        <<sigc1(scomp,sigmapprime,m)<<'\n'
        <<sigc2(scomp,sigmapprime,m)<<endl;
    */
    /*for(double d=-0.01;d<=0.01;d=d+0.001)
    {
        scomp = 8.6 + ii*d;
        sigmapprime = 3.984 + ii*d;
        cout<<real(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpplus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<real(sigpminus_witheps(scomp,sigmapprime,m,eps))<<'\t'
        <<imag(sigpminus_witheps(scomp,sigmapprime,m,eps))<<endl;
    }*/

    /* some testing functions */

    /*
    
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;
    double eps = 0.001;
    double s = 8.90;
    comp sigb = sigmab(a,m);

    comp q = pmom(s,sigb-ii*eps,m);
    comp phib = phibthreshold(a,m);
    cout<<"phibth:"<<phib<<endl;

    comp somex1 = x_for_qc1_minus(s,q,a,m,eps);
    comp somex2 = x_for_qc1_plus(s,q,a,m,eps);

    comp qc1 = q_c1(s,a,m,eps);
    comp qc2 = q_c2(s,a,m,eps);

    cout<<"qc1:"<<qc1<<endl;
    cout<<"qc2:"<<qc2<<endl;

    comp x1 = x_for_pk(s,qc1,q,m,eps);
    comp x2 = x_for_pk(s,qc2,q,m,eps);
    cout<<"x1:"<<x1<<endl;
    cout<<"x2:"<<x2<<endl;

    //double realx = real(x);

    //cout<<x_for_qc1_minus(s,q,a,m,eps)<<'\t'<<x_for_qc1_plus(s,q,a,m,eps)<<endl;
    comp pcutplusx1 = pcut_plus_comp(x1,s,q,m,eps);
    comp pcutplusx2 = pcut_plus_comp(x2,s,q,m,eps);
    comp pcutminusx1 = pcut_minus_comp(real(x1),s,q,m,eps);
    comp pcutminusx2 = pcut_minus_comp(real(x2),s,q,m,eps);

    cout<<"pcutplusx1:"<<pcutplusx1<<endl;
    cout<<"pcutplusx2:"<<pcutplusx2<<endl;
    cout<<"pcutminusx1:"<<pcutminusx1<<endl;
    cout<<"pcutminusx2:"<<pcutminusx2<<endl;
    //cout<<q_c2(s,a,m,eps)<<'\t'<<pcut_minus_comp(x,s,q,m,eps)<<'\t'<<pcut_plus_comp(x,s,q,m,eps)<<endl; 
    //cout<<pcut_plus()

    */
    /*for(double s3=7.0;s3<=real(phib);s3=s3+0.001)
    {
        cout<<"s3:"<<s3<<'\t'
            <<"qc1:"<<q_c1(s3,a,m,eps)<<'\t'
            <<"qc2:"<<q_c2(s3,a,m,eps)<<'\t'
            <<"qplus:"<<q_plus(s3,a,m,eps)<<'\t'
            <<"qminus:"<<q_minus(s3,a,m,eps)<<endl;
    }*/

    /*
    vector<comp> qvec;
    mom_vector_maker_2(qvec,s,0.0,1.3,a,m,0.001,1000,0.001);
    for(int i=0;i<qvec.size();++i)
    {
        cout<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
    }
    */
   /*comp ii = {0.0,1.0};
   cout<<mysqrt(8.72+ii*0.1)<<endl;
   comp a = {-1.0,0.0};
   double b = -1.0;
   cout<<c<<endl;*/
   
   return 0;
}