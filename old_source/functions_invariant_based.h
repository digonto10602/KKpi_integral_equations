#ifndef FUNCTIONS_INVARIANT_BASED
#define FUNCTIONS_INVARIANT_BASED

#include<bits/stdc++.h>
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

double Hfunc(   comp sigmap,
                comp sigmapprime,
                double m    )
{
    return Jfunc(real(sigmap/(4.0*m*m)))*Jfunc(real(sigmapprime/(4.0*m*m)));
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