#include<bits/stdc++.h>
#include "functions_invariant_based.h"
#include "invariant_vector_maker.h"

using namespace std;

typedef complex<double> comp;

void kernelplotcontour()
{
    ofstream fout;

    double a = 16.0;
    double m = 1.0;
    double s = 8.15;
    double eps = 1.0e-3;
    comp scomp = (comp) s;
    comp sigmapprime = sigmab(a,m);
    double sigpRinitial = real(sigc1(s,sigmab(a,m),m)) - 0.5;
    double sigpRfinal = 4.1;
    double delsigpR = abs(sigpRinitial-sigpRfinal)/2000.0;
    double sigpIinitial = -0.15;
    double sigpIfinal = 0.15;
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

int main()
{
    kernelplotcontour();
    ofstream fout;

    double a = 16.0;
    double m = 1.0;
    double s = 8.15;
    double eps = 1.0e-3;
    comp scomp = (comp) s;
    comp sigmapprime = sigmab(a,m);
    double sigpRinitial = 3.55;
    double sigpRfinal = 5.0;
    double delsigpR = abs(sigpRinitial-sigpRfinal)/2000.0;
    double sigpIinitial = -0.1;
    double sigpIfinal = 0.1;
    double delsigpI = abs(sigpIinitial - sigpIfinal)/2000.0;

    double points = 2000.0;

    vector<comp> sigvec;

    sigma_vector_maker(sigvec, s, 0.0, sigmax(s,m), a, m, points, 1.0e-5);


    string filename="sigvec.dat";
    fout.open(filename.c_str());

    for(int i=0;i<sigvec.size();++i)
    {
        fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;
        fout<<endl;
    }
    fout.close();
    //cout<<sigvec.size()<<endl;

    return 0;

}