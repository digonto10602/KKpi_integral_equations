#include<bits/stdc++.h>
#include "functions_invariant_based.h"
#include "invariant_vector_maker.h"

using namespace std;

typedef complex<double> comp;

comp GS1(    comp s,
            comp sigmap,
            comp sigmapprime,
            double m,
            double epsilon      )
{
    comp p = mysqrt(kallentriangle(s,sigmap,m*m))/(2.0*sqrt(s));
    comp pprime = mysqrt(kallentriangle(s,sigmapprime,m*m))/(2.0*sqrt(s));

    comp zpp = zppprime(s,sigmap,sigmapprime,m,epsilon);

    return -1.0/(4.0*p*pprime)*(log( (zpp - 1.0)/(zpp + 1.0) ));

}

void kernelplotcontour()
{
    ofstream fout;

    double a = 16.0;
    double m = 1.0;
    //double s = 8.65;
    double eps = 1.0e-3;
    
    comp sigmapprime = sigmab(a,m);
    double sigpRinitial = 3.55;
    double sigpRfinal = 5.0;
    double delsigpR = abs(sigpRinitial-sigpRfinal)/500.0;
    double sigpIinitial = -0.1;
    double sigpIfinal = 0.1;
    double delsigpI = abs(sigpIinitial - sigpIfinal)/500.0;
    
    double sinitial = 8.72;
    double sfinal = real(phibthreshold(a,m));
    double dels = abs(sinitial - sfinal)/50.0;

    cout<<setprecision(10)<<dels<<endl;

    for(int count=0;count<50;++count)
    { 

        double s = sinitial + (double)count*dels;

        comp scomp = (comp) s;

    string filename = "kernelcontour_a_" + to_string((int)a) + "_count_" + to_string(count) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename.c_str());
    for(double sigmapR=sigpRinitial;sigmapR<=sigpRfinal;sigmapR=sigmapR + delsigpR)
    {
        
        for(double sigmapI=sigpIinitial;sigmapI<=sigpIfinal;sigmapI=sigmapI + delsigpI)
        //double someeps = 1.0e-3;
        {
            comp sigmapcomp = (comp) sigmapR + ii*sigmapI;

            //comp ope = kernel(scomp,sigmapcomp,sigmapprime,a,m,eps);
            comp ope = kernel_1(scomp,sigmapprime,sigmapcomp,a,m,eps);

            fout<<s<<'\t'<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
            cout<<s<<'\t'<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
        }
        fout<<endl;
    }   

    fout.close();

    string filename1 = "threshkernelcontour_a_" + to_string((int)a) + "_count_" + to_string(count) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename1.c_str());

    fout<<real(sigpplus(s,sigmab(a,m),m))<<'\t'<<real(sigpminus(s,sigmab(a,m),m))<<'\t'<<real(sigc1(s,sigmab(a,m),m))<<'\t'<<real(sigc2(s,sigmab(a,m),m))<<'\t'<<pow(sqrt(s)-m,2.0)<<'\t'<<0<<'\t'<<0<<endl;

    fout.close();

    string filename2 = "sigveccontour_a_" + to_string((int)a) + "_count_" + to_string(count) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename2.c_str());

    vector<comp> sigvec;
    sigma_vector_maker(sigvec,s,0.0,sigmax(s,m),a,m,1000.0,1.0e-5);
    for(int i=0;i<sigvec.size();++i)
    fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;

    fout.close();
    }
}


void kernelplotcontour_1()
{
    ofstream fout;

    double a = 16.0;
    double m = 1.0;
    double s = 8.65;
    double eps = 1.0e-3;
    
    comp sigmapprime = sigmab(a,m);
    double sigpRinitial = 3.55;
    double sigpRfinal = 4.10;
    double delsigpR = abs(sigpRinitial-sigpRfinal)/500.0;
    double sigpIinitial = -0.1;
    double sigpIfinal = 0.1;
    double delsigpI = abs(sigpIinitial - sigpIfinal)/500.0;
    
    double sinitial = 8.72;
    double sfinal = real(phibthreshold(a,m));
    double dels = abs(sinitial - sfinal)/50.0;

    cout<<setprecision(10)<<dels<<endl;

    //for(int count=0;count<50;++count)
     

        //double s = sinitial + (double)count*dels;

        comp scomp = (comp) s;

    string filename = "kernelcontour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename.c_str());
    for(double sigmapR=sigpRinitial;sigmapR<=sigpRfinal;sigmapR=sigmapR + delsigpR)
    {
        
        for(double sigmapI=sigpIinitial;sigmapI<=sigpIfinal;sigmapI=sigmapI + delsigpI)
        //double someeps = 1.0e-3;
        {
            comp sigmapcomp = (comp) sigmapR + ii*sigmapI;

            //comp ope = kernel(scomp,sigmapcomp,sigmapprime,a,m,eps);
            comp ope = kernel(scomp,sigmapprime,sigmapcomp,a,m,eps);

            fout<<s<<'\t'<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
            cout<<s<<'\t'<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
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
    sigma_vector_maker_1(sigvec,s,0.0,sigmax(s,m),a,m,1000.0,1.0e-6);
    for(int i=0;i<sigvec.size();++i)
    fout<<real(sigvec[i])<<'\t'<<imag(sigvec[i])<<'\t'<<0<<endl;

    fout.close();
    
}



void tauplotcontour()
{
    ofstream fout;

    double a = 16.0;
    double m = 1.0;
    //double s = 8.65;
    double eps = 1.0e-3;
    
    comp sigmapprime = sigmab(a,m);
    double sigpRinitial = 3.55;
    double sigpRfinal = 5.0;
    double delsigpR = abs(sigpRinitial-sigpRfinal)/500.0;
    double sigpIinitial = -0.1;
    double sigpIfinal = 0.1;
    double delsigpI = abs(sigpIinitial - sigpIfinal)/500.0;
    
    double sinitial = 8.8;
    double sfinal = real(phibthreshold(a,m));
    double dels = abs(sinitial - sfinal)/50.0;

    cout<<setprecision(10)<<dels<<endl;

    
        double s = sinitial ;

        comp scomp = (comp) s;

    string filename = "tau1contour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename.c_str());
    for(double sigmapR=sigpRinitial;sigmapR<=sigpRfinal;sigmapR=sigmapR + delsigpR)
    {
        
        for(double sigmapI=sigpIinitial;sigmapI<=sigpIfinal;sigmapI=sigmapI + delsigpI)
        //double someeps = 1.0e-3;
        {
            comp sigmapcomp = (comp) sigmapR + ii*sigmapI;

            //comp ope = kernel(scomp,sigmapcomp,sigmapprime,a,m,eps);
            comp ope = tau1(scomp,sigmapcomp,m);

            fout<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
            cout<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
        }
        fout<<endl;
    }   

    fout.close();

    
}



void M2plotcontour()
{
    ofstream fout;

    double a = 16.0;
    double m = 1.0;
    //double s = 8.65;
    double eps = 1.0e-3;
    
    comp sigmapprime = sigmab(a,m);
    double sigpRinitial = 3.55;
    double sigpRfinal = 5.0;
    double delsigpR = abs(sigpRinitial-sigpRfinal)/500.0;
    double sigpIinitial = -0.1;
    double sigpIfinal = 0.1;
    double delsigpI = abs(sigpIinitial - sigpIfinal)/500.0;
    
    double sinitial = 8.72;
    double sfinal = real(phibthreshold(a,m));
    double dels = abs(sinitial - sfinal)/50.0;

    cout<<setprecision(10)<<dels<<endl;

    
        double s = sinitial ;

        comp scomp = (comp) s;

    string filename = "M2_2contour_a_" + to_string((int)a) + "_s_" + to_string(s) + "_sigpsb_" + to_string(real(sigmab(a,m))) + ".dat";

    fout.open(filename.c_str());
    for(double sigmapR=sigpRinitial;sigmapR<=sigpRfinal;sigmapR=sigmapR + delsigpR)
    {
        
        for(double sigmapI=sigpIinitial;sigmapI<=sigpIfinal;sigmapI=sigmapI + delsigpI)
        //double someeps = 1.0e-3;
        {
            comp sigmapcomp = (comp) sigmapR + ii*sigmapI;

            //comp ope = kernel(scomp,sigmapcomp,sigmapprime,a,m,eps);
            comp ope = M2kfunc(a,sigmapcomp,m,1.0e-5);

            fout<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
            cout<<sigmapR<<'\t'<<sigmapI<<'\t'<<real(ope)<<'\t'<<imag(ope)<<endl;
        }
        fout<<endl;
    }   

    fout.close();

    
}



int main()
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    double sinitial = 8.72;
    double dels = 0.0256555;

    //comp s = {8.65,0.0};
    double m = 1.0;
    double a = 16.0;
    double eps = 1.0e-5;
    comp sigb = sigmab(a,m);
    comp sigp = sigb;//{2.0,0.001};
    comp sigk = sigb;//{2.0,0.0};
    vector<double> svec(10);
    svec[0] = 8.72;
    svec[1] = 8.74565548557;
    svec[2] = 8.77131097113;
    svec[3] = 8.7969664567;
    svec[4] = 8.82262194227;
    svec[5] = 8.84827742783;
    svec[6] = 8.8739329134;
    svec[7] = 8.89958839897;
    svec[8] = 8.92524388453;
    svec[9] = 8.9508993701;

    /*
    OPE.dat Sebastian;
    8.72 10.7918485699 -0.00187969839005
    8.74565548557 12.1920083172 -0.00243732248005
    8.77131097113 14.0104730497 -0.00328631214243
    8.7969664567 16.4715768577 -0.00467549260778
    8.82262194227 20.000559326 -0.00719856383388
    8.84827742783 25.5280650758 -0.012616024249
    8.8739329134 35.6776870728 -0.0287427058856
    8.89958839897 65.6572737154 -0.186324618997
    8.92524388453 66.4503788425 -68.7747846353
    8.9508993701 56.1307912796 -137.606104325
*/

    //cout<<GS1(s,sigp,sigk,m,eps)<<endl;
    //cout<<tau(s,sigp,m)<<endl;
    /*for(int i=0;i<10;++i)
    {
        //comp s = (comp) sinitial + (comp)i*dels;
        //double sdoub = sinitial + i*dels;


        comp Gs = GS((comp)svec[i],sigb,sigb,m,1.0e-5);
        cout<<setprecision(10)<<svec[i]<<'\t'<<setprecision(10)<<real(Gs)<<'\t'<<setprecision(10)<<imag(Gs)<<endl;
    }*/
    /*for(double s2k=3.9;s2k<=4.2;s2k=s2k+0.0005)
    {
        comp s2kcomp = (comp) s2k ;//+ (comp) ii*1.0e-5;
        comp M2 = M2kfunc(a,s2kcomp,m,1.0e-3);
        cout<<setprecision(10)<<s2k<<'\t'<<real(M2)<<'\t'<<imag(M2)<<endl;
    }*/

    /*for(double s2k=0;s2k<=6;s2k=s2k+0.0005)
    {
        comp s2kcomp = (comp) s2k ;//+ (comp) ii*1.0e-5;
        comp M2 = Jfunc(s2k);
        cout<<setprecision(10)<<s2k<<'\t'<<real(M2)<<'\t'<<imag(M2)<<endl;
    }*/
    /*vector<comp> sigvec;

    sigma_vector_maker(sigvec,8.65,0.0,sigmax(8.65,1.0),16.0,1.0,100.0,1.0e-5);

    for(int i=0;i<sigvec.size();++i)
    {
        cout<<sigvec[i]<<'\t'<<tau(8.65,sigvec[i],m)<<endl;
    }*/


    //for(int i=0;i<10;++i)
    /*for(double s2k=0.0;s2k<=5.0;s2k=s2k+0.5)
    {
        double s = 8.6;
        comp s2kcomp = s2k + ii*1.0e-6;
        cout<<setprecision(10)<<s2k<<'\t'<<real(tau1(s,s2kcomp,m))<<'\t'<<imag(tau1(s,s2kcomp,m))<<endl;
    
        cout<<setprecision(10)<<s2k<<'\t'<<real(tau(s,s2kcomp,m))<<'\t'<<imag(tau(s,s2kcomp,m))<<endl;
    }*/

    /*for(int i=0;i<10;++i)
    {
        comp K = kernel_1(svec[i],sigb,sigb,a,m,1.0e-5);
        cout<<setprecision(10)<<svec[i]<<'\t'<<real(K)<<'\t'<<imag(K)<<endl;
    }*/

    //double somenum = -1.0;
    //cout<<mysqrt(somenum)<<'\t'<<mysqrt1(somenum)<<endl;

    //kernelplotcontour_1();
    //tauplotcontour();
    //M2plotcontour();
    double s = 8.65;
    /*vector<comp> sigvec;
    sigma_vector_maker_1(sigvec,s,0.0,sigmax(s,m),a,m,100.0,1.0e-6);

    for(int i=0;i<sigvec.size();++i)
    {
        cout<<"sigvec["<<i<<"]:"<<sigvec[i]<<endl;
        cout<<"tau:"<<tau(s,sigvec[i],1.0)<<'\t'<<tau1(s,sigvec[i],1.0)<<endl;
        //cout<<"M2k:"<<M2kfunc(a,sigvec[i],m,1.0e-3)<<endl;
    }*/

    //cout<<M2kfunc(a,3.72534 + ii*0.0423318,m,1.0e-3)<<endl;
    
    /*comp s2k = 3.8 + ii*1.0e-4;
    comp s2p = sigb;
    cout<<"M2:"<<M2kfunc(16.0,s2k,1.0,1.0e-5)<<endl;
    cout<<"tau:"<<tau(8.65,s2k,1.0)<<endl;
    cout<<"Gs:"<<GS(8.65,sigb,s2k,1.0,1.0e-5)<<endl;

    cout<<"kernel:"<<kernel(8.65,sigb,s2k,16.0,1.0,1.0e-5)<<endl;*/

    vector<comp> sigvec;

    sigma_vector_maker_2(sigvec,8.6,0.0,sigmax(8.6,1.0),16,1.0,100,20,1.0e-5);
    for(int i=0;i<sigvec.size();++i)
    {
        cout<<sigvec[i]<<endl;
    }

    cout<<sigc1(8.6,sigb,1.0)<<endl;
    cout<<sigc2(8.6,sigb,1.0)<<endl;
    cout<<radius(8.6,sigb,1.0)<<endl;
    return 0;    
}