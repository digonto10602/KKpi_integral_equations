//macos:g++-11 printer.cpp -o printer -O3 -std=c++14 -I ~/homebrew/include/eigen3/ -L ~/homebrew/lib/
//wahab:g++ printer.cpp -o printer -O3 -std=c++14 -I /cm/shared/applications/eigen/3.3.7/include/eigen3/
//wahab_gpu:nvcc -I /cm/shared/applications/cuda-toolkit/10.0.130/include -I /cm/shared/applications/eigen/3.3.7/include/eigen3/ printer.cpp -o printer_gpu -O3 -std=c++14 -L /cm/shared/applications/cuda-toolkit/10.0.130/lib64 -lcublas -lcusolver -lcudart
#include<bits/stdc++.h>
#include "functions_invariant_based.h"
#include "invariant_vector_maker.h"
#include "integralequation_invariant_based.h"
#include "interpolator.h"
#include "solvers_gpu.h"

using namespace std;

typedef complex<double> comp;

void Mphib_belowthreshold_vs_s3()
{
    double pi = acos(-1.0);
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.7;
    double sfinal = 8.95;
    double dels = abs(sinitial-sfinal)/50.0;
    double points = 1000.0;
    double eps = 1.0e-5;//2.0*pi/points;

    string filename="Mphib_a_" + to_string(a) + "_N_" + to_string(points) + "_eps_" + to_string(eps) + ".dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;


    for(double s=sinitial;s<=sfinal+0.05;s=s+dels)
    {
        if(s>=real(phibthreshold(a,m))) eps = 2.0*pi/points;
        if(s>=sfinal) dels = 0.05/100.0;
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);

        vector<comp> sigvec;

        sigma_vector_maker(sigvec,s,sigmamin,sigmamax,a,m,points,1.0e-8);
    
        int size = sigvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker(Gvec,s,sigvec,sigb,m,eps);

        Bmat_maker(Bmat,s,sigvec,sigvec,a,m,eps);

        double relerror;

        

        //LinearSolver(Bmat,dsol,Gvec,relerror);
        cusolverComplex(Bmat,Gvec,dsol,size);

        comp result;

        interpolator_ds_integraleq(dsol,sigvec,s,sigb,sigb,a,m,eps,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        result = gsq*result;

        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        count = count + 1;
    }

    fout.close();

}


void Mphib_belowthreshold_vs_s3_vs_L()
{
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.6;
    double sfinal = 9.0;
    double dels = abs(sinitial-sfinal)/100.0;
    double points; 
    double eps = 1.0e-3;

    


    for(double s=sinitial;s<=sfinal;s=s+dels)
    {
        string filename="Mphib_a_" + to_string(a) + "_s3_" + to_string(s) +  ".dat";
        ofstream fout;
        fout.open(filename.c_str());

        for(double L=1000.0;L<=5000.0;L=L+500.0)
        {
            double pi = acos(-1.0);
            points = L;
            eps = 2.0*pi/L;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = sigmab(a,m);

            vector<comp> sigvec;

            sigma_vector_maker(sigvec,s,sigmamin,sigmamax,a,m,points,eps);
    
            int size = sigvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            Gvec_maker(Gvec,s,sigvec,sigb,m,eps);

            Bmat_maker(Bmat,s,sigvec,sigvec,a,m,eps);

            double relerror;

            //LinearSolver(Bmat,dsol,Gvec,relerror);
            cusolverComplex(Bmat,Gvec,dsol,size);

            comp result;

            interpolator_ds_integraleq(dsol,sigvec,s,sigb,sigb,a,m,eps,result);

            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            result = gsq*result;

            //cout<<"res:"<<result<<endl;
            cout<<"s3:"<<s<<'\t'
                <<"dels:"<<dels<<'\t'
                <<"L:"<<L<<'\t'
                <<"res:"<<result<<endl;

            cout<<"------------------------"<<endl;

            fout<<s<<'\t'<<points<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        }

        fout.close();
    }

    

}

void Mphib_belowthreshold_vs_eps()
{
    double a = 16.0;
    double m = 1.0;

    double s = 8.89;
    double sinitial = 8.7;
    double sfinal = 8.95;
    double dels = abs(sinitial-sfinal)/100.0;
    double points = 5000.0;
    //double eps = 1.0e-5;

    string filename="MphibEpstest_a_" + to_string(a) + "_N_" + to_string(points)  + ".dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;


    //for(double s=sinitial;s<=sfinal+0.05;s=s+dels)
    for(double eps=1.0e-5;eps<=1.0e-2;eps=eps+100)
    {
        if(s>=sfinal) dels = 0.05/100.0;
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);

        vector<comp> sigvec;

        sigma_vector_maker(sigvec,s,sigmamin,sigmamax,a,m,points,eps*10);
    
        int size = sigvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker(Gvec,s,sigvec,sigb,m,eps);

        Bmat_maker(Bmat,s,sigvec,sigvec,a,m,eps);

        double relerror;

        

        LinearSolver(Bmat,dsol,Gvec,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

        comp result;

        interpolator_ds_integraleq(dsol,sigvec,s,sigb,sigb,a,m,eps,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        result = gsq*result;

        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"eps:"<<eps<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        fout<<eps<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        count = count + 1;
    }

    fout.close();

}

int main()
{
    //Mphib_belowthreshold_vs_s3();
    Mphib_belowthreshold_vs_eps();

    return 0;
}
