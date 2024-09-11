//(not valied anymore) macos:g++-11 printer.cpp -o printer -O3 -std=c++14 -I ~/homebrew/include/eigen3/ -L ~/homebrew/lib/

//macos new:g++-11 printer_cpu.cpp -o printer -O3 -std=c++11 -I /usr/local/include/eigen3/
//wahab:g++ printer.cpp -o printer -O3 -std=c++14 -I /cm/shared/applications/eigen/3.3.7/include/eigen3/
//wahab_gpu:nvcc -I /cm/shared/applications/cuda-toolkit/10.0.130/include -I /cm/shared/applications/eigen/3.3.7/include/eigen3/ printer.cpp -o printer_gpu -O3 -std=c++14 -L /cm/shared/applications/cuda-toolkit/10.0.130/lib64 -lcublas -lcusolver -lcudart
#include<bits/stdc++.h>
#include "functions_invariant_based.h"
#include "invariant_vector_maker.h"
#include "integralequation_invariant_based.h"
#include "interpolator.h"
#include "solvers.h"

using namespace std;

typedef complex<double> comp;

void Mphib_belowthreshold_vs_s3()
{
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.77;//8.72;//real(phibthreshold(a,m));//
    double sfinal = 8.79;real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/50.0;
    double points1 = 1000.0;
    double points2 = 500.0;
    double eps = 1.0e-7;
    double box = 10.0;

    string filename="justtestMphib_a_" + to_string((int)a) + "_N1_" + to_string((int)points1) + "_N2_" + to_string((int)points2) + "_epspow_" + to_string((int)abs(log10(eps))) + "_box_" + to_string((int)box) + ".dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    double eta = 20.0;


    for(double s=sinitial;s<=sfinal;s=s+dels)
    { 
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        //cout<<"s:"<<s<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> sigvec;

        //sigma_vector_maker_2(sigvec,s,sigmamin,sigmamax,a,m,points1,points2,1.0e-10);
        //sigma_vector_maker_abovethreshold(sigvec,s,sigmamin,sigmamax,points);
        double eps1 = abs((1.0e-2)*eps*eps_energy_factor_minus(s,sigb,m));
        cout<<"eps given to sigvec_maker = "<<eps1<<endl;
        sigma_vector_maker_6(sigvec,s,sigmamin,sigmamax,a,m,eps1,points1,points2,box);
    
        int size = sigvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker(Gvec,s,sigvec,sigb,m,eps);
        Gvec = -1.0*Gvec;

        Bmat_maker(Bmat,s,sigvec,sigvec,a,m,eps);

        double relerror;

        

        LinearSolver(Bmat,dsol,Gvec,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

        comp result;

        interpolator_ds_integraleq(dsol,sigvec,s,sigb,sigb,a,m,eps,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        result = gsq*result;
        //result = rhophib(s,a,m)*result;

        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}


void dSqqs2k2msq_belowthreshold_vs_s3()
{
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.72;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/50.0;
    double points1 = 1000.0;
    double points2 = 500.0;
    double eps = 0.01;
    double box = 10.0;

    string filename="dSqqs2q2msq_a_" + to_string((int)a) + "_N_" + to_string((int)points1) + "_eps_" + to_string(eps) + ".dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    double eta = 20.0;


    for(double s=sinitial;s<=sfinal;s=s+dels)
    { 
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        //cout<<"s:"<<s<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> sigvec;

        //sigma_vector_maker_2(sigvec,s,sigmamin,sigmamax,a,m,points1,points2,1.0e-10);
        //sigma_vector_maker_abovethreshold(sigvec,s,sigmamin,sigmamax,points);
        double eps1 = abs((1.0e-2)*eps*eps_energy_factor_minus(s,sigb,m));
        cout<<"eps given to sigvec_maker = "<<eps1<<endl;
        sigma_vector_maker_6(sigvec,s,sigmamin,sigmamax,a,m,eps1,points1,points2,box);
    
        int size = sigvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker(Gvec,s,sigvec,sigb,m,eps);
        Gvec = -1.0*Gvec;

        Bmat_maker(Bmat,s,sigvec,sigvec,a,m,eps);

        double relerror;

        

        LinearSolver(Bmat,dsol,Gvec,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

        comp result;

        interpolator_ds_integraleq(dsol,sigvec,s,sigb,sigb,a,m,eps,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        result = gsq*result;
        //result = rhophib(s,a,m)*result;

        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}


void Mphib_belowthreshold_vs_N1()
{
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.72;//8.72;//real(phibthreshold(a,m));//
    double sfinal = 8.79;real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/50.0;
    double points1 = 1000.0;
    double points2 = 200.0;
    double eps = 1.0e-5;

    //string filename="justtestMphib_a_" + to_string(a) + "_N1_" + to_string(points1) + "_N2_" + to_string(points2) + "_epspow_" + to_string((int)abs(log10(eps))) + "_nobuffer.dat";
    
    string filename="N1test_Mphib_a_" + to_string(a) + "_N2_" + to_string(points2) + "_epspow_" + to_string((int)abs(log10(eps))) + ".dat";
    
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    double eta = 20.0;

    double s = sinitial;

    //for(double s=sinitial;s<=sfinal;s=s+dels)
    for(int matsizeN1=100;matsizeN1<=2000;matsizeN1=matsizeN1+100)
    { 
        points1 = matsizeN1;
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        //cout<<"s:"<<s<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> sigvec;

        //sigma_vector_maker_2(sigvec,s,sigmamin,sigmamax,a,m,points1,points2,1.0e-10);
        //sigma_vector_maker_abovethreshold(sigvec,s,sigmamin,sigmamax,points);
        double eps1 = abs((1.0e-2)*eps*eps_energy_factor_minus(s,sigb,m));
        cout<<"eps given to sigvec_maker = "<<eps1<<endl;
        sigma_vector_maker_5(sigvec,s,sigmamin,sigmamax,a,m,eps1,points1,points2,10.0);
    
        int size = sigvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker(Gvec,s,sigvec,sigb,m,eps);
        Gvec = -1.0*Gvec;

        Bmat_maker(Bmat,s,sigvec,sigvec,a,m,eps);

        double relerror;

        

        LinearSolver(Bmat,dsol,Gvec,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

        comp result;

        interpolator_ds_integraleq(dsol,sigvec,s,sigb,sigb,a,m,eps,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        result = gsq*result;
        //result = rhophib(s,a,m)*result;

        cout<<"N1:"<<points1<<endl;
        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        //fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        fout<<s<<'\t'<<points1<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}

void Mphib_belowthreshold_vs_N2()
{
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.77;//8.72;//real(phibthreshold(a,m));//
    double sfinal = 8.79;real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/50.0;
    double points1 = 1000.0;
    double points2 = 200.0;
    double eps = 1.0e-4;
    double box = 100;

    //string filename="justtestMphib_a_" + to_string(a) + "_N1_" + to_string(points1) + "_N2_" + to_string(points2) + "_epspow_" + to_string((int)abs(log10(eps))) + "_nobuffer.dat";
    

    string filename =   "N2test_Mphib_a_" + to_string(a) 
                        + "_N1_" + to_string((int)points1) 
                        + "_epspow_" + to_string((int)abs(log10(eps))) 
                        + "_box_" + to_string((int) box) 
                        + ".dat";
    
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    double eta = 20.0;

    double s = sinitial;

    //for(double s=sinitial;s<=sfinal;s=s+dels)
    for(int matsizeN1=10;matsizeN1<=500;matsizeN1=matsizeN1+10)
    { 
        points2 = matsizeN1;
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        //cout<<"s:"<<s<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> sigvec;

        //sigma_vector_maker_2(sigvec,s,sigmamin,sigmamax,a,m,points1,points2,1.0e-10);
        //sigma_vector_maker_abovethreshold(sigvec,s,sigmamin,sigmamax,points);
        double eps1 = abs((1.0e-3)*eps*eps_energy_factor_minus(s,sigb,m));
        cout<<"eps given to sigvec_maker = "<<eps1<<endl;
        sigma_vector_maker_6(sigvec,s,sigmamin,sigmamax,a,m,eps1,points1,points2,box);
    
        int size = sigvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker(Gvec,s,sigvec,sigb,m,eps);
        Gvec = -1.0*Gvec;

        Bmat_maker(Bmat,s,sigvec,sigvec,a,m,eps);

        double relerror;

        

        LinearSolver(Bmat,dsol,Gvec,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

        comp result;

        interpolator_ds_integraleq(dsol,sigvec,s,sigb,sigb,a,m,eps,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        result = gsq*result;
        //result = rhophib(s,a,m)*result;

        cout<<"N2:"<<points2<<endl;
        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        //fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        fout<<s<<'\t'<<6.0*points2<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}

void Mphib_belowthreshold_vs_N2_vs_box()
{
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.72;//8.85;//real(phibthreshold(a,m));//
    double sfinal = 8.79;real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/50.0;
    double points1 = 1000.0;
    double points2 = 200.0;
    double eps = 1.0e-7;
    double s = sinitial;
    //double box = 1000000;

    //string filename="justtestMphib_a_" + to_string(a) + "_N1_" + to_string(points1) + "_N2_" + to_string(points2) + "_epspow_" + to_string((int)abs(log10(eps))) + "_nobuffer.dat";
    
    
    for(double box=1.0;box<=1000.0;box=box*10.0)
    {
    string filename =   "N2test_Mphib_a_" + to_string(a) 
                        + "_s_" + to_string((float)s)
                        + "_N1_" + to_string((int)points1) 
                        + "_epspow_" + to_string((int)abs(log10(eps))) 
                        + "_box_" + to_string((int) box) 
                        + "_pow2eneps.dat";
    
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    double eta = 20.0;

    

    //for(double s=sinitial;s<=sfinal;s=s+dels)
    for(int matsizeN1=300;matsizeN1<=500;matsizeN1=matsizeN1+10)
    { 
        points2 = matsizeN1;
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        //cout<<"s:"<<s<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> sigvec;

        //sigma_vector_maker_2(sigvec,s,sigmamin,sigmamax,a,m,points1,points2,1.0e-10);
        //sigma_vector_maker_abovethreshold(sigvec,s,sigmamin,sigmamax,points);
        
        /* this is 10^{-2} times smaller than the
           imaginary part of the sigma_+ */
        double eps1 = abs((1.0e-2)*eps*eps_energy_factor_minus(s,sigb,m)); 
        cout<<"delta1:"<<eps*eps_energy_factor_minus(s,sigb,m)<<'\t'<<"imagSigPlus:"<<imag(sigpplus_witheps(s,sigb,m,eps))<<endl;
        //double eps1 = abs(0.5*eps*eps_energy_factor_minus(s,sigb,m)); 
        
        cout<<"eps given to sigvec_maker = "<<eps1<<endl;
        sigma_vector_maker_6(sigvec,s,sigmamin,sigmamax,a,m,eps1,points1,points2,box);
    
        int size = sigvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker(Gvec,s,sigvec,sigb,m,eps);
        Gvec = -1.0*Gvec;

        Bmat_maker(Bmat,s,sigvec,sigvec,a,m,eps);

        double relerror;

        

        LinearSolver(Bmat,dsol,Gvec,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

        comp result;

        interpolator_ds_integraleq(dsol,sigvec,s,sigb,sigb,a,m,eps,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        result = gsq*result;
        //result = rhophib(s,a,m)*result;

        cout<<"N2:"<<points2<<endl;
        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        //fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        fout<<setprecision(10)<<s<<'\t'<<5.0*points2<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();
    cout<<"save file:"<<filename<<endl;
    }

}

void Mphib_belowthreshold_vs_eps()
{
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.72;//8.72;//real(phibthreshold(a,m));//
    double sfinal = 8.79;real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/50.0;
    double points1 = 1000.0;
    double points2 = 200.0;
    //double eps = 1.0e-5;

    //string filename="justtestMphib_a_" + to_string(a) + "_N1_" + to_string(points1) + "_N2_" + to_string(points2) + "_epspow_" + to_string((int)abs(log10(eps))) + "_nobuffer.dat";
    
    string filename="epstest_Mphib_a_" + to_string(a) + "_N1_" + to_string(points1) + "_N2_" + to_string(points2)  + ".dat";
    
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    double eta = 20.0;

    double s = sinitial;

    //for(double s=sinitial;s<=sfinal;s=s+dels)
    //for(int matsizeN1=10;matsizeN1<=300;matsizeN1=matsizeN1+10)]
    for(double eps=0.01;eps>=1.0e-8;eps=eps/2.0)
    { 
        //points2 = matsizeN1;
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        //cout<<"s:"<<s<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> sigvec;

        //sigma_vector_maker_2(sigvec,s,sigmamin,sigmamax,a,m,points1,points2,1.0e-10);
        //sigma_vector_maker_abovethreshold(sigvec,s,sigmamin,sigmamax,points);
        double eps1 = abs((1.0e-2)*eps*eps_energy_factor_minus(s,sigb,m));
        cout<<"eps given to sigvec_maker = "<<eps1<<endl;
        sigma_vector_maker_5(sigvec,s,sigmamin,sigmamax,a,m,eps1,points1,points2,10.0);
    
        int size = sigvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker(Gvec,s,sigvec,sigb,m,eps);
        Gvec = -1.0*Gvec;

        Bmat_maker(Bmat,s,sigvec,sigvec,a,m,eps);

        double relerror;

        

        LinearSolver(Bmat,dsol,Gvec,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

        comp result;

        interpolator_ds_integraleq(dsol,sigvec,s,sigb,sigb,a,m,eps,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        result = gsq*result;
        //result = rhophib(s,a,m)*result;

        cout<<"eps:"<<eps<<endl;
        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        //fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        fout<<s<<'\t'<<eps<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}

void Mphib_belowthreshold_vs_contourbox()
{
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.72;//8.72;//real(phibthreshold(a,m));//
    double sfinal = 8.79;real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/50.0;
    double points1 = 1000.0;
    double points2 = 200.0;
    double eps = 1.0e-5;

    //string filename="justtestMphib_a_" + to_string(a) + "_N1_" + to_string(points1) + "_N2_" + to_string(points2) + "_epspow_" + to_string((int)abs(log10(eps))) + "_nobuffer.dat";
    
    string filename="contourboxtest_Mphib_a_" + to_string(a) + "_N1_" + to_string(points1) + "_N2_" + to_string(points2) + "_epspow_" + to_string((int)abs(log10(eps))) + ".dat";
    
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    double eta = 20.0;

    double s = sinitial;

    //for(double s=sinitial;s<=sfinal;s=s+dels)
    //for(int matsizeN1=10;matsizeN1<=300;matsizeN1=matsizeN1+10)]
    //for(double eps=0.01;eps>=1.0e-8;eps=eps/2.0)
    for(double box=0.0;box<=100000.0;box=box+2000.0)
    { 
        //points2 = matsizeN1;
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        //cout<<"s:"<<s<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> sigvec;

        //sigma_vector_maker_2(sigvec,s,sigmamin,sigmamax,a,m,points1,points2,1.0e-10);
        //sigma_vector_maker_abovethreshold(sigvec,s,sigmamin,sigmamax,points);
        double eps1 = abs((1.0e-2)*eps*eps_energy_factor_minus(s,sigb,m));
        cout<<"eps given to sigvec_maker = "<<eps1<<endl;
        sigma_vector_maker_5(sigvec,s,sigmamin,sigmamax,a,m,eps1,points1,points2,box);
    
        int size = sigvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker(Gvec,s,sigvec,sigb,m,eps);
        Gvec = -1.0*Gvec;

        Bmat_maker(Bmat,s,sigvec,sigvec,a,m,eps);

        double relerror;

        

        LinearSolver(Bmat,dsol,Gvec,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

        comp result;

        interpolator_ds_integraleq(dsol,sigvec,s,sigb,sigb,a,m,eps,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        result = gsq*result;
        //result = rhophib(s,a,m)*result;

        cout<<"box:"<<box<<endl;
        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        //fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        fout<<setprecision(10)<<s<<'\t'<<box<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        
        count = count + 1;
        cout<<"-------------------------"<<endl;
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

            LinearSolver(Bmat,dsol,Gvec,relerror);
            //cusolverComplex(Bmat,Gvec,dsol,size);

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


/*void Mphib_belowthreshold_vs_eps()
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
*/

void Mphib_belowthreshold_vs_complex_s3()
{
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.77;//8.72;//real(phibthreshold(a,m));//
    double sfinal = 8.79;//real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/50.0;
    double points1 = 1000.0;
    double points2 = 500.0;
    double eps = 0.0;
    double box = 0.0;
    double imS = -0.01;

    string filename=    "MphibComplexS3_a_" + to_string((int)a) 
                        + "_N1_" + to_string((int)points1) 
                        + "_N2_" + to_string((int)points2) 
                        + "_imS_" + to_string(imS) 
                        + "_eps_" + to_string((int)eps)
                        + "_box_" + to_string((int)box) + ".dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    double eta = 20.0;


    for(double s=sinitial;s<=sfinal;s=s+dels)
    { 
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        comp scomp = s + ii*imS;
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(scomp,m);
        comp sigb = sigmab(a,m);
        //cout<<"s:"<<s<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> sigvec;

        //sigma_vector_maker_2(sigvec,s,sigmamin,sigmamax,a,m,points1,points2,1.0e-10);
        //sigma_vector_maker_abovethreshold(sigvec,s,sigmamin,sigmamax,points);
        double eps1 = 0.0;//abs((1.0e-2)*eps*eps_energy_factor_minus(s,sigb,m));
        cout<<"eps given to sigvec_maker = "<<eps1<<endl;
        sigma_vector_maker_8(sigvec,scomp,sigmamin,sigmamax,a,m,eps1,points1,points2,box);
    
        int size = sigvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker(Gvec,scomp,sigvec,sigb,m,eps);
        Gvec = -1.0*Gvec;

        Bmat_maker(Bmat,scomp,sigvec,sigvec,a,m,eps);

        double relerror;

        

        LinearSolver(Bmat,dsol,Gvec,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

        comp result;

        interpolator_ds_integraleq(dsol,sigvec,scomp,sigb,sigb,a,m,eps,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        result = gsq*result;
        //result = rhophib(s,a,m)*result;

        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        fout<<real(scomp)<<'\t'imag(scomp)<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}



int main()
{
    double a = 16.0;
    double m = 1.0;

    cout<<phibthreshold(a,m)<<endl;
    cout<<Gleftbranchpoint(a,m)<<endl;
    cout<<Grightbranchpoint(a,m)<<endl;

    comp ii = {0.0,1.0};
    comp s = 8.72 ;
    comp sigb = sigmab(a,m);

    cout<<"sigp:"<<sigpplus_witheps(s,sigb,m,0.0)<<endl;
    cout<<"sigm:"<<sigpminus_witheps(s,sigb,m,0.0)<<endl;
    cout<<"maxInt:"<<pow(sqrt(s)-m,2)<<endl;

    //Mphib_belowthreshold_vs_eps();
    //Mphib_belowthreshold_vs_N1();
    //Mphib_belowthreshold_vs_N2();
    //Mphib_belowthreshold_vs_N2_vs_box();
    //Mphib_belowthreshold_vs_eps();
    //Mphib_belowthreshold_vs_contourbox();
    //Mphib_belowthreshold_vs_s3();
    Mphib_belowthreshold_vs_complex_s3();

    return 0;
}
