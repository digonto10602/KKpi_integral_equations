
//macos:g++-11 printer.cpp -o printer -O3 -std=c++14 -I ~/homebrew/include/eigen3/ -L ~/homebrew/lib/
//wahab:g++ printer.cpp -o printer -O3 -std=c++14 -I /cm/shared/applications/eigen/3.3.7/include/eigen3/
//wahab_gpu:nvcc -I /cm/shared/applications/cuda-toolkit/10.0.130/include -I /cm/shared/applications/eigen/3.3.7/include/eigen3/ printer_gpu_momrep.cpp  -O3 -std=c++14 -L /cm/shared/applications/cuda-toolkit/10.0.130/lib64 -lcublas -lcusolver -lcudart -o printer_gpu
//wahab_gpu_omp:nvcc -I /cm/shared/applications/cuda-toolkit/10.0.130/include -I /cm/shared/applications/eigen/3.3.7/include/eigen3/ printer_gpu_momrep.cpp  -O3 -std=c++14 -L /cm/shared/applications/cuda-toolkit/10.0.130/lib64 -lcublas -lcusolver -lcudart -lgomp -Xcompiler -fopenmp -o printer_gpu
//ubuntupc_omp:nvcc117 -I/usr/include/eigen3/ -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/include -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/lib64 printer_gpu_momrep.cpp  -O3 -std=c++14 -lcublas -lcusolver -lcudart -lgomp -Xcompiler -fopenmp -o printer_gpu
//ubuntupc_omp1:nvcc -I/usr/include/eigen3/ -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/include -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/lib64 ../printer_gpu_momrep_from_ac_final.cpp -O3 -std=c++14 -lcublas -lcusolver -lcudart -lgomp -Xcompiler -fopenmp -o printer_gpu
#include<bits/stdc++.h>
#include "functions_momrep_based.h"
#include "momentum_vector_maker.h"
#include "integralequation_momrep_based.h"
#include "interpolator_momrep.h"
#include "solvers_gpu.h"
#include<omp.h>
#include <Eigen/Eigenvalues> //for vertex factor calculations
#include "LUSolver_modified.h"
#include "fastGL_based_functions.h"
#include "SA_method_functions.h"
#include <sys/time.h>

using namespace std;

typedef complex<double> comp;

void Mphib_belowthreshold_vs_s3()
{
    double a = 500.0;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;
    double sinitial = 8.998;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/500.0;//abs(8.70275246887-8.70115128532);//
    double points1 = 2000.0;
    //double points2 = 500.0;
    double eps = 1.0e-4;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.01;
    //double box = 10.0;

    string filename="Mphib_momrep_lastpollsection_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_epspow_" + to_string((int)abs(log10(eps))) 
                    + ".dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    //double eta = 20.0;


    for(double s=sinitial;s<=sfinal;s=s+dels)
    { 
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        cout<<"am = "<<a<<endl;
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp qval = pmom(s,sigb,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        
        //cout<<"s:"<<s<<endl;
        cout<<"run: "<<count+1<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;

        mom_vector_maker_4(qvec,s,kmin,kmax,a,m,eps,eps_for_m2k,points1,qvec_r);
        cout<<"qvec size = "<<qvec.size()<<endl;
    
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec = -1.0*Gvec;

        Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);

        double relerror;

        

        //LinearSolver_2(Bmat,dsol,Gvec,relerror);
        cusolverComplex(Bmat,Gvec,dsol,size);

        comp result;

        interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        result = gsq*result;

        
        comp rhopb = rhophib(s,a,m);
        comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
        comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
        comp GSval = GS_pk(s,qval,qval,m,eps);
        //result = rhophib(s,a,m)*result;

        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        fout<<setprecision(16)<<s<<'\t'
            <<real(result)<<'\t'
            <<imag(result)<<'\t'
            <<real(1.0/result)<<'\t'
            <<imag(1.0/result)<<'\t'
            <<real(result + GSval)<<'\t'
            <<imag(result + GSval)<<'\t'
            <<real(1.0/(result + GSval))<<'\t'
            <<imag(1.0/(result + GSval))<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}


void Mphib_belowthreshold_vs_s3_fixed_s3imag_contour45()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    double s3imag = 0.05;
    double sinitial = 8.72;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;

    double delspoints = 100.0;
    double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
    double points1 = 5000.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.01;
    //double box = 10.0;

    string filename="fixeds3imag_Mphib_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_contour44.dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    //double eta = 20.0;

    for(int i=0;i<delspoints+1;++i)
    //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    { 
        double s3real = sinitial + i*dels;
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);

        comp s = s3real + ii*s3imag;
        cout<<"am = "<<a<<endl;
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp qval = pmom(s,sigb,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        
        //cout<<"s:"<<s<<endl;
        cout<<"run: "<<count+1<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;

        int tag1=0,tag2=0;
            
            
        if(s3imag<=0.0)
        //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        else if(s3imag>0.0)
        {
            //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            //mom_vector_maker_45(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
            cout<<"tag1 = "<<tag1<<endl;
            cout<<"tag2 = "<<tag2<<endl;
                
                
        }

        cout<<"qvec created with size = "<<qvec.size()<<endl;
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

            
        if(s3imag<=0.0)
        {
            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        }
        else if(s3imag>0.0)
        {
            Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
        }
        //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
        Gvec = -1.0*Gvec;

        //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
        Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
        //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);

            
        //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
        double relerror;
            
        
        //cout<<"did this = "<<count<<" times"<<endl;
        //cout<<"i val = "<<i<<'\t';
        //cout<<"j val = "<<j<<endl;
            
        //LinearSolver_2(Bmat,dsol,Gvec,relerror);

        time(&time_start);
        cusolverComplex(Bmat,Gvec,dsol,size);
        time(&time_end);

        double time_taken = double((double)time_end - (double)time_start);
        cout<<"Time taken to solve : "<<fixed 
            <<setprecision(5)<<time_taken;
        cout<<" sec"<<endl;


        comp result;

        
        if(s3imag<=0.0)
        {
            interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
        }
        else if(s3imag>0.0)
        {
            interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
        }
        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        comp ds = result;
        result = gsq*result;

        
        comp rhopb = rhophib(s,a,m);
        comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
        comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
        comp GSval = GS_pk(s,qval,qval,m,eps);
        //result = rhophib(s,a,m)*result;

        cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        fout<<setprecision(16)<<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<real(result)<<'\t'
            <<imag(result)<<'\t'
            <<real(ds)<<'\t'
            <<imag(ds)<<'\t'
            <<real(mphib2)<<'\t'
            <<imag(mphib2)<<'\t'
            <<real(mphib2denom)<<'\t'
            <<imag(mphib2denom)<<'\t'
            <<a<<'\t'
            <<N<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}

void Mphib_belowthreshold_vs_s3_fixed_s3imag_contour46()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    double s3imag = 0.05;
    double sinitial = 8.72;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;

    double delspoints = 100.0;
    double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
    double points1 = 5000.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.01;
    //double box = 10.0;

    string filename="fixeds3imag_Mphib_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_contour46.dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    //double eta = 20.0;

    for(int i=0;i<delspoints+1;++i)
    //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    { 
        double s3real = sinitial + i*dels;
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);

        comp s = s3real + ii*s3imag;
        cout<<"am = "<<a<<endl;
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp qval = pmom(s,sigb,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        
        //cout<<"s:"<<s<<endl;
        cout<<"run: "<<count+1<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;

        int tag1=0,tag2=0;
            
            
        if(s3imag<=0.0)
        //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        else if(s3imag>0.0)
        {
            //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            //mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            mom_vector_maker_46(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
            cout<<"tag1 = "<<tag1<<endl;
            cout<<"tag2 = "<<tag2<<endl;
                
                
        }

        cout<<"qvec created with size = "<<qvec.size()<<endl;
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

            
        if(s3imag<=0.0)
        {
            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        }
        else if(s3imag>0.0)
        {
            Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
        }
        //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
        Gvec = -1.0*Gvec;

        //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
        Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
        //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);

            
        //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
        double relerror;
            
        
        //cout<<"did this = "<<count<<" times"<<endl;
        //cout<<"i val = "<<i<<'\t';
        //cout<<"j val = "<<j<<endl;
            
        //LinearSolver_2(Bmat,dsol,Gvec,relerror);

        time(&time_start);
        cusolverComplex(Bmat,Gvec,dsol,size);
        time(&time_end);

        double time_taken = double(time_end - time_start);
        cout<<"Time taken to solve : "<<fixed 
            <<time_taken<<setprecision(5);
        cout<<" sec"<<endl;


        comp result;

        
        if(s3imag<=0.0)
        {
            interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
        }
        else if(s3imag>0.0)
        {
            interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
        }
        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        comp ds = result;
        result = gsq*result;

        
        comp rhopb = rhophib(s,a,m);
        comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
        comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
        comp GSval = GS_pk(s,qval,qval,m,eps);
        //result = rhophib(s,a,m)*result;

        cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        fout<<setprecision(16)<<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<real(result)<<'\t'
            <<imag(result)<<'\t'
            <<real(ds)<<'\t'
            <<imag(ds)<<'\t'
            <<real(mphib2)<<'\t'
            <<imag(mphib2)<<'\t'
            <<real(mphib2denom)<<'\t'
            <<imag(mphib2denom)<<'\t'
            <<a<<'\t'
            <<N<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}

void Mphib_belowthreshold_vs_s3_fixed_s3imag_contour47()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    double s3imag = 0.05;
    double sinitial = 8.72;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;

    double delspoints = 100.0;
    double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
    double points1 = 5000.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.01;
    //double box = 10.0;

    string filename="fixeds3imag_Mphib_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_contour47.dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    //double eta = 20.0;

    for(int i=0;i<delspoints+1;++i)
    //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    { 
        double s3real = sinitial + i*dels;
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);

        comp s = s3real + ii*s3imag;
        cout<<"am = "<<a<<endl;
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp qval = pmom(s,sigb,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        
        //cout<<"s:"<<s<<endl;
        cout<<"run: "<<count+1<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;

        int tag1=0,tag2=0;
            
            
        if(s3imag<=0.0)
        //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        else if(s3imag>0.0)
        {
            //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            //mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
            cout<<"tag1 = "<<tag1<<endl;
            cout<<"tag2 = "<<tag2<<endl;
                
                
        }

        cout<<"qvec created with size = "<<qvec.size()<<endl;
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

            
        if(s3imag<=0.0)
        {
            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        }
        else if(s3imag>0.0)
        {
            Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
        }
        //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
        Gvec = -1.0*Gvec;

        //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
        Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
        //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);

            
        //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
        double relerror;
            
        
        //cout<<"did this = "<<count<<" times"<<endl;
        //cout<<"i val = "<<i<<'\t';
        //cout<<"j val = "<<j<<endl;
            
        //LinearSolver_2(Bmat,dsol,Gvec,relerror);

        time(&time_start);
        cusolverComplex(Bmat,Gvec,dsol,size);
        time(&time_end);

        double time_taken = double(time_end - time_start);
        cout<<"Time taken to solve : "<<fixed 
            <<time_taken<<setprecision(5);
        cout<<" sec"<<endl;


        comp result;

        
        if(s3imag<=0.0)
        {
            interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
        }
        else if(s3imag>0.0)
        {
            interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
        }
        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        comp ds = result;
        result = gsq*result;

        
        comp rhopb = rhophib(s,a,m);
        comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
        comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
        comp GSval = GS_pk(s,qval,qval,m,eps);
        //result = rhophib(s,a,m)*result;

        cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        fout<<setprecision(16)<<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<real(result)<<'\t'
            <<imag(result)<<'\t'
            <<real(ds)<<'\t'
            <<imag(ds)<<'\t'
            <<real(mphib2)<<'\t'
            <<imag(mphib2)<<'\t'
            <<real(mphib2denom)<<'\t'
            <<imag(mphib2denom)<<'\t'
            <<a<<'\t'
            <<N<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}

void Mphib_belowthreshold_vs_s3_fixed_s3imag_contour43()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    double s3imag = -0.0001;
    double sinitial = 8.5;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;

    double delspoints = 100.0;
    double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
    double points1 = 500.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.01;
    //double box = 10.0;

    string filename="fixeds3imag_Mphib_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_contour43_1.dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    //double eta = 20.0;

    for(int i=0;i<delspoints+1;++i)
    //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    { 
        double s3real = sinitial + i*dels;
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);

        comp s = s3real + ii*s3imag;
        cout<<"am = "<<a<<endl;
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp qval = pmom(s,sigb,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        
        //cout<<"s:"<<s<<endl;
        cout<<"run: "<<count+1<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;

        int tag1=0,tag2=0;
            
            
        if(s3imag<=0.0)
        {//mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            cout<<"contour_43 chosen"<<endl;
            mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        }
        else if(s3imag>0.0)
        {
            //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            //mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            cout<<"contour_47 chosen"<<endl;
            mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
            cout<<"tag1 = "<<tag1<<endl;
            cout<<"tag2 = "<<tag2<<endl;
                
                
        }

        cout<<"qvec created with size = "<<qvec.size()<<endl;
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

            
        if(s3imag<=0.0)
        {
            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        }
        else if(s3imag>0.0)
        {
            Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
        }
        //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
        Gvec = -1.0*Gvec;

        //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
        Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
        //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);

            
        //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
        double relerror;
            
        
        //cout<<"did this = "<<count<<" times"<<endl;
        //cout<<"i val = "<<i<<'\t';
        //cout<<"j val = "<<j<<endl;
            
        //LinearSolver_2(Bmat,dsol,Gvec,relerror);

        time(&time_start);
        cusolverComplex(Bmat,Gvec,dsol,size);
        time(&time_end);

        double time_taken = double(time_end - time_start);
        cout<<"Time taken to solve : "<<fixed 
            <<time_taken<<setprecision(5);
        cout<<" sec"<<endl;


        comp result;

        
        if(s3imag<=0.0)
        {
            interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
        }
        else if(s3imag>0.0)
        {
            interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
        }
        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        comp ds = result;
        result = gsq*result;

        
        comp rhopb = rhophib(s,a,m);
        comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
        comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
        comp GSval = GS_pk(s,qval,qval,m,eps);
        //result = rhophib(s,a,m)*result;

        cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        fout<<setprecision(16)<<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<real(result)<<'\t'
            <<imag(result)<<'\t'
            <<real(ds)<<'\t'
            <<imag(ds)<<'\t'
            <<real(mphib2)<<'\t'
            <<imag(mphib2)<<'\t'
            <<real(mphib2denom)<<'\t'
            <<imag(mphib2denom)<<'\t'
            <<a<<'\t'
            <<N<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}

void Mphib_secondsheet_belowthreshold_vs_s3_fixed_s3imag_contour43_with_weights()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    double s3imag = -1.0e-7;
    double gvalleft = real(Gleftbranchpoint(a,m));
    double gvalright = real(Grightbranchpoint(a,m));
    
    double sinitial = 8.7;//gvalright;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;

    double delspoints = 500.0;
    double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
    double points1 = 900.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;//1.0e-5;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.01;
    //double box = 10.0;

    string filename="Mphib_HS_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_contour43.dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    //double eta = 20.0;

    for(int i=0;i<delspoints+1;++i)
    //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    { 
        double s3real = sinitial + i*dels;
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);

        comp s = s3real + ii*s3imag;
        cout<<"am = "<<a<<endl;
        //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp qval = pmom(s,sigb,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        
        //cout<<"s:"<<s<<endl;
        cout<<"run: "<<count+1<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;
        vector<comp> weights;

        int tag1=0,tag2=0;
            
            
        if(s3imag<=0.0)
        {//mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            cout<<"contour_43 chosen"<<endl;
            //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            mom_vector_maker_43_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        }
        else if(s3imag>0.0)
        {
            //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            //mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            cout<<"contour_47 chosen"<<endl;
            mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
            cout<<"tag1 = "<<tag1<<endl;
            cout<<"tag2 = "<<tag2<<endl;
                
                
        }

        cout<<"qvec created with size = "<<qvec.size()<<endl;
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

            
        if(s3imag<=0.0)
        {
            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        }
        else if(s3imag>0.0)
        {
            Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
        }
        //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
        Gvec = -1.0*Gvec;

        //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
        //Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
        //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);
        Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

            
        //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
        double relerror;
            
        
        //cout<<"did this = "<<count<<" times"<<endl;
        //cout<<"i val = "<<i<<'\t';
        //cout<<"j val = "<<j<<endl;
            
        //LinearSolver_2(Bmat,dsol,Gvec,relerror);

        time(&time_start);
        //LinearSolver_2(Bmat,dsol,Gvec,relerror);

        cusolverComplex(Bmat,Gvec,dsol,size);
        time(&time_end);

        double time_taken = double(time_end - time_start);
        cout<<"Time taken to solve : "<<fixed 
            <<time_taken<<setprecision(5);
        cout<<" sec"<<endl;


        comp result;

        
        if(s3imag<=0.0)
        {
            //interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
            interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
        }
        else if(s3imag>0.0)
        {
            interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
        }
        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        comp ds = result;
        result = gsq*result;

        
        comp rhopb = rhophib(s,a,m);
        comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
        comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
        comp GSval = GS_pk(s,qval,qval,m,eps);
        //result = rhophib(s,a,m)*result;

        cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        fout<<setprecision(16)<<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<real(result)<<'\t'
            <<imag(result)<<'\t'
            <<real(ds)<<'\t'
            <<imag(ds)<<'\t'
            <<real(mphib2)<<'\t'
            <<imag(mphib2)<<'\t'
            <<real(mphib2denom)<<'\t'
            <<imag(mphib2denom)<<'\t'
            <<a<<'\t'
            <<N<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}

void Mphib_secondsheet_belowthreshold_vs_s3_fixed_s3imag_contour47_with_weights()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = 10.0;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    double s3imag = 0.0;
    double gvalleft = real(Gleftbranchpoint(a,m));
    double gvalright = real(Grightbranchpoint(a,m));
    
    double sinitial = gvalright;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;

    double delspoints = 500.0;
    double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
    double points1 = 200.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 1.0e-5;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.01;
    //double box = 10.0;

    string filename="Mphib_secondsheet_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_contour43_1.dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    //double eta = 20.0;

    for(int i=0;i<delspoints+1;++i)
    //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    { 
        double s3real = sinitial + i*dels;
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);

        comp s = s3real + ii*s3imag;
        cout<<"am = "<<a<<endl;
        //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp qval = pmom(s,sigb,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        
        //cout<<"s:"<<s<<endl;
        cout<<"run: "<<count+1<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;
        vector<comp> weights;

        int tag1=0,tag2=0;
            
            
        if(s3imag<=0.0)
        {//mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            cout<<"contour_43 chosen"<<endl;
            //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            mom_vector_maker_43_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        }
        else if(s3imag>0.0)
        {
            //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            //mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            cout<<"contour_47 chosen"<<endl;
            mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
            cout<<"tag1 = "<<tag1<<endl;
            cout<<"tag2 = "<<tag2<<endl;
                
                
        }

        cout<<"qvec created with size = "<<qvec.size()<<endl;
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

            
        if(s3imag<=0.0)
        {
            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        }
        else if(s3imag>0.0)
        {
            Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
        }
        //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
        Gvec = -1.0*Gvec;

        //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
        //Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
        //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);
        Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

            
        //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
        double relerror;
            
        
        //cout<<"did this = "<<count<<" times"<<endl;
        //cout<<"i val = "<<i<<'\t';
        //cout<<"j val = "<<j<<endl;
            
        //LinearSolver_2(Bmat,dsol,Gvec,relerror);

        time(&time_start);
        LinearSolver_2(Bmat,dsol,Gvec,relerror);

        //cusolverComplex(Bmat,Gvec,dsol,size);
        time(&time_end);

        double time_taken = double(time_end - time_start);
        cout<<"Time taken to solve : "<<fixed 
            <<time_taken<<setprecision(5);
        cout<<" sec"<<endl;


        comp result;

        
        if(s3imag<=0.0)
        {
            //interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
            interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
        }
        else if(s3imag>0.0)
        {
            interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
        }
        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        comp ds = result;
        result = gsq*result;

        
        comp rhopb = rhophib(s,a,m);
        comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
        comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
        comp GSval = GS_pk(s,qval,qval,m,eps);
        //result = rhophib(s,a,m)*result;

        cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        fout<<setprecision(16)<<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<real(result)<<'\t'
            <<imag(result)<<'\t'
            <<real(ds)<<'\t'
            <<imag(ds)<<'\t'
            <<real(mphib2)<<'\t'
            <<imag(mphib2)<<'\t'
            <<real(mphib2denom)<<'\t'
            <<imag(mphib2denom)<<'\t'
            <<a<<'\t'
            <<N<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}


void Mphib_secondsheet_belowthreshold_vs_s3_a_fixed_s3imag_with_weights()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = 10.0;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    double s3imag = 0.0;//-1.0e-5;//this was 0 before
    

    double delspoints = 500.0;
    double points1 = 500.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 1.0e-5;//0.0;//this was 1.0e-5 before
    double eps_for_m2k = 0.0;
    double qvec_r = 0.101;
    //double box = 10.0;
    int acount = 250;

    double astart = 11.0;
    double afinal = 20.0;
    double dela = abs(astart - afinal)/250;
    
    for(double a=astart;a<=afinal;a=a+dela)
    {
        //double a = 3.692;
            //if(a>=3.7) dela = 0.101;
        //else if(a>=17.0) dela = 0.51;
        //else if(a>=1000.0) dela = 300.01;

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));

        double pth = real(phibthreshold(a,m));
        double firstdistance = abs(gvalright - pth);
    
        double sinitial = gvalright - firstdistance;//(gvalleft + (gvalleft+gvalright)/2.0)/2.0;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
        double sfinal = real(phibthreshold(a,m));//8.82;
        double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
    
        string filename="Mphib_secondsheet_acount_" + to_string(acount) + ".dat";
        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        //double eta = 20.0;

        for(int i=0;i<delspoints+1;++i)
        //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
        { 
            double s3real = sinitial + i*dels;
            //double s = 8.78;
            //if(s>=sfinal) dels = 0.05/100.0;
            //eps = eps_above_threshold(eta,s,a,m,points);

            comp s = s3real + ii*s3imag;
            cout<<"am = "<<a<<endl;
            //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = sigmab(a,m);
            comp qval = pmom(s,sigb,m);

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        
            //cout<<"s:"<<s<<endl;
            cout<<"run: "<<count+1<<endl;
            cout<<"a = "<<a<<'\t'<<" dela = "<<dela<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            int switch_for_gvec_fixer = 0;//0 means use the gvec_fixer_1
                                          //1 means dont use it
            
            if(s3imag<=0.0)
            {
                //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                cout<<"contour_43 chosen"<<endl;
                //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                mom_vector_maker_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1);
                //mom_vector_maker_43_with_weights_with_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                
                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
            }
            else if(s3imag>0.0)
            {
                //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                cout<<"contour_47 chosen"<<endl;
                //mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
                mom_vector_maker_seba_imspos(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                cout<<"tag1 = "<<tag1<<endl;
                cout<<"tag2 = "<<tag2<<endl;
                
                
            }

            cout<<"qvec created with size = "<<qvec.size()<<endl;
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            
            if(s3imag<=0.0)
            {
                Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            }
            else if(s3imag>0.0)
            {
                Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            }
            //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            Gvec = -1.0*Gvec;

            //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
            //Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
            //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);
            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

            
            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        
            //cout<<"did this = "<<count<<" times"<<endl;
            //cout<<"i val = "<<i<<'\t';
            //cout<<"j val = "<<j<<endl;
            
            //LinearSolver_2(Bmat,dsol,Gvec,relerror);

            time(&time_start);
            LinearSolver_2(Bmat,dsol,Gvec,relerror);

            //cusolverComplex(Bmat,Gvec,dsol,size);
            time(&time_end);

            double time_taken = double(time_end - time_start);
            cout<<"Time taken to solve : "<<fixed 
                <<time_taken<<setprecision(5);
            cout<<" sec"<<endl;


            comp result;

        
            if(s3imag<=0.0)
            {
                //interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
                interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            }
            else if(s3imag>0.0)
            {
                //interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
                interpolator_ds_integraleq_momrep_2eps_withtags_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
            }
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            comp ds = result;
            result = gsq*result;

        
            comp rhopb = rhophib(s,a,m);
            cout<<"rho phib with ims= "<<rhopb<<endl;
            cout<<"rho phib with reals="<<rhophib(real(s),a,m)<<endl;
            cout<<"rho phib with ims pos="<<rhophib(real(s) + ii*abs(imag(s)),a,m)<<endl;
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;

            cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            fout<<setprecision(16)<<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(ds)<<'\t'
                <<imag(ds)<<'\t'
                <<real(mphib2)<<'\t'
                <<imag(mphib2)<<'\t'
                <<real(mphib2denom)<<'\t'
                <<imag(mphib2denom)<<'\t'
                <<a<<'\t'
                <<N<<'\t'
                <<gvalleft<<'\t'
                <<gvalright<<'\t'
                <<real(phibthreshold(a,m))<<endl;
            count = count + 1;
            
            cout<<"-------------------------"<<endl;
        }
        acount = acount + 1;
        fout.close();

    }

}

void Mphib_secondsheet_belowthreshold_vs_s3_N_fixed_s3imag_with_weights()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = 16;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    double s3imag = 0.1;//1.0e-5;
    

    double delspoints = 2000.0;
    double points1 = 500.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    //double box = 10.0;
    int acount = 0;

    double dela = 0.051;
    //for(double a=2.401;a<=16.03;a=a+dela)
    vector<double> epsvec(5);
    epsvec[0] = 0.000010;
    epsvec[1] = 0.000007;
    epsvec[2] = 0.000005;
    epsvec[3] = 0.000003;
    epsvec[4] = 0.000001;

    //for(int epsi=0;epsi<epsvec.size();++epsi)
    //{

    //    eps = epsvec[epsi];

    //int N = 500;
    //for(int N=900;N<=1000;N=N+20)
    {
        //if(a>=3.7) dela = 0.101;
        //else if(a>=17.0) dela = 0.51;
        //else if(a>=1000.0) dela = 300.01;

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        double sinitial = 8.7;//6.5;//8.3;//5.0;//8.974;//8.6;//8.975;//8.96815;//7.5;//8.7;//7.5;//8.75;//1.001;//gvalright - 0.5;//(gvalleft + (gvalleft+gvalright)/2.0)/2.0;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
        double sfinal = real(phibthreshold(a,m));//8.82;
        double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
        points1 = N;
        //string filename="Mphib_secondsheet_acount_" + to_string(acount) + ".dat";
        string filename="Mphib_a_" + to_string(a) + "_N_" + to_string((int)N) + "_temp1.dat";
        
        //string filename="testGL_Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_2.dat";
        //string filename="tmp.dat";
        

        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        //double eta = 20.0;

        for(int i=0;i<delspoints+1;++i) 
        //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
        { 
            double s3real = sinitial + i*dels;
            //comp s = 8.65 + ii*0.05;
            //if(s>=sfinal) dels = 0.05/100.0;
            //eps = eps_above_threshold(eta,s,a,m,points);
            //s3real = 8.8; //this was added
            //s3imag = 0.5; //this was added 

            comp s = s3real + ii*s3imag;
            cout<<"am = "<<a<<endl;
            //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = sigmab(a,m);
            comp qval = pmom(s,sigb,m);
            cout<<"q = "<<qval<<endl;

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        

            //cout<<"s:"<<s<<endl;
            cout<<"run: "<<count+1<<endl;
            cout<<"a = "<<a<<'\t'<<" dela = "<<dela<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            
            int switch_for_gvec_fixer = 0; //0 means we have used contour and turn on gvec fixer
                                           //1 means we have used straight line and turn off gvec fixer
            
            if(s3imag<0.0)
            {
                //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                cout<<"contour_43 chosen"<<endl;
                //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                mom_vector_maker_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1);
                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
            }
            else if(s3imag>=0.0)
            {
                //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                cout<<"contour_47 chosen"<<endl;
                //mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
                int changed_points = points1 + 5;
                vector<comp> tmp_qvec;
                vector<comp> tmp_weights;
                //mom_vector_maker_seba_imspos(tmp_qvec,tmp_weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)changed_points,tag1,tag2,switch_for_gvec_fixer);
                //mom_vector_maker_seba_imspos(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                //mom_vector_maker_47_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
                //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                
                //mom_vector_maker_seba_imspos_2(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer); //change this
            
                //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                
                /*
                mom_vector_maker_seba_imsneg(tmp_qvec,tmp_weights,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1);
                for(int i=0;i<tmp_qvec.size();++i)
                {
                    comp qv = real(tmp_qvec[i]) - ii*imag(tmp_qvec[i]);
                    comp w = real(tmp_weights[i]) - ii*imag(tmp_weights[i]);
                    qvec.push_back(qv);
                    weights.push_back(w);
                }
                cout<<"tag1 = "<<tag1<<endl;
                cout<<"tag2 = "<<tag2<<endl;
                */

                mom_vector_maker_seba_imspos_5(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)points1,tag1,tag2,switch_for_gvec_fixer);
        
                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                
                //this portion is for loading ready made contours:
                /*ifstream fin;
                string contour_file = "for_digonto.txt";
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
                qvec = contour;*/
                
                
            }

            cout<<"qvec created with size = "<<qvec.size()<<endl;
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            
            if(s3imag<0.0)
            {
                Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            }
            else if(s3imag>=0.0)
            {
                Gvec_maker_momrep_withtags_1(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            }

            switch_for_gvec_fixer=0;
            if(switch_for_gvec_fixer==0)
            {
                //Gvec_fixer_1(Gvec,Gvec,qvec,s,qval,m,eps,tag1,tag2); 
                //Gvec_fixer_2(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
                
                //this two were used for the paper calculations
                //Gvec_fixer_3(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
                //tag_fixer_for_gvec_fixer3(Gvec,qvec,s,qval,m,eps,tag1);
                //--------------------------------------------//


                //Gvec_fixer_5(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
                
                Gvec_fixer_6(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
               
                
            }

            cout<<"tag fixer = "<<tag1<<'\t'<<tag2<<endl;
            
            
            //sebatest part ---------
            ofstream fout2;//("opefile.dat");
            string opefile = "opefile.dat";
            fout2.open(opefile.c_str());
            ofstream fout3; 
            string opefile1 = "opefile1.dat";

            fout3.open(opefile1.c_str());

            for(int i=0;i<qvec.size();++i)
            {
                fout2<<i<<'\t'<<real(Gvec[i])<<'\t'<<imag(Gvec[i])<<endl;
                fout3<<i<<'\t'<<real(GS_pk(s,qvec[i],qval,m,eps))<<'\t'<<imag(GS_pk(s,qvec[i],qval,m,eps))<<'\t'<<real(GS_pk_secondsheet(s,qvec[i],qval,m,eps))<<'\t'<<imag(GS_pk_secondsheet(s,qvec[i],qval,m,eps))<<endl;
            }
            fout2.close();
            fout3.close();
            
            //-----------------------
            
            //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            Eigen::VectorXcd interpolater_Gvec = Gvec;
            Gvec = -1.0*Gvec;

            //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
            //Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
            //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);
            //Bmat_maker_momrep_2eps_withtags_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,tag1,tag2);
            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

            cout<<"determinant of B = "<<Bmat.determinant()<<endl;
            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        
            //cout<<"did this = "<<count<<" times"<<endl;
            //cout<<"i val = "<<i<<'\t';
            //cout<<"j val = "<<j<<endl;
            
            //LinearSolver_2(Bmat,dsol,Gvec,relerror);

            time(&time_start);
            LinearSolver_2(Bmat,dsol,Gvec,relerror);

            //cusolverComplex(Bmat,Gvec,dsol,size);
            time(&time_end);

            double time_taken = double(time_end - time_start);
            cout<<"Time taken to solve : "<<fixed 
                <<time_taken<<setprecision(5);
            cout<<" sec"<<endl;


            comp result;

        
            if(s3imag<0.0)
            {
                //interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
                interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            }
            else if(s3imag>=0.0)
            {
                //interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
                //interpolator_ds_integraleq_momrep_2eps_withtags_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
                interpolator_ds_integraleq_momrep_2eps_withtags_with_weights_usingGvec(dsol,qvec,weights,interpolater_Gvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
            
            }
            //result = dsol[10];
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            comp ds = result;
            result = gsq*result;

        
            comp rhopb = rhophib(s,a,m);
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;

            cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            fout<<setprecision(16)<<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(ds)<<'\t'
                <<imag(ds)<<'\t'
                <<real(mphib2)<<'\t'
                <<imag(mphib2)<<'\t'
                <<real(mphib2denom)<<'\t'
                <<imag(mphib2denom)<<'\t'
                <<a<<'\t'
                <<N<<'\t'
                <<gvalleft<<'\t'
                <<gvalright<<'\t'
                <<real(phibthreshold(a,m))<<'\t'
                <<eps<<endl;
            count = count + 1;
            
            cout<<"-------------------------"<<endl;
        }
        acount = acount + 1;
        fout.close();

    }

    //}

}

void Mphib_secondsheet_belowthreshold_vs_s3_N_3d_with_weights()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = -6.4;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    //double s3imag = 0.0;//1.0e-5;
    

    double delspoints = 100.0;
    double points1 = 250.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    //double box = 10.0;
    int acount = 0;

    double dela = 0.051;
    //for(double a=2.401;a<=16.03;a=a+dela)
    vector<double> epsvec(5);
    epsvec[0] = 0.000010;
    epsvec[1] = 0.000007;
    epsvec[2] = 0.000005;
    epsvec[3] = 0.000003;
    epsvec[4] = 0.000001;

    //for(int epsi=0;epsi<epsvec.size();++epsi)
    //{

    //    eps = epsvec[epsi];

    //int N = 500;
    //for(int N=900;N<=1000;N=N+20)
    {
        //if(a>=3.7) dela = 0.101;
        //else if(a>=17.0) dela = 0.51;
        //else if(a>=1000.0) dela = 300.01;

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        double sinitial = 8.85;//6.5;//8.3;//5.0;//8.974;//8.6;//8.975;//8.96815;//7.5;//8.7;//7.5;//8.75;//1.001;//gvalright - 0.5;//(gvalleft + (gvalleft+gvalright)/2.0)/2.0;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
        double sfinal = 9.05;//real(phibthreshold(a,m));//8.82;
        double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
        points1 = N;
        //string filename="Mphib_secondsheet_acount_" + to_string(acount) + ".dat";
        string filename="Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_3d.dat";
        
        //string filename="testGL_Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_2.dat";
        //string filename="tmp.dat";
        
        double simaginitial = -0.15;
        double simagfinal = 0.151;
        double delsimagpoints = 50;
        double delsimag = abs(simaginitial - simagfinal)/delsimagpoints;

        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        //double eta = 20.0;
        //test = 8.90465,0.01713

        for(int j=0;j<delsimagpoints+1;++j)
        {
            double s3imag = simaginitial + j*delsimag;
       
        for(int i=0;i<delspoints+1;++i)
        //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
        { 
            double s3real = sinitial + i*dels;
            //double s = 8.78;
            //if(s>=sfinal) dels = 0.05/100.0;
            //eps = eps_above_threshold(eta,s,a,m,points);

            comp s =  s3real + ii*s3imag;
            cout<<"am = "<<a<<endl;
            //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = 2.0*m*m;//sigmab(a,m);
            comp qval = pmom(s,sigb,m);
            cout<<"q = "<<qval<<endl;

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        

            //cout<<"s:"<<s<<endl;
            cout<<"run: "<<count+1<<endl;
            cout<<"a = "<<a<<'\t'<<" dela = "<<dela<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            
            int switch_for_gvec_fixer = 0; //0 means we have used contour and turn on gvec fixer
                                           //1 means we have used straight line and turn off gvec fixer
            
            if(s3imag<0.0)
            {
                //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                cout<<"contour_43 chosen"<<endl;
                //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                
                //this one was using last
                //mom_vector_maker_43_with_weights_with_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                
                //mom_vector_maker_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1);
                
                line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                switch_for_gvec_fixer = 1;
            }
            else if(s3imag>=0.0)
            {
                //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                cout<<"contour_47 chosen"<<endl;
                //mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
                int changed_points = points1 + 5;
                vector<comp> tmp_qvec;
                vector<comp> tmp_weights;
                //mom_vector_maker_seba_imspos(tmp_qvec,tmp_weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)changed_points,tag1,tag2,switch_for_gvec_fixer);
                //mom_vector_maker_seba_imspos(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                
                //mom_vector_maker_47_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
                
                //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                //mom_vector_maker_seba_imspos_2(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                

                //mom_vector_maker_43_with_weights_with_seba_imsneg(qvec,weights,real(s)-ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_seba_imspos_2_with_contour47(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                /*for(int ll=0;ll<qvec.size();++ll)
                {
                    comp tempq = qvec[ll];
                    qvec[ll] = real(tempq) - ii*imag(tempq); 
                    comp tempweights = weights[ll];
                    weights[ll] = real(weights[ll]) - ii*imag(weights[ll]);
                }*/

                //mom_vector_maker_seba_imspos_5(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                
                line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                switch_for_gvec_fixer = 1;
                tag1 = 0;
                tag2 = 0;

                /*if(s3imag>1.0e-3)
                mom_vector_maker_seba_imspos_5(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                else
                mom_vector_maker_seba_imspos_2(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                */
                /*
                if(tag1!=0)
                    switch_for_gvec_fixer = 0;
                else 
                    switch_for_gvec_fixer = 1;
                */
                //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                
                /*for(int i=0;i<tmp_qvec.size()-5;++i)
                {
                    qvec.push_back(tmp_qvec[i]);
                    weights.push_back(tmp_weights[i]);
                }*/
                cout<<"tag1 = "<<tag1<<endl;
                cout<<"tag2 = "<<tag2<<endl;
                
                
            }

            cout<<"qvec created with size = "<<qvec.size()<<endl;
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            
            if(s3imag<0.0)
            {
                Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            }
            else if(s3imag>=0.0)
            {
                Gvec_maker_momrep_withtags_1(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            }

            if(switch_for_gvec_fixer==0)
                cout<<"gvec fixer switch = "<<switch_for_gvec_fixer<<",deform contour"<<endl;
            else 
                cout<<"gvec fixer switch = "<<switch_for_gvec_fixer<<",take straight line"<<endl;

            if(switch_for_gvec_fixer==0)
            {   
                //cout<<"here the problem lies 1"<<endl;
                //Gvec_fixer_3(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
                
                //we used this before 
                //Gvec_fixer_4(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);//Gvec fixer has a contour maker inside of it
                
                //cout<<"here the problem lies 2"<<endl;
                //tag_fixer_for_gvec_fixer3(Gvec,qvec,s,qval,m,eps,tag1); //this works for gvec_fixer_4 as well
                //--------------------------------//

                Gvec_fixer_6(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);//Gvec fixer has a contour maker inside of it
                
                /*if(s3imag>1.0e-3)
                Gvec_fixer_6(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);//Gvec fixer has a contour maker inside of it
                else 
                {
                    //we used this before 
                    Gvec_fixer_3(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);//Gvec fixer has a contour maker inside of it
                
                    //cout<<"here the problem lies 2"<<endl;
                    tag_fixer_for_gvec_fixer3(Gvec,qvec,s,qval,m,eps,tag1); //this works for gvec_fixer_4 as well
                    //--------------------------------//
 
                }*/
                //cout<<"here the problem lies 3"<<endl;
               
            }

            cout<<"tag fixer = "<<tag1<<'\t'<<tag2<<endl;
            //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            Eigen::VectorXcd interpolater_Gvec = Gvec;
            Gvec = -1.0*Gvec;

            //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
            //Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
            //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);
            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

            
            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        
            //cout<<"did this = "<<count<<" times"<<endl;
            //cout<<"i val = "<<i<<'\t';
            //cout<<"j val = "<<j<<endl;
            
            //LinearSolver_2(Bmat,dsol,Gvec,relerror);

            time(&time_start);
            LinearSolver_2(Bmat,dsol,Gvec,relerror);

            //cusolverComplex(Bmat,Gvec,dsol,size);
            time(&time_end);

            double time_taken = double(time_end - time_start);
            cout<<"Time taken to solve : "<<fixed 
                <<time_taken<<setprecision(5);
            cout<<" sec"<<endl;


            comp result;

        
            if(s3imag<0.0)
            {
                //interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
                interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            }
            else if(s3imag>=0.0)
            {
                //interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
                //interpolator_ds_integraleq_momrep_2eps_withtags_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
                interpolator_ds_integraleq_momrep_2eps_withtags_with_weights_usingGvec(dsol,qvec,weights,interpolater_Gvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
            
            }
            //result = dsol[10];
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            comp ds = result;
            result = gsq*result;

            //result = Bmat.determinant();

        
            comp rhopb = rhophib(s,a,m);
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;

            cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            fout<<setprecision(16)<<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(ds)<<'\t'
                <<imag(ds)<<'\t'
                <<real(mphib2)<<'\t'
                <<imag(mphib2)<<'\t'
                <<real(mphib2denom)<<'\t'
                <<imag(mphib2denom)<<'\t'
                <<a<<'\t'
                <<N<<'\t'
                <<gvalleft<<'\t'
                <<gvalright<<'\t'
                <<real(phibthreshold(a,m))<<'\t'
                <<eps<<endl;
            count = count + 1;
            
            cout<<"-------------------------"<<endl;
        }
        acount = acount + 1;
        
        }
        fout.close();

    }

    //}

}

void ds_belowthreshold_vs_s3_N_3d_with_weights()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = -6.40;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    //double s3imag = 0.0;//1.0e-5;
    

    double delspoints = 100.0;
    double points1 = 250.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    //double box = 10.0;
    int acount = 0;

    double dela = 0.051;
    //for(double a=2.401;a<=16.03;a=a+dela)
    vector<double> epsvec(5);
    epsvec[0] = 0.000010;
    epsvec[1] = 0.000007;
    epsvec[2] = 0.000005;
    epsvec[3] = 0.000003;
    epsvec[4] = 0.000001;

    //for(int epsi=0;epsi<epsvec.size();++epsi)
    //{

    //    eps = epsvec[epsi];

    //int N = 500;
    //for(int N=900;N<=1000;N=N+20)
    {
        //if(a>=3.7) dela = 0.101;
        //else if(a>=17.0) dela = 0.51;
        //else if(a>=1000.0) dela = 300.01;

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        double sinitial = 8.99;//6.5;//8.3;//5.0;//8.974;//8.6;//8.975;//8.96815;//7.5;//8.7;//7.5;//8.75;//1.001;//gvalright - 0.5;//(gvalleft + (gvalleft+gvalright)/2.0)/2.0;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
        double sfinal = 9.06;//real(phibthreshold(a,m));//8.82;
        double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
        points1 = N;
        //string filename="Mphib_secondsheet_acount_" + to_string(acount) + ".dat";
        string filename="Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_3d.dat";
        
        //string filename="testGL_Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_2.dat";
        //string filename="tmp.dat";
        
        double simaginitial = -0.0005;
        double simagfinal = 0.0000001;
        double delsimagpoints = 100;
        double delsimag = abs(simaginitial - simagfinal)/delsimagpoints;

        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        //double eta = 20.0;
        //test = 8.90465,0.01713

        for(int j=0;j<delsimagpoints+1;++j)
        {
            double s3imag = simaginitial + j*delsimag;
       
        for(int i=0;i<delspoints+1;++i)
        //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
        { 
            double s3real = sinitial + i*dels;
            //double s = 8.78;
            //if(s>=sfinal) dels = 0.05/100.0;
            //eps = eps_above_threshold(eta,s,a,m,points);

            comp s =  s3real + ii*s3imag;
            cout<<"am = "<<a<<endl;
            //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = 2.0*m*m;//sigmab(a,m);
            comp qval = pmom(s,sigb,m);
            cout<<"q = "<<qval<<endl;

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        

            //cout<<"s:"<<s<<endl;
            cout<<"run: "<<count+1<<endl;
            cout<<"a = "<<a<<'\t'<<" dela = "<<dela<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            
            int switch_for_gvec_fixer = 0; //0 means we have used contour and turn on gvec fixer
                                           //1 means we have used straight line and turn off gvec fixer
            
            if(s3imag<=0.0)
            {
                //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                cout<<"contour_43 chosen"<<endl;
                //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                
                //this one was using last
                //mom_vector_maker_43_with_weights_with_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                
                //mom_vector_maker_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1);
                comp m2kbc = M2kbranchcut_right_momrep_plus_eps(s,m,eps_for_m2k);
                double shift = real(m2kbc)/10;
                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                //contour_for_resonance(qvec,weights,kmin,kmax,shift,points1);
                //contour_for_resonance_1(qvec,weights,m2kbc,kmin,kmax,shift,points1);
                cout<<"m2k bc ="<<m2kbc<<endl;
                double reshift = 0.2;
                double imshift = -0.3;
                contour_for_resonance_2(qvec,weights,kmin,reshift,imshift,kmax,points1);
                
                switch_for_gvec_fixer = 1;
            }
            else if(s3imag>0.0)
            {
                //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                cout<<"contour_47 chosen"<<endl;
                //mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
                int changed_points = points1 + 5;
                vector<comp> tmp_qvec;
                vector<comp> tmp_weights;
                //mom_vector_maker_seba_imspos(tmp_qvec,tmp_weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)changed_points,tag1,tag2,switch_for_gvec_fixer);
                //mom_vector_maker_seba_imspos(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                
                //mom_vector_maker_47_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
                
                //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                //mom_vector_maker_seba_imspos_2(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                

                //mom_vector_maker_43_with_weights_with_seba_imsneg(qvec,weights,real(s)-ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_seba_imspos_2_with_contour47(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                /*for(int ll=0;ll<qvec.size();++ll)
                {
                    comp tempq = qvec[ll];
                    qvec[ll] = real(tempq) - ii*imag(tempq); 
                    comp tempweights = weights[ll];
                    weights[ll] = real(weights[ll]) - ii*imag(weights[ll]);
                }*/

                //mom_vector_maker_seba_imspos_5(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                
                comp m2kbc = M2kbranchcut_right_momrep_plus_eps(s,m,eps_for_m2k);
                cout<<"m2k bc ="<<m2kbc<<endl;
                double shift = real(m2kbc)/10;
                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                //contour_for_resonance(qvec,weights,kmin,kmax,shift,points1);
                //contour_for_resonance_1(qvec,weights,m2kbc,kmin,kmax,shift,points1);

                double reshift = 0.2;
                double imshift = -0.3;
                contour_for_resonance_2(qvec,weights,kmin,reshift,imshift,kmax,points1);
                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                

                switch_for_gvec_fixer = 1;
                tag1 = 0;
                tag2 = 0;

                
                cout<<"tag1 = "<<tag1<<endl;
                cout<<"tag2 = "<<tag2<<endl;
                
                
            }

            cout<<"qvec created with size = "<<qvec.size()<<endl;
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
          

            cout<<"tag fixer = "<<tag1<<'\t'<<tag2<<endl;
            //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            Eigen::VectorXcd interpolater_Gvec = Gvec;
            Gvec = -1.0*Gvec;

            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

            
            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        
            //cout<<"did this = "<<count<<" times"<<endl;
            //cout<<"i val = "<<i<<'\t';
            //cout<<"j val = "<<j<<endl;
            
            //LinearSolver_2(Bmat,dsol,Gvec,relerror);

            time(&time_start);
            LinearSolver_2(Bmat,dsol,Gvec,relerror);

            //cusolverComplex(Bmat,Gvec,dsol,size);
            time(&time_end);

            double time_taken = double(time_end - time_start);
            cout<<"Time taken to solve : "<<fixed 
                <<time_taken<<setprecision(5);
            cout<<" sec"<<endl;


            comp result;

            interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            
        
            
            //result = dsol[10];
            comp m2k = M2kfunc(a,sigb,m, eps_for_m2k);
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            comp ds = result;
            result = m2k*result*m2k;

            //result = Bmat.determinant();

        
            comp rhopb = rhophib(s,a,m);
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;

            cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            fout<<setprecision(16)<<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(ds)<<'\t'
                <<imag(ds)<<'\t'
                <<real(mphib2)<<'\t'
                <<imag(mphib2)<<'\t'
                <<real(mphib2denom)<<'\t'
                <<imag(mphib2denom)<<'\t'
                <<a<<'\t'
                <<N<<'\t'
                <<gvalleft<<'\t'
                <<gvalright<<'\t'
                <<real(phibthreshold(a,m))<<'\t'
                <<eps<<endl;
            count = count + 1;
            
            cout<<"-------------------------"<<endl;
        }
        acount = acount + 1;
        
        }
        fout.close();

    }

    //}

}



void m2k_cut_tagger(    vector<comp> &qvec,
                        comp s,
                        double m,
                        int &tag    )
{
    //here we figure out the tag for m2k_cut that 
    //crosses the real q axis. 
    vector<comp> m2k_cut_vec;
    int m2k_sign1 = 0;
    int m2k_sign2 = 0;
    comp m2kcut_in_p_prev;
    int tag_for_m2k = 0;
    for(int i=0;i<5000;++i)
    {
        comp ii = {0.0,1.0};
        double pi = acos(-1.0);
        double theta = -90*pi/180;
        //since the cut is rotated -90 degree
        //we can just write it straightforward
        double delsigk = abs(5.0)/5000;
        double delsig = -i*delsigk;
        comp sigk = 4.0 + (comp)ii*delsig;
        //comp m2kcut = m2k_cut_tracer(sigk, m, theta);
        comp m2kcut_in_p = pmom(s,sigk,m);
        m2k_cut_vec.push_back(m2kcut_in_p);
        m2k_sign1 = sign_checker(imag(m2kcut_in_p));
        if(m2k_sign1==0 || m2k_sign2==0)
        {
            m2k_sign2 = m2k_sign1;
            m2kcut_in_p_prev = m2kcut_in_p;
            continue;
        }
        else if(m2k_sign1==m2k_sign2)
        {
            m2k_sign2 = m2k_sign1;
            m2kcut_in_p_prev = m2kcut_in_p;
            continue;

        } 
        else 
        {
            comp avg_m2k_cut = (m2kcut_in_p + m2kcut_in_p_prev)/2.0;
            for(int ind=0;ind<qvec.size();++ind)
            {
                if(real(qvec[ind])<=real(avg_m2k_cut))
                {
                    tag_for_m2k = ind;
                }
            }

        }
                //cout<<sigk<<'\t'<<m2kcut_in_p<<endl;
                //fout<<real(m2kcut_in_p)<<'\t'<<imag(m2kcut_in_p)<<endl;
    }

    tag = tag_for_m2k;
}



//this checks the transition between sheets for 
//m2k or the kernel depending on the set output
//the contour is straight, and we rotate the cut 
//structure of m2k to be in the negative imaginary
//s2k axis. 
void check_m2k_cut_transition()
{
    comp ii = {0.0,1.0};
    double a = -6.4;
    double m = 1.0;
    double sreal = 9.10;
    double simag = 0.0;
    comp s = sreal + ii*simag;
    double epsilon = 0.0;

    double pi = acos(-1.0);
    double theta = -90.0*pi/180.0;

    double kmin = 0.0;
    comp kmax = pmom(s,0.0,m);
    double N = 500.0;
    vector<comp> qvec;
    vector<comp> weights;

    line_maker_with_weights(qvec,weights,kmin,kmax,N);

    string m2kfile = "m2k_transition_test.dat";

    ofstream fout;
    fout.open(m2kfile.c_str());

    for(int i=0;i<qvec.size();++i)
    {
        comp p = pmom(s,2.0,m);
        comp k = qvec[i];
        comp sigk = sigma_p(s,k,m);

        comp m2k = M2kfunc_rotatedcut(a,sigk,m,epsilon,theta);
        comp m2k_secondsheet = M2kfunc_secondsheet_rotatedcut(a,sigk,m,epsilon,theta);

        comp kern = kernel_pk_2eps_for_resonances(s,p,k,a,m,epsilon,epsilon,theta,0);
        fout<<i<<'\t'
            <<real(k)<<'\t'
            <<imag(k)<<'\t'
            <<real(sigk)<<'\t'
            <<imag(sigk)<<'\t'
            <<real(m2k)<<'\t'
            <<imag(m2k)<<'\t'
            <<real(m2k_secondsheet)<<'\t'
            <<imag(m2k_secondsheet)<<'\t'
            <<real(kern)<<'\t'
            <<imag(kern)<<endl;

    }
    fout.close();
}

//this section prints outs the amplitude for 3d plotting 
//the purpose is to look for three body resonances
void ds_belowthreshold_vs_s3_N_3d_with_weights_for_resonances()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = -1.50;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    //double s3imag = 0.0;//1.0e-5;
    

    double delspoints = 100.0;
    double points1 = 250.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    //double box = 10.0;
    int acount = 0;

    double dela = 0.051;
    //for(double a=2.401;a<=16.03;a=a+dela)
    vector<double> epsvec(5);
    epsvec[0] = 0.000010;
    epsvec[1] = 0.000007;
    epsvec[2] = 0.000005;
    epsvec[3] = 0.000003;
    epsvec[4] = 0.000001;

    //for(int epsi=0;epsi<epsvec.size();++epsi)
    //{

    //    eps = epsvec[epsi];

    //int N = 500;
    //for(int N=900;N<=1000;N=N+20)
    {
        //if(a>=3.7) dela = 0.101;
        //else if(a>=17.0) dela = 0.51;
        //else if(a>=1000.0) dela = 300.01;

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        double sinitial = 7.1;//6.5;//8.3;//5.0;//8.974;//8.6;//8.975;//8.96815;//7.5;//8.7;//7.5;//8.75;//1.001;//gvalright - 0.5;//(gvalleft + (gvalleft+gvalright)/2.0)/2.0;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
        double sfinal = 9.3;//real(phibthreshold(a,m));//8.82;
        double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
        points1 = N;
        //string filename="Mphib_secondsheet_acount_" + to_string(acount) + ".dat";
        string filename="Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_3d.dat";
        
        //string filename="testGL_Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_2.dat";
        //string filename="tmp.dat";
        
        double simaginitial = -1.505;
        double simagfinal = 0.1011;
        double delsimagpoints = 100;
        double delsimag = abs(simaginitial - simagfinal)/delsimagpoints;

        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        //double eta = 20.0;
        //test = 8.90465,0.01713

        for(int j=0;j<delsimagpoints+1;++j)
        {
            double s3imag = simaginitial + j*delsimag;
       
        for(int i=0;i<delspoints+1;++i)
        //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
        { 
            double s3real = sinitial + i*dels;
            //double s = 8.78;
            //if(s>=sfinal) dels = 0.05/100.0;
            //eps = eps_above_threshold(eta,s,a,m,points);

            comp s =  s3real + ii*s3imag;
            cout<<"am = "<<a<<endl;
            //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = 2.0*m*m;//sigmab(a,m);
            comp qval = pmom(s,sigb,m);
            cout<<"q = "<<qval<<endl;

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        

            //cout<<"s:"<<s<<endl;
            cout<<"run: "<<count+1<<endl;
            cout<<"a = "<<a<<'\t'<<" dela = "<<dela<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            
            int switch_for_gvec_fixer = 0; //0 means we have used contour and turn on gvec fixer
                                           //1 means we have used straight line and turn off gvec fixer
            
            if(s3imag<0.0)
            {
                //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                cout<<"contour_43 chosen"<<endl;
                //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                
                //this one was using last
                //mom_vector_maker_43_with_weights_with_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                
                //mom_vector_maker_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1);
                comp m2kbc = M2kbranchcut_right_momrep_plus_eps(s,m,eps_for_m2k);
                double shift = real(m2kbc)/10;
                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                //contour_for_resonance(qvec,weights,kmin,kmax,shift,points1);
                //contour_for_resonance_1(qvec,weights,m2kbc,kmin,kmax,shift,points1);
                cout<<"m2k bc ="<<m2kbc<<endl;
                double reshift = 0.2;
                double imshift = -0.3;
                //contour_for_resonance_2(qvec,weights,kmin,reshift,imshift,kmax,points1);

                if(real(s)<9.0)
                line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                else 
                contour_for_resonance_3(qvec,weights,s,m,kmin,kmax,eps_for_m2k,-0.1,points1);

                //contour_for_resonance_3(qvec,weights,s,m,kmin,kmax,eps_for_m2k,-0.1,points1);
                
                switch_for_gvec_fixer = 1;

                /*for(int i=0;i<qvec.size();++i)
                {
                    qvec[i] = qvec[i] - ii*0.000001;
                }*/
            }
            else if(s3imag>=0.0)
            {
                //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                cout<<"contour_47 chosen"<<endl;
                //mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
                int changed_points = points1 + 5;
                vector<comp> tmp_qvec;
                vector<comp> tmp_weights;
                //mom_vector_maker_seba_imspos(tmp_qvec,tmp_weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)changed_points,tag1,tag2,switch_for_gvec_fixer);
                //mom_vector_maker_seba_imspos(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                
                //mom_vector_maker_47_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
                
                //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                //mom_vector_maker_seba_imspos_2(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                

                //mom_vector_maker_43_with_weights_with_seba_imsneg(qvec,weights,real(s)-ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_seba_imspos_2_with_contour47(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                /*for(int ll=0;ll<qvec.size();++ll)
                {
                    comp tempq = qvec[ll];
                    qvec[ll] = real(tempq) - ii*imag(tempq); 
                    comp tempweights = weights[ll];
                    weights[ll] = real(weights[ll]) - ii*imag(weights[ll]);
                }*/

                //mom_vector_maker_seba_imspos_5(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                
                comp m2kbc = M2kbranchcut_right_momrep_plus_eps(s,m,eps_for_m2k);
                cout<<"m2k bc ="<<m2kbc<<endl;
                double shift = real(m2kbc)/10;
                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                //contour_for_resonance(qvec,weights,kmin,kmax,shift,points1);
                //contour_for_resonance_1(qvec,weights,m2kbc,kmin,kmax,shift,points1);

                double reshift = 0.2;
                double imshift = -0.3;
                //contour_for_resonance_2(qvec,weights,kmin,reshift,imshift,kmax,points1);
                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);

                if(real(s)<9.0)
                line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                else                
                contour_for_resonance_3(qvec,weights,s,m,kmin,kmax,eps_for_m2k,-0.1,points1);
            
                switch_for_gvec_fixer = 1;
                tag1 = 0;
                tag2 = 0;

                
                cout<<"tag1 = "<<tag1<<endl;
                cout<<"tag2 = "<<tag2<<endl;
                
                
            }

            cout<<"qvec created with size = "<<qvec.size()<<endl;
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
          

            cout<<"tag fixer = "<<tag1<<'\t'<<tag2<<endl;
            //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            Eigen::VectorXcd interpolater_Gvec = Gvec;
            Gvec = -1.0*Gvec;

            int tag_for_m2k = 0;
            

            //if(real(s)>=9.0 && imag(s)<=0.0)
            if(real(s)>=9.0)
            
            {
                double pi = acos(-1.0);
                double theta = -90.0*pi/180.0;
                //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
                //cout<<"m2k cut tagger = "<<tag_for_m2k<<endl;
                Bmat_maker_momrep_2eps_with_weights_rotated_m2k_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta);
                //Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }
            else
            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

            
            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        
            //cout<<"did this = "<<count<<" times"<<endl;
            //cout<<"i val = "<<i<<'\t';
            //cout<<"j val = "<<j<<endl;
            
            //LinearSolver_2(Bmat,dsol,Gvec,relerror);

            time(&time_start);
            LinearSolver_2(Bmat,dsol,Gvec,relerror);

            //cusolverComplex(Bmat,Gvec,dsol,size);
            time(&time_end);

            double time_taken = double(time_end - time_start);
            cout<<"Time taken to solve : "<<fixed 
                <<time_taken<<setprecision(5);
            cout<<" sec"<<endl;


            comp result;

            //if(real(s)>=9.0 && imag(s)<=0.0)
            if(real(s)>=9.0)
            {
                double pi = acos(-1.0);
                double theta = -90.0*pi/180.0;
                //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
                interpolator_ds_integraleq_momrep_2eps_with_weights_rotated_m2k_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,result);
            }
            else 
            interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            
            //result = Bmat.determinant();
            
            //result = dsol[10];
            comp m2k = M2kfunc(a,sigb,m, eps_for_m2k);
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            comp ds = result;
            result = m2k*result*m2k;

            //result = Bmat.determinant();

        
            comp rhopb = rhophib(s,a,m);
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;

            cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            fout<<setprecision(16)<<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(ds)<<'\t'
                <<imag(ds)<<'\t'
                <<real(mphib2)<<'\t'
                <<imag(mphib2)<<'\t'
                <<real(mphib2denom)<<'\t'
                <<imag(mphib2denom)<<'\t'
                <<a<<'\t'
                <<N<<'\t'
                <<gvalleft<<'\t'
                <<gvalright<<'\t'
                <<real(phibthreshold(a,m))<<'\t'
                <<eps<<endl;
            count = count + 1;
            
            cout<<"-------------------------"<<endl;
        }
        acount = acount + 1;
        
        }
        fout.close();

    }

    //}

}

void ds_belowthreshold_vs_s3_N_3d_with_weights_for_resonances_1()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = -8.1;//-6.40;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    //double s3imag = 0.0;//1.0e-5;
    

    double delspoints = 100.0;
    double points1 = 250.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    //double box = 10.0;
    int acount = 0;

    double dela = 0.051;
    //for(double a=2.401;a<=16.03;a=a+dela)
    vector<double> epsvec(5);

    {

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        double sinitial = 8.99;//6.5;//8.3;//5.0;//8.974;//8.6;//8.975;//8.96815;//7.5;//8.7;//7.5;//8.75;//1.001;//gvalright - 0.5;//(gvalleft + (gvalleft+gvalright)/2.0)/2.0;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
        double sfinal = 9.01;//real(phibthreshold(a,m));//8.82;
        double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
        points1 = N;
        //string filename="Mphib_secondsheet_acount_" + to_string(acount) + ".dat";
        string filename="Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_3d.dat";
        
        //string filename="testGL_Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_2.dat";
        //string filename="tmp.dat";
        
        double simaginitial = -0.011;
        double simagfinal = 0.01011;
        double delsimagpoints = 100;
        double delsimag = abs(simaginitial - simagfinal)/delsimagpoints;

        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        //double eta = 20.0;
        //test = 8.90465,0.01713

        for(int j=0;j<delsimagpoints+1;++j)
        {
            double s3imag = simaginitial + j*delsimag;
       
        for(int i=0;i<delspoints+1;++i)
        //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
        { 
            double s3real = sinitial + i*dels;

            comp s =  s3real + ii*s3imag;
            cout<<"am = "<<a<<endl;
            //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = 2.0*m*m;//sigmab(a,m);
            comp qval = pmom(s,sigb,m);
            cout<<"q = "<<qval<<endl;

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        

            //cout<<"s:"<<s<<endl;
            cout<<"run: "<<count+1<<endl;
            cout<<"a = "<<a<<'\t'<<" dela = "<<dela<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            
            int switch_for_gvec_fixer = 0; //0 means we have used contour and turn on gvec fixer
                                           //1 means we have used straight line and turn off gvec fixer
            
            
            int tag_for_m2k = 0;

            if(s3imag<0.0)
            {
                //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                cout<<"contour_43 chosen"<<endl;
                //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                
                //this one was using last
                //mom_vector_maker_43_with_weights_with_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                
                //mom_vector_maker_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1);
                comp m2kbc = M2kbranchcut_right_momrep_plus_eps(s,m,eps_for_m2k);
                double shift = real(m2kbc)/10;
                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                //contour_for_resonance(qvec,weights,kmin,kmax,shift,points1);
                //contour_for_resonance_1(qvec,weights,m2kbc,kmin,kmax,shift,points1);
                cout<<"m2k bc ="<<m2kbc<<endl;
                double reshift = 0.2;
                double imshift = -0.3;
                //contour_for_resonance_2(qvec,weights,kmin,reshift,imshift,kmax,points1);

                if(real(s)<9.0)
                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                contour_for_resonance_4(qvec,weights,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k);

                else 
                contour_for_resonance_3(qvec,weights,s,m,kmin,kmax,eps_for_m2k,-0.1,points1);

                //contour_for_resonance_3(qvec,weights,s,m,kmin,kmax,eps_for_m2k,-0.1,points1);
                
                switch_for_gvec_fixer = 1;

                /*for(int i=0;i<qvec.size();++i)
                {
                    qvec[i] = qvec[i] - ii*0.000001;
                }*/
            }
            else if(s3imag>=0.0)
            {
                //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                cout<<"contour_47 chosen"<<endl;
                //mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
                int changed_points = points1 + 5;
                vector<comp> tmp_qvec;
                vector<comp> tmp_weights;
                //mom_vector_maker_seba_imspos(tmp_qvec,tmp_weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)changed_points,tag1,tag2,switch_for_gvec_fixer);
                //mom_vector_maker_seba_imspos(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                
                //mom_vector_maker_47_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
                
                //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                //mom_vector_maker_seba_imspos_2(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                

                //mom_vector_maker_43_with_weights_with_seba_imsneg(qvec,weights,real(s)-ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_seba_imspos_2_with_contour47(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                /*for(int ll=0;ll<qvec.size();++ll)
                {
                    comp tempq = qvec[ll];
                    qvec[ll] = real(tempq) - ii*imag(tempq); 
                    comp tempweights = weights[ll];
                    weights[ll] = real(weights[ll]) - ii*imag(weights[ll]);
                }*/

                //mom_vector_maker_seba_imspos_5(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                
                comp m2kbc = M2kbranchcut_right_momrep_plus_eps(s,m,eps_for_m2k);
                cout<<"m2k bc ="<<m2kbc<<endl;
                double shift = real(m2kbc)/10;
                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                //contour_for_resonance(qvec,weights,kmin,kmax,shift,points1);
                //contour_for_resonance_1(qvec,weights,m2kbc,kmin,kmax,shift,points1);

                double reshift = 0.2;
                double imshift = -0.3;
                //contour_for_resonance_2(qvec,weights,kmin,reshift,imshift,kmax,points1);
                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);

                if(real(s)<9.0)
                line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                else                
                contour_for_resonance_3(qvec,weights,s,m,kmin,kmax,eps_for_m2k,-0.1,points1);
            
                switch_for_gvec_fixer = 1;
                tag1 = 0;
                tag2 = 0;

                
                cout<<"tag1 = "<<tag1<<endl;
                cout<<"tag2 = "<<tag2<<endl;
                
                
            }

            cout<<"qvec created with size = "<<qvec.size()<<endl;
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
          

            //cout<<"tag fixer = "<<tag1<<'\t'<<tag2<<endl;
            //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            Eigen::VectorXcd interpolater_Gvec = Gvec;
            Gvec = -1.0*Gvec;

            

            //if(real(s)>=9.0 && imag(s)<=0.0)
            if(real(s)>=9.0)
            
            {
                double pi = acos(-1.0);
                double theta = -90.0*pi/180.0;
                //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
                //cout<<"m2k cut tagger = "<<tag_for_m2k<<endl;
                Bmat_maker_momrep_2eps_with_weights_rotated_m2k_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta);
                //Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }
            else if(real(s)<9.0 && imag(s)<=0.0)
            {
                double theta = 0.0;
                Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }
            else
            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

            
            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        
            //cout<<"did this = "<<count<<" times"<<endl;
            //cout<<"i val = "<<i<<'\t';
            //cout<<"j val = "<<j<<endl;
            
            //LinearSolver_2(Bmat,dsol,Gvec,relerror);

            time(&time_start);
            LinearSolver_2(Bmat,dsol,Gvec,relerror);

            //cusolverComplex(Bmat,Gvec,dsol,size);
            time(&time_end);

            double time_taken = double(time_end - time_start);
            cout<<"Time taken to solve : "<<fixed 
                <<time_taken<<setprecision(5);
            cout<<" sec"<<endl;


            comp result;

            //if(real(s)>=9.0 && imag(s)<=0.0)
            if(real(s)>=9.0)
            {
                double pi = acos(-1.0);
                double theta = -90.0*pi/180.0;
                //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
                interpolator_ds_integraleq_momrep_2eps_with_weights_rotated_m2k_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,result);
            }
            else if(real(s)<9.0 && imag(s)<=0.0)
            {
                double theta = 0.0;
                interpolator_ds_integraleq_momrep_2eps_with_weights_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,tag_for_m2k,result);
            
            }
            else 
            interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            
            //result = Bmat.determinant();
            
            //result = dsol[10];
            comp m2k = M2kfunc(a,sigb,m, eps_for_m2k);
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            comp ds = result;
            result = m2k*result*m2k;

            //result = Bmat.determinant();

        
            comp rhopb = rhophib(s,a,m);
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;

            cout<<setprecision(16)<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            fout<<setprecision(16)<<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(ds)<<'\t'
                <<imag(ds)<<'\t'
                <<real(mphib2)<<'\t'
                <<imag(mphib2)<<'\t'
                <<real(mphib2denom)<<'\t'
                <<imag(mphib2denom)<<'\t'
                <<a<<'\t'
                <<N<<'\t'
                <<gvalleft<<'\t'
                <<gvalright<<'\t'
                <<real(phibthreshold(a,m))<<'\t'
                <<eps<<endl;
            count = count + 1;
            
            cout<<"-------------------------"<<endl;
        }
        acount = acount + 1;
        
        }
        fout.close();

    }

    //}

}

void ds_belowthreshold_vs_s3_N_3d_with_weights_for_resonances_secondsheet_bottom_half()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = -8.10;//-6.40;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    //double s3imag = 0.0;//1.0e-5;
    

    double delspoints = 100.0;
    double points1 = 250.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    //double box = 10.0;
    int acount = 0;

    double dela = 0.051;
    //for(double a=2.401;a<=16.03;a=a+dela)
    vector<double> epsvec(5);

    {

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        double sinitial = 8.99;//6.5;//8.3;//5.0;//8.974;//8.6;//8.975;//8.96815;//7.5;//8.7;//7.5;//8.75;//1.001;//gvalright - 0.5;//(gvalleft + (gvalleft+gvalright)/2.0)/2.0;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
        double sfinal = 9.01;//real(phibthreshold(a,m));//8.82;
        double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
        points1 = N;
        //string filename="Mphib_secondsheet_acount_" + to_string(acount) + ".dat";
        string filename="dS_sheet2_bottomhalf_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_3d.dat";
        
        //string filename="testGL_Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_2.dat";
        //string filename="tmp.dat";
        
        double simaginitial = -0.0101;
        double simagfinal = 0.01011;
        double delsimagpoints = 100;
        double delsimag = abs(simaginitial - simagfinal)/delsimagpoints;

        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        //double eta = 20.0;
        //test = 8.90465,0.01713

        for(int j=0;j<delsimagpoints+1;++j)
        {
            double s3imag = simaginitial + j*delsimag;
       
        for(int i=0;i<delspoints+1;++i)
        //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
        { 
            double s3real = sinitial + i*dels;

            comp s =  s3real + ii*s3imag;
            cout<<"am = "<<a<<endl;
            //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = 2.0*m*m;//sigmab(a,m);
            comp qval = pmom(s,sigb,m);
            cout<<"q = "<<qval<<endl;

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        

            //cout<<"s:"<<s<<endl;
            cout<<"run: "<<count+1<<endl;
            cout<<"a = "<<a<<'\t'<<" dela = "<<dela<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            
            int switch_for_gvec_fixer = 0; //0 means we have used contour and turn on gvec fixer
                                           //1 means we have used straight line and turn off gvec fixer
            
            
            int tag_for_m2k = 0;

            if(s3imag<0.0)
            {
                contour_for_resonance_4(qvec,weights,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k);

            }
            else if(s3imag>=0.0)
            {
                line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                //contour_for_resonance_3(qvec,weights,s,m,kmin,kmax,eps_for_m2k,-0.1,points1);
                
            }

            cout<<"tag at i = "<<tag_for_m2k<<endl;
            cout<<"qvec created with size = "<<qvec.size()<<endl;
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
          
            Eigen::VectorXcd interpolater_Gvec = Gvec;
            Gvec = -1.0*Gvec;

            

            if(imag(s)<0.0)
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }
            else if(imag(s)>=0.0)
            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);


            //if(real(s)>=9.0 && imag(s)<=0.0)
            /*if(real(s)>=9.0)
            
            {
                double pi = acos(-1.0);
                double theta = -90.0*pi/180.0;
                //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
                //cout<<"m2k cut tagger = "<<tag_for_m2k<<endl;
                Bmat_maker_momrep_2eps_with_weights_rotated_m2k_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta);
                //Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }
            else if(real(s)<9.0 && imag(s)<=0.0)
            {
                //double theta = 0.0;
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                //Bmat_maker_momrep_2eps_with_weights_rotated_m2k_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta);
                
                Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }
            else
            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);
            */
            
            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        

            time(&time_start);
            LinearSolver_2(Bmat,dsol,Gvec,relerror);

            //cusolverComplex(Bmat,Gvec,dsol,size);
            time(&time_end);

            double time_taken = double(time_end - time_start);
            cout<<"Time taken to solve : "<<fixed 
                <<time_taken<<setprecision(5);
            cout<<" sec"<<endl;


            comp result;


            if(imag(s)<0.0)
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                interpolator_ds_integraleq_momrep_2eps_with_weights_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,tag_for_m2k,result);
            }
            else 
            interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            

            //if(real(s)>=9.0 && imag(s)<=0.0)
            /*if(real(s)>=9.0)
            {
                double pi = acos(-1.0);
                double theta = -90.0*pi/180.0;
                //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
                interpolator_ds_integraleq_momrep_2eps_with_weights_rotated_m2k_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,result);
            }
            else if(real(s)<9.0 && imag(s)<=0.0)
            {
                double theta = -90.0*pi/180.0;
                //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
                interpolator_ds_integraleq_momrep_2eps_with_weights_rotated_m2k_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,result);
            
            }
            else 
            interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            */

            //result = Bmat.determinant();
            
            //result = dsol[10];
            comp m2k = M2kfunc(a,sigb,m, eps_for_m2k);
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            comp ds = result;
            result = m2k*result*m2k;

            //result = Bmat.determinant();

        
            comp rhopb = rhophib(s,a,m);
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;

            cout<<setprecision(16)<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            fout<<setprecision(16)<<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(ds)<<'\t'
                <<imag(ds)<<'\t'
                <<real(mphib2)<<'\t'
                <<imag(mphib2)<<'\t'
                <<real(mphib2denom)<<'\t'
                <<imag(mphib2denom)<<'\t'
                <<a<<'\t'
                <<N<<'\t'
                <<gvalleft<<'\t'
                <<gvalright<<'\t'
                <<real(phibthreshold(a,m))<<'\t'
                <<eps<<endl;
            count = count + 1;
            
            cout<<"-------------------------"<<endl;
        }
        acount = acount + 1;
        
        }
        fout.close();

    }

    //}

}

void ds_belowthreshold_vs_s3_N_3d_with_weights_for_resonances_secondsheet_full()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = -8.10;//-6.40;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    //double s3imag = 0.0;//1.0e-5;
    

    double delspoints = 100.0;
    double points1 = 250.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    //double box = 10.0;
    int acount = 0;

    double dela = 0.051;
    //for(double a=2.401;a<=16.03;a=a+dela)
    vector<double> epsvec(5);

    {

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        double sinitial = 8.99;//6.5;//8.3;//5.0;//8.974;//8.6;//8.975;//8.96815;//7.5;//8.7;//7.5;//8.75;//1.001;//gvalright - 0.5;//(gvalleft + (gvalleft+gvalright)/2.0)/2.0;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
        double sfinal = 9.01;//real(phibthreshold(a,m));//8.82;
        double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
        points1 = N;
        //string filename="Mphib_secondsheet_acount_" + to_string(acount) + ".dat";
        string filename="dS_sheet2_full_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_3d.dat";
        
        //string filename="testGL_Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_2.dat";
        //string filename="tmp.dat";
        
        double simaginitial = -0.0101;
        double simagfinal = 0.01011;
        double delsimagpoints = 100;
        double delsimag = abs(simaginitial - simagfinal)/delsimagpoints;

        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        //double eta = 20.0;
        //test = 8.90465,0.01713

        for(int j=0;j<delsimagpoints+1;++j)
        {
            double s3imag = simaginitial + j*delsimag;
       
        for(int i=0;i<delspoints+1;++i)
        //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
        { 
            double s3real = sinitial + i*dels;

            comp s =  s3real + ii*s3imag;
            cout<<"am = "<<a<<endl;
            //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = 2.0*m*m;//sigmab(a,m);
            comp qval = pmom(s,sigb,m);
            cout<<"q = "<<qval<<endl;

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        

            //cout<<"s:"<<s<<endl;
            cout<<"run: "<<count+1<<endl;
            cout<<"a = "<<a<<'\t'<<" dela = "<<dela<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            
            int switch_for_gvec_fixer = 0; //0 means we have used contour and turn on gvec fixer
                                           //1 means we have used straight line and turn off gvec fixer
            
            
            int tag_for_m2k = 0;

            if(s3imag<0.0)
            {
                contour_for_resonance_4(qvec,weights,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k);

            }
            else if(s3imag>=0.0)
            {
                contour_for_resonance_6(qvec,weights,a,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k);
                
            }

            cout<<"tag at i = "<<tag_for_m2k<<endl;
            cout<<"qvec created with size = "<<qvec.size()<<endl;
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
          
            Eigen::VectorXcd interpolater_Gvec = Gvec;
            Gvec = -1.0*Gvec;

            

            if(imag(s)<0.0)
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }
            else if(imag(s)>=0.0)
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }

            //if(real(s)>=9.0 && imag(s)<=0.0)
            /*if(real(s)>=9.0)
            
            {
                double pi = acos(-1.0);
                double theta = -90.0*pi/180.0;
                //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
                //cout<<"m2k cut tagger = "<<tag_for_m2k<<endl;
                Bmat_maker_momrep_2eps_with_weights_rotated_m2k_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta);
                //Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }
            else if(real(s)<9.0 && imag(s)<=0.0)
            {
                //double theta = 0.0;
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                //Bmat_maker_momrep_2eps_with_weights_rotated_m2k_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta);
                
                Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }
            else
            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);
            */
            
            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        

            time(&time_start);
            LinearSolver_2(Bmat,dsol,Gvec,relerror);

            //cusolverComplex(Bmat,Gvec,dsol,size);
            time(&time_end);

            double time_taken = double(time_end - time_start);
            cout<<"Time taken to solve : "<<fixed 
                <<time_taken<<setprecision(5);
            cout<<" sec"<<endl;


            comp result;


            if(imag(s)<0.0)
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                interpolator_ds_integraleq_momrep_2eps_with_weights_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,tag_for_m2k,result);
            }
            else 
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                interpolator_ds_integraleq_momrep_2eps_with_weights_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,tag_for_m2k,result);
            }

            //if(real(s)>=9.0 && imag(s)<=0.0)
            /*if(real(s)>=9.0)
            {
                double pi = acos(-1.0);
                double theta = -90.0*pi/180.0;
                //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
                interpolator_ds_integraleq_momrep_2eps_with_weights_rotated_m2k_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,result);
            }
            else if(real(s)<9.0 && imag(s)<=0.0)
            {
                double theta = -90.0*pi/180.0;
                //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
                interpolator_ds_integraleq_momrep_2eps_with_weights_rotated_m2k_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,result);
            
            }
            else 
            interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            */

            //result = Bmat.determinant();
            
            //result = dsol[10];
            comp m2k = M2kfunc(a,sigb,m, eps_for_m2k);
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            comp ds = result;
            result = m2k*result*m2k;

            //result = Bmat.determinant();

        
            comp rhopb = rhophib(s,a,m);
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;

            cout<<setprecision(16)<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            fout<<setprecision(16)<<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(ds)<<'\t'
                <<imag(ds)<<'\t'
                <<real(mphib2)<<'\t'
                <<imag(mphib2)<<'\t'
                <<real(mphib2denom)<<'\t'
                <<imag(mphib2denom)<<'\t'
                <<a<<'\t'
                <<N<<'\t'
                <<gvalleft<<'\t'
                <<gvalright<<'\t'
                <<real(phibthreshold(a,m))<<'\t'
                <<eps<<endl;
            count = count + 1;
            
            cout<<"-------------------------"<<endl;
        }
        acount = acount + 1;
        
        }
        fout.close();

    }

    //}

}

void ds_belowthreshold_vs_s3_N_3d_with_weights_for_resonances_sheet_minus1_full()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = -9.10;//-6.40;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    //double s3imag = 0.0;//1.0e-5;
    

    double delspoints = 100.0;
    double points1 = 250.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    //double box = 10.0;
    int acount = 0;

    double dela = 0.051;
    //for(double a=2.401;a<=16.03;a=a+dela)
    vector<double> epsvec(5);

    {

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        double sinitial = 8.98;//6.5;//8.3;//5.0;//8.974;//8.6;//8.975;//8.96815;//7.5;//8.7;//7.5;//8.75;//1.001;//gvalright - 0.5;//(gvalleft + (gvalleft+gvalright)/2.0)/2.0;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
        double sfinal = 9.02;//real(phibthreshold(a,m));//8.82;
        double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
        points1 = N;
        //string filename="Mphib_secondsheet_acount_" + to_string(acount) + ".dat";
        string filename="dS_sheet-1_full_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_3d.dat";
        
        //string filename="testGL_Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_2.dat";
        //string filename="tmp.dat";
        
        double simaginitial = -0.01501;
        double simagfinal = 0.015011;
        double delsimagpoints = 100;
        double delsimag = abs(simaginitial - simagfinal)/delsimagpoints;

        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        //double eta = 20.0;
        //test = 8.90465,0.01713

        for(int j=0;j<delsimagpoints+1;++j)
        {
            double s3imag = simaginitial + j*delsimag;
       
        for(int i=0;i<delspoints+1;++i)
        //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
        { 
            double s3real = sinitial + i*dels;

            comp s =  s3real + ii*s3imag;
            cout<<"am = "<<a<<endl;
            //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = 2.0*m*m;//sigmab(a,m);
            comp qval = pmom(s,sigb,m);
            cout<<"q = "<<qval<<endl;

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        

            //cout<<"s:"<<s<<endl;
            cout<<"run: "<<count+1<<endl;
            cout<<"a = "<<a<<'\t'<<" dela = "<<dela<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            
            int switch_for_gvec_fixer = 0; //0 means we have used contour and turn on gvec fixer
                                           //1 means we have used straight line and turn off gvec fixer
            
            
            int tag_for_m2k = 0;

            if(s3imag<0.0)
            {
                contour_for_resonance_for_sheet_minus1_imsneg(qvec,weights,a,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k);
                
                for(int qind=0;qind<qvec.size();++qind)
                {
                    comp reqvec = real(qvec[qind]);
                    comp imqvec = -imag(qvec[qind]);
                    qvec[qind] = reqvec + ii*imqvec;

                    comp reweight = real(weights[qind]);
                    comp imweight = -imag(weights[qind]);
                    weights[qind] = reweight + ii*imweight;
                }
                //contour_for_resonance_4(qvec,weights,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k);

            }
            else if(s3imag>=0.0)
            {
                contour_for_resonance_for_sheet_minus1_imspos(qvec,weights,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k);
            
                for(int qind=0;qind<qvec.size();++qind)
                {
                    comp reqvec = real(qvec[qind]);
                    comp imqvec = -imag(qvec[qind]);
                    qvec[qind] = reqvec + ii*imqvec;

                    comp reweight = real(weights[qind]);
                    comp imweight = -imag(weights[qind]);
                    weights[qind] = reweight + ii*imweight;
                }
                //contour_for_resonance_6(qvec,weights,a,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k);
                
            }

            cout<<"tag at i = "<<tag_for_m2k<<endl;
            cout<<"qvec created with size = "<<qvec.size()<<endl;
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
          
            Eigen::VectorXcd interpolater_Gvec = Gvec;
            Gvec = -1.0*Gvec;

            

            if(imag(s)<0.0)
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }
            else if(imag(s)>=0.0)
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }

            //if(real(s)>=9.0 && imag(s)<=0.0)
            /*if(real(s)>=9.0)
            
            {
                double pi = acos(-1.0);
                double theta = -90.0*pi/180.0;
                //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
                //cout<<"m2k cut tagger = "<<tag_for_m2k<<endl;
                Bmat_maker_momrep_2eps_with_weights_rotated_m2k_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta);
                //Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }
            else if(real(s)<9.0 && imag(s)<=0.0)
            {
                //double theta = 0.0;
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                //Bmat_maker_momrep_2eps_with_weights_rotated_m2k_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta);
                
                Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }
            else
            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);
            */
            
            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        

            time(&time_start);
            LinearSolver_2(Bmat,dsol,Gvec,relerror);

            //cusolverComplex(Bmat,Gvec,dsol,size);
            time(&time_end);

            double time_taken = double(time_end - time_start);
            cout<<"Time taken to solve : "<<fixed 
                <<time_taken<<setprecision(5);
            cout<<" sec"<<endl;


            comp result;


            if(imag(s)<0.0)
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                interpolator_ds_integraleq_momrep_2eps_with_weights_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,tag_for_m2k,result);
            }
            else 
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                interpolator_ds_integraleq_momrep_2eps_with_weights_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,tag_for_m2k,result);
            }

            //if(real(s)>=9.0 && imag(s)<=0.0)
            /*if(real(s)>=9.0)
            {
                double pi = acos(-1.0);
                double theta = -90.0*pi/180.0;
                //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
                interpolator_ds_integraleq_momrep_2eps_with_weights_rotated_m2k_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,result);
            }
            else if(real(s)<9.0 && imag(s)<=0.0)
            {
                double theta = -90.0*pi/180.0;
                //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
                interpolator_ds_integraleq_momrep_2eps_with_weights_rotated_m2k_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,result);
            
            }
            else 
            interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            */

            //result = Bmat.determinant();
            
            //result = dsol[10];
            comp m2k = M2kfunc(a,sigb,m, eps_for_m2k);
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            comp ds = result;
            result = m2k*result*m2k;

            //result = Bmat.determinant();

        
            comp rhopb = rhophib(s,a,m);
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;

            cout<<setprecision(16)<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            fout<<setprecision(16)<<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(ds)<<'\t'
                <<imag(ds)<<'\t'
                <<real(mphib2)<<'\t'
                <<imag(mphib2)<<'\t'
                <<real(mphib2denom)<<'\t'
                <<imag(mphib2denom)<<'\t'
                <<a<<'\t'
                <<N<<'\t'
                <<gvalleft<<'\t'
                <<gvalright<<'\t'
                <<real(phibthreshold(a,m))<<'\t'
                <<eps<<endl;
            count = count + 1;
            
            cout<<"-------------------------"<<endl;
        }
        acount = acount + 1;
        
        }
        fout.close();

    }

    //}

}


void ds_belowthreshold_vs_s3_N_3d_with_weights_for_resonances_secondsheet_full_a_loop()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = -8.10;//-6.40;
    double m = 1.0;

    //cout<<"am = "<<a<<endl;
    //double s = 8.65;

    //double s3imag = 0.0;//1.0e-5;
    

    double delspoints = 100.0;
    double points1 = 250.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    //double box = 10.0;
    //int acount = 0;

    //for(double a=2.401;a<=16.03;a=a+dela)
    vector<double> epsvec(5);

    double ainitial = -25.1;
    double afinal = -102.3;
    double apoints = 40.0;
    double dela = abs(ainitial - afinal)/apoints;

    //for(int acount=0;acount<(int)apoints;++acount)
    {
        //a = ainitial - acount*dela;
        int acount = 2;
        a = -81.07;
        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        double sinitial = 8.9975;//6.5;//8.3;//5.0;//8.974;//8.6;//8.975;//8.96815;//7.5;//8.7;//7.5;//8.75;//1.001;//gvalright - 0.5;//(gvalleft + (gvalleft+gvalright)/2.0)/2.0;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
        double sfinal = 8.9990;//9.001;//25;//real(phibthreshold(a,m));//8.82;
        double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
        points1 = N;
        string filename="dS_sheet2_acount_" + to_string(acount) + ".dat";
        //string filename="dS_sheet2_full_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_3d.dat";
        
        //string filename="testGL_Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_2.dat";
        //string filename="tmp.dat";
        
        double simaginitial = -0.000025106;//-0.010051;
        double simagfinal = 0.000075151;//0.0100511;
        double delsimagpoints = 100;
        double delsimag = abs(simaginitial - simagfinal)/delsimagpoints;

        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        //double eta = 20.0;
        //test = 8.90465,0.01713

        for(int j=0;j<delsimagpoints+1;++j)
        {
            double s3imag = simaginitial + j*delsimag;
       
        for(int i=0;i<delspoints+1;++i)
        //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
        { 
            double s3real = sinitial + i*dels;

            comp s =  s3real + ii*s3imag;
            //cout<<"am = "<<a<<endl;
            //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = 2.0*m*m;//sigmab(a,m);
            comp qval = pmom(s,sigb,m);
            //cout<<"q = "<<qval<<endl;

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        

            //cout<<"s:"<<s<<endl;
            //cout<<"run: "<<count+1<<endl;
            //cout<<"a = "<<a<<'\t'<<" dela = "<<dela<<endl;
            //cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            //cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            
            int switch_for_gvec_fixer = 0; //0 means we have used contour and turn on gvec fixer
                                           //1 means we have used straight line and turn off gvec fixer
            
            
            int tag_for_m2k = 0;

            if(s3imag<0.0)
            {
                contour_for_resonance_4(qvec,weights,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k);

            }
            else if(s3imag>=0.0)
            {
                contour_for_resonance_6(qvec,weights,a,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k);
                
            }

            //cout<<"tag at i = "<<tag_for_m2k<<endl;
            //cout<<"qvec created with size = "<<qvec.size()<<endl;
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
          
            Eigen::VectorXcd interpolater_Gvec = Gvec;
            Gvec = -1.0*Gvec;

            

            if(imag(s)<0.0)
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }
            else if(imag(s)>=0.0)
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }

            //if(real(s)>=9.0 && imag(s)<=0.0)
            /*if(real(s)>=9.0)
            
            {
                double pi = acos(-1.0);
                double theta = -90.0*pi/180.0;
                //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
                //cout<<"m2k cut tagger = "<<tag_for_m2k<<endl;
                Bmat_maker_momrep_2eps_with_weights_rotated_m2k_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta);
                //Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }
            else if(real(s)<9.0 && imag(s)<=0.0)
            {
                //double theta = 0.0;
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                //Bmat_maker_momrep_2eps_with_weights_rotated_m2k_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta);
                
                Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }
            else
            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);
            */
            
            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        

            time(&time_start);
            LinearSolver_2(Bmat,dsol,Gvec,relerror);

            //cusolverComplex(Bmat,Gvec,dsol,size);
            time(&time_end);

            double time_taken = double(time_end - time_start);
            //cout<<"Time taken to solve : "<<fixed 
            //    <<time_taken<<setprecision(5);
            //cout<<" sec"<<endl;


            comp result;


            if(imag(s)<0.0)
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                interpolator_ds_integraleq_momrep_2eps_with_weights_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,tag_for_m2k,result);
            }
            else 
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                interpolator_ds_integraleq_momrep_2eps_with_weights_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,tag_for_m2k,result);
            }

            //if(real(s)>=9.0 && imag(s)<=0.0)
            /*if(real(s)>=9.0)
            {
                double pi = acos(-1.0);
                double theta = -90.0*pi/180.0;
                //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
                interpolator_ds_integraleq_momrep_2eps_with_weights_rotated_m2k_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,result);
            }
            else if(real(s)<9.0 && imag(s)<=0.0)
            {
                double theta = -90.0*pi/180.0;
                //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
                interpolator_ds_integraleq_momrep_2eps_with_weights_rotated_m2k_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,result);
            
            }
            else 
            interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            */

            //result = Bmat.determinant();
            
            //result = dsol[10];
            comp m2k = M2kfunc(a,sigb,m, eps_for_m2k);
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            comp ds = result;
            result = m2k*result*m2k;

            //result = Bmat.determinant();

        
            comp rhopb = rhophib(s,a,m);
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;

            //cout<<setprecision(16)<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            fout<<setprecision(16)<<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(ds)<<'\t'
                <<imag(ds)<<'\t'
                <<real(mphib2)<<'\t'
                <<imag(mphib2)<<'\t'
                <<real(mphib2denom)<<'\t'
                <<imag(mphib2denom)<<'\t'
                <<a<<'\t'
                <<N<<'\t'
                <<gvalleft<<'\t'
                <<gvalright<<'\t'
                <<real(phibthreshold(a,m))<<'\t'
                <<eps<<endl;
            count = count + 1;
            
            //cout<<"-------------------------"<<endl;
        }
        //acount = acount + 1;
        
        }
        cout<<"a = "<<a<<", file = "<<filename<<" done"<<endl;
        fout.close();

    }

    //}

}


void ds_belowthreshold_vs_s3_N_3d_with_weights_for_resonances_firstsheet()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = -8.10;//-6.40;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    //double s3imag = 0.0;//1.0e-5;
    

    double delspoints = 100.0;
    double points1 = 250.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    //double box = 10.0;
    int acount = 0;

    double dela = 0.051;
    //for(double a=2.401;a<=16.03;a=a+dela)
    vector<double> epsvec(5);

    {

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        double sinitial = 8.99;//6.5;//8.3;//5.0;//8.974;//8.6;//8.975;//8.96815;//7.5;//8.7;//7.5;//8.75;//1.001;//gvalright - 0.5;//(gvalleft + (gvalleft+gvalright)/2.0)/2.0;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
        double sfinal = 9.01;//real(phibthreshold(a,m));//8.82;
        double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
        points1 = N;
        //string filename="Mphib_secondsheet_acount_" + to_string(acount) + ".dat";
        string filename="ds_sheetI_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_3d.dat";
        
        //string filename="testGL_Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_2.dat";
        //string filename="tmp.dat";
        
        double simaginitial = -0.0101;
        double simagfinal = 0.01011;
        double delsimagpoints = 100;
        double delsimag = abs(simaginitial - simagfinal)/delsimagpoints;

        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        //double eta = 20.0;
        //test = 8.90465,0.01713

        for(int j=0;j<delsimagpoints+1;++j)
        {
            double s3imag = simaginitial + j*delsimag;
       
        for(int i=0;i<delspoints+1;++i)
        //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
        { 
            double s3real = sinitial + i*dels;

            comp s =  s3real + ii*s3imag;
            cout<<"am = "<<a<<endl;
            //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = 2.0*m*m;//sigmab(a,m);
            comp qval = pmom(s,sigb,m);
            cout<<"q = "<<qval<<endl;

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        

            //cout<<"s:"<<s<<endl;
            cout<<"run: "<<count+1<<endl;
            cout<<"a = "<<a<<'\t'<<" dela = "<<dela<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            
            int switch_for_gvec_fixer = 0; //0 means we have used contour and turn on gvec fixer
                                           //1 means we have used straight line and turn off gvec fixer
            
            
            int tag_for_m2k = 0;

            if(s3imag<0.0)
            {
                //contour_for_resonance_4(qvec,weights,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k);
                line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                
            }
            else if(s3imag>=0.0)
            {
                line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                //contour_for_resonance_3(qvec,weights,s,m,kmin,kmax,eps_for_m2k,-0.1,points1);
                
            }

            cout<<"tag at i = "<<tag_for_m2k<<endl;
            cout<<"qvec created with size = "<<qvec.size()<<endl;
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
          
            Eigen::VectorXcd interpolater_Gvec = Gvec;
            Gvec = -1.0*Gvec;

            

            if(imag(s)<0.0)
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

                //Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }
            else if(imag(s)>=0.0)
            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);


            //if(real(s)>=9.0 && imag(s)<=0.0)
            /*if(real(s)>=9.0)
            
            {
                double pi = acos(-1.0);
                double theta = -90.0*pi/180.0;
                //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
                //cout<<"m2k cut tagger = "<<tag_for_m2k<<endl;
                Bmat_maker_momrep_2eps_with_weights_rotated_m2k_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta);
                //Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }
            else if(real(s)<9.0 && imag(s)<=0.0)
            {
                //double theta = 0.0;
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                //Bmat_maker_momrep_2eps_with_weights_rotated_m2k_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta);
                
                Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }
            else
            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);
            */
            
            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        

            time(&time_start);
            LinearSolver_2(Bmat,dsol,Gvec,relerror);

            //cusolverComplex(Bmat,Gvec,dsol,size);
            time(&time_end);

            double time_taken = double(time_end - time_start);
            cout<<"Time taken to solve : "<<fixed 
                <<time_taken<<setprecision(5);
            cout<<" sec"<<endl;


            comp result;


            if(imag(s)<0.0)
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                interpolator_ds_integraleq_momrep_2eps_with_weights_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,tag_for_m2k,result);
            }
            else 
            interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            

            //if(real(s)>=9.0 && imag(s)<=0.0)
            /*if(real(s)>=9.0)
            {
                double pi = acos(-1.0);
                double theta = -90.0*pi/180.0;
                //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
                interpolator_ds_integraleq_momrep_2eps_with_weights_rotated_m2k_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,result);
            }
            else if(real(s)<9.0 && imag(s)<=0.0)
            {
                double theta = -90.0*pi/180.0;
                //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
                interpolator_ds_integraleq_momrep_2eps_with_weights_rotated_m2k_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,result);
            
            }
            else 
            interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            */

            //result = Bmat.determinant();
            
            //result = dsol[10];
            comp m2k = M2kfunc(a,sigb,m, eps_for_m2k);
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            comp ds = result;
            result = m2k*result*m2k;

            //result = Bmat.determinant();

        
            comp rhopb = rhophib(s,a,m);
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;

            cout<<setprecision(16)<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            fout<<setprecision(16)<<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(ds)<<'\t'
                <<imag(ds)<<'\t'
                <<real(mphib2)<<'\t'
                <<imag(mphib2)<<'\t'
                <<real(mphib2denom)<<'\t'
                <<imag(mphib2denom)<<'\t'
                <<a<<'\t'
                <<N<<'\t'
                <<gvalleft<<'\t'
                <<gvalright<<'\t'
                <<real(phibthreshold(a,m))<<'\t'
                <<eps<<endl;
            count = count + 1;
            
            cout<<"-------------------------"<<endl;
        }
        acount = acount + 1;
        
        }
        fout.close();

    }

    //}

}




void Mphib_firstsheet_abovethreshold_vs_s3_N_fixed_s3imag_with_weights()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = 16;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    double s3imag = 0.0;//1.0e-5;
    double eta = 25.0;
    

    double delspoints = 2000.0;
    double points1 = 500.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    //double box = 10.0;
    int acount = 0;

    double dela = 0.051;
    //for(double a=2.401;a<=16.03;a=a+dela)
    vector<double> epsvec(5);
    epsvec[0] = 0.000010;
    epsvec[1] = 0.000007;
    epsvec[2] = 0.000005;
    epsvec[3] = 0.000003;
    epsvec[4] = 0.000001;

    //for(int epsi=0;epsi<epsvec.size();++epsi)
    //{

    //    eps = epsvec[epsi];

    //int N = 500;
    //for(int N=900;N<=1000;N=N+20)
    {
        //if(a>=3.7) dela = 0.101;
        //else if(a>=17.0) dela = 0.51;
        //else if(a>=1000.0) dela = 300.01;

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        double sinitial = real(phibthreshold(a,m));//8.75;//5.0;//8.974;//8.6;//8.975;//8.96815;//7.5;//8.7;//7.5;//8.75;//1.001;//gvalright - 0.5;//(gvalleft + (gvalleft+gvalright)/2.0)/2.0;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
        double sfinal = 9.0*m*m;//real(phibthreshold(a,m));//8.82;
        double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
        points1 = N;
        //string filename="Mphib_secondsheet_acount_" + to_string(acount) + ".dat";
        string filename="Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_abovethreshold.dat";
        
        //string filename="testGL_Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_2.dat";
        //string filename="tmp.dat";
        

        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        //double eta = 20.0;

        for(int i=0;i<delspoints+1;++i)
        //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
        { 
            double s3real = sinitial + i*dels;
            //double s = 8.78;
            //if(s>=sfinal) dels = 0.05/100.0;
            //eps = eps_above_threshold(eta,s,a,m,N);
            //eps_for_m2k = eps;

            comp s = s3real + ii*s3imag;
            eps = eps_above_threshold(eta,s,a,m,N);
            eps_for_m2k = eps;

            cout<<"am = "<<a<<endl;
            //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = sigmab(a,m);
            comp qval = pmom(s,sigb,m);
            cout<<"q = "<<qval<<endl;

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        

            //cout<<"s:"<<s<<endl;
            cout<<"run: "<<count+1<<endl;
            cout<<"a = "<<a<<'\t'<<" dela = "<<dela<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            
            int switch_for_gvec_fixer = 0; //0 means we have used contour and turn on gvec fixer
                                           //1 means we have used straight line and turn off gvec fixer
            
            line_maker_with_weights(qvec,weights,kmin,kmax,points1);
            
            if(s3imag<=0.0)
            {
                //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                cout<<"contour_43 chosen"<<endl;
                //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1);
                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
            }
            else if(s3imag>0.0)
            {
                //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                cout<<"contour_47 chosen"<<endl;
                //mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
                int changed_points = points1 + 5;
                vector<comp> tmp_qvec;
                vector<comp> tmp_weights;
                //mom_vector_maker_seba_imspos(tmp_qvec,tmp_weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)changed_points,tag1,tag2,switch_for_gvec_fixer);
                //mom_vector_maker_seba_imspos(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                //mom_vector_maker_47_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
                //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                //mom_vector_maker_seba_imspos_2(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
            
                //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                
                /*for(int i=0;i<tmp_qvec.size()-5;++i)
                {
                    qvec.push_back(tmp_qvec[i]);
                    weights.push_back(tmp_weights[i]);
                }*/
                cout<<"tag1 = "<<tag1<<endl;
                cout<<"tag2 = "<<tag2<<endl;
                
                
            }

            cout<<"qvec created with size = "<<qvec.size()<<endl;
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            
            if(s3imag<=0.0)
            {
                Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            }
            else if(s3imag>0.0)
            {
                Gvec_maker_momrep_withtags_1(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            }

            switch_for_gvec_fixer=1;
            if(switch_for_gvec_fixer==0)
            {
                //Gvec_fixer_1(Gvec,Gvec,qvec,s,qval,m,eps,tag1,tag2); 
                //Gvec_fixer_2(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
                //Gvec_fixer_3(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
                //tag_fixer_for_gvec_fixer3(Gvec,qvec,s,qval,m,eps,tag1);
               
            }

            cout<<"tag fixer = "<<tag1<<'\t'<<tag2<<endl;
            //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            Eigen::VectorXcd interpolater_Gvec = Gvec;
            Gvec = -1.0*Gvec;

            //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
            //Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
            //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);
            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

            
            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        
            //cout<<"did this = "<<count<<" times"<<endl;
            //cout<<"i val = "<<i<<'\t';
            //cout<<"j val = "<<j<<endl;
            
            //LinearSolver_2(Bmat,dsol,Gvec,relerror);

            time(&time_start);
            LinearSolver_2(Bmat,dsol,Gvec,relerror);

            //cusolverComplex(Bmat,Gvec,dsol,size);
            time(&time_end);

            double time_taken = double(time_end - time_start);
            cout<<"Time taken to solve : "<<fixed 
                <<time_taken<<setprecision(5);
            cout<<" sec"<<endl;


            comp result;

        
            if(s3imag<=0.0)
            {
                //interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
                interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            }
            else if(s3imag>0.0)
            {
                //interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
                //interpolator_ds_integraleq_momrep_2eps_withtags_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
                interpolator_ds_integraleq_momrep_2eps_withtags_with_weights_usingGvec(dsol,qvec,weights,interpolater_Gvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
            
            }
            //result = dsol[10];
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            comp ds = result;
            result = gsq*result;
            comp mphib = result;

        
            comp rhopb = rhophib(s,a,m);
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;
            comp qcotdel = 8.0*pi*sqrt(s)*real(1.0/mphib); 
            double error = abs((imag(1.0/mphib) + rhopb)/rhopb)*100.0;

            cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            cout<<"ERROR = "<<error<<endl;
            fout<<setprecision(16)<<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(qcotdel)<<'\t'
                <<imag(qcotdel)<<'\t'
                <<real(mphib2)<<'\t'
                <<imag(mphib2)<<'\t'
                <<real(qval*qval)<<'\t'
                <<imag(qval*qval)<<'\t'
                <<a<<'\t'
                <<N<<'\t'
                <<gvalleft<<'\t'
                <<gvalright<<'\t'
                <<real(phibthreshold(a,m))<<'\t'
                <<eps<<'\t'
                <<error<<endl;
            count = count + 1;
            
            cout<<"-------------------------"<<endl;
        }
        acount = acount + 1;
        fout.close();

    }

    //}

}

void Mphib_firstsheet_abovethreshold_SAmethod_vs_s3_N_fixed_s3imag_with_weights()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = 2;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    double s3imag = 1.0e-5;
    double eta = 15.0;
    

    double delspoints = 200.0;
    double points1 = 500.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    //double box = 10.0;
    int acount = 0;

    double dela = 0.051;
    //for(double a=2.401;a<=16.03;a=a+dela)
    vector<double> epsvec(5);
    epsvec[0] = 0.000010;
    epsvec[1] = 0.000007;
    epsvec[2] = 0.000005;
    epsvec[3] = 0.000003;
    epsvec[4] = 0.000001;

    //for(int epsi=0;epsi<epsvec.size();++epsi)
    //{

    //    eps = epsvec[epsi];

    //int N = 5000;
    //for(int N=1000;N<=6000;N=N+1000)
    {
        //if(a>=3.7) dela = 0.101;
        //else if(a>=17.0) dela = 0.51;
        //else if(a>=1000.0) dela = 300.01;

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        double sinitial = real(phibthreshold(a,m));//8.75;//5.0;//8.974;//8.6;//8.975;//8.96815;//7.5;//8.7;//7.5;//8.75;//1.001;//gvalright - 0.5;//(gvalleft + (gvalleft+gvalright)/2.0)/2.0;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
        double sfinal = 9.0*m*m;//real(phibthreshold(a,m));//8.82;
        double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
        points1 = N;
        //string filename="Mphib_secondsheet_acount_" + to_string(acount) + ".dat";
        string filename="Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_abovethreshold_SA_eps0.dat";
        
        //string filename="testGL_Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_2.dat";
        //string filename="tmp.dat";
        

        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        //double eta = 20.0;

        for(int i=0;i<delspoints+1;++i)
        //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
        { 
            double s3real = sinitial + i*dels;
            //double s = 8.78;
            //if(s>=sfinal) dels = 0.05/100.0;
            //eps = eps_above_threshold(eta,s,a,m,N);
            //eps_for_m2k = eps;

            comp s = s3real + ii*s3imag;
            eps = eps_momrep_above_threshold(eta,s,a,m,N);
            eps_for_m2k = eps;

            cout<<"am = "<<a<<endl;
            //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = sigmab(a,m);
            comp qval = pmom(s,sigb,m);
            cout<<"q = "<<qval<<endl;

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        

            //cout<<"s:"<<s<<endl;
            cout<<"run: "<<count+1<<endl;
            cout<<"a = "<<a<<'\t'<<" dela = "<<dela<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            
            int switch_for_gvec_fixer = 0; //0 means we have used contour and turn on gvec fixer
                                           //1 means we have used straight line and turn off gvec fixer
            
            line_maker_with_weights(qvec,weights,kmin,kmax,points1);
            
            
            cout<<"qvec created with size = "<<qvec.size()<<endl;
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);
            Eigen::VectorXcd Kphib(size);

            
           
            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);

            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            
            Gvec = -1.0*gsq*Gvec;

            deltaBmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);
            double relerror;
            

            time(&time_start);
            LinearSolver_2(Bmat,Kphib,Gvec,relerror);

            //cusolverComplex(Bmat,Gvec,Kphib,size);
            time(&time_end);

            double time_taken = double(time_end - time_start);
            cout<<"Time taken to solve : "<<fixed 
                <<time_taken<<setprecision(5);
            cout<<" sec"<<endl;


            comp result;

        
            interpolator_Kphib_integraleq_momrep_2eps_with_weights(Kphib,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);

            //result = dsol[10];
            
            comp ds = result;
            //result = gsq*result;
            

        
            comp rhopb = rhophib(s,a,m);
            comp mphib = result/(1.0 - ii*result*rhopb);
            comp mphib2 = mphib/(1.0 + 2.0*ii*rhopb*mphib);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*mphib);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;
            comp qcotdel = 8.0*pi*sqrt(s)*real(1.0/mphib); 
            double error = abs((imag(1.0/mphib) + rhopb)/rhopb)*100.0;

            cout<<setprecision(16);
            cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            cout<<"rhoM = "<<rhopb*mphib<<endl;
            cout<<"imMinv = "<<imag(1.0/mphib)<<'\t'<<" rho = "<<rhopb<<endl;
            cout<<"ERROR = "<<error<<endl;
            fout<<setprecision(16)<<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(mphib)<<'\t'
                <<imag(mphib)<<'\t'
                <<real(qcotdel)<<'\t'
                <<imag(qcotdel)<<'\t'
                <<real(rhopb*mphib)<<'\t'
                <<imag(rhopb*mphib)<<'\t'
                <<real(qval*qval)<<'\t'
                <<imag(qval*qval)<<'\t'
                <<a<<'\t'
                <<N<<'\t'
                <<gvalleft<<'\t'
                <<gvalright<<'\t'
                <<real(phibthreshold(a,m))<<'\t'
                <<eps<<'\t'
                <<error<<endl;
            count = count + 1;
            
            cout<<"-------------------------"<<endl;
        }
        acount = acount + 1;
        fout.close();

    }

    //}

}


void dSqqsigma2msq_belowthreshold_vs_s3_N_fixed_s3imag_with_weights()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    double s3imag = 1.0e-5;
    

    double delspoints = 2000.0;
    double points1 = 500.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    //double box = 10.0;
    int acount = 0;

    double dela = 0.051;
    //for(double a=2.401;a<=16.03;a=a+dela)
    vector<double> epsvec(5);
    epsvec[0] = 0.000010;
    epsvec[1] = 0.000007;
    epsvec[2] = 0.000005;
    epsvec[3] = 0.000003;
    epsvec[4] = 0.000001;

    //for(int epsi=0;epsi<epsvec.size();++epsi)
    //{

    //    eps = epsvec[epsi];

    //for(int N=900;N<=1000;N=N+20)
    {
        //if(a>=3.7) dela = 0.101;
        //else if(a>=17.0) dela = 0.51;
        //else if(a>=1000.0) dela = 300.01;

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        double sinitial = 8.975;//8.7;//7.5;//3.5;//7.5;//8.7;//gvalright - 0.5;//(gvalleft + (gvalleft+gvalright)/2.0)/2.0;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
        double sfinal = real(phibthreshold(a,m));//8.82;
        double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
        points1 = N;
        //string filename="Mphib_secondsheet_acount_" + to_string(acount) + ".dat";
        //string filename="Mphib_a_" + to_string((int)a) + "_N_" + to_string(N) + "_epscount_" + to_string(epsi) + ".dat";
        
        string filename="dsqq2msq_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_2.dat";
        

        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        //double eta = 20.0;

        for(int i=0;i<delspoints+1;++i)
        //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
        { 
            double s3real = sinitial + i*dels;
            //double s = 8.78;
            //if(s>=sfinal) dels = 0.05/100.0;
            //eps = eps_above_threshold(eta,s,a,m,points);

            comp s = s3real + ii*s3imag;
            cout<<"am = "<<a<<endl;
            //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = 2.0*m*m;//sigmab(a,m);
            comp qval = pmom(s,sigb,m);

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        

            //cout<<"s:"<<s<<endl;
            cout<<"run: "<<count+1<<endl;
            cout<<"a = "<<a<<'\t'<<" dela = "<<dela<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            
            line_maker_with_weights(qvec,weights,kmin,kmax,points1);
            
            
            //if(s3imag<=0.0)
            //{
                //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //cout<<"contour_43 chosen"<<endl;
                //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1);
                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
            //}
            //else if(s3imag>0.0)
            //{
                //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //cout<<"contour_47 chosen"<<endl;
                //mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
                //cout<<"tag1 = "<<tag1<<endl;
                //cout<<"tag2 = "<<tag2<<endl;
                
                
            //}

            cout<<"qvec created with size = "<<qvec.size()<<endl;
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            
            if(s3imag<=0.0)
            {
                Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            }
            else if(s3imag>0.0)
            {
                Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            }
            //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            Gvec = -1.0*Gvec;

            //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
            //Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
            //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);
            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

            
            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        
            //cout<<"did this = "<<count<<" times"<<endl;
            //cout<<"i val = "<<i<<'\t';
            //cout<<"j val = "<<j<<endl;
            
            //LinearSolver_2(Bmat,dsol,Gvec,relerror);

            time(&time_start);
            LinearSolver_2(Bmat,dsol,Gvec,relerror);

            //cusolverComplex(Bmat,Gvec,dsol,size);
            time(&time_end);

            double time_taken = double(time_end - time_start);
            cout<<"Time taken to solve : "<<fixed 
                <<time_taken<<setprecision(5);
            cout<<" sec"<<endl;


            comp result;

        
            if(s3imag<=0.0)
            {
                //interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
                interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            }
            else if(s3imag>0.0)
            {
                interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
            }
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            comp ds = result;
            result = gsq*result;

        
            comp rhopb = rhophib(s,a,m);
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;

            cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            fout<<setprecision(16)<<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(ds)<<'\t'
                <<imag(ds)<<'\t'
                <<real(mphib2)<<'\t'
                <<imag(mphib2)<<'\t'
                <<real(mphib2denom)<<'\t'
                <<imag(mphib2denom)<<'\t'
                <<a<<'\t'
                <<N<<'\t'
                <<gvalleft<<'\t'
                <<gvalright<<'\t'
                <<real(phibthreshold(a,m))<<'\t'
                <<eps<<endl;
            count = count + 1;
            
            cout<<"-------------------------"<<endl;
        }
        acount = acount + 1;
        fout.close();

    }

    //}

}


void Mphib_belowthreshold_vs_s3_fixed_s3imag_contourlinear()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    double s3imag = 0.05;
    double sinitial = 8.72;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;

    double delspoints = 100.0;
    double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
    double points1 = 5000.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.01;
    //double box = 10.0;

    string filename="fixeds3imag_Mphib_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_linear.dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    //double eta = 20.0;

    for(int i=0;i<delspoints+1;++i)
    //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    { 
        double s3real = sinitial + i*dels;
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);

        comp s = s3real + ii*s3imag;
        cout<<"am = "<<a<<endl;
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp qval = pmom(s,sigb,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        
        //cout<<"s:"<<s<<endl;
        cout<<"run: "<<count+1<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;

        comp delk = (kmax-kmin)/points1;

        int tag1,tag2;
            
        mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);    
        
        //if(s3imag<=0.0)
        //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        //else if(s3imag>0.0)
        //{
            //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        //    mom_vector_maker_45(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
        //    cout<<"tag1 = "<<tag1<<endl;
        //    cout<<"tag2 = "<<tag2<<endl;
                
                
        //}

        cout<<"qvec created with size = "<<qvec.size()<<endl;
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            
        
        //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
        Gvec = -1.0*Gvec;

        //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
        Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
        //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);

            
        //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
        double relerror;
            
        
        //cout<<"did this = "<<count<<" times"<<endl;
        //cout<<"i val = "<<i<<'\t';
        //cout<<"j val = "<<j<<endl;
            
        //LinearSolver_2(Bmat,dsol,Gvec,relerror);

        time(&time_start);
        cusolverComplex(Bmat,Gvec,dsol,size);
        time(&time_end);

        double time_taken = double(time_end - time_start);
        cout<<"Time taken to solve : "<<fixed 
            <<time_taken<<setprecision(5);
        cout<<" sec"<<endl;
        
            
        comp result;

        interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
        
        
        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        comp ds = result;
        result = gsq*result;

        
        comp rhopb = rhophib(s,a,m);
        comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
        comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
        comp GSval = GS_pk(s,qval,qval,m,eps);
        //result = rhophib(s,a,m)*result;

        cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        fout<<setprecision(16)<<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<real(result)<<'\t'
            <<imag(result)<<'\t'
            <<real(ds)<<'\t'
            <<imag(ds)<<'\t'
            <<real(mphib2)<<'\t'
            <<imag(mphib2)<<'\t'
            <<real(mphib2denom)<<'\t'
            <<imag(mphib2denom)<<'\t'
            <<a<<'\t'
            <<N<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}


void Mphib_belowthreshold_vs_s3_for_different_am()
{
    //double a = 500.0;
    double m = 1.0;

    //cout<<"am = "<<a<<endl;
    //double s = 8.65;
    
    double points1 = 2000.0;
    //double points2 = 500.0;
    double eps = 0.0;//1.0e-4;
    double eps_for_m2k = eps;
    double qvec_r = 0.01;
    //double box = 10.0;

    int totalapoints = 50;
    double ainitial = 1.38;
    double afinal = 1.45;
    double dela = abs(ainitial - afinal)/totalapoints;

    for(int i=0;i<totalapoints;++i)
    {
    
        string filename="Mphib_momrep_acount_" + to_string(i) 
                        + "_N_" + to_string((int)points1)  
                        + ".dat";
        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        double a = ainitial + i*dela;


        double delspoints = 200.0;
        double sinitial = 4.5;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
        double sfinal = real(phibthreshold(a,m));//8.82;
        double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//

        //double eta = 20.0;

        for(int j=0;j<delspoints+1;++j)
        //for(double s=sinitial;s<=sfinal;s=s+dels)
        { 
            double s = sinitial + j*dels;
            //double s = 8.78;
            //if(s>=sfinal) dels = 0.05/100.0;
            //eps = eps_above_threshold(eta,s,a,m,points);
            cout<<"am = "<<a<<endl;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = sigmab(a,m);
            comp qval = pmom(s,sigb-ii*eps_for_m2k,m);

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        
            //cout<<"s:"<<s<<endl;
            cout<<"run: "<<count+1<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;
            cout<<"eps_for_m2k"<<eps_for_m2k<<endl;

            vector<comp> qvec;
            
            mom_vector_maker_42(qvec,s,kmin,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            //mom_vector_maker_4(qvec,s,kmin,kmax,a,m,eps,eps_for_m2k,points1,qvec_r);
            cout<<"qvec size = "<<qvec.size()<<endl;
    
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            Gvec = -1.0*Gvec;

            Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);

            double relerror;

        

            //LinearSolver_2(Bmat,dsol,Gvec,relerror);
            cusolverComplex(Bmat,Gvec,dsol,size);

            comp result;

            interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);

            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            result = gsq*result;

        
            comp rhopb = rhophib(s,a,m);
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;

            cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            fout<<setprecision(16)<<s<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                //<<real(1.0/result)<<'\t'
                //<<imag(1.0/result)<<'\t'
                <<real(mphib2)<<'\t'
                <<imag(mphib2)<<'\t'
                //<<real(result + GSval)<<'\t'
                //<<imag(result + GSval)<<'\t'
                <<real(mphib2denom)<<'\t'
                <<imag(mphib2denom)<<'\t'
                //<<real(1.0/(result + GSval))<<'\t'
                //<<imag(1.0/(result + GSval))<<endl;
                <<a<<'\t'
                <<eps<<'\t'
                <<eps_for_m2k<<endl;
            count = count + 1;
            cout<<"-------------------------"<<endl;
        }

        fout.close();
    }
}


void Mphib_belowthreshold_vs_s3_3d_omp()
{
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    //double sinitial = 8.88;//8.72;//real(phibthreshold(a,m));//
    //double sfinal = real(phibthreshold(a,m));//8.82;
    //double dels = abs(sinitial-sfinal)/400.0;
    double points1 = 1000.0;
    double s3realpoints = 5.0;
    double s3imagpoints = 5.0;
    double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.01;
    //double box = 10.0;

    string filename="Mphib_momrep_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_3d"
                    //+ "_epspow_" + to_string((int)abs(log10(eps))) 
                    + "_withoutomp.dat";
    ofstream fout;
    

    int count = 0;

    double eta = 20.0;

    double s3realinitial = 8.72;
    double s3realfinal = real(phibthreshold(a,m));
    double dels3real = abs(s3realinitial-s3realfinal)/s3realpoints;

    double s3imaginitial = -0.1;
    double s3imagfinal = 0.1;
    double dels3imag = abs(s3imaginitial - s3imagfinal)/s3imagpoints;
    //cout<<"dels3real = "<<dels3real<<'\t'<<"dels3imag = "<<dels3imag<<endl;

    //vector<vector<comp> > resultvec(points1,points1);
    comp **resultvec = new comp*[(int)s3imagpoints+1];
    for(int i=0; i<(int)s3imagpoints+1; ++i)
    {
        resultvec[i] = new comp[(int)s3realpoints+1];
    }

    //comp resultvec[(int)s3imagpoints+1][(int)s3realpoints+1];

    int i;
    //#pragma omp parallel for shared(resultvec)
    omp_set_num_threads(10);
    //#pragma omp distribute parallel for shared(resultvec,count)

    
    cudaStream_t nstreams[(int)s3imagpoints+1];
    
    #pragma omp parallel for shared(resultvec)
    for(i=0;i<(int)s3imagpoints+1;++i)
    {
        //for(double s=sinitial;s<=sfinal;s=s+dels)
        //std::cout<<"from Thread = "<<omp_get_thread_num<<endl;
        //std::cout<<"total Thread = "<<omp_get_num_threads<<endl;
        int printcount = 0;
        
        for(int j=0;j<(int)s3realpoints+1;++j)
        {
            double s3imag = s3imaginitial + i*dels3imag;
            double s3real = s3realinitial + j*dels3real;

            
            comp s = s3real + ii*s3imag;
            //double s = 8.78;
            //if(s>=sfinal) dels = 0.05/100.0;
            //eps = eps_above_threshold(eta,s,a,m,points);
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = sigmab(a,m);
            comp qval = pmom(s,sigb,m);

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        
            //cout<<"s:"<<s<<endl;
            //cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            //cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;

            comp delk = abs(kmin - kmax)/points1;
            mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
            //mom_vector_maker_4(qvec,s,kmin,kmax,a,m,eps,eps_for_m2k,points1,qvec_r);
    
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            
            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            Gvec = -1.0*Gvec;

            //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
            Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);

            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        
            //cout<<"did this = "<<count<<" times"<<endl;
            //cout<<"i val = "<<i<<'\t';
            //cout<<"j val = "<<j<<endl;
            
            //LinearSolver_2(Bmat,dsol,Gvec,relerror);
            

            
            cudaStream_t somethingstream;
            cudaStreamCreate(&somethingstream);

            
            //{
            //    cusolverComplex(Bmat,Gvec,dsol,size);
            //}
            //cout<<"stream no. for i = "<<i<<'\t'<<" is streamid = "<<&somethingstream<<endl;
            //cout<<"stream no. for i = "<<i<<'\t'<<" is streamid = "<<*nstreams[i]<<endl;

            //#pragma omp critical
            //{
                cusolverComplexAsync(Bmat,Gvec,dsol,size,somethingstream);
            //}
            comp result;

            cudaStreamDestroy(somethingstream);

            //std::cout<<"cuda work finished from i ="<<i<<'\t'<<" j = "<<j<<endl;

            interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);

            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            result = gsq*result;
            //result = rhophib(s,a,m)*result;
            
            resultvec[i][j] = result;

            #pragma omp critical
            {
                //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
                //fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<'\t'<<real(1.0/result)<<'\t'<<imag(1.0/result)<<endl;
                
                
                count = count + 1;
                cout<<count<<" out of "<<((s3realpoints+1)*(s3imagpoints+1))<<" runs finished"<<endl;
                //cout<<"run = "<<count<<endl;
                //cout<<"-------------------------"<<endl;
            }

            //if(printcount==0)
            //{
                
                //printcount = 1;
            //}
        }
    }

    cudaDeviceReset();

    //fout.open(filename.c_str());
    for(int k=0;k<s3imagpoints+1;++k)
    {
        for(int l=0;l<s3realpoints+1;++l)
        {
            //cout<<"k="<<k<<'\t'<<"l="<<l<<endl;
            //cout<<"dels3imag="<<dels3imag<<'\t'<<"dels3real="<<dels3real<<endl;
            double s3imag = s3imaginitial + k*dels3imag;
            double s3real = s3realinitial + l*dels3real;
            //cout<<"s3imag = "<<s3imag<<'\t'<<"s3real = "<<s3real<<endl;
            comp s = s3real + ii*s3imag;

            cout<<s3real<<'\t'<<s3imag<<'\t'<<real(resultvec[k][l])<<'\t'<<imag(resultvec[k][l])<<endl;
            //fout<<s3real<<'\t'<<s3imag<<'\t'<<real(resultvec[k][l])<<'\t'<<imag(resultvec[k][l])<<endl;
        }
    }

    //fout.close();

    for(int i=0; i<s3imagpoints+1; ++i)
    {
        delete [] resultvec[i];
    }

    delete [] resultvec;
}


void Mphib_belowthreshold_vs_s3_3d_omp_qvecsecondsheet()
{
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    //double sinitial = 8.88;//8.72;//real(phibthreshold(a,m));//
    //double sfinal = real(phibthreshold(a,m));//8.82;
    //double dels = abs(sinitial-sfinal)/400.0;
    double points1 = 1000.0;
    double s3realpoints = 100.0;
    double s3imagpoints = 100.0;
    double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.01;
    //double box = 10.0;

    string filename="Mphib_qtranspose_momrep_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_3d"
                    //+ "_epspow_" + to_string((int)abs(log10(eps))) 
                    + "_withHpk.dat";//_followq_poss.dat";
    ofstream fout;
    

    int count = 0;

    double eta = 20.0;

    double s3realinitial = 8.72;
    double s3realfinal = real(phibthreshold(a,m));
    double dels3real = abs(s3realinitial-s3realfinal)/s3realpoints;

    double s3imaginitial = -0.1;
    double s3imagfinal = 0.1;
    double dels3imag = abs(s3imaginitial - s3imagfinal)/s3imagpoints;
    //cout<<"dels3real = "<<dels3real<<'\t'<<"dels3imag = "<<dels3imag<<endl;

    //vector<vector<comp> > resultvec(points1,points1);
    comp **resultvec = new comp*[(int)s3imagpoints+1];
    for(int i=0; i<(int)s3imagpoints+1; ++i)
    {
        resultvec[i] = new comp[(int)s3realpoints+1];
    }

    comp **dsvec = new comp*[(int)s3imagpoints+1];
    for(int i=0; i<(int)s3imagpoints+1; ++i)
    {
        dsvec[i] = new comp[(int)s3realpoints+1];
    }

    //comp resultvec[(int)s3imagpoints+1][(int)s3realpoints+1];

    int i;
    //#pragma omp parallel for shared(resultvec)
    //omp_set_num_threads(10);
    //#pragma omp distribute parallel for shared(resultvec,count)

    
    //cudaStream_t nstreams[(int)s3imagpoints+1];
    
    //#pragma omp parallel for shared(resultvec)
    for(i=0;i<(int)s3imagpoints+1;++i)
    {
        //for(double s=sinitial;s<=sfinal;s=s+dels)
        //std::cout<<"from Thread = "<<omp_get_thread_num<<endl;
        //std::cout<<"total Thread = "<<omp_get_num_threads<<endl;
        int printcount = 0;
        
        for(int j=0;j<(int)s3realpoints+1;++j)
        {
            double s3imag = s3imaginitial + i*dels3imag;
            double s3real = s3realinitial + j*dels3real;

            comp sigk = 3.80*m*m;
            comp s = s3real + ii*s3imag;
            //double s = 8.78;
            //if(s>=sfinal) dels = 0.05/100.0;
            //eps = eps_above_threshold(eta,s,a,m,points);
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = sigmab(a,m);
            //comp qval = pmom(s,sigb,m);
            comp qval = pmom(s,sigb,m);
            
            //if(s3imag<=0.0) qval = -qval;
            cout<<"qval = "<<qval<<endl;

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        
            //cout<<"s:"<<s<<endl;
            //cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            //cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;

            comp delk = abs(kmin - kmax)/points1;
            int tag1,tag2;
            
            //mom_vector_maker_6(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,points1,qvec_r,tag1,tag2);
            //mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
            //mom_vector_maker_4(qvec,s,kmin,kmax,a,m,eps,eps_for_m2k,points1,qvec_r);
            if(s3imag<=0.0)
            //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            else if(s3imag>0.0)
            {
                //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                
                for(int qind=0;qind<qvec.size();++qind)
                {
                    comp reqvec = real(qvec[qind]);
                    comp imqvec = -imag(qvec[qind]);
                    qvec[qind] = reqvec + ii*imqvec;
                }
            }

            cout<<"qvec created with size = "<<qvec.size()<<endl;
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            
            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            Gvec = -1.0*Gvec;

            //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
            Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
            //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);

            
            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        
            //cout<<"did this = "<<count<<" times"<<endl;
            //cout<<"i val = "<<i<<'\t';
            //cout<<"j val = "<<j<<endl;
            
            //LinearSolver_2(Bmat,dsol,Gvec,relerror);
            

            
            //cudaStream_t somethingstream;
            //cudaStreamCreate(&somethingstream);

            //#pragma omp critical
            //{
            
            cusolverComplex(Bmat,Gvec,dsol,size);
            
            //}
            //cout<<"stream no. for i = "<<i<<'\t'<<" is streamid = "<<&somethingstream<<endl;
            //cout<<"stream no. for i = "<<i<<'\t'<<" is streamid = "<<*nstreams[i]<<endl;

            //#pragma omp critical
            //{
                //cusolverComplexAsync(Bmat,Gvec,dsol,size,somethingstream);
            //}
            comp result;

            //cudaStreamDestroy(somethingstream);

            //std::cout<<"cuda work finished from i ="<<i<<'\t'<<" j = "<<j<<endl;

            interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
            //interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
            comp m2kfunc = M2kfunc(a,sigk,m,eps_for_m2k);
            comp m2ksq = m2kfunc*m2kfunc;
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            dsvec[i][j] = result;
            comp dsres = result;
            result = gsq*result;
            //result = m2ksq*result;
            //result = rhophib(s,a,m)*result;
            
            resultvec[i][j] = result;

            //#pragma omp critical
            {
                cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<(count + 1)/((s3realpoints+1)*(s3imagpoints+1))*100.0<<"%"<<endl;
                cout<<"ds:"<<dsres<<endl;
                //fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<'\t'<<real(1.0/result)<<'\t'<<imag(1.0/result)<<endl;
                
                
                count = count + 1;
                cout<<count<<" out of "<<((s3realpoints+1)*(s3imagpoints+1))<<" runs finished"<<endl;
                cout<<"=============================="<<endl;
                cout<<endl;
                //if(count==1) break;
                //cout<<"run = "<<count<<endl;
                //cout<<"-------------------------"<<endl;
            }

            //if(printcount==0)
            //{
                
                //printcount = 1;
            //}
        }
    }

    //cudaDeviceReset();

    fout.open(filename.c_str());
    for(int k=0;k<s3imagpoints+1;++k)
    {
        for(int l=0;l<s3realpoints+1;++l)
        {
            //cout<<"k="<<k<<'\t'<<"l="<<l<<endl;
            //cout<<"dels3imag="<<dels3imag<<'\t'<<"dels3real="<<dels3real<<endl;
            double s3imag = s3imaginitial + k*dels3imag;
            double s3real = s3realinitial + l*dels3real;
            //cout<<"s3imag = "<<s3imag<<'\t'<<"s3real = "<<s3real<<endl;
            comp s = s3real + ii*s3imag;

            cout<<s3real<<'\t'<<s3imag<<'\t'<<real(resultvec[k][l])<<'\t'<<imag(resultvec[k][l])<<'\t'<<real(dsvec[k][l])<<'\t'<<imag(dsvec[k][l])<<endl;
            fout<<s3real<<'\t'<<s3imag<<'\t'<<real(resultvec[k][l])<<'\t'<<imag(resultvec[k][l])<<'\t'<<real(dsvec[k][l])<<'\t'<<imag(dsvec[k][l])<<endl;
        }
    }

    fout.close();

    for(int i=0; i<s3imagpoints+1; ++i)
    {
        delete [] resultvec[i];
    }

    delete [] resultvec;
}

void Mphib_belowthreshold_vs_s3_3d_omp_qvec_undersheet()
{
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    //double sinitial = 8.88;//8.72;//real(phibthreshold(a,m));//
    //double sfinal = real(phibthreshold(a,m));//8.82;
    //double dels = abs(sinitial-sfinal)/400.0;
    double points1 = 2000.0;
    double s3realpoints = 100.0;
    double s3imagpoints = 100.0;
    double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.01;
    //double box = 10.0;

    string filename="Mphib_momrep_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_3d"
                    //+ "_epspow_" + to_string((int)abs(log10(eps))) 
                    + "_opesecondsheet_smooth.dat";//_followq_poss.dat";
    ofstream fout;
    

    int count = 0;

    double eta = 20.0;

    double s3realinitial = 8.72;
    double s3realfinal = real(phibthreshold(a,m));
    double dels3real = abs(s3realinitial-s3realfinal)/s3realpoints;

    double s3imaginitial = -0.1;
    double s3imagfinal = 0.1;
    double dels3imag = abs(s3imaginitial - s3imagfinal)/s3imagpoints;
    //cout<<"dels3real = "<<dels3real<<'\t'<<"dels3imag = "<<dels3imag<<endl;

    //vector<vector<comp> > resultvec(points1,points1);
    comp **resultvec = new comp*[(int)s3imagpoints+1];
    for(int i=0; i<(int)s3imagpoints+1; ++i)
    {
        resultvec[i] = new comp[(int)s3realpoints+1];
    }

    comp **dsvec = new comp*[(int)s3imagpoints+1];
    for(int i=0; i<(int)s3imagpoints+1; ++i)
    {
        dsvec[i] = new comp[(int)s3realpoints+1];
    }

    //comp resultvec[(int)s3imagpoints+1][(int)s3realpoints+1];

    int i;
    //#pragma omp parallel for shared(resultvec)
    //omp_set_num_threads(10);
    //#pragma omp distribute parallel for shared(resultvec,count)

    
    //cudaStream_t nstreams[(int)s3imagpoints+1];
    
    //#pragma omp parallel for shared(resultvec)
    for(i=0;i<(int)s3imagpoints+1;++i)
    {
        //for(double s=sinitial;s<=sfinal;s=s+dels)
        //std::cout<<"from Thread = "<<omp_get_thread_num<<endl;
        //std::cout<<"total Thread = "<<omp_get_num_threads<<endl;
        int printcount = 0;
        
        for(int j=0;j<(int)s3realpoints+1;++j)
        {
            double s3imag = s3imaginitial + i*dels3imag;
            double s3real = s3realinitial + j*dels3real;

            comp sigk = 3.80*m*m;
            comp s = s3real + ii*s3imag;
            //double s = 8.78;
            //if(s>=sfinal) dels = 0.05/100.0;
            //eps = eps_above_threshold(eta,s,a,m,points);
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = sigmab(a,m);
            //comp qval = pmom(s,sigb,m);
            comp qval = pmom(s,sigb,m);
            
            //if(s3imag<=0.0) qval = -qval;
            cout<<"qval = "<<qval<<endl;

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        
            //cout<<"s:"<<s<<endl;
            //cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            //cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;

            comp delk = abs(kmin - kmax)/points1;
            int tag1,tag2;
            
            //mom_vector_maker_6(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,points1,qvec_r,tag1,tag2);
            //mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
            //mom_vector_maker_4(qvec,s,kmin,kmax,a,m,eps,eps_for_m2k,points1,qvec_r);
            if(s3imag<=0.0)
            //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            else if(s3imag>0.0)
            {
                //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
                cout<<"tag1 = "<<tag1<<endl;
                cout<<"tag2 = "<<tag2<<endl;
                
                
            }

            cout<<"qvec created with size = "<<qvec.size()<<endl;
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            
            if(s3imag<=0.0)
            {
                Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            }
            else if(s3imag>0.0)
            {
                Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            }
            //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            Gvec = -1.0*Gvec;

            //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
            Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
            //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);

            
            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        
            //cout<<"did this = "<<count<<" times"<<endl;
            //cout<<"i val = "<<i<<'\t';
            //cout<<"j val = "<<j<<endl;
            
            //LinearSolver_2(Bmat,dsol,Gvec,relerror);
            

            
            //cudaStream_t somethingstream;
            //cudaStreamCreate(&somethingstream);

            //#pragma omp critical
            //{
            
            cusolverComplex(Bmat,Gvec,dsol,size);
            
            //}
            //cout<<"stream no. for i = "<<i<<'\t'<<" is streamid = "<<&somethingstream<<endl;
            //cout<<"stream no. for i = "<<i<<'\t'<<" is streamid = "<<*nstreams[i]<<endl;

            //#pragma omp critical
            //{
                //cusolverComplexAsync(Bmat,Gvec,dsol,size,somethingstream);
            //}
            comp result;

            //cudaStreamDestinterpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);roy(somethingstream);

            //std::cout<<"cuda work finished from i ="<<i<<'\t'<<" j = "<<j<<endl;

            if(s3imag<=0.0)
            {
                interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
            }
            else if(s3imag>0.0)
            {
                interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
            }
            //interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
            comp m2kfunc = M2kfunc(a,sigk,m,eps_for_m2k);
            comp m2ksq = m2kfunc*m2kfunc;
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            dsvec[i][j] = result;
            comp dsres = result;
            result = gsq*result;
            //result = m2ksq*result;
            //result = rhophib(s,a,m)*result;
            
            resultvec[i][j] = result;

            //#pragma omp critical
            {
                cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<(count + 1)/((s3realpoints+1)*(s3imagpoints+1))*100.0<<"%"<<endl;
                cout<<"ds:"<<dsres<<endl;
                //fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<'\t'<<real(1.0/result)<<'\t'<<imag(1.0/result)<<endl;
                
                
                count = count + 1;
                cout<<count<<" out of "<<((s3realpoints+1)*(s3imagpoints+1))<<" runs finished"<<endl;
                cout<<"=============================="<<endl;
                cout<<endl;
                //if(count==1) break;
                //cout<<"run = "<<count<<endl;
                //cout<<"-------------------------"<<endl;
            }

            //if(printcount==0)
            //{
                
                //printcount = 1;
            //}
        }
    }

    //cudaDeviceReset();

    fout.open(filename.c_str());
    for(int k=0;k<s3imagpoints+1;++k)
    {
        for(int l=0;l<s3realpoints+1;++l)
        {
            //cout<<"k="<<k<<'\t'<<"l="<<l<<endl;
            //cout<<"dels3imag="<<dels3imag<<'\t'<<"dels3real="<<dels3real<<endl;
            double s3imag = s3imaginitial + k*dels3imag;
            double s3real = s3realinitial + l*dels3real;
            //cout<<"s3imag = "<<s3imag<<'\t'<<"s3real = "<<s3real<<endl;
            comp s = s3real + ii*s3imag;

            cout<<s3real<<'\t'<<s3imag<<'\t'<<real(resultvec[k][l])<<'\t'<<imag(resultvec[k][l])<<'\t'<<real(dsvec[k][l])<<'\t'<<imag(dsvec[k][l])<<endl;
            fout<<s3real<<'\t'<<s3imag<<'\t'<<real(resultvec[k][l])<<'\t'<<imag(resultvec[k][l])<<'\t'<<real(dsvec[k][l])<<'\t'<<imag(dsvec[k][l])<<endl;
        }
    }

    fout.close();

    for(int i=0; i<s3imagpoints+1; ++i)
    {
        delete [] resultvec[i];
    }

    delete [] resultvec;
}


void M_belowthreshold_vs_s3_3d_omp_sigkloop()
{
    comp ii = {0.0,1.0};
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    //double sinitial = 8.88;//8.72;//real(phibthreshold(a,m));//
    //double sfinal = real(phibthreshold(a,m));//8.82;
    //double dels = abs(sinitial-sfinal)/400.0;
    double points1 = 1000.0;
    double s3realpoints = 100.0;
    double s3imagpoints = 100.0;
    double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.01;
    //double box = 10.0;

    double sigkinitial = 3.975*m*m;
    double sigkfinal = 4.0*m*m;
    double sigkpoints = 5.0;
    double delsigk = abs(sigkinitial - sigkfinal)/sigkpoints;
    double sigkcount = 15;

    for(int kk=0;kk<sigkpoints;++kk)
    {
        double sigkval = sigkinitial + kk*delsigk;
    
    string filename="M_momrep_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_3d"
                    //+ "_epspow_" + to_string((int)abs(log10(eps))) 
                    + "_sigkcount_" + to_string((int)sigkcount) + "_sigk3.9t04.dat";//_followq_poss.dat";
    ofstream fout;
    

    int count = 0;

    double eta = 20.0;

    double s3realinitial = 8.72;
    double s3realfinal = real(phibthreshold(a,m));
    double dels3real = abs(s3realinitial-s3realfinal)/s3realpoints;

    double s3imaginitial = -0.1;
    double s3imagfinal = 0.1;
    double dels3imag = abs(s3imaginitial - s3imagfinal)/s3imagpoints;
    //cout<<"dels3real = "<<dels3real<<'\t'<<"dels3imag = "<<dels3imag<<endl;

    //vector<vector<comp> > resultvec(points1,points1);
    comp **resultvec = new comp*[(int)s3imagpoints+1];
    for(int i=0; i<(int)s3imagpoints+1; ++i)
    {
        resultvec[i] = new comp[(int)s3realpoints+1];
    }

    comp **dsvec = new comp*[(int)s3imagpoints+1];
    for(int i=0; i<(int)s3imagpoints+1; ++i)
    {
        dsvec[i] = new comp[(int)s3realpoints+1];
    }

    //comp resultvec[(int)s3imagpoints+1][(int)s3realpoints+1];

    int i;
    //#pragma omp parallel for shared(resultvec)
    //omp_set_num_threads(10);
    //#pragma omp distribute parallel for shared(resultvec,count)

    
    //cudaStream_t nstreams[(int)s3imagpoints+1];
    
    //#pragma omp parallel for shared(resultvec)
    for(i=0;i<(int)s3imagpoints+1;++i)
    {
        //for(double s=sinitial;s<=sfinal;s=s+dels)
        //std::cout<<"from Thread = "<<omp_get_thread_num<<endl;
        //std::cout<<"total Thread = "<<omp_get_num_threads<<endl;
        int printcount = 0;
        
        for(int j=0;j<(int)s3realpoints+1;++j)
        {
            double s3imag = s3imaginitial + i*dels3imag;
            double s3real = s3realinitial + j*dels3real;

            comp sigk = sigkval;//3.80*m*m;
            comp s = s3real + ii*s3imag;
            //double s = 8.78;
            //if(s>=sfinal) dels = 0.05/100.0;
            //eps = eps_above_threshold(eta,s,a,m,points);
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = sigmab(a,m);
            //comp qval = pmom(s,sigb,m);
            comp qval = pmom(s,sigk,m);
            
            //if(s3imag<=0.0) qval = -qval;
            cout<<"qval = "<<qval<<endl;

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        
            //cout<<"s:"<<s<<endl;
            //cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            //cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;

            comp delk = abs(kmin - kmax)/points1;
            int tag1,tag2;
            
            //mom_vector_maker_6(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,points1,qvec_r,tag1,tag2);
            mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
            //mom_vector_maker_4(qvec,s,kmin,kmax,a,m,eps,eps_for_m2k,points1,qvec_r);
    
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            
            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            Gvec = -1.0*Gvec;

            //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
            Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
            //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);


            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        
            //cout<<"did this = "<<count<<" times"<<endl;
            //cout<<"i val = "<<i<<'\t';
            //cout<<"j val = "<<j<<endl;
            
            //LinearSolver_2(Bmat,dsol,Gvec,relerror);
            

            
            //cudaStream_t somethingstream;
            //cudaStreamCreate(&somethingstream);

            //#pragma omp critical
            //{
                cusolverComplex(Bmat,Gvec,dsol,size);
            //}
            //cout<<"stream no. for i = "<<i<<'\t'<<" is streamid = "<<&somethingstream<<endl;
            //cout<<"stream no. for i = "<<i<<'\t'<<" is streamid = "<<*nstreams[i]<<endl;

            //#pragma omp critical
            //{
                //cusolverComplexAsync(Bmat,Gvec,dsol,size,somethingstream);
            //}
            comp result;

            //cudaStreamDestroy(somethingstream);

            //std::cout<<"cuda work finished from i ="<<i<<'\t'<<" j = "<<j<<endl;

            interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
            //interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
            comp m2kfunc = M2kfunc(a,sigk,m,eps_for_m2k);
            comp m2ksq = m2kfunc*m2kfunc;
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            dsvec[i][j] = result;
            comp dsres = result;
            //result = gsq*result;
            result = m2ksq*result;
            //result = rhophib(s,a,m)*result;
            
            resultvec[i][j] = result;

            //#pragma omp critical
            {
                cout<<"sigk:"<<sigk<<endl;
                cout<<"s:"<<s<<'\t'<<"dsres:"<<dsres<<'\n'<<"Dsres:"<<result<<'\t'<<"run:"<<(count + 1)/((s3realpoints+1)*(s3imagpoints+1))*100.0<<"%"<<endl;
                //fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<'\t'<<real(1.0/result)<<'\t'<<imag(1.0/result)<<endl;
                
                
                count = count + 1;
                cout<<count<<" out of "<<((s3realpoints+1)*(s3imagpoints+1))<<" runs finished"<<endl;

                //if(count==1) break;
                //cout<<"run = "<<count<<endl;
                //cout<<"-------------------------"<<endl;
            }

            //if(printcount==0)
            //{
                
                //printcount = 1;
            //}
        }
    }

    //cudaDeviceReset();

    fout.open(filename.c_str());
    for(int k=0;k<s3imagpoints+1;++k)
    {
        for(int l=0;l<s3realpoints+1;++l)
        {
            //cout<<"k="<<k<<'\t'<<"l="<<l<<endl;
            //cout<<"dels3imag="<<dels3imag<<'\t'<<"dels3real="<<dels3real<<endl;
            double s3imag = s3imaginitial + k*dels3imag;
            double s3real = s3realinitial + l*dels3real;
            //cout<<"s3imag = "<<s3imag<<'\t'<<"s3real = "<<s3real<<endl;
            comp s = s3real + ii*s3imag;

            cout<<s3real<<'\t'<<s3imag<<'\t'<<real(resultvec[k][l])<<'\t'<<imag(resultvec[k][l])<<'\t'<<sigkval<<endl;
            fout<<s3real<<'\t'<<s3imag<<'\t'<<real(resultvec[k][l])<<'\t'<<imag(resultvec[k][l])<<'\t'<<real(dsvec[k][l])<<'\t'<<imag(dsvec[k][l])<<'\t'<<sigkval<<endl;
        }
    }

    fout.close();

    for(int i=0; i<s3imagpoints+1; ++i)
    {
        delete [] resultvec[i];
    }

    delete [] resultvec;

    sigkcount = sigkcount + 1;

    }
}


void Mphib_belowthreshold_vs_eps()
{
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.70;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/100.0;
    double points1 = 5000.0;
    double points2 = 500.0;
    //double eps = 1.0e-5;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.01;
    //double box = 10.0;

    string filename="Mphib_momrep_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_contour4.dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    double eta = 20.0;

    double s = 8.72;

    double epsinitial = 0.1;
    double epsfinal = 1.0e-5;
    double epspoint = 200;
    double deleps = abs(epsinitial - epsfinal)/epspoint;

    //for(double s=sinitial;s<=sfinal;s=s+dels)
    for(double eps=epsfinal;eps<=epsinitial;eps=eps+deleps)
    { 
        //double eps = epsfinal + i*deleps;
        if(eps<1.0e-4) deleps = 1.0e-5;
        else if(eps>=1.0e-4 && eps<1.0e-3) deleps = 1.0e-4;
        else if(eps>=1.0e-3 && eps<1.0e-2) deleps = 1.0e-3;
        else if(eps>=1.0e-2 && eps<1.0e-1) deleps = 1.0e-2;
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp qval = pmom(s,sigb,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        
        //cout<<"s:"<<s<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;

        mom_vector_maker_4(qvec,s,kmin,kmax,a,m,eps,eps_for_m2k,points1/1.0,qvec_r);
    
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec = -1.0*Gvec;

        Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);

        double relerror;

        

        //LinearSolver_2(Bmat,dsol,Gvec,relerror);
        cusolverComplex(Bmat,Gvec,dsol,size);

        comp result;

        interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        result = gsq*result;
        //result = rhophib(s,a,m)*result;

        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        cout<<"eps:"<<eps<<'\t'<<"size:"<<size<<endl;
        fout<<setprecision(16)<<s<<'\t'<<eps<<'\t'<<size<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}

void Mphib_belowthreshold_vs_N()
{
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.70;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/100.0;
    double points1 = 3000.0;
    double points2 = 500.0;
    //double eps = 1.0e-5;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.01;
    //double box = 10.0;
    double eps = 1.0e-2;

    string filename="Mphib_momrep_a_" + to_string((int)a) 
                    + "_eps_" + to_string(eps) 
                    + "_contour4.dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    double eta = 20.0;

    double s = 8.72;

    double epsinitial = 0.01;
    double epsfinal = 1.0e-5;
    double epspoint = 50;
    double deleps = abs(epsinitial - epsfinal)/epspoint;

    double Ninitial = 500.0;
    double Nfinal = 5000.0;
    double delN = 500.0;

    //for(double s=sinitial;s<=sfinal;s=s+dels)
    //for(double eps=epsfinal;eps<=epsinitial;eps=eps+deleps)
    for(double N=Ninitial;N<=Nfinal;N=N+delN)
    { 
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        points1 = N;
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp qval = pmom(s,sigb,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        
        //cout<<"s:"<<s<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;

        mom_vector_maker_4(qvec,s,kmin,kmax,a,m,eps,eps_for_m2k,points1/1.0,qvec_r);
    
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec = -1.0*Gvec;

        Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);

        double relerror;

        

        LinearSolver_2(Bmat,dsol,Gvec,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

        comp result;

        interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        result = gsq*result;
        //result = rhophib(s,a,m)*result;

        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        cout<<"eps:"<<eps<<'\t'<<"size:"<<size<<endl;
        fout<<setprecision(16)<<s<<'\t'<<eps<<'\t'<<size<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}


void Mphib_belowthreshold_vs_s3_opposite()
{
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.77;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/100.0;
    double points1 = 500.0;
    double points2 = 500.0;
    double eps = 0.01;
    //double box = 10.0;

    string filename="Mphib_momrep_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_epspow_" + to_string((int)abs(log10(eps))) 
                    + ".dat";
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
        comp qval = pmom(s,sigb,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        
        //cout<<"s:"<<s<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;

        mom_vector_maker_opposite_1(qvec,s,kmin,kmax,a,m,eps,points1);
    
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec = -1.0*Gvec;

        Bmat_maker_momrep(Bmat,s,qvec,qvec,a,m,eps);

        double relerror;

        

        //LinearSolver(Bmat,dsol,Gvec,relerror);
        cusolverComplex(Bmat,Gvec,dsol,size);

        comp result;

        interpolator_ds_integraleq_momrep(dsol,qvec,s,qval,qval,a,m,eps,result);

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

void Mphib_belowthreshold_vs_s3_N_eps()
{
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double srealinitial = 8.72;//real(phibthreshold(a,m));//
    double srealfinal = real(phibthreshold(a,m));//8.82;
    double delsreal = abs(srealinitial-srealfinal)/100.0;
    double points1 = 500.0;
    double points2 = 500.0;
    //double eps = 0.0;
    //double box = 10.0;

   

    int count = 0;

    double eta = 20.0;

    double simaginitial = -0.1;
    double simagfinal = 0.1;
    double delsimag = abs(simaginitial - simagfinal)/20.0;

    double epsinitial = 0.02;
    double epsfinal = 0.06;
    double deleps = 0.01;

    double Ninitial = 250.0;
    double Nfinal = 1000.0;
    double delN = 150.0;

    for(int si=0;si<=20;++si)
    {
        double s3imag = simaginitial + si*delsimag;

        for(int sr=0;sr<=100;++sr)
        {
            double s3real = srealinitial + sr*delsreal;
            comp s = s3real + ii*s3imag;

            for(int eps1=0;eps1<=4;++eps1)
            {
                string filename=    "Mphib_momrep_a_" + to_string((int)a) 
                                    + "_sr_" + to_string(sr)
                                    + "_si_" + to_string(si)
                                    + "_eps_" + to_string(eps1) 
                                    + ".dat";
                ofstream fout;
                fout.open(filename.c_str());
                
                double eps = epsinitial + eps1*deleps;

                for(double N=Ninitial;N<=Nfinal;N=N+delN)
                {
                    points1 = N;
                    //double s = 8.78;
                    //if(s>=sfinal) dels = 0.05/100.0;
                    //eps = eps_above_threshold(eta,s,a,m,points);
                    comp sigmamin = {0.0,0.0};
                    comp sigmamax = sigmax(s,m);
                    comp sigb = sigmab(a,m);
                    comp qval = pmom(s,sigb,m);

                    comp kmin = 0.0;
                    comp kmax = pmom(s,0.0,m);
        
                    //cout<<"s:"<<s<<endl;
                    cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
                    cout<<"eps:"<<eps<<endl;

                    vector<comp> qvec;

                    mom_vector_maker_3(qvec,s,kmin,kmax,a,m,eps,0.0,points1,0.01);
    
                    int size = qvec.size();
                    Eigen::MatrixXcd Bmat(size,size);
                    Eigen::VectorXcd Gvec(size);
                    Eigen::VectorXcd dsol(size);

                    Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
                    Gvec = -1.0*Gvec;

                    Bmat_maker_momrep(Bmat,s,qvec,qvec,a,m,eps);

                    double relerror;

        

                    LinearSolver(Bmat,dsol,Gvec,relerror);
                    //cusolverComplex(Bmat,Gvec,dsol,size);

                    comp result;

                    interpolator_ds_integraleq_momrep(dsol,qvec,s,qval,qval,a,m,eps,result);

                    comp gfunc = gfuncConst(a,0.0,m);
                    comp gsq = gfunc*gfunc;
                    result = gsq*result;
                    //result = rhophib(s,a,m)*result;

                    cout<<"s:"<<s<<'\t'<<" N:"<<qvec.size()<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
                    fout<<real(s)<<'\t'
                        <<imag(s)<<'\t'
                        <<eps<<'\t'
                        <<size<<'\t'
                        <<real(result)<<'\t'
                        <<imag(result)<<'\t'
                        <<real(1.0/result)<<'\t'
                        <<imag(1.0/result)<<endl;
                    count = count + 1;
                    cout<<"-------------------------"<<endl;
                }

            fout.close();

            }
        }
    }

    

}

comp dSqqs2q2msq_func(  double s,
                        double a  )
{
    //double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.72;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/100.0;
    double points1 = 1000.0;
    double points2 = 500.0;
    double eps = 0.0;
    //double box = 10.0;

    //string filename="dSqqs2q2msq_momrep_a_" + to_string((int)a) 
    //                + "_N_" + to_string((int)points1) 
    //                + "_eps_" + to_string(eps) 
    //                + ".dat";
    //ofstream fout;
    //fout.open(filename.c_str());

    //int count = 0;

    //double eta = 20.0;


    //for(double s=sinitial;s<=sfinal;s=s+dels)
    //{ 
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp sigq = 2.0*m*m;
        comp qval = pmom(s,sigq,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        double delk = abs(kmax-kmin)/points1;
        
        //cout<<"s:"<<s<<endl;
        //cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        //cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;

        mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
        //mom_vector_maker_1(qvec,s,kmin,kmax,a,m,eps,points1);
    
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec = -1.0*Gvec;

        Bmat_maker_momrep(Bmat,s,qvec,qvec,a,m,eps);

        double relerror;

        

        LinearSolver(Bmat,dsol,Gvec,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

        comp result;

        interpolator_ds_integraleq_momrep(dsol,qvec,s,qval,qval,a,m,eps,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;

        //result = result + GS_pk(s,qval,qval,m,eps);
        //result = gsq*result;
        //result = rhophib(s,a,m)*result;

        //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        //fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        //count = count + 1;
        //cout<<"-------------------------"<<endl;
    //}

    //fout.close();

    return 1.0/result;

}

comp dSqqs2q2msq_func_with_N(   double s,
                                double a,
                                double N  )
{
    //double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.72;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/100.0;
    double points1 = N;
    double points2 = 500.0;
    double eps = 0.0;
    //double box = 10.0;

    //string filename="dSqqs2q2msq_momrep_a_" + to_string((int)a) 
    //                + "_N_" + to_string((int)points1) 
    //                + "_eps_" + to_string(eps) 
    //                + ".dat";
    //ofstream fout;
    //fout.open(filename.c_str());

    //int count = 0;

    //double eta = 20.0;


    //for(double s=sinitial;s<=sfinal;s=s+dels)
    //{ 
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp sigq = 2.0*m*m;
        comp qval = pmom(s,sigq,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        double delk = abs(kmax-kmin)/points1;
        
        //cout<<"s:"<<s<<endl;
        //cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        //cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;

        int firstpart_for_nonunimesh = 1;
        int secondpart_for_nonunimesh = 4;
        //mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
        //mom_vector_maker_linear_2(qvec,kmin,kmax,points1,1,firstpart_for_nonunimesh,secondpart_for_nonunimesh);
        //mom_vector_maker_1(qvec,s,kmin,kmax,a,m,eps,points1);
        line_maker(qvec,kmin,kmax,points1);
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        //Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec = -1.0*Gvec;

        Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,0.0);

        double relerror;

        

        //LinearSolver(Bmat,dsol,Gvec,relerror);
        cusolverComplex(Bmat,Gvec,dsol,size);

        comp result;

        interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,0.0,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;

        //result = result + GS_pk(s,qval,qval,m,eps);
        //result = gsq*result;
        //result = rhophib(s,a,m)*result;

        //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        //fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        //count = count + 1;
        //cout<<"-------------------------"<<endl;
    //}

    //fout.close();

    return 1.0/result;

}

comp dSqqs2q2msq_func_with_N_using_determinant(     double s,
                                                    double a,
                                                    double N  )
{
    //double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.72;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/100.0;
    double points1 = N;
    double points2 = 500.0;
    double eps = 0.0;
    //double box = 10.0;

    //string filename="dSqqs2q2msq_momrep_a_" + to_string((int)a) 
    //                + "_N_" + to_string((int)points1) 
    //                + "_eps_" + to_string(eps) 
    //                + ".dat";
    //ofstream fout;
    //fout.open(filename.c_str());

    //int count = 0;

    //double eta = 20.0;


    //for(double s=sinitial;s<=sfinal;s=s+dels)
    //{ 
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp sigq = 2.0*m*m;
        comp qval = pmom(s,sigq,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        double delk = abs(kmax-kmin)/points1;
        
        //cout<<"s:"<<s<<endl;
        //cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        //cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;

        int firstpart_for_nonunimesh = 1;
        int secondpart_for_nonunimesh = 4;
        //mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
        //mom_vector_maker_linear_2(qvec,kmin,kmax,points1,1,firstpart_for_nonunimesh,secondpart_for_nonunimesh);
        //mom_vector_maker_1(qvec,s,kmin,kmax,a,m,eps,points1);
        line_maker(qvec,kmin,kmax,points1);
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        //Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        //Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        //Gvec = -1.0*Gvec;

        Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,0.0);

        double relerror;

        

        //LinearSolver(Bmat,dsol,Gvec,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

        comp result = Bmat.determinant();

        //interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,0.0,result);

        //comp gfunc = gfuncConst(a,0.0,m);
        //comp gsq = gfunc*gfunc;

        //result = result + GS_pk(s,qval,qval,m,eps);
        //result = gsq*result;
        //result = rhophib(s,a,m)*result;

        //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        //fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        //count = count + 1;
        //cout<<"-------------------------"<<endl;
    //}

    //fout.close();

    return result;

}

comp dSqqs2q2msq_func_with_N_using_determinant_with_LUcuda(     double s,
                                                                double a,
                                                                double N  )
{
    //double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.72;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/100.0;
    double points1 = N;
    double points2 = 500.0;
    double eps = 0.0;
    //double box = 10.0;

    //string filename="dSqqs2q2msq_momrep_a_" + to_string((int)a) 
    //                + "_N_" + to_string((int)points1) 
    //                + "_eps_" + to_string(eps) 
    //                + ".dat";
    //ofstream fout;
    //fout.open(filename.c_str());

    //int count = 0;

    //double eta = 20.0;


    //for(double s=sinitial;s<=sfinal;s=s+dels)
    //{ 
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp sigq = 2.0*m*m;
        comp qval = pmom(s,sigq,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        double delk = abs(kmax-kmin)/points1;
        
        //cout<<"s:"<<s<<endl;
        //cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        //cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;

        int firstpart_for_nonunimesh = 1;
        int secondpart_for_nonunimesh = 4;
        //mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
        //mom_vector_maker_linear_2(qvec,kmin,kmax,points1,1,firstpart_for_nonunimesh,secondpart_for_nonunimesh);
        //mom_vector_maker_1(qvec,s,kmin,kmax,a,m,eps,points1);
        line_maker(qvec,kmin,kmax,points1);
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        //Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec = -1.0*Gvec;

        Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,0.0);

        double relerror;

        comp result;

        struct timeval stop, start;
        gettimeofday(&start, NULL);
        LUSolver_complex_cuda(Bmat,Gvec,dsol,size,result);
        //result = Bmat.determinant();
        gettimeofday(&stop, NULL);
        double actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/1000000.0;
        cout<<"time taken to solve = "<<actual_sol_time<<endl;
    
        //LinearSolver(Bmat,dsol,Gvec,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

        //comp result = Bmat.determinant();

        //interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,0.0,result);

        //comp gfunc = gfuncConst(a,0.0,m);
        //comp gsq = gfunc*gfunc;

        //result = result + GS_pk(s,qval,qval,m,eps);
        //result = gsq*result;
        //result = rhophib(s,a,m)*result;

        //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        //fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        //count = count + 1;
        //cout<<"-------------------------"<<endl;
    //}

    //fout.close();

    return result;

}

comp Mphib_secondsheet_denom_func_with_N_using_weights(     double s,
                                                            double a,
                                                            double N  )
{
    //double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.72;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/100.0;
    double points1 = N;
    double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0;
    //double box = 10.0;

    //string filename="dSqqs2q2msq_momrep_a_" + to_string((int)a) 
    //                + "_N_" + to_string((int)points1) 
    //                + "_eps_" + to_string(eps) 
    //                + ".dat";
    //ofstream fout;
    //fout.open(filename.c_str());

    //int count = 0;

    //double eta = 20.0;


    //for(double s=sinitial;s<=sfinal;s=s+dels)
    //{ 
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp sigq = 2.0*m*m;
        comp qval = pmom(s,sigb,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        double delk = abs(kmax-kmin)/points1;
        
        //cout<<"s:"<<s<<endl;
        //cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        //cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;
        vector<comp> weights;

        int firstpart_for_nonunimesh = 1;
        int secondpart_for_nonunimesh = 4;
        //mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
        //mom_vector_maker_linear_2(qvec,kmin,kmax,points1,1,firstpart_for_nonunimesh,secondpart_for_nonunimesh);
        //mom_vector_maker_1(qvec,s,kmin,kmax,a,m,eps,points1);
        //mom_vector_maker_43_with_weights(qvec,weights,s,kmin,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        line_maker_with_weights(qvec,weights,kmin,kmax,points1);
        //line_maker(qvec,kmin,kmax,points1);
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        //Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec = -1.0*Gvec;

        //Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,0.0);
        Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);
        double relerror;

        comp result;

        struct timeval stop, start;
        gettimeofday(&start, NULL);
        //LUSolver_complex_cuda(Bmat,Gvec,dsol,size,result);
        //result = Bmat.determinant();
        LinearSolver(Bmat,dsol,Gvec,relerror);
        gettimeofday(&stop, NULL);
        double actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/1000000.0;
        cout<<"time taken to solve = "<<actual_sol_time<<endl;
    
        //LinearSolver(Bmat,dsol,Gvec,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

        //comp result = Bmat.determinant();

        //interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,0.0,result);
        interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        comp ds = result;
        result = gsq*result;

        
        comp rhopb = rhophib(s,a,m);
        comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
        comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
        comp GSval = GS_pk(s,qval,qval,m,eps);

        //result = result + GS_pk(s,qval,qval,m,eps);
        //result = gsq*result;
        //result = rhophib(s,a,m)*result;

        //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        //fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        //count = count + 1;
        //cout<<"-------------------------"<<endl;
    //}

    //fout.close();

    return mphib2denom;

}

comp dSqqs2q2msq_func_with_N_using_weights(     double s,
                                                            double a,
                                                            double N  )
{
    //double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.72;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/100.0;
    double points1 = N;
    double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0;
    //double box = 10.0;

    //string filename="dSqqs2q2msq_momrep_a_" + to_string((int)a) 
    //                + "_N_" + to_string((int)points1) 
    //                + "_eps_" + to_string(eps) 
    //                + ".dat";
    //ofstream fout;
    //fout.open(filename.c_str());

    //int count = 0;

    //double eta = 20.0;


    //for(double s=sinitial;s<=sfinal;s=s+dels)
    //{ 
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp sigq = 2.0*m*m;
        comp qval = pmom(s,sigq,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        double delk = abs(kmax-kmin)/points1;
        
        //cout<<"s:"<<s<<endl;
        //cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        //cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;
        vector<comp> weights;

        int firstpart_for_nonunimesh = 1;
        int secondpart_for_nonunimesh = 4;
        //mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
        //mom_vector_maker_linear_2(qvec,kmin,kmax,points1,1,firstpart_for_nonunimesh,secondpart_for_nonunimesh);
        //mom_vector_maker_1(qvec,s,kmin,kmax,a,m,eps,points1);
        //mom_vector_maker_43_with_weights(qvec,weights,s,kmin,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        //line_maker(qvec,kmin,kmax,points1);
        line_maker_with_weights(qvec,weights,kmin,kmax,points1);
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        //Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec = -1.0*Gvec;

        //Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,0.0);
        Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);
        double relerror;

        comp result;

        struct timeval stop, start;
        gettimeofday(&start, NULL);
        //LUSolver_complex_cuda(Bmat,Gvec,dsol,size,result);
        result = Bmat.determinant();
        //LinearSolver(Bmat,dsol,Gvec,relerror);
        gettimeofday(&stop, NULL);
        double actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/1000000.0;
        cout<<"time taken to solve = "<<actual_sol_time<<endl;
    
        //LinearSolver(Bmat,dsol,Gvec,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

        //comp result = Bmat.determinant();

        //interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,0.0,result);
        //interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        comp ds = result;
        //result = gsq*result;

        
        comp rhopb = rhophib(s,a,m);
        comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
        comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
        comp GSval = GS_pk(s,qval,qval,m,eps);

        //result = result + GS_pk(s,qval,qval,m,eps);
        //result = gsq*result;
        //result = rhophib(s,a,m)*result;

        //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        //fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        //count = count + 1;
        //cout<<"-------------------------"<<endl;
    //}

    //fout.close();

    return result;//mphib2denom;

}

comp KDF_func_with_N_using_weights(     comp KDF,
                                        double s,
                                        double a,
                                        double N  )
{
    //double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.72;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/100.0;
    double points1 = N;
    double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0;
    //double box = 10.0;

    //string filename="dSqqs2q2msq_momrep_a_" + to_string((int)a) 
    //                + "_N_" + to_string((int)points1) 
    //                + "_eps_" + to_string(eps) 
    //                + ".dat";
    //ofstream fout;
    //fout.open(filename.c_str());

    //int count = 0;

    //double eta = 20.0;


    //for(double s=sinitial;s<=sfinal;s=s+dels)
    //{ 
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp sigq = 2.0*m*m;
        comp qval = pmom(s,sigq,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        double delk = abs(kmax-kmin)/points1;
        
        //cout<<"s:"<<s<<endl;
        //cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        //cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;
        vector<comp> weights;

        int firstpart_for_nonunimesh = 1;
        int secondpart_for_nonunimesh = 4;
        //mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
        //mom_vector_maker_linear_2(qvec,kmin,kmax,points1,1,firstpart_for_nonunimesh,secondpart_for_nonunimesh);
        //mom_vector_maker_1(qvec,s,kmin,kmax,a,m,eps,points1);
        //mom_vector_maker_43_with_weights(qvec,weights,s,kmin,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        //line_maker(qvec,kmin,kmax,points1);
        line_maker_with_weights(qvec,weights,kmin,kmax,points1);
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::MatrixXcd Gmat(size,size);
        Eigen::MatrixXcd dsol(size,size);

        //Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gmat_maker_momrep(Gmat,s,qvec,qvec,m,eps);
        Gmat = -1.0*Gmat;

        //Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,0.0);
        Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);
        double relerror;

        comp result;

        struct timeval stop, start;
        gettimeofday(&start, NULL);
        //LUSolver_complex_cuda(Bmat,Gvec,dsol,size,result);
        //result = Bmat.determinant();
        LinearSolver_3(Bmat,dsol,Gmat,relerror);
        gettimeofday(&stop, NULL);
        double actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/1000000.0;
        cout<<"time taken to solve = "<<actual_sol_time<<endl;
    
        //LinearSolver(Bmat,dsol,Gvec,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

        //comp result = Bmat.determinant();

        //interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,0.0,result);
        //interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        comp ds = result;
        //result = gsq*result;

        
        comp rhopb = rhophib(s,a,m);
        comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
        comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
        comp GSval = GS_pk(s,qval,qval,m,eps);

        F3infvol(dsol,qvec,weights,qvec,weights,s,a,m,eps,result);
        //result = result + GS_pk(s,qval,qval,m,eps);
        //result = gsq*result;
        //result = rhophib(s,a,m)*result;

        //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        //fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        //count = count + 1;
        //cout<<"-------------------------"<<endl;
    //}

    //fout.close();

    return result + 1.0/KDF;//mphib2denom;

}



void dSqqs2k2msq_pole_searching(    double init_guess,
                                    double a,
                                    double &pole    )
{
    double initialguess_s = init_guess;
    int max_iteration = 100;
    double dels = 0.001;
    //double a = 16.0;
    double s = initialguess_s;
    
    double epsilon = 1.0e-10;
    double tolerance = 1.0e-10;

    for(int i=0;i<max_iteration;++i)
    {
        double f = (double) real(dSqqs2q2msq_func(s,a));
        double fprime = (double) real(dSqqs2q2msq_func(s+dels,a) - dSqqs2q2msq_func(s,a))/dels;

        cout<<"a = "<<a<<endl;
        cout<<"iteration:"<<i+1<<'\t'<<" s:"<<s<<'\t'<<" f:"<<f<<'\t'<<" fprime:"<<fprime<<endl;

        if(abs(fprime)<epsilon)
        {
            cout<<"fprime smaller than epsilon = "<<epsilon<<endl; 
            break;
        }
        double temp_s = s;
        s = s - f/fprime;

        if(abs(s-temp_s)<=tolerance)
        {
            cout<<"s-temp_s is smaller than tolerance"<<endl;
            break;
        }
    }
    cout<<"pole found at s = "<<s<<endl;
    pole = s;
}


void dSqqs2k2msq_pole_searching_bissection(     double init_guess1,
                                                double init_guess2,
                                                double a,
                                                double &pole    )
{
    double apoint = init_guess1;
    double bpoint = init_guess2;
    int max_iteration = 100;
    double dels = 0.001;
    //double a = 16.0;
    //double s = apoint;
    
    double epsilon = 1.0e-10;
    double tolerance = 1.0e-15;
    double sa = apoint;
    double sb = bpoint;

    for(int i=0;i<max_iteration;++i)
    {
        
        double sc = (sa + sb)/2.0;
        double fc = (double) real(dSqqs2q2msq_func(sc,a));
        double fa = (double) real(dSqqs2q2msq_func(sa,a));
        //double fprime = (double) real(dSqqs2q2msq_func(s+dels,a) - dSqqs2q2msq_func(s,a))/dels;

        cout<<"a = "<<a<<endl;
        cout<<"iteration:"<<i+1<<'\t'<<" s:"<<sc<<'\t'<<" f:"<<fc<<endl;

        if(abs(fc)<epsilon || (sb-sa)/2.0 < tolerance )
        {
            cout<<"---------------------------"<<endl;
            cout<<"fc:"<<fc<<'\t'<<"eps:"<<epsilon<<endl;
            cout<<"b-a/2 : "<<(sb-sa)/2.0<<endl;
            pole = sc;
            cout<<"pole found at s = "<<setprecision(16)<<sc<<endl;
            cout<<"---------------------------"<<endl;
            cout<<endl;
            break;
        }
        if(fc>0.0 && fa>0.0) 
        {
            sa = sc;
        }
        else 
        {
            sb = sc;
        }

        if(i==max_iteration-1)
        {
            cout<<"method failed, max iteration reached"<<endl;
        }
    }
    //cout<<"pole found at s = "<<s<<endl;
    //pole = s;
}

void dSqqs2k2msq_pole_searching_bissection_with_N(      double init_guess1,
                                                        double init_guess2,
                                                        double a,
                                                        double &pole,
                                                        double N    )
{
    double apoint = init_guess1;
    double bpoint = init_guess2;
    int max_iteration = 100;
    double dels = 0.001;
    //double a = 16.0;
    //double s = apoint;
    
    double epsilon = 1.0e-15;
    double tolerance = 1.0e-15;
    double sa = apoint;
    double sb = bpoint;

    for(int i=0;i<max_iteration;++i)
    {
        
        double sc = (sa + sb)/2.0;
        double fc = (double) real(dSqqs2q2msq_func_with_N(sc,a,N));
        double fa = (double) real(dSqqs2q2msq_func_with_N(sa,a,N));
        //double fprime = (double) real(dSqqs2q2msq_func(s+dels,a) - dSqqs2q2msq_func(s,a))/dels;

        cout<<"a = "<<a<<'\t'<<" N = "<<N<<endl;
        cout<<"init_guess1 = "<<setprecision(16)<<init_guess1<<'\n'
            <<"init_guess2 = "<<setprecision(16)<<init_guess2<<endl;
        cout<<" s_a = "<<setprecision(16)<<sa<<'\n'
            <<" s_b = "<<setprecision(16)<<sb<<'\n'
            <<" s_c = "<<setprecision(16)<<sc<<endl;
        cout<<"iteration:"<<i+1<<'\t'<<" s:"<<sc<<'\t'<<" f:"<<fc<<endl;

        if(abs(fc)<epsilon || abs(sb-sa)/2.0 < tolerance )
        {
            cout<<"---------------------------"<<endl;
            cout<<"fc:"<<fc<<'\t'<<"eps:"<<epsilon<<endl;
            cout<<"b-a/2 : "<<(sb-sa)/2.0<<endl;
            pole = sc;
            cout<<"pole found at s = "<<setprecision(16)<<sc<<endl;
            cout<<"---------------------------"<<endl;
            cout<<endl;
            break;
        }
        if(fc>0.0 && fa>0.0) 
        {
            sa = sc;
        }
        else 
        {
            sb = sc;
        }

        if(i==max_iteration-1)
        {
            cout<<"method failed, max iteration reached"<<endl;
        }
    }
    //cout<<"pole found at s = "<<s<<endl;
    //pole = s;
}

void dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(    double init_guess1,
                                                                        double init_guess2,
                                                                        double a,
                                                                        double &pole,
                                                                        double N    )
{
    double apoint = init_guess1;
    double bpoint = init_guess2;
    int max_iteration = 100;
    double dels = 0.001;
    //double a = 16.0;
    //double s = apoint;
    
    double epsilon = 1.0e-15;
    double tolerance = 1.0e-15;
    double sa = apoint;
    double sb = bpoint;

    for(int i=0;i<max_iteration;++i)
    {
        
        double sc = (sa + sb)/2.0;
        //double fc = (double) real(dSqqs2q2msq_func_with_N_using_determinant(sc,a,N));
        //double fa = (double) real(dSqqs2q2msq_func_with_N_using_determinant(sa,a,N));
        double fc = (double) real(dSqqs2q2msq_func_with_N_using_determinant_with_LUcuda(sc,a,N));
        double fa = (double) real(dSqqs2q2msq_func_with_N_using_determinant_with_LUcuda(sa,a,N));
        
        //double fprime = (double) real(dSqqs2q2msq_func(s+dels,a) - dSqqs2q2msq_func(s,a))/dels;
        int sign_b = 0;
        double fb = 0;
        int sign_c = sign_checker(fc);
        int sign_a = sign_checker(fa);
        if(i==0)
        {
            //fb = (double) real(dSqqs2q2msq_func_with_N_using_determinant(sb,a,N));
            fb = (double) real(dSqqs2q2msq_func_with_N_using_determinant_with_LUcuda(sb,a,N));
            sign_b = sign_checker(fb);
            cout<<"sign a = "<<sign_a<<endl;
            cout<<"sign_b = "<<sign_b<<endl;
            cout<<"sign_c = "<<sign_c<<endl;
          
            if(sign_a==sign_b)
            {
                cout<<"two opposite limits are of the same sign"<<endl;
                cout<<setprecision(16);
                cout<<"s_a = "<<sa<<'\t'<<" s_b = "<<sb<<endl;
                cout<<"f_a = "<<fa<<'\t'<<" f_b = "<<fb<<endl;
                cout<<"please choose different limits"<<endl;
                exit(128);
            }
        }
        cout<<"a = "<<a<<'\t'<<" N = "<<N<<endl;
        cout<<"init_guess1 = "<<setprecision(16)<<init_guess1<<'\n'
            <<"init_guess2 = "<<setprecision(16)<<init_guess2<<endl;
        cout<<" s_a = "<<setprecision(16)<<sa<<'\n'
            <<" s_b = "<<setprecision(16)<<sb<<'\n'
            <<" s_c = "<<setprecision(16)<<sc<<endl;
        cout<<"iteration:"<<i+1<<'\t'<<" s:"<<sc<<'\t'<<" f:"<<fc<<endl;

        if(abs(fc)<epsilon || abs(sb-sa)/2.0 < tolerance )
        {
            cout<<"---------------------------"<<endl;
            cout<<"fc:"<<fc<<'\t'<<"eps:"<<epsilon<<endl;
            cout<<"b-a/2 : "<<(sb-sa)/2.0<<endl;
            pole = sc;
            cout<<"pole found at s = "<<setprecision(16)<<sc<<endl;
            cout<<"---------------------------"<<endl;
            cout<<endl;
            break;
        }
        //if(fc>0.0 && fa>0.0) 
        if(sign_c==sign_a)
        {
            sa = sc;
        }
        else 
        {
            sb = sc;
        }

        if(i==max_iteration-1)
        {
            cout<<"method failed, max iteration reached"<<endl;
        }
    }
    //cout<<"pole found at s = "<<s<<endl;
    //pole = s;
}

void Mphib_secondsheet_denom_pole_searching_bissection_with_N_using_weights(    double init_guess1,
                                                                        double init_guess2,
                                                                        double a,
                                                                        double &pole,
                                                                        double N    )
{
    double apoint = init_guess1;
    double bpoint = init_guess2;
    int max_iteration = 100;
    double dels = 0.001;
    //double a = 16.0;
    //double s = apoint;
    
    double epsilon = 1.0e-15;
    double tolerance = 1.0e-15;
    double sa = apoint;
    double sb = bpoint;

    for(int i=0;i<max_iteration;++i)
    {
        
        double sc = (sa + sb)/2.0;
        //double fc = (double) real(dSqqs2q2msq_func_with_N_using_determinant(sc,a,N));
        //double fa = (double) real(dSqqs2q2msq_func_with_N_using_determinant(sa,a,N));
        double fc = (double) real(Mphib_secondsheet_denom_func_with_N_using_weights(sc,a,N));
        double fa = (double) real(Mphib_secondsheet_denom_func_with_N_using_weights(sa,a,N));
        
        //double fprime = (double) real(dSqqs2q2msq_func(s+dels,a) - dSqqs2q2msq_func(s,a))/dels;
        int sign_b = 0;
        double fb = 0;
        int sign_c = sign_checker(fc);
        int sign_a = sign_checker(fa);
        if(i==0)
        {
            //fb = (double) real(dSqqs2q2msq_func_with_N_using_determinant(sb,a,N));
            fb = (double) real(Mphib_secondsheet_denom_func_with_N_using_weights(sb,a,N));
            sign_b = sign_checker(fb);
            cout<<"sign a = "<<sign_a<<endl;
            cout<<"sign_b = "<<sign_b<<endl;
            cout<<"sign_c = "<<sign_c<<endl;
          
            //if(sign_a==sign_b)
            if(sign_b==-1 && sign_a==+1)
            {
                cout<<"two opposite limits are of the opposite sign"<<endl;
                cout<<setprecision(16);
                cout<<"s_a = "<<sa<<'\t'<<" s_b = "<<sb<<endl;
                cout<<"f_a = "<<fa<<'\t'<<" f_b = "<<fb<<endl;
                cout<<"pole cant be found"<<endl;
                pole = 0;
                break;
                //exit(128);
            }
            else if(sign_a==-1 && sign_b==-1)
            {
                cout<<"two opposite limits are negative signed"<<endl;
                cout<<setprecision(16);
                cout<<"s_a = "<<sa<<'\t'<<" s_b = "<<sb<<endl;
                cout<<"f_a = "<<fa<<'\t'<<" f_b = "<<fb<<endl;
                cout<<"pole cant be found"<<endl;
                pole = 0;
                break;
            }
        }
        cout<<"a = "<<a<<'\t'<<" N = "<<N<<endl;
        cout<<"init_guess1 = "<<setprecision(16)<<init_guess1<<'\n'
            <<"init_guess2 = "<<setprecision(16)<<init_guess2<<endl;
        cout<<" s_a = "<<setprecision(16)<<sa<<'\n'
            <<" s_b = "<<setprecision(16)<<sb<<'\n'
            <<" s_c = "<<setprecision(16)<<sc<<endl;
        cout<<"iteration:"<<i+1<<'\t'<<" s:"<<sc<<'\t'<<" f:"<<fc<<endl;

        if(abs(fc)<epsilon || abs(sb-sa)/2.0 < tolerance )
        {
            cout<<"---------------------------"<<endl;
            cout<<"fc:"<<fc<<'\t'<<"eps:"<<epsilon<<endl;
            cout<<"b-a/2 : "<<(sb-sa)/2.0<<endl;
            pole = sc;
            cout<<"pole found at s = "<<setprecision(16)<<sc<<endl;
            cout<<"---------------------------"<<endl;
            cout<<endl;
            break;
        }
        //if(fc>0.0 && fa>0.0) 
        if(sign_c==sign_a)
        {
            sa = sc;
        }
        else 
        {
            sb = sc;
        }

        if(i==max_iteration-1)
        {
            cout<<"method failed, max iteration reached"<<endl;
        }
    }
    //cout<<"pole found at s = "<<s<<endl;
    //pole = s;
}

void dSqqs2k2msq_pole_searching_bissection_with_N_using_weights(    double init_guess1,
                                                                        double init_guess2,
                                                                        double a,
                                                                        double &pole,
                                                                        double N    )
{
    double apoint = init_guess1;
    double bpoint = init_guess2;
    int max_iteration = 100;
    double dels = 0.001;
    //double a = 16.0;
    //double s = apoint;
    
    double epsilon = 1.0e-15;
    double tolerance = 1.0e-15;
    double sa = apoint;
    double sb = bpoint;

    for(int i=0;i<max_iteration;++i)
    {
        
        double sc = (sa + sb)/2.0;
        //double fc = (double) real(dSqqs2q2msq_func_with_N_using_determinant(sc,a,N));
        //double fa = (double) real(dSqqs2q2msq_func_with_N_using_determinant(sa,a,N));
        double fc = (double) real(dSqqs2q2msq_func_with_N_using_weights(sc,a,N));
        double fa = (double) real(dSqqs2q2msq_func_with_N_using_weights(sa,a,N));
        
        //double fprime = (double) real(dSqqs2q2msq_func(s+dels,a) - dSqqs2q2msq_func(s,a))/dels;
        int sign_b = 0;
        double fb = 0;
        int sign_c = sign_checker(fc);
        int sign_a = sign_checker(fa);
        if(i==0)
        {
            //fb = (double) real(dSqqs2q2msq_func_with_N_using_determinant(sb,a,N));
            fb = (double) real(dSqqs2q2msq_func_with_N_using_weights(sb,a,N));
            sign_b = sign_checker(fb);
            cout<<"sign a = "<<sign_a<<endl;
            cout<<"sign_b = "<<sign_b<<endl;
            cout<<"sign_c = "<<sign_c<<endl;
          
            if(sign_a==sign_b)
            {
                cout<<"two opposite limits are of the same sign"<<endl;
                cout<<setprecision(16);
                cout<<"s_a = "<<sa<<'\t'<<" s_b = "<<sb<<endl;
                cout<<"f_a = "<<fa<<'\t'<<" f_b = "<<fb<<endl;
                cout<<"please choose different limits"<<endl;
                exit(128);
            }
        }
        cout<<"a = "<<a<<'\t'<<" N = "<<N<<endl;
        cout<<"init_guess1 = "<<setprecision(16)<<init_guess1<<'\n'
            <<"init_guess2 = "<<setprecision(16)<<init_guess2<<endl;
        cout<<" s_a = "<<setprecision(16)<<sa<<'\n'
            <<" s_b = "<<setprecision(16)<<sb<<'\n'
            <<" s_c = "<<setprecision(16)<<sc<<endl;
        cout<<"iteration:"<<i+1<<'\t'<<" s:"<<sc<<'\t'<<" f:"<<fc<<endl;

        if(abs(fc)<epsilon || abs(sb-sa)/2.0 < tolerance )
        {
            cout<<"---------------------------"<<endl;
            cout<<"fc:"<<fc<<'\t'<<"eps:"<<epsilon<<endl;
            cout<<"b-a/2 : "<<(sb-sa)/2.0<<endl;
            pole = sc;
            cout<<"pole found at s = "<<setprecision(16)<<sc<<endl;
            cout<<"---------------------------"<<endl;
            cout<<endl;
            break;
        }
        //if(fc>0.0 && fa>0.0) 
        if(sign_c==sign_a)
        {
            sa = sc;
        }
        else 
        {
            sb = sc;
        }

        if(i==max_iteration-1)
        {
            cout<<"method failed, max iteration reached"<<endl;
        }
    }
    //cout<<"pole found at s = "<<s<<endl;
    //pole = s;
}

void KDF_F3_pole_searching_bissection_with_N_using_weights(     comp KDF,        
                                                                double init_guess1,
                                                                double init_guess2,
                                                                double a,
                                                                double &pole,
                                                                double N    )
{
    double apoint = init_guess1;
    double bpoint = init_guess2;
    int max_iteration = 100;
    double dels = 0.001;
    //double a = 16.0;
    //double s = apoint;
    
    double epsilon = 1.0e-15;
    double tolerance = 1.0e-15;
    double sa = apoint;
    double sb = bpoint;

    for(int i=0;i<max_iteration;++i)
    {
        
        double sc = (sa + sb)/2.0;
        //double fc = (double) real(dSqqs2q2msq_func_with_N_using_determinant(sc,a,N));
        //double fa = (double) real(dSqqs2q2msq_func_with_N_using_determinant(sa,a,N));
        double fc = (double) real(KDF_func_with_N_using_weights(KDF,sc,a,N));
        double fa = (double) real(KDF_func_with_N_using_weights(KDF,sa,a,N));
        
        //double fprime = (double) real(dSqqs2q2msq_func(s+dels,a) - dSqqs2q2msq_func(s,a))/dels;
        int sign_b = 0;
        double fb = 0.0;//(double) real(KDF_func_with_N_using_weights(KDF,sb,a,N));
        //sign_b = sign_checker(fb);
        int sign_c = sign_checker(fc);
        int sign_a = sign_checker(fa);
        cout<<"sign a = "<<sign_a<<endl;
        if(i==0)
        cout<<"sign_b = "<<sign_b<<endl;
        cout<<"sign_c = "<<sign_c<<endl;
        if(i==0)
        {
            //fb = (double) real(dSqqs2q2msq_func_with_N_using_determinant(sb,a,N));
            fb = (double) real(KDF_func_with_N_using_weights(KDF,sb,a,N));
            sign_b = sign_checker(fb);
            cout<<"sign a = "<<sign_a<<endl;
            cout<<"sign_b = "<<sign_b<<endl;
            cout<<"sign_c = "<<sign_c<<endl;
          
            if(sign_a==sign_b)
            {
                cout<<"two opposite limits are of the same sign"<<endl;
                cout<<setprecision(16);
                cout<<"s_a = "<<sa<<'\t'<<" s_b = "<<sb<<endl;
                cout<<"f_a = "<<fa<<'\t'<<" f_b = "<<fb<<endl;
                /*cout<<"incrementing sb"<<endl;
                for(int i=0;;++i)
                {
                    double dels = 0.000005;
                    sb = sb + dels;
                    cout<<"incrementing sb from = "<<init_guess2<<" to "<<sb<<endl;
                    fb = (double) real(KDF_func_with_N_using_weights(KDF,sb,a,N));
                    sign_b = sign_checker(fb);
                    if(sign_a==sign_b) continue;
                    else break;
                }*/
                cout<<"pole cant be found"<<endl;
                pole = 0;
                break;
                //exit(128);
            }
        }
        cout<<"a = "<<a<<'\t'<<" N = "<<N<<endl;
        cout<<"init_guess1 = "<<setprecision(16)<<init_guess1<<'\n'
            <<"init_guess2 = "<<setprecision(16)<<init_guess2<<endl;
        cout<<" s_a = "<<setprecision(16)<<sa<<'\n'
            <<" s_b = "<<setprecision(16)<<sb<<'\n'
            <<" s_c = "<<setprecision(16)<<sc<<endl;
        cout<<"iteration:"<<i+1<<'\t'<<" s:"<<sc<<'\t'<<" f:"<<fc<<endl;

        if(abs(fc)<epsilon || abs(sb-sa)/2.0 < tolerance )
        //if(abs(fc)<epsilon)
        {
            cout<<"---------------------------"<<endl;
            cout<<"fc:"<<fc<<'\t'<<"eps:"<<epsilon<<endl;
            cout<<"b-a/2 : "<<(sb-sa)/2.0<<endl;
            pole = sc;
            cout<<"pole found at s = "<<setprecision(18)<<sc<<endl;
            cout<<"---------------------------"<<endl;
            cout<<endl;
            break;
        }
        //if(fc>0.0 && fa>0.0) 
        if(sign_c==sign_a)
        {
            sa = sc;
        }
        else 
        {
            sb = sc;
        }

        if(i==max_iteration-1)
        {
            cout<<"method failed, max iteration reached"<<endl;
        }
    }
    //cout<<"pole found at s = "<<s<<endl;
    //pole = s;
}



void poletest()
{
    double m = 1.0;
    double pole = 0.0;
    vector<double> polevec;
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.99;
    double init_guess2 = 0.0;//(double) real(phibthreshold(a,m));
    for(double a = 50.0; a<1000.0; a = a + 50.0)
    {
        double dels = 0.00001;
        double tempval = 0.0;
        double set_s = (double) real(phibthreshold(a,m));
        init_guess2 = (double) real(phibthreshold(a,m));
    
       //(double) real(phibthreshold(a,m));
        for(int i=0;i<100;++i)
        {
            double s = set_s - i*dels;
            tempval = real(dSqqs2q2msq_func(s,a));
            cout<<"tempval = "<<tempval<<endl;
            if(tempval<0.0)
            {
                init_guess2 = s;
                
                cout<<"phibth = "<<set_s<<'\t'
                    <<"inti_guess2 = "<<init_guess2<<endl;
                break;
            }
            if(i==100-1)
            {
                cout<<"init_guess2 not found"<<endl;
                cout<<"tempval = "<<tempval<<endl;
            }
        }
        
        
        dSqqs2k2msq_pole_searching_bissection(init_guess1,init_guess2,a,pole);
        polevec.push_back(pole);
        //init_guess2 = pole;
    }

    for(int i=0;i<polevec.size();++i)
    {
        double a1 = 50.0 + i*50.0;
        double sb = polevec[i];
        cout<<setprecision(16)<<a1<<'\t'<<sb<<'\t'
            <<sqrt(sb) - 3.0<<endl;
    }
}



void poletest_bs1_vs_N()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.60;
    double init_guess2 = 8.99;//(double) real(phibthreshold(a,m));
    
    for(int N=4000; N<=7000; N=N+500)
    {
        vector<double> polevec;
        vector<double> avec;
        double dela = 50.0;

        string filename = "bs1st_N_" + to_string(N) + "_00.dat";
        
        for(double a = 500.0; a<=100000.0; a = a + dela)
        {
            if(a>=1000.0) dela = 1000.0;
            if(a>=10000.0) dela = 2500.0;
            //double dels = 0.00001;
            //double tempval = 0.0;
            //double set_s = (double) real(phibthreshold(a,m));
            //init_guess2 = (double) real(phibthreshold(a,m));
    
            //(double) real(phibthreshold(a,m));
            
        
            cout<<"BS1 test"<<endl;
            dSqqs2k2msq_pole_searching_bissection_with_N(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
            double sb = polevec[i];
            double am = avec[i];
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
        }
        fout.close();
        cout<<"==============================="<<endl;
    }

    
}

void poletest_bs1_vs_singleN()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.60;
    double init_guess2 = 8.99;//(double) real(phibthreshold(a,m));
    
    //for(int N=4000; N<=7000; N=N+500)
    int N = 2000;
    {
        vector<double> polevec;
        vector<double> avec;
        double dela = 5000.0;

        string filename = "bs1st_N_" + to_string(N) + "_withanalyticH.dat";
        
        for(double a = 5000.0; a<=100000.0; a = a + dela)
        {
            //if(a>=1000.0) dela = 1000.0;
            //if(a>=10000.0) dela = 2500.0;
            //double dels = 0.00001;
            //double tempval = 0.0;
            //double set_s = (double) real(phibthreshold(a,m));
            //init_guess2 = (double) real(phibthreshold(a,m));
    
            //(double) real(phibthreshold(a,m));
            
        
            cout<<"BS1 test"<<endl;
            dSqqs2k2msq_pole_searching_bissection_with_N(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
            double sb = polevec[i];
            double am = avec[i];
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
        }
        fout.close();
        cout<<"==============================="<<endl;
    }

    
}

void poletest_bs1_vs_singleN_singlea(   double a,
                                        int N,
                                        double &spole )
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.60;
    double init_guess2 = 8.99;//(double) real(phibthreshold(a,m));
    
    //for(int N=4000; N<=7000; N=N+500)
    //int N = 1000;
    {
        vector<double> polevec;
        vector<double> avec;
        double dela = 5000.0;

        //string filename = "bs1st_N_" + to_string(N) + "_withanalyticH.dat";
        
        //for(double a = 5000.0; a<=100000.0; a = a + dela)
        //double a = 16.0;
        string filename = "bs1st_a_" + to_string((int)a) + "_N_" + to_string(N) + "_smoothcutoff.dat";
        
        {
            //if(a>=1000.0) dela = 1000.0;
            //if(a>=10000.0) dela = 2500.0;
            //double dels = 0.00001;
            //double tempval = 0.0;
            //double set_s = (double) real(phibthreshold(a,m));
            //init_guess2 = (double) real(phibthreshold(a,m));
    
            //(double) real(phibthreshold(a,m));
            
        
            cout<<"BS1 test"<<endl;
            dSqqs2k2msq_pole_searching_bissection_with_N(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
            double sb = polevec[i];
            double am = avec[i];
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
            spole = polevec[i];
        }
        fout.close();
        cout<<"==============================="<<endl;
    }

    
}

void poletest_bs1_vs_singleN_singlea_using_determinant(     double a,
                                                            int N,
                                                            double &spole )
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.60;
    double init_guess2 = 8.99;//(double) real(phibthreshold(a,m));
    
    //for(int N=4000; N<=7000; N=N+500)
    //int N = 1000;
    {
        vector<double> polevec;
        vector<double> avec;
        double dela = 5000.0;

        //string filename = "bs1st_N_" + to_string(N) + "_withanalyticH.dat";
        
        //for(double a = 5000.0; a<=100000.0; a = a + dela)
        //double a = 16.0;
        string filename = "bs1st_a_" + to_string((int)a) + "_N_" + to_string(N) + "_hardcutoff.dat";
        
        {
            //if(a>=1000.0) dela = 1000.0;
            //if(a>=10000.0) dela = 2500.0;
            //double dels = 0.00001;
            //double tempval = 0.0;
            //double set_s = (double) real(phibthreshold(a,m));
            //init_guess2 = (double) real(phibthreshold(a,m));
    
            //(double) real(phibthreshold(a,m));
            
        
            cout<<"BS1 test"<<endl;
            dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
            double sb = polevec[i];
            double am = avec[i];
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
            
        }
        fout.close();
        cout<<"==============================="<<endl;
        spole = polevec[0];
    }

    
}

void poletest_bs1_vs_N_using_determinant()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.60;
    double init_guess2 = 8.99;//(double) real(phibthreshold(a,m));
    
    for(int N=3000; N<=6000; N=N+500)
    //int N = 1000;
    {
        vector<double> polevec;
        vector<double> avec;
        double dela = 100.0;

        //string filename = "bs1st_N_" + to_string(N) + "_withanalyticH.dat";
        string filename = "bs1st_N_" + to_string(N) + "_smoothcutoff_LU.dat";
        
        for(double a = 300.0; a<=100000.0; a = a + dela)
        //double a = 16.0;
        
        {
            if(a>=1000.0) dela = 1000.0;
            if(a>=10000.0) dela = 5000.0;
            //double dels = 0.00001;
            //double tempval = 0.0;
            //double set_s = (double) real(phibthreshold(a,m));
            //init_guess2 = (double) real(phibthreshold(a,m));
    
            //(double) real(phibthreshold(a,m));
            
        
            cout<<"BS1 test"<<endl;
            dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
            double sb = polevec[i];
            double am = avec[i];
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
            //spole = polevec[i];
        }
        fout.close();
        cout<<"==============================="<<endl;
    }

    
}

void poletest_bs1_vs_N_using_weights()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    double a = 16.0;//1.45;//6.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.5;//7.1;//8.4;//5.97;////8.4;//7.0;//8.60;
    double init_guess2 = 8.8;//(double) real(phibthreshold(a,m));////8.6;//7.5;//8.99;//(double) real(phibthreshold(a,m));
    //string filename = "bs1st_a_" +to_string((int)a) +"_N_hardcutoff.dat";
        
    //for(int N=500; N<=1000; N=N+100)
    int N = 500;
    {
        vector<double> polevec;
        vector<double> avec;
        double dela = 100.0;

        //string filename = "bs1st_N_" + to_string(N) + "_withanalyticH.dat";
        //string filename = "bs1st_N_" + to_string(N) + "_smoothcutoff.dat";
        string filename = "bs1st_a_" +to_string((int)a) +"_N_" + to_string(N) + ".dat";
        

        for(double a = 2.0; a<=100000.0; a = a + dela)
        //double a = 16.0;
        
        {
            if(a>=1000.0) dela = 1000.0;
            if(a>=10000.0) dela = 5000.0;
            //double dels = 0.00001;
            //double tempval = 0.0;
            //double set_s = (double) real(phibthreshold(a,m));
            //init_guess2 = (double) real(phibthreshold(a,m));
    
            //(double) real(phibthreshold(a,m));
            
        
            cout<<"BS1 test"<<endl;
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            dSqqs2k2msq_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            
            polevec.push_back(pole);
            avec.push_back(a);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
            double sb = polevec[i];
            double am = avec[i];
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt((double) real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt((double) real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<endl;
            //spole = polevec[i];
        }
        fout.close();
        cout<<"==============================="<<endl;
    }

    
}

void poletest_bs1_vs_N_using_weights1()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    //double a = 16.0;//1.45;//6.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.575160206957158 - 0.01;//6.5;//7.1;//8.4;//5.97;////8.4;//7.0;//8.60;
    ////8.6;//7.5;//8.99;//(double) real(phibthreshold(a,m));
    //string filename = "bs1st_a_" +to_string((int)a) +"_N_hardcutoff.dat";
        
    //for(int N=500; N<=1000; N=N+100)
    int N = 500;
    {
        vector<double> polevec;
        vector<double> avec;
        double dela = 0.5;

        //string filename = "bs1st_N_" + to_string(N) + "_withanalyticH.dat";
        //string filename = "bs1st_N_" + to_string(N) + "_smoothcutoff.dat";
        string filename = "bs1st_N_" + to_string(N) + "_hardcutoff_1.dat";
        fout.open(filename.c_str());

        double init_guess2 = 8.933456276189833;

        for(double a = 10.0; a<=100000.0; a = a + dela)
        //double a = 16.0;
        
        {
            //(double) real(phibthreshold(a,m));
            if(a>=10.0) dela = 10.0;
            if(a>=100.0) dela = 100.0;
            if(a>=1000.0) dela = 1000.0;
            if(a>=10000.0) dela = 10000.0;
            if(a>=100000.0) dela = 100000.0;
            //double dels = 0.00001;
            //double tempval = 0.0;
            //double set_s = (double) real(phibthreshold(a,m));
            //init_guess2 = (double) real(phibthreshold(a,m));
    
            //(double) real(phibthreshold(a,m));
            
        
            cout<<"BS1 test"<<endl;
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            dSqqs2k2msq_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            
            polevec.push_back(pole);
            avec.push_back(a);
            double am = a;
            double sb = pole ;
            
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt((double) real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt((double) real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<(double) real(phibthreshold(am,m))<<endl;
            init_guess1 = pole - 0.03;
            //init_guess2 = (double) real(phibthreshold(am,m));
        }

        
        
        fout.close();
        cout<<"==============================="<<endl;
    }

    
}



void poletest_bs1_vs_singleN_singlea_using_weights(double a,
                                                            int N,
                                                            double &spole)
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.60;
    double init_guess2 = 8.99;//(double) real(phibthreshold(a,m));
    
    //for(int N=3000; N<=6000; N=N+500)
    //int N = 500;
    {
        vector<double> polevec;
        vector<double> avec;
        double dela = 100.0;

        //string filename = "bs1st_N_" + to_string(N) + "_withanalyticH.dat";
        //string filename = "bs1st_N_" + to_string(N) + "_hardcutoff.dat";
        
        //for(double a = 500.0; a<=100000.0; a = a + dela)
        //double a = 16.0;
        
        {
            if(a>=1000.0) dela = 1000.0;
            if(a>=10000.0) dela = 5000.0;
            //double dels = 0.00001;
            //double tempval = 0.0;
            //double set_s = (double) real(phibthreshold(a,m));
            //init_guess2 = (double) real(phibthreshold(a,m));
    
            //(double) real(phibthreshold(a,m));
            
        
            cout<<"BS1 test"<<endl;
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            dSqqs2k2msq_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            
            polevec.push_back(pole);
            avec.push_back(a);
            spole = pole;
            //init_guess2 = pole;
        }

        //fout.open(filename.c_str());
        /*for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
            double sb = polevec[i];
            double am = avec[i];
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt((double) real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt((double) real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<endl;
            //spole = polevec[i];
        }
        fout.close();*/
        cout<<"==============================="<<endl;
    }

    
}

void poletest_bs1_vs_singleN_singlea_using_weights_for_vertexfactor(double a,
                                                            int N,
                                                            double &spole)
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.70;
    double init_guess2 = 8.90;//(double) real(phibthreshold(a,m));
    
    //for(int N=3000; N<=6000; N=N+500)
    //int N = 500;
    {
        vector<double> polevec;
        vector<double> avec;
        double dela = 100.0;

        //string filename = "bs1st_N_" + to_string(N) + "_withanalyticH.dat";
        //string filename = "bs1st_N_" + to_string(N) + "_hardcutoff.dat";
        
        //for(double a = 500.0; a<=100000.0; a = a + dela)
        //double a = 16.0;
        
        {
            if(a>=1000.0) dela = 1000.0;
            if(a>=10000.0) dela = 5000.0;
            //double dels = 0.00001;
            //double tempval = 0.0;
            //double set_s = (double) real(phibthreshold(a,m));
            //init_guess2 = (double) real(phibthreshold(a,m));
    
            //(double) real(phibthreshold(a,m));
            
        
            cout<<"BS1 test"<<endl;
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            dSqqs2k2msq_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            
            polevec.push_back(pole);
            avec.push_back(a);
            spole = pole;
            //init_guess2 = pole;
        }

        //fout.open(filename.c_str());
        /*for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
            double sb = polevec[i];
            double am = avec[i];
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt((double) real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt((double) real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<endl;
            //spole = polevec[i];
        }
        fout.close();*/
        cout<<"==============================="<<endl;
    }

    
}



void poletest_bs2_vs_N()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.99;
    double init_guess2 = 0.0;//(double) real(phibthreshold(a,m));
    
    for(int N=4000; N<=7000; N=N+500)
    {
        vector<double> polevec;
        vector<double> avec;
        string filename = "bs2nd_N_" + to_string(N) + "_00.dat";
        double dela = 50.0;
        
        for(double a = 500.0; a<=100000.0; a = a + dela)
        {
            if(a>=1000.0) dela = 1000.0;
            if(a>=10000.0) dela = 2500.0;
            double dels = 0.00001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m));
            init_guess2 = (double) real(phibthreshold(a,m));
            double prevtempval = 0.0;
            double prevs = 0.0;
    
            //(double) real(phibthreshold(a,m));
            for(int i=0;i<100;++i)
            {
                double s = set_s - i*dels;
                tempval = real(dSqqs2q2msq_func_with_N(s,a,N));
                cout<<"am = "<<a<<'\t'<<"N = "<<N<<endl;
                cout<<"value of s = "<<setprecision(16)<<s<<endl;
                cout<<"value of dS = "<<setprecision(16)<<tempval<<endl;
                if(i==0)
                {
                    prevtempval = tempval;
                    prevs = s;
                    continue;
                } 
                if(tempval<0.0)
                {
                    init_guess2 = s;
                
                    cout<<"phibth = "<<setprecision(16)<<set_s<<'\t'
                        <<"init_guess1 = "<<setprecision(16)<<init_guess1<<endl;
                    break;
                }
                prevtempval = tempval;
                prevs = s;
                if(i==100-1)
                {
                    cout<<"init_guess2 not found"<<endl;
                    cout<<"tempval = "<<tempval<<endl;
                }
            }
        
        
            cout<<"BS2 test"<<endl;
            dSqqs2k2msq_pole_searching_bissection_with_N(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
            double sb = polevec[i];
            double am = avec[i];
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
        }
        fout.close();
        cout<<"==============================="<<endl;
    }

    
}


void poletest_bs2_vs_singleN()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.995;
    double init_guess2 = 0.0;//(double) real(phibthreshold(a,m));
    
    //for(int N=4000; N<=7000; N=N+500)
    int N = 2000;
    {
        vector<double> polevec;
        vector<double> avec;
        string filename = "bs2nd_N_" + to_string(N) + "_withanalyticH.dat";
        double dela = 5000.0;
        
        for(double a = 5000.0; a<=100000.0; a = a + dela)
        {
            //if(a>=1000.0) dela = 1000.0;
            //if(a>=10000.0) dela = 2500.0;
            double dels = 0.00001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m));
            init_guess2 = (double) real(phibthreshold(a,m));
            double prevtempval = 0.0;
            double prevs = 0.0;
    
            //(double) real(phibthreshold(a,m));
            for(int i=0;i<100;++i)
            {
                double s = set_s - i*dels;
                tempval = real(dSqqs2q2msq_func_with_N(s,a,N));
                cout<<"am = "<<a<<'\t'<<"N = "<<N<<endl;
                cout<<"value of s = "<<setprecision(16)<<s<<endl;
                cout<<"value of dS = "<<setprecision(16)<<tempval<<endl;
                if(i==0)
                {
                    prevtempval = tempval;
                    prevs = s;
                    continue;
                } 
                if(tempval<0.0)
                {
                    init_guess2 = s;
                
                    cout<<"phibth = "<<setprecision(16)<<set_s<<'\t'
                        <<"init_guess1 = "<<setprecision(16)<<init_guess1<<endl;
                    break;
                }
                prevtempval = tempval;
                prevs = s;
                if(i==100-1)
                {
                    cout<<"init_guess2 not found"<<endl;
                    cout<<"tempval = "<<tempval<<endl;
                }
            }
        
        
            cout<<"BS2 test"<<endl;
            dSqqs2k2msq_pole_searching_bissection_with_N(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
            double sb = polevec[i];
            double am = avec[i];
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
        }
        fout.close();
        cout<<"==============================="<<endl;
    }

    
}

void poletest_bs2_vs_singleN_singlea(   double a,
                                        int N,
                                        double &spole)
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.995;
    double init_guess2 = 0.0;//(double) real(phibthreshold(a,m));
    
    //for(int N=4000; N<=7000; N=N+500)
    //int N = 2000;
    {
        vector<double> polevec;
        vector<double> avec;
        string filename = "bs2nd_a_" + to_string((int)a) + "_N_" + to_string(N) + "_smoothcutoff.dat";
        double dela = 5000.0;
        
        //for(double a = 5000.0; a<=100000.0; a = a + dela)
        {
            //if(a>=1000.0) dela = 1000.0;
            //if(a>=10000.0) dela = 2500.0;
            double dels = 0.00001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m));
            init_guess2 = (double) real(phibthreshold(a,m));
            double prevtempval = 0.0;
            double prevs = 0.0;
    
            //(double) real(phibthreshold(a,m));
            for(int i=0;i<100;++i)
            {
                double s = set_s - i*dels;
                tempval = real(dSqqs2q2msq_func_with_N(s,a,N));
                cout<<"am = "<<a<<'\t'<<"N = "<<N<<endl;
                cout<<"value of s = "<<setprecision(16)<<s<<endl;
                cout<<"value of dS = "<<setprecision(16)<<tempval<<endl;
                if(i==0)
                {
                    prevtempval = tempval;
                    prevs = s;
                    continue;
                } 
                if(tempval<0.0)
                {
                    init_guess2 = s;
                
                    cout<<"phibth = "<<setprecision(16)<<set_s<<'\t'
                        <<"init_guess1 = "<<setprecision(16)<<init_guess1<<endl;
                    break;
                }
                prevtempval = tempval;
                prevs = s;
                if(i==100-1)
                {
                    cout<<"init_guess2 not found"<<endl;
                    cout<<"tempval = "<<tempval<<endl;
                }
            }
        
        
            cout<<"BS2 test"<<endl;
            dSqqs2k2msq_pole_searching_bissection_with_N(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
            double sb = polevec[i];
            double am = avec[i];
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
        }
        fout.close();
        cout<<"==============================="<<endl;
        spole = polevec[0];
    }

    
}

void poletest_bs2_vs_singleN_singlea_using_determinant(   double a,
                                        int N,
                                        double &spole)
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.995;
    double init_guess2 = 8.9999;//.0;//(double) real(phibthreshold(a,m));
    
    //for(int N=4000; N<=7000; N=N+500)
    //int N = 2000;
    {
        vector<double> polevec;
        vector<double> avec;
        string filename = "bs2nd_a_" + to_string((int)a) + "_N_" + to_string(N) + "_hardcutoff_LU.dat";
        double dela = 5000.0;
        
        //for(double a = 5000.0; a<=100000.0; a = a + dela)
        {
            //if(a>=1000.0) dela = 1000.0;
            //if(a>=10000.0) dela = 2500.0;
            double dels = 0.00001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m));
            //init_guess2 = (double) real(phibthreshold(a,m));
            double prevtempval = 0.0;
            double prevs = 0.0;
    
            //(double) real(phibthreshold(a,m));
            /*for(int i=0;i<100;++i)
            {
                double s = set_s - i*dels;
                tempval = real(dSqqs2q2msq_func_with_N(s,a,N));
                cout<<"am = "<<a<<'\t'<<"N = "<<N<<endl;
                cout<<"value of s = "<<setprecision(16)<<s<<endl;
                cout<<"value of dS = "<<setprecision(16)<<tempval<<endl;
                if(i==0)
                {
                    prevtempval = tempval;
                    prevs = s;
                    continue;
                } 
                if(tempval<0.0)
                {
                    init_guess2 = s;
                
                    cout<<"phibth = "<<setprecision(16)<<set_s<<'\t'
                        <<"init_guess1 = "<<setprecision(16)<<init_guess1<<endl;
                    break;
                }
                prevtempval = tempval;
                prevs = s;
                if(i==100-1)
                {
                    cout<<"init_guess2 not found"<<endl;
                    cout<<"tempval = "<<tempval<<endl;
                }
            }*/
        
        
            cout<<"BS2 test"<<endl;
            dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
            double sb = polevec[i];
            double am = avec[i];
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
        }
        fout.close();
        cout<<"==============================="<<endl;
        spole = polevec[0];
    }

    
}

void poletest_bs2_vs_N_using_determinant()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.995;
    double init_guess2 = 8.9999;//.0;//(double) real(phibthreshold(a,m));
    
    //for(int N=3000; N<=6000; N=N+500)
    int N = 7000;
    {
        vector<double> polevec;
        vector<double> avec;
        string filename = "bs2nd_N_" + to_string(N) + "_smoothcutoff_LU.dat";
        double dela = 100.0;
        
        fout.open(filename.c_str());
        for(double a = 300.0; a<=100000.0; a = a + dela)
        {
            if(a>=1000.0) dela = 1000.0;
            if(a>=10000.0) dela = 5000.0;
            double dels = 0.00001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m));
            //init_guess2 = (double) real(phibthreshold(a,m));
            double prevtempval = 0.0;
            double prevs = 0.0;
    
            //(double) real(phibthreshold(a,m));
            /*for(int i=0;i<100;++i)
            {
                double s = set_s - i*dels;
                tempval = real(dSqqs2q2msq_func_with_N(s,a,N));
                cout<<"am = "<<a<<'\t'<<"N = "<<N<<endl;
                cout<<"value of s = "<<setprecision(16)<<s<<endl;
                cout<<"value of dS = "<<setprecision(16)<<tempval<<endl;
                if(i==0)
                {
                    prevtempval = tempval;
                    prevs = s;
                    continue;
                } 
                if(tempval<0.0)
                {
                    init_guess2 = s;
                
                    cout<<"phibth = "<<setprecision(16)<<set_s<<'\t'
                        <<"init_guess1 = "<<setprecision(16)<<init_guess1<<endl;
                    break;
                }
                prevtempval = tempval;
                prevs = s;
                if(i==100-1)
                {
                    cout<<"init_guess2 not found"<<endl;
                    cout<<"tempval = "<<tempval<<endl;
                }
            }*/
        
        
            cout<<"BS2 test"<<endl;
            dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);

            double sb = pole;
            double am = a;
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<endl;
            //init_guess2 = pole;
        }
        fout.close();
        cout<<"==============================="<<endl;
        //spole = polevec[0];
    }

    
}

void poletest_bs2_vs_N_using_weights()
{
    ofstream fout;
    double a = 16.0;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.974;//8.99;//8.995;
    double init_guess2 = (double) real(phibthreshold(a,m));//8.9999;//.0;//(double) real(phibthreshold(a,m));
    
    string filename = "bs2nd_a_"+to_string((int)a)+"_N_smoothcutoff.dat";
    fout.open(filename.c_str());
    for(int N=500; N<=1000; N=N+100)
    //int N = 500;
    {
        
        vector<double> polevec;
        vector<double> avec;
        //string filename = "bs2nd_N_" + to_string(N) + "_smoothcutoff.dat";
        //string filename = "bs2nd_a_"+to_string((int)a)+"_N_" + to_string(N) + "_smoothcutoff.dat";
        double dela = 100.0;
        
        //fout.open(filename.c_str());
        //for(double a = 500.0; a<=100000.0; a = a + dela)
        {
            if(a>=1000.0) dela = 1000.0;
            if(a>=10000.0) dela = 5000.0;
            double dels = 0.00001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m));
            //init_guess2 = (double) real(phibthreshold(a,m));
            double prevtempval = 0.0;
            double prevs = 0.0;
    
            //(double) real(phibthreshold(a,m));
            /*for(int i=0;i<100;++i)
            {
                double s = set_s - i*dels;
                tempval = real(dSqqs2q2msq_func_with_N(s,a,N));
                cout<<"am = "<<a<<'\t'<<"N = "<<N<<endl;
                cout<<"value of s = "<<setprecision(16)<<s<<endl;
                cout<<"value of dS = "<<setprecision(16)<<tempval<<endl;
                if(i==0)
                {
                    prevtempval = tempval;
                    prevs = s;
                    continue;
                } 
                if(tempval<0.0)
                {
                    init_guess2 = s;
                
                    cout<<"phibth = "<<setprecision(16)<<set_s<<'\t'
                        <<"init_guess1 = "<<setprecision(16)<<init_guess1<<endl;
                    break;
                }
                prevtempval = tempval;
                prevs = s;
                if(i==100-1)
                {
                    cout<<"init_guess2 not found"<<endl;
                    cout<<"tempval = "<<tempval<<endl;
                }
            }*/
        
        
            cout<<"BS2 test"<<endl;
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            dSqqs2k2msq_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            
            polevec.push_back(pole);
            avec.push_back(a);

            double sb = pole;
            double am = a;
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt((double) real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt((double) real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<endl;
            //init_guess2 = pole;
        }
        //fout.close();
        cout<<"==============================="<<endl;
        //spole = polevec[0];
    }
    fout.close();

    
}

void poletest_bs2_vs_N_using_weights1()
{
    ofstream fout;
    //double a = 16.0;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.590045949378075+0.15;//8.999098205302705;//8.99;//8.995;
    
    string filename = "bs2nd_N_hardcutoff_3.dat";
    fout.open(filename.c_str());
    //for(int N=500; N<=1000; N=N+100)
    int N = 500;
    {
        
        vector<double> polevec;
        vector<double> avec;
        //string filename = "bs2nd_N_" + to_string(N) + "_smoothcutoff.dat";
        //string filename = "bs2nd_a_"+to_string((int)a)+"_N_" + to_string(N) + "_smoothcutoff.dat";
        double dela = 0.5;//100.0;
        
        //fout.open(filename.c_str());
        for(double a = 13.0; a<=15.0; a = a + dela)
        {
            double init_guess2 = (double) real(phibthreshold(a,m));//8.9999;//.0;//(double) real(phibthreshold(a,m));
            
            /*if(a>=10.0) dela = 10.0;
            if(a>=100.0) dela = 100.0;
            if(a>=1000.0) dela = 1000.0;
            if(a>=10000.0) dela = 10000.0;
            if(a>=100000.0) dela = 100000.0;
            if(a>=1000000.0) dela = 1000000.0;*/
            dela = 0.5;
            double dels = 0.00001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m));
            //init_guess2 = (double) real(phibthreshold(a,m));
            double prevtempval = 0.0;
            double prevs = 0.0;
    
            //(double) real(phibthreshold(a,m));
            /*for(int i=0;i<100;++i)
            {
                double s = set_s - i*dels;
                tempval = real(dSqqs2q2msq_func_with_N(s,a,N));
                cout<<"am = "<<a<<'\t'<<"N = "<<N<<endl;
                cout<<"value of s = "<<setprecision(16)<<s<<endl;
                cout<<"value of dS = "<<setprecision(16)<<tempval<<endl;
                if(i==0)
                {
                    prevtempval = tempval;
                    prevs = s;
                    continue;
                } 
                if(tempval<0.0)
                {
                    init_guess2 = s;
                
                    cout<<"phibth = "<<setprecision(16)<<set_s<<'\t'
                        <<"init_guess1 = "<<setprecision(16)<<init_guess1<<endl;
                    break;
                }
                prevtempval = tempval;
                prevs = s;
                if(i==100-1)
                {
                    cout<<"init_guess2 not found"<<endl;
                    cout<<"tempval = "<<tempval<<endl;
                }
            }*/
        
        
            cout<<"BS2 test"<<endl;
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            dSqqs2k2msq_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            
            polevec.push_back(pole);
            avec.push_back(a);

            double sb = pole;
            double am = a;
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt((double) real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt((double) real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<(double) real(phibthreshold(am,m))<<endl;
            //init_guess2 = pole;
        }
        //fout.close();
        cout<<"==============================="<<endl;
        //spole = polevec[0];
    }
    fout.close();

    
}


void poletest_bs3_vs_N()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.99;
    double init_guess2 = 0.0;//(double) real(phibthreshold(a,m));
    
    
    for(int N=3000; N<=7000; N=N+500)
    {
        vector<double> polevec;
        string filename = "bs3rd_N_" + to_string(N) + "_00.dat";
        vector<double> avec;
        double dela = 50.0;
        
        for(double a = 1000.0; a<=4000.0; a = a + dela)
        {
            if(a>=1000) dela = 1000;
            if(a>=10000) dela = 2500;
            double dels = 0.000001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m));
            double prevtempval=INT_MAX;
            double prevs=0.0;
    
            //(double) real(phibthreshold(a,m));
            for(int i=0;i<100;++i)
            {
                double s = set_s - i*dels;
                tempval = real(dSqqs2q2msq_func_with_N(s,a,N));
                cout<<"am = "<<setprecision(16)<<a<<'\t'<<"N = "<<N<<endl;
                cout<<"value of s = "<<setprecision(16)<<s<<endl;
                cout<<"value of dS = "<<setprecision(16)<<tempval<<endl;
                if(i==0)
                {
                    prevtempval = tempval;
                    prevs = s;
                    continue;
                } 
                if(tempval<0.0 && prevtempval>0.0)
                {
                    init_guess1 = prevs;
                
                    cout<<"phibth = "<<setprecision(16)<<set_s<<'\t'
                        <<"init_guess1 = "<<setprecision(16)<<init_guess1<<endl;
                    break;
                }
                
                prevtempval = tempval;
                prevs = s;
                if(i==100-1)
                {
                    cout<<"init_guess2 not found"<<endl;
                    cout<<"tempval = "<<tempval<<endl;
                    cout<<"setting to current s = "<<s<<endl;
                    init_guess1 = s;
                }
            }
            //init_guess1 = 8.99990000001;
            init_guess2 = (double) real(phibthreshold(a,m)); 
        
            
        
            dSqqs2k2msq_pole_searching_bissection_with_N(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
            double am = avec[i];
            double sb = polevec[i];
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
        }
        fout.close();
        cout<<"==============================="<<endl;
    }

    
}

void poletest_bs3_vs_singleN()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.99;
    double init_guess2 = 0.0;//(double) real(phibthreshold(a,m));
    
    int N = 2000;
    //for(int N=3000; N<=7000; N=N+500)
    {
        vector<double> polevec;
        string filename = "bs3rd_N_" + to_string(N) + "_1.dat";
        vector<double> avec;
        double dela = 5000.0;
        
        for(double a = 5000.0; a<=100000.0; a = a + dela)
        {
            //if(a>=1000) dela = 1000;
            //if(a>=10000) dela = 2500;
            double dels = 0.000001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m));
            double prevtempval=INT_MAX;
            double prevs=0.0;
    
            //(double) real(phibthreshold(a,m));
            for(int i=0;i<100;++i)
            {
                double s = set_s - i*dels;
                tempval = real(dSqqs2q2msq_func_with_N(s,a,N));
                cout<<"am = "<<setprecision(16)<<a<<'\t'<<"N = "<<N<<endl;
                cout<<"value of s = "<<setprecision(16)<<s<<endl;
                cout<<"value of dS = "<<setprecision(16)<<tempval<<endl;
                if(i==0)
                {
                    prevtempval = tempval;
                    prevs = s;
                    continue;
                } 
                if(tempval<0.0 && prevtempval>0.0)
                {
                    init_guess1 = prevs;
                
                    cout<<"phibth = "<<setprecision(16)<<set_s<<'\t'
                        <<"init_guess1 = "<<setprecision(16)<<init_guess1<<endl;
                    break;
                }
                
                prevtempval = tempval;
                prevs = s;
                if(i==100-1)
                {
                    cout<<"init_guess2 not found"<<endl;
                    cout<<"tempval = "<<tempval<<endl;
                    cout<<"setting to current s = "<<s<<endl;
                    init_guess1 = s;
                }
            }
            //init_guess1 = 8.99990000001;
            init_guess2 = (double) real(phibthreshold(a,m)); 
        
            
        
            dSqqs2k2msq_pole_searching_bissection_with_N(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
            double am = avec[i];
            double sb = polevec[i];
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
        }
        fout.close();
        cout<<"==============================="<<endl;
    }

    
}

void poletest_bs3_vs_singleN_singlea(   double a,
                                        int N,
                                        double &spole)
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.99;
    double init_guess2 = 0.0;//(double) real(phibthreshold(a,m));
    
    //int N = 2000;
    //for(int N=3000; N<=7000; N=N+500)
    {
        vector<double> polevec;
        string filename = "bs3rd_a_" + to_string((int)a) + "_N_" + to_string(N) + "_smoothcutoff.dat";
        vector<double> avec;
        double dela = 5000.0;
        
        //for(double a = 5000.0; a<=100000.0; a = a + dela)
        {
            //if(a>=1000) dela = 1000;
            //if(a>=10000) dela = 2500;
            double dels = 0.000001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m));
            double prevtempval=INT_MAX;
            double prevs=0.0;
    
            //(double) real(phibthreshold(a,m));
            for(int i=0;i<100;++i)
            {
                double s = set_s - i*dels;
                tempval = real(dSqqs2q2msq_func_with_N(s,a,N));
                cout<<"am = "<<setprecision(16)<<a<<'\t'<<"N = "<<N<<endl;
                cout<<"value of s = "<<setprecision(16)<<s<<endl;
                cout<<"value of dS = "<<setprecision(16)<<tempval<<endl;
                if(i==0)
                {
                    prevtempval = tempval;
                    prevs = s;
                    continue;
                } 
                if(tempval<0.0 && prevtempval>0.0)
                {
                    init_guess1 = prevs;
                
                    cout<<"phibth = "<<setprecision(16)<<set_s<<'\t'
                        <<"init_guess1 = "<<setprecision(16)<<init_guess1<<endl;
                    break;
                }
                
                prevtempval = tempval;
                prevs = s;
                if(i==100-1)
                {
                    cout<<"init_guess2 not found"<<endl;
                    cout<<"tempval = "<<tempval<<endl;
                    cout<<"setting to current s = "<<s<<endl;
                    init_guess1 = s;
                }
            }
            //init_guess1 = 8.99990000001;
            init_guess2 = (double) real(phibthreshold(a,m)); 
        
            
        
            dSqqs2k2msq_pole_searching_bissection_with_N(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
            double am = avec[i];
            double sb = polevec[i];
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
        }
        fout.close();
        cout<<"==============================="<<endl;
        spole = polevec[0];
    }

    
}

void poletest_bs3_vs_singleN_singlea_using_determinant(     double a,
                                                            int N,
                                                            double &spole)
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.99995;
    double init_guess2 = (double) real(phibthreshold(a,m));
    
    //int N = 2000;
    //for(int N=3000; N<=7000; N=N+500)
    {
        vector<double> polevec;
        string filename = "bs3rd_a_" + to_string((int)a) + "_N_" + to_string(N) + "_smoothcutoff_LU.dat";
        vector<double> avec;
        double dela = 5000.0;
        
        //for(double a = 5000.0; a<=100000.0; a = a + dela)
        {
            //if(a>=1000) dela = 1000;
            //if(a>=10000) dela = 2500;
            double dels = 0.000001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m));
            double prevtempval=INT_MAX;
            double prevs=0.0;
    
            //(double) real(phibthreshold(a,m));
            /*for(int i=0;i<100;++i)
            {
                double s = set_s - i*dels;
                tempval = real(dSqqs2q2msq_func_with_N(s,a,N));
                cout<<"am = "<<setprecision(16)<<a<<'\t'<<"N = "<<N<<endl;
                cout<<"value of s = "<<setprecision(16)<<s<<endl;
                cout<<"value of dS = "<<setprecision(16)<<tempval<<endl;
                if(i==0)
                {
                    prevtempval = tempval;
                    prevs = s;
                    continue;
                } 
                if(tempval<0.0 && prevtempval>0.0)
                {
                    init_guess1 = prevs;
                
                    cout<<"phibth = "<<setprecision(16)<<set_s<<'\t'
                        <<"init_guess1 = "<<setprecision(16)<<init_guess1<<endl;
                    break;
                }
                
                prevtempval = tempval;
                prevs = s;
                if(i==100-1)
                {
                    cout<<"init_guess2 not found"<<endl;
                    cout<<"tempval = "<<tempval<<endl;
                    cout<<"setting to current s = "<<s<<endl;
                    init_guess1 = s;
                }
            }*/
            //init_guess1 = 8.99990000001;
            //init_guess2 = (double) real(phibthreshold(a,m)); 
        
            
            cout<<"BS3 test"<<endl;
            dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
            double am = avec[i];
            double sb = polevec[i];
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
        }
        fout.close();
        cout<<"==============================="<<endl;
        spole = polevec[0];
    }

    
}


void poletest_bs3_vs_N_using_determinant()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.99995;
    double init_guess2 = 0.0;//(double) real(phibthreshold(a,m));
    
    int N = 7000;
    //for(int N=3000; N<=6000; N=N+500)
    {
        vector<double> polevec;
        string filename = "bs3rd_N_" + to_string(N) + "_smoothcutoff_LU.dat";
        vector<double> avec;
        double dela = 100.0;
        fout.open(filename.c_str());
        for(double a = 500.0; a<=100000.0; a = a + dela)
        {
            if(a==500.0) init_guess1 = 8.9997;
            else init_guess1 = pole;

            if(a>=1000) dela = 1000;
            if(a>=10000) dela = 5000;
            double dels = 0.000001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m));
            double prevtempval=INT_MAX;
            double prevs=0.0;
    
            //(double) real(phibthreshold(a,m));
            /*for(int i=0;i<100;++i)
            {
                double s = set_s - i*dels;
                tempval = real(dSqqs2q2msq_func_with_N(s,a,N));
                cout<<"am = "<<setprecision(16)<<a<<'\t'<<"N = "<<N<<endl;
                cout<<"value of s = "<<setprecision(16)<<s<<endl;
                cout<<"value of dS = "<<setprecision(16)<<tempval<<endl;
                if(i==0)
                {
                    prevtempval = tempval;
                    prevs = s;
                    continue;
                } 
                if(tempval<0.0 && prevtempval>0.0)
                {
                    init_guess1 = prevs;
                
                    cout<<"phibth = "<<setprecision(16)<<set_s<<'\t'
                        <<"init_guess1 = "<<setprecision(16)<<init_guess1<<endl;
                    break;
                }
                
                prevtempval = tempval;
                prevs = s;
                if(i==100-1)
                {
                    cout<<"init_guess2 not found"<<endl;
                    cout<<"tempval = "<<tempval<<endl;
                    cout<<"setting to current s = "<<s<<endl;
                    init_guess1 = s;
                }
            }*/
            //init_guess1 = 8.99990000001;
            init_guess2 = (double) real(phibthreshold(a,m)); 
        
            
            cout<<"BS3 test"<<endl;
            dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);

            double am = a;
            double sb = pole;
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            //init_guess2 = pole;
        }

        fout.close();
        cout<<"==============================="<<endl;
        
        
        //spole = polevec[0];
    }

    
}

void poletest_bs3_vs_N_using_weights()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.99997;//8.99995;
    double init_guess2 = 0.0;//(double) real(phibthreshold(a,m));
    
    int N = 500;
    //for(int N=3000; N<=6000; N=N+500)
    {
        vector<double> polevec;
        string filename = "bs3rd_N_" + to_string(N) + "_hardcutoff.dat";
        vector<double> avec;
        double dela = 100.0;
        fout.open(filename.c_str());
        for(double a = 500.0; a<=100000.0; a = a + dela)
        {
            //if(a==500.0) init_guess1 = 8.99997;
            //else init_guess1 = (pole + init_guess1)/2.0;

            if(a>=1000) dela = 1000;
            if(a>=10000) dela = 5000;
            double dels = 0.00000025;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m));
            double prevtempval=INT_MAX;
            double prevs=0.0;
            init_guess2 = (double) real(phibthreshold(a,m));
            for(int i=0;;++i)
            {
                //double check_s = (double) real(phibthreshold(a,m)) - i*dels;
                double check_s = init_guess1 + i*dels;
                
                double check_func = real(dSqqs2q2msq_func_with_N_using_weights(check_s,a,N));
                cout<<" guessing init_guess1 at s="<<check_s<<'\t'<<" f = "<<check_func<<endl;
                //if(sign_checker(check_func)!=sign_checker(real(dSqqs2q2msq_func_with_N_using_weights(init_guess2,a,N))))
                if(sign_checker(check_func)!=sign_checker(init_guess1))
                
                {
                    //init_guess1 = check_s;
                    init_guess2 = check_s;
                    break;
                }
            }
             
        
            
            cout<<"BS3 test"<<endl;
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            dSqqs2k2msq_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            
            polevec.push_back(pole);
            avec.push_back(a);

            double am = a;
            double sb = pole;
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt((double) real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt((double) real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            //init_guess2 = pole;
        }

        fout.close();
        cout<<"==============================="<<endl;
        
        
        //spole = polevec[0];
    }

    
}

void poletest_bs3_vs_N_using_weights1()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.99997;//8.99995;
    double init_guess2 = 0.0;//(double) real(phibthreshold(a,m));
    
    int N = 500;
    //for(int N=3000; N<=6000; N=N+500)
    {
        vector<double> polevec;
        string filename = "bs3rd_N_" + to_string(N) + "_hardcutoff_1.dat";
        vector<double> avec;
        double dela = 1000.0;
        fout.open(filename.c_str());
        for(double a = 10000.0; a<=20000.0; a = a + dela)
        {
            //if(a==500.0) init_guess1 = 8.99997;
            //else init_guess1 = (pole + init_guess1)/2.0;

            if(a>=1000) dela = 100;
            if(a>=2000) dela = 1000;
            if(a>=10000) dela = 1000;
            if(a>=100000) dela = 100000;
            double dels = 0.00000025;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m));
            double prevtempval=INT_MAX;
            double prevs=0.0;
            init_guess2 = (double) real(phibthreshold(a,m));
            for(int i=0;;++i)
            {
                //double check_s = (double) real(phibthreshold(a,m)) - i*dels;
                double check_s = init_guess1 + i*dels;
                
                double check_func = real(dSqqs2q2msq_func_with_N_using_weights(check_s,a,N));
                cout<<" guessing init_guess1 at s="<<check_s<<'\t'<<" f = "<<check_func<<endl;
                //if(sign_checker(check_func)!=sign_checker(real(dSqqs2q2msq_func_with_N_using_weights(init_guess2,a,N))))
                if(sign_checker(check_func)!=sign_checker(init_guess1))
                
                {
                    //init_guess1 = check_s;
                    init_guess2 = check_s;
                    break;
                }
            }
             
        
            
            cout<<"BS3 test"<<endl;
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            dSqqs2k2msq_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            
            polevec.push_back(pole);
            avec.push_back(a);

            double am = a;
            double sb = pole;
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt((double) real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt((double) real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            //init_guess2 = pole;
        }

        fout.close();
        cout<<"==============================="<<endl;
        
        
        //spole = polevec[0];
    }

    
}


void poletest_bs4_vs_N_using_weights()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.99997;//8.99995;
    double init_guess2 = 0.0;//(double) real(phibthreshold(a,m));
    
    int N = 500;
    //for(int N=3000; N<=6000; N=N+500)
    {
        vector<double> polevec;
        string filename = "bs4th_N_" + to_string(N) + "_hardcutoff.dat";
        vector<double> avec;
        double dela = 100.0;
        fout.open(filename.c_str());
        for(double a = 7000.0; a<=100000.0; a = a + dela)
        {
            //if(a==500.0) init_guess1 = 8.99997;
            //else init_guess1 = (pole + init_guess1)/2.0;

            if(a>=1000) dela = 1000;
            if(a>=10000) dela = 5000;
            double dels = 0.0000000025;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m));
            double prevtempval=INT_MAX;
            double prevs=0.0;
            init_guess2 = (double) real(phibthreshold(a,m));
            for(int i=0;;++i)
            {
                double check_s = (double) real(phibthreshold(a,m)) - i*dels;
                //double check_s = init_guess1 + i*dels;
                
                double check_func = real(dSqqs2q2msq_func_with_N_using_weights(check_s,a,N));
                cout<<" guessing init_guess1 at s="<<check_s<<'\t'<<" f = "<<check_func<<endl;
                if(sign_checker(check_func)!=sign_checker(real(dSqqs2q2msq_func_with_N_using_weights(init_guess2,a,N))))
                //if(sign_checker(check_func)!=sign_checker(init_guess1))
                
                {
                    init_guess1 = check_s;
                    //init_guess2 = check_s;
                    break;
                }
            }
             
        
            
            cout<<"BS4 test"<<endl;
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            dSqqs2k2msq_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            
            polevec.push_back(pole);
            avec.push_back(a);

            double am = a;
            double sb = pole;
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt((double) real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt((double) real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            //init_guess2 = pole;
        }

        fout.close();
        cout<<"==============================="<<endl;
        
        
        //spole = polevec[0];
    }

    
}

void poletest_virtual_vs_N_using_weights()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.99995;
    double init_guess2 = 0.0;//(double) real(phibthreshold(a,m));
    
    int N = 500;
    //for(int N=3000; N<=6000; N=N+500)
    {
        double a = 6.0;
        vector<double> polevec;
        //string filename = "bs_virtual_N_" + to_string(N) + "_0_smoothcutoff.dat";
        string filename = "bs_virtual_a_"+to_string(a)+"_N_" + to_string(N) + "_smoothcutoff.dat";
        
        vector<double> avec;
        double dela = 0.0005;
        fout.open(filename.c_str());
        
        //for(double a = 4194.015; a<=5873.320000; a = a + dela)
        
        //for(double a = 1.001; a<=10000.00; a = a + dela)
        {
            //if(a==500.0) init_guess1 = 8.9997;
            //else init_guess1 = pole;
            if(a>=1.02) dela = 0.3;
            if(a>=10) dela = 1.01;
            if(a>=20) dela = 3.01;
            if(a>=100) dela = 30.01;
            if(a>=1000) dela = 100.01;
            double dels = 0.000001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m)) - 1.0e-9;
            double prevtempval=INT_MAX;
            double prevs=0.0;
    
            //(double) real(phibthreshold(a,m));
            /*for(int i=0;i<100;++i)
            {
                double s = set_s - i*dels;
                tempval = real(dSqqs2q2msq_func_with_N(s,a,N));
                cout<<"am = "<<setprecision(16)<<a<<'\t'<<"N = "<<N<<endl;
                cout<<"value of s = "<<setprecision(16)<<s<<endl;
                cout<<"value of dS = "<<setprecision(16)<<tempval<<endl;
                if(i==0)
                {
                    prevtempval = tempval;
                    prevs = s;
                    continue;
                } 
                if(tempval<0.0 && prevtempval>0.0)
                {
                    init_guess1 = prevs;
                
                    cout<<"phibth = "<<setprecision(16)<<set_s<<'\t'
                        <<"init_guess1 = "<<setprecision(16)<<init_guess1<<endl;
                    break;
                }
                
                prevtempval = tempval;
                prevs = s;
                if(i==100-1)
                {
                    cout<<"init_guess2 not found"<<endl;
                    cout<<"tempval = "<<tempval<<endl;
                    cout<<"setting to current s = "<<s<<endl;
                    init_guess1 = s;
                }
            }*/
            //init_guess1 = 8.99990000001;
            init_guess1 = (double) real(Grightbranchpoint(a,m)) + 1.0e-12;
            init_guess2 = (double) real(phibthreshold(a,m)) - 1.0e-12; 
        
            
            //cout<<"Virtual Poles test"<<endl;
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            Mphib_secondsheet_denom_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);

            double am = a;
            double sb = pole;
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<(double) sqrt(real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<setprecision(16)<<(double) sqrt(real(phibthreshold(am,m))) - sqrt(real(Grightbranchpoint(am,m)))<<'\t'<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<(double) sqrt(real(phibthreshold(am,m))) - sqrt(sb)<<'\t'
                <<N<<'\t'
                <<setprecision(16)<<(double) sqrt(real(phibthreshold(am,m))) - sqrt(real(Grightbranchpoint(am,m)))<<'\t'
                <<(double) real(phibthreshold(am,m))<<endl;
            //init_guess2 = pole;

            //cout<<setprecision(16)<<am<<'\t'<<sqrt(real(phibthreshold(am,m))) - sqrt(real(Grightbranchpoint(a,m)))<<endl;
        }

        fout.close();
        //cout<<"==============================="<<endl;
        
        
        //spole = polevec[0];
    }

    
}

void poletest_BS1_vs_N_using_weights_KDF()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    comp KDF = {100.0,0.0};
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.92;
    double init_guess2 = 8.995;//(double) real(phibthreshold(a,m));
    
    int N = 500;
    //for(int N=3000; N<=6000; N=N+500)
    {
        vector<double> polevec;
        string filename = "bs1_KDF_N_" + to_string(N) + "_1.dat";
        vector<double> avec;
        double dela = 100.0;
        fout.open(filename.c_str());
        //double a = 16.0;
        for(double a = 500.0; a<=100000.0; a = a + dela)
        {
            if(a==500.0) init_guess1 = init_guess1;
            else init_guess1 = pole;

            if(a>=1000) dela = 1000;
            if(a>=10000) dela = 5000;
            double dels = 0.000001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m)) - 1.0e-9;
            double prevtempval=INT_MAX;
            double prevs=0.0;
    
            
            //init_guess1 = 8.99990000001;
            //init_guess1 = (double) real(Grightbranchpoint(a,m)) + 1.0e-8;
            //init_guess2 = (double) real(phibthreshold(a,m)) - 1.0e-9; 
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            //init_guess1 = pole;
        
            
            cout<<"KDF Poles test"<<endl;
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            //Mphib_secondsheet_denom_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            KDF_F3_pole_searching_bissection_with_N_using_weights(KDF,init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);

            double am = a;
            double sb = pole;
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<(double) sqrt(real(phibthreshold(am,m))) - sqrt(sb) <<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<(double) sqrt(real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            //init_guess2 = pole;
        }

        fout.close();
        cout<<"==============================="<<endl;
        
        
        //spole = polevec[0];
    }

    
}

void poletest_BS1_vs_N_singlea_using_weights_KDF()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    comp KDF = {100.0,0.0};
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.92;
    double init_guess2 = 8.995;//(double) real(phibthreshold(a,m));
    
    int N = 500;
    //for(int N=3000; N<=6000; N=N+500)
    {
        vector<double> polevec;
        string filename = "bs1_KDF_N_" + to_string(N) + "_tmp.dat";
        vector<double> avec;
        double dela = 100.0;
        fout.open(filename.c_str());
        double a = 100000.0;
        //for(double a = 500.0; a<=100000.0; a = a + dela)
        {
            //if(a==500.0) init_guess1 = init_guess1;
            //else init_guess1 = pole;

            if(a>=1000) dela = 1000;
            if(a>=10000) dela = 5000;
            double dels = 0.000001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m)) - 1.0e-9;
            double prevtempval=INT_MAX;
            double prevs=0.0;
    
            
            //init_guess1 = 8.99990000001;
            //init_guess1 = (double) real(Grightbranchpoint(a,m)) + 1.0e-8;
            //init_guess2 = (double) real(phibthreshold(a,m)) - 1.0e-9; 
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            //init_guess1 = pole;
        
            
            cout<<"KDF Poles test"<<endl;
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            //Mphib_secondsheet_denom_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            KDF_F3_pole_searching_bissection_with_N_using_weights(KDF,init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);

            double am = a;
            double sb = pole;
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<(double) sqrt(real(phibthreshold(am,m))) - sqrt(sb) <<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<(double) sqrt(real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            //init_guess2 = pole;
        }

        fout.close();
        cout<<"==============================="<<endl;
        
        
        //spole = polevec[0];
    }

    
}

void poletest_BS2_vs_N_using_weights_KDF()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    comp KDF = {100.0,0.0};
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.99966;
    double init_guess2 = 8.99994;//(double) real(phibthreshold(a,m));
    
    int N = 500;
    //for(int N=3000; N<=6000; N=N+500)
    {
        vector<double> polevec;
        string filename = "bs2_KDF_N_" + to_string(N) + "_1.dat";
        vector<double> avec;
        double dela = 100.0;
        fout.open(filename.c_str());
        //double a = 16.0;
        for(double a = 500.0; a<=100000.0; a = a + dela)
        {
            if(a==500.0) init_guess1 = init_guess1;
            else init_guess1 = (pole + init_guess1)/2.0;

            if(a>=1000) dela = 1000;
            if(a>=10000) dela = 5000;
            double dels = 0.000001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m)) - 1.0e-9;
            double prevtempval=INT_MAX;
            double prevs=0.0;
    
            
            //init_guess1 = 8.99990000001;
            //init_guess1 = (double) real(Grightbranchpoint(a,m)) + 1.0e-8;
            //init_guess2 = (double) real(phibthreshold(a,m)) - 1.0e-9; 
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            //init_guess1 = pole;
        
            
            cout<<"KDF Poles test"<<endl;
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            //Mphib_secondsheet_denom_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            KDF_F3_pole_searching_bissection_with_N_using_weights(KDF,init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);

            double am = a;
            double sb = pole;
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<(double) sqrt(real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<(double) sqrt(real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            //init_guess2 = pole;
        }

        fout.close();
        cout<<"==============================="<<endl;
        
        
        //spole = polevec[0];
    }

    
}

void poletest_BS2_vs_N_singlea_using_weights_KDF()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    comp KDF = {100.0,0.0};
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.99966;//8.999756491345767;//8.999780000771434;//8.9997;
    double init_guess2 = 8.99999;//(double) real(phibthreshold(a,m));
    
    int N = 500;
    //for(int N=3000; N<=6000; N=N+500)
    {
        vector<double> polevec;
        string filename = "bs2_KDF_N_" + to_string(N) + "_tmp.dat";
        vector<double> avec;
        double dela = 100.0;
        fout.open(filename.c_str());
        double a = 10000.0;
        //for(double a = 500.0; a<=100000.0; a = a + dela)
        {
            //if(a==500.0) init_guess1 = init_guess1;
            //else init_guess1 = pole;

            if(a>=1000) dela = 1000;
            if(a>=10000) dela = 5000;
            double dels = 0.000001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m)) - 1.0e-9;
            double prevtempval=INT_MAX;
            double prevs=0.0;
    
            
            //init_guess1 = 8.99990000001;
            //init_guess1 = (double) real(Grightbranchpoint(a,m)) + 1.0e-8;
            //init_guess2 = (double) real(phibthreshold(a,m)) - 1.0e-9; 
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            //init_guess1 = pole;
        
            
            cout<<"KDF Poles test"<<endl;
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            //Mphib_secondsheet_denom_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            KDF_F3_pole_searching_bissection_with_N_using_weights(KDF,init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);

            double am = a;
            double sb = pole;
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<(double) sqrt(real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<(double) sqrt(real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            //init_guess2 = pole;
        }

        fout.close();
        cout<<"==============================="<<endl;
        
        
        //spole = polevec[0];
    }

    
}


void poletest_BS3_vs_N_using_weights_KDF()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    comp KDF = {100.0,0.0};
    double temp_pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.999975;//8.999975;
    double init_guess2 = 0.0;//(double) real(phibthreshold(a,m));
    
    int N = 500;
    //for(int N=3000; N<=6000; N=N+500)
    {
        vector<double> polevec;
        string filename = "bs3_KDF_N_" + to_string(N) + "_2.dat";
        vector<double> avec;
        double dela = 100.0;
        fout.open(filename.c_str());
        //double a = 16.0;
        for(double a = 1000.0; a<=10000.0; a = a + dela)
        {   
            init_guess1 = 8.999975;
            //if(a==500.0) init_guess1 = init_guess1;
            //else init_guess1 = (pole + init_guess1)/2.0;
            double dels = 0.000000235;;
            if(a>=1000){
                dela = 1000;
                dels = 0.000000235;
            } 
            if(a>=10000){
                dela = 5000;
                dels = 0.0000000235;
            }

            
            //double dels = 0.0000002;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m)) - 1.0e-9;
            double prevtempval=INT_MAX;
            double prevs=0.0;
    
            init_guess2 = (double) real(phibthreshold(a,m));

            comp check_func_at_initguess1 = KDF_func_with_N_using_weights(KDF,init_guess1,a,N);
            if((double) real(check_func_at_initguess1)>0.0)
            {
                for(int i=0; ;++i)
                {
                    double check_s = init_guess1 + i*dels; 
                    if(check_s>9.0)
                    {
                        cout<<"init guess not found"<<endl;
                        break;
                    }
                    comp check_func = KDF_func_with_N_using_weights(KDF,check_s,a,N);
                    cout<<setprecision(16);
                    cout<<"guessing init_guess2 with s = "<<check_s<<'\t'<<" func = "<<real(check_func)<<endl;
                
                    //if(sign_checker(init_guess1)!=sign_checker((double) real(check_func)))
                    if((double)real(check_func)<0.0)
                    {
                        cout<<"init_guess set at s = "<<check_s<<endl;
                        init_guess2 = check_s;
                        break;
                    }
                }
            }
            else 
            {
                double temp_s1;
                for(int i=0; ;++i)
                {

                    double check_s = init_guess1 + i*dels; 
                    if(check_s>9.0)
                    {
                        cout<<"init guess not found"<<endl;
                        break;
                    }
                    comp check_func = KDF_func_with_N_using_weights(KDF,check_s,a,N);
                    cout<<setprecision(16);
                    cout<<"guessing init_guess1 with s = "<<check_s<<'\t'<<" func = "<<real(check_func)<<endl;
                
                    //if(sign_checker(init_guess1)!=sign_checker((double) real(check_func)))
                    if((double)real(check_func)>0.0)
                    {
                        cout<<"init_guess set at s = "<<check_s<<endl;
                        init_guess1 = check_s;
                        temp_s1 = init_guess1;
                        break;
                    }
                }
                for(int i=0; ;++i)
                {
                    double check_s = temp_s1 + i*dels; 
                    if(check_s>9.0)
                    {
                        cout<<"init guess not found"<<endl;
                        break;
                    }
                    comp check_func = KDF_func_with_N_using_weights(KDF,check_s,a,N);
                    cout<<setprecision(16);
                    cout<<"guessing init_guess2 with s = "<<check_s<<'\t'<<" func = "<<real(check_func)<<endl;
                
                    //if(sign_checker(init_guess1)!=sign_checker((double) real(check_func)))
                    if((double)real(check_func)<0.0)
                    {
                        cout<<"init_guess set at s = "<<check_s<<endl;
                        init_guess2 = check_s;
                        //temp_s1 = init_guess2;
                        break;
                    }
                }
            }
            
            //init_guess1 = 8.99990000001;
            //init_guess1 = (double) real(Grightbranchpoint(a,m)) + 1.0e-8;
            //init_guess2 = (double) real(phibthreshold(a,m)) - 1.0e-9; 
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            //init_guess1 = pole;
            //temp_pole = pole;
            
            cout<<"KDF Poles test"<<endl;
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            //Mphib_secondsheet_denom_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            KDF_F3_pole_searching_bissection_with_N_using_weights(KDF,init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);

            double am = a;
            double sb = pole;
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<(double) sqrt(real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<(double) sqrt(real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            //init_guess2 = pole;
        }

        fout.close();
        cout<<"==============================="<<endl;
        
        
        //spole = polevec[0];
    }

    
}

void poletest_BS3_vs_N_singlea_using_weights_KDF()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    comp KDF = {100.0,0.0};
    double temp_pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.999975;
    double init_guess2 = 0.0;//(double) real(phibthreshold(a,m));
    
    int N = 1000;
    //for(int N=3000; N<=6000; N=N+500)
    {
        vector<double> polevec;
        string filename = "bs3_KDF_N_" + to_string(N) + "_1.dat";
        vector<double> avec;
        double dela = 100.0;
        fout.open(filename.c_str());
        double a = 500.0;
        //for(double a = 500.0; a<=100000.0; a = a + dela)
        {
            //if(a==500.0) init_guess1 = init_guess1;
            //else init_guess1 = (pole + init_guess1)/2.0;

            if(a>=1000) dela = 1000;
            if(a>=10000) dela = 5000;
            double dels = 0.000002;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m)) - 1.0e-9;
            double prevtempval=INT_MAX;
            double prevs=0.0;
    
            init_guess2 = (double) real(phibthreshold(a,m));

            for(int i=0; ;++i)
            {
                double check_s = init_guess2 - i*dels; 
                comp check_func = KDF_func_with_N_using_weights(KDF,check_s,a,N);
                if(real(check_func)>0.0)
                {
                    init_guess1 = check_s;
                    break;
                }
            }
            
            cout<<"KDF Poles test"<<endl;
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            //Mphib_secondsheet_denom_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            KDF_F3_pole_searching_bissection_with_N_using_weights(KDF,init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);

            double am = a;
            double sb = pole;
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<(double) sqrt(real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<(double) sqrt(real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            //init_guess2 = pole;
        }

        fout.close();
        cout<<"==============================="<<endl;
        
        
        //spole = polevec[0];
    }

    
}

void poletest_master_vs_N_using_weights_KDF()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    comp KDF = {100.0,0.0};
    double temp_pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.9999753;
    double init_guess2 = 0.0;//(double) real(phibthreshold(a,m));
    
    int N = 500;
    //for(int N=3000; N<=6000; N=N+500)
    {
        vector<double> polevec;
        string filename = "bs3_KDF_N_" + to_string(N) + ".dat";
        vector<double> avec;
        double dela = 100.0;
        fout.open(filename.c_str());
        //double a = 16.0;
        for(double a = 500.0; a<=100000.0; a = a + dela)
        {
            //if(a==500.0) init_guess1 = 8.9997;
            //else init_guess1 = temp_pole;

            if(a>=1000) dela = 1000;
            if(a>=10000) dela = 5000;
            double dels = 0.000001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m)) - 1.0e-9;
            double prevtempval=INT_MAX;
            double prevs=0.0;
    
            init_guess2 = (double) real(phibthreshold(a,m));
            //init_guess1 = 8.99990000001;
            //init_guess1 = (double) real(Grightbranchpoint(a,m)) + 1.0e-8;
            //init_guess2 = (double) real(phibthreshold(a,m)) - 1.0e-9; 
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            //init_guess1 = pole;
            //temp_pole = pole;
            
            cout<<"KDF Poles test"<<endl;
            //dSqqs2k2msq_pole_searching_bissection_with_N_using_determinant(init_guess1,init_guess2,a,pole,N);
            //Mphib_secondsheet_denom_pole_searching_bissection_with_N_using_weights(init_guess1,init_guess2,a,pole,N);
            KDF_F3_pole_searching_bissection_with_N_using_weights(KDF,init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);

            double am = a;
            double sb = pole;
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<(double) sqrt(real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<(double) sqrt(real(phibthreshold(am,m))) - sqrt(sb)<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            //init_guess2 = pole;
        }

        fout.close();
        cout<<"==============================="<<endl;
        
        
        //spole = polevec[0];
    }

    
}



//small fixes for anomalous BS3 poles, where
//pole matches phib threshold
void poletest_bs3_vs_N_fix()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.99998;
    double init_guess2 = 0.0;//(double) real(phibthreshold(a,m));
    
    
    for(int N=4000; N<=7000; N=N+500)
    {
        //int N = 4000;
        vector<double> polevec;
        string filename = "bs3rd_N_" + to_string(N) + "_00_fix.dat";
        vector<double> avec;
        double dela = 50.0;
        
        for(double a = 750.0; a<=6000.0; a = a + dela)
        {
            if(a>=1000) dela = 1000;
            if(a>=10000) dela = 2500;
            double dels = 0.000001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m));
            double prevtempval=INT_MAX;
            double prevs=0.0;
    
            //(double) real(phibthreshold(a,m));
            for(int i=0;i<100;++i)
            {
                double s = set_s - i*dels;
                tempval = real(dSqqs2q2msq_func_with_N(s,a,N));
                cout<<"am = "<<setprecision(16)<<a<<'\t'<<"N = "<<N<<endl;
                cout<<"value of s = "<<setprecision(16)<<s<<endl;
                cout<<"value of dS = "<<setprecision(16)<<tempval<<endl;
                if(i==0)
                {
                    prevtempval = tempval;
                    prevs = s;
                    continue;
                } 
                if(tempval<0.0)
                {
                    init_guess2 = s;
                
                    cout<<"phibth = "<<setprecision(16)<<set_s<<'\t'
                        <<"init_guess2 = "<<setprecision(16)<<init_guess2<<endl;
                    break;
                }
                
                prevtempval = tempval;
                prevs = s;
                if(i==100-1)
                {
                    cout<<"init_guess2 not found"<<endl;
                    cout<<"tempval = "<<tempval<<endl;
                    cout<<"setting to current s = "<<s<<endl;
                    init_guess2 = s;
                }
            }
            //init_guess1 = 8.99990000001;
            //init_guess2 = (double) real(phibthreshold(a,m)); 
        
            
        
            dSqqs2k2msq_pole_searching_bissection_with_N(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
            double am = avec[i];
            double sb = polevec[i];
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
        }
        fout.close();
        cout<<"==============================="<<endl;
    }

    
}

void poletest_bs3_vs_singleN_fix()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.99998;
    double init_guess2 = 0.0;//(double) real(phibthreshold(a,m));
    
    
    //for(int N=4000; N<=7000; N=N+500)
    int N = 2000;
    {
        //int N = 4000;
        vector<double> polevec;
        string filename = "bs3rd_N_" + to_string(N) + ".dat";
        vector<double> avec;
        double dela = 5000.0;
        
        for(double a = 5000.0; a<=100000.0; a = a + dela)
        {
            //if(a>=1000) dela = 1000;
            //if(a>=10000) dela = 2500;
            double dels = 0.000001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m));
            double prevtempval=INT_MAX;
            double prevs=0.0;
    
            //(double) real(phibthreshold(a,m));
            for(int i=0;i<100;++i)
            {
                double s = set_s - i*dels;
                tempval = real(dSqqs2q2msq_func_with_N(s,a,N));
                cout<<"am = "<<setprecision(16)<<a<<'\t'<<"N = "<<N<<endl;
                cout<<"value of s = "<<setprecision(16)<<s<<endl;
                cout<<"value of dS = "<<setprecision(16)<<tempval<<endl;
                if(i==0)
                {
                    prevtempval = tempval;
                    prevs = s;
                    continue;
                } 
                if(tempval<0.0)
                {
                    init_guess2 = s;
                
                    cout<<"phibth = "<<setprecision(16)<<set_s<<'\t'
                        <<"init_guess2 = "<<setprecision(16)<<init_guess2<<endl;
                    break;
                }
                
                prevtempval = tempval;
                prevs = s;
                if(i==100-1)
                {
                    cout<<"init_guess2 not found"<<endl;
                    cout<<"tempval = "<<tempval<<endl;
                    cout<<"setting to current s = "<<s<<endl;
                    init_guess2 = s;
                }
            }
            //init_guess1 = 8.99990000001;
            //init_guess2 = (double) real(phibthreshold(a,m)); 
        
            
        
            dSqqs2k2msq_pole_searching_bissection_with_N(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            avec.push_back(a);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
            double am = avec[i];
            double sb = polevec[i];
            cout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
            fout<<setprecision(16)<<am<<'\t'<<setprecision(16)<<sb<<'\t'
                <<setprecision(16)<<sqrt(sb) - 3.0<<'\t'<<N<<'\t'<<setprecision(16)<<(double) real(phibthreshold(am,m))<<endl;
        }
        fout.close();
        cout<<"==============================="<<endl;
    }

    
}


void testing_Mphib_belowthreshold_vs_s3()
{
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.70;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/10.0;
    double points1 = 2000.0;
    double points2 = 500.0;
    double eps = 1.0e-7;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.01;
    //double box = 10.0;

    string filename="Mphib_momrep_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_epspow_" + to_string((int)abs(log10(eps))) 
                    + ".dat";
    ofstream fout;
    //fout.open(filename.c_str());

    int count = 0;

    double eta = 20.0;


    //for(double s=sinitial;s<=sfinal;s=s+dels)
    //{ 
        double s = 8.72;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp qval = pmom(s,sigb,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        
        //cout<<"s:"<<s<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;

        mom_vector_maker_4(qvec,s,kmin,kmax,a,m,eps,eps_for_m2k,points1,qvec_r);
    
        int size = qvec.size();
        cout<<"size = "<<size<<endl;
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec = -1.0*Gvec;

        Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);

        double relerror;

        

        //LinearSolver(Bmat,dsol,Gvec,relerror);
        cusolverComplex(Bmat,Gvec,dsol,size);

        comp result;

        interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        result = gsq*result;
        //result = rhophib(s,a,m)*result;

        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    //}

    //fout.close();

}

comp Mphib_func_with_N(   double s,
                                double a,
                                double N,
                                double eps  )
{
    //double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.72;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/100.0;
    double points1 = N;
    double points2 = 500.0;
    //double eps = 0.0;
    //double box = 10.0;

    //string filename="dSqqs2q2msq_momrep_a_" + to_string((int)a) 
    //                + "_N_" + to_string((int)points1) 
    //                + "_eps_" + to_string(eps) 
    //                + ".dat";
    //ofstream fout;
    //fout.open(filename.c_str());

    //int count = 0;

    //double eta = 20.0;


    //for(double s=sinitial;s<=sfinal;s=s+dels)
    //{ 
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp sigq = sigb;//2.0*m*m;
        comp qval = pmom(s,sigq,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        double delk = abs(kmax-kmin)/points1;
        
        //cout<<"s:"<<s<<endl;
        //cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        //cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;

        //mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
        double eps_for_m2k = 0.0;
        double qvec_r = 0.01;
        //mom_vector_maker_4(qvec,s,kmin,kmax,a,m,eps,points1);
        mom_vector_maker_4(qvec,s,kmin,kmax,a,m,eps,eps_for_m2k,points1,qvec_r);

        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec = -1.0*Gvec;

        Bmat_maker_momrep(Bmat,s,qvec,qvec,a,m,eps);

        double relerror;

        

        //LinearSolver(Bmat,dsol,Gvec,relerror);
        cusolverComplex(Bmat,Gvec,dsol,size);

        comp result;

        interpolator_ds_integraleq_momrep(dsol,qvec,s,qval,qval,a,m,eps,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;

        //result = gsq*result;
        result = result + GS_pk(s,qval,qval,m,eps);
        result = gsq*result;
        //result = rhophib(s,a,m)*result;

        //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        //fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        //count = count + 1;
        //cout<<"-------------------------"<<endl;
    //}

    //fout.close();

    return 1.0/real(result);

}

void Mphib_pole_searching_bissection_with_N(      double init_guess1,
                                                        double init_guess2,
                                                        double a,
                                                        double &pole,
                                                        double N,
                                                        double eps    )
{
    double apoint = init_guess1;
    double bpoint = init_guess2;
    int max_iteration = 100;
    double dels = 0.001;
    //double a = 16.0;
    //double s = apoint;
    
    double epsilon = 1.0e-10;
    double tolerance = 1.0e-15;
    double sa = apoint;
    double sb = bpoint;

    for(int i=0;i<max_iteration;++i)
    {
        
        double sc = (sa + sb)/2.0;
        double fc = (double) real(Mphib_func_with_N(sc,a,N,eps));
        double fa = (double) real(Mphib_func_with_N(sa,a,N,eps));
        //double fprime = (double) real(dSqqs2q2msq_func(s+dels,a) - dSqqs2q2msq_func(s,a))/dels;

        cout<<"a = "<<a<<endl;
        cout<<"iteration:"<<i+1<<'\t'<<" s:"<<sc<<'\t'<<" f:"<<fc<<endl;

        if(abs(fc)<epsilon || (sb-sa)/2.0 < tolerance )
        {
            cout<<"---------------------------"<<endl;
            cout<<"fc:"<<fc<<'\t'<<"eps:"<<epsilon<<endl;
            cout<<"b-a/2 : "<<(sb-sa)/2.0<<endl;
            pole = sc;
            cout<<"pole found at s = "<<setprecision(16)<<sc<<endl;
            cout<<"---------------------------"<<endl;
            cout<<endl;
            break;
        }
        if(fc>0.0 && fa>0.0) 
        {
            sa = sc;
        }
        else 
        {
            sb = sc;
        }

        if(i==max_iteration-1)
        {
            cout<<"method failed, max iteration reached"<<endl;
        }
    }
    //cout<<"pole found at s = "<<s<<endl;
    //pole = s;
}

void poletest_with_Mphib_bs1_vs_N()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    double eps = 1.0e-4;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.70;
    double init_guess2 = 8.99;//(double) real(phibthreshold(a,m));
    
    for(int N=500; N<=4000; N=N+500)
    {
        vector<double> polevec;
        string filename = "bs1st_N_" + to_string(N) + ".dat";
        
        double dela = 500.0;
        for(double a = 500.0; a<=3000.0; a = a + dela)
        {
            //if(a==500.0) dela = 500.0;
            //double dels = 0.00001;
            //double tempval = 0.0;
            //double set_s = (double) real(phibthreshold(a,m));
            //init_guess2 = (double) real(phibthreshold(a,m));
    
            //(double) real(phibthreshold(a,m));
            
        
        
            Mphib_pole_searching_bissection_with_N(init_guess1,init_guess2,a,pole,N,eps);
            polevec.push_back(pole);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 500.0 + i*500.0;
            double sb = polevec[i];
            cout<<setprecision(16)<<a1<<'\t'<<sb<<'\t'
                <<sqrt(sb) - 3.0<<'\t'<<N<<endl;
            fout<<setprecision(16)<<a1<<'\t'<<sb<<'\t'
                <<sqrt(sb) - 3.0<<'\t'<<N<<endl;
        }
        fout.close();
        cout<<"==============================="<<endl;
    }

    
}

void poletest_with_Mphib_bs2_vs_N()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    double eps = 1.0e-4;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.99;
    double init_guess2 = 0.0;//(double) real(phibthreshold(a,m));
    
    for(int N=500; N<=4000; N=N+500)
    {
        vector<double> polevec;
        string filename = "bs2nd_N" + to_string(N) + ".dat";
        
        for(double a = 500.0; a<=3000.0; a = a + 500.0)
        {
            double dels = 0.00001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m));
            init_guess2 = (double) real(phibthreshold(a,m));
    
            //(double) real(phibthreshold(a,m));
            /*for(int i=0;i<100;++i)
            {
                double s = set_s - i*dels;
                tempval = real(Mphib_func_with_N(s,a,N,eps));
                cout<<"value of dS = "<<tempval<<endl;
                if(tempval<0.0)
                {
                    init_guess2 = s;
                
                    cout<<"phibth = "<<set_s<<'\t'
                        <<"init_guess2 = "<<init_guess2<<endl;
                    break;
                }
                if(i==100-1)
                {
                    cout<<"init_guess2 not found"<<endl;
                    cout<<"tempval = "<<tempval<<endl;
                }
            }*/
        
        
            Mphib_pole_searching_bissection_with_N(init_guess1,init_guess2,a,pole,N,eps);
            polevec.push_back(pole);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 500.0 + i*500.0;
            double sb = polevec[i];
            cout<<setprecision(16)<<a1<<'\t'<<sb<<'\t'
                <<sqrt(sb) - 3.0<<'\t'<<N<<endl;
            fout<<setprecision(16)<<a1<<'\t'<<sb<<'\t'
                <<sqrt(sb) - 3.0<<'\t'<<N<<endl;
        }
        fout.close();
        cout<<"==============================="<<endl;
    }

    
}


void dSqqs2q2msq_belowthreshold_vs_s3_withmeshchoice( int someA, int someB, int N)
{
    double a = 500.0;
    double m = 1.0;

    //double s = 8.65;
    double spole = 8.99999999916955;//8.999999440008509;//8.999983333329475;
    double delspole = real(phibthreshold(a,m)) - spole;
    
    double sinitial = 8.999974;//spole - 0.0000001;//8.999000;//spole;// - delspole;//8.95;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.99;
    double spoints = 500.0;
    
    double dels = abs(sinitial-sfinal)/spoints;
    double points1 = N;
    double points2 = 500.0;
    double eps = 0.0;
    //double box = 10.0;
    int firstpart_for_nonunimesh = someA;
    int secondpart_for_nonunimesh = someB;

    string filename="non-uniform_dSqqs2q2msq_BS3_momrep_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_eps_" + to_string(eps) 
                    + "linvecmak_" + to_string(firstpart_for_nonunimesh) 
                    + to_string(secondpart_for_nonunimesh) 
                    + ".dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    double eta = 20.0;


    for(int i=0;i<=spoints;++i)
    //for(double s=sinitial;s<=sfinal;s=s+dels)
    { 
        double s = sinitial + (double)i*dels;
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        cout<<phibthreshold(a,m)<<endl;
        cout<<"am = "<<a<<endl;
        //cout<<spole<<endl;
        //cout<<spole-delspole<<endl;
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp sigq = 2.0*m*m;
        comp qval = pmom(s,sigq,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        double delk = abs(kmax-kmin)/points1;
        
        //cout<<"s:"<<s<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;

        if(firstpart_for_nonunimesh==0)
        mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
        else
        mom_vector_maker_linear_2(qvec,kmin,kmax,points1,1,firstpart_for_nonunimesh,secondpart_for_nonunimesh);
        //mom_vector_maker_1(qvec,s,kmin,kmax,a,m,eps,points1);
    
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec = -1.0*Gvec;

        Bmat_maker_momrep(Bmat,s,qvec,qvec,a,m,eps);

        comp detB = Bmat.determinant();

        double relerror;

        

        //LinearSolver_2(Bmat,dsol,Gvec,relerror);
        cusolverComplex(Bmat,Gvec,dsol,size);

        relerror = abs((Bmat*dsol - Gvec).norm()/Gvec.norm());
        comp result;

        interpolator_ds_integraleq_momrep(dsol,qvec,s,qval,qval,a,m,eps,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        //comp result = (s-spole)*result;
        //result = gsq*result;
        //result = rhophib(s,a,m)*result;

        cout<<"s = "<<setprecision(16)<<s<<'\t'<<"dels ="<<dels<<endl;
        cout<<"res = "<<result<<'\t'<<" run = "<<count + 1<<endl;
        cout<<"detB = "<<detB<<'\t'<<" error = "<<relerror<<endl;
        cout<<"A = "<<someA<<'\t'<<"B = "<<someB<<endl;
        cout<<"N = "<<N<<endl;
        fout<<setprecision(16)<<s<<'\t'
            <<real(result)<<'\t'
            <<imag(result)<<'\t'
            <<real(1.0/result)<<'\t'
            <<imag(1.0/result)<<'\t'
            <<real(detB)<<'\t'
            <<imag(detB)<<'\t'
            <<relerror<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}

void dSqqs2q2msq_belowthreshold_vs_s3()
{
    double a = 30000.0;
    double m = 1.0;

    //double s = 8.65;
    double spole = 8.999999483027445;//8.999999440008509;//8.999983333329475;
    double delspole = real(phibthreshold(a,m)) - spole;
    
    double sinitial = 8.99999948;//8.99998;//spole - 0.0000001;//8.999000;//spole;// - delspole;//8.95;//8.72;//real(phibthreshold(a,m));//
    double sfinal = 8.9999995;//real(phibthreshold(a,m));//8.99;
    double spoints = 200.0;
    
    double dels = abs(sinitial-sfinal)/spoints;
    double points1 = 500;//N;
    double points2 = 500.0;
    double eps = 0.0;
    //double box = 10.0;
    //int firstpart_for_nonunimesh = someA;
    //int secondpart_for_nonunimesh = someB;

    string filename="dSqqs2q2msq_BS3_momrep_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_eps_" + to_string(eps) 
                    + "_weights.dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    double eta = 20.0;


    for(int i=0;i<=spoints;++i)
    //for(double s=sinitial;s<=sfinal;s=s+dels)
    { 
        //if(i!=0) break;
        double s = sinitial + (double)i*dels;
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        cout<<"phibth = "<<phibthreshold(a,m)<<endl;
        cout<<"am = "<<a<<endl;
        //cout<<spole<<endl;
        //cout<<spole-delspole<<endl;
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp sigq = 2.0*m*m;
        comp qval = pmom(s,sigq,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        double delk = abs(kmax-kmin)/points1;
        
        //cout<<"s:"<<s<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;
        vector<comp> weights;

        //if(firstpart_for_nonunimesh==0)
        //mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
        line_maker_with_weights(qvec,weights,kmin,kmax,points1);
        //for(int i=0;i<qvec.size();++i)
        //{
            //cout<<"i="<<i<<'\t'<<"q="<<qvec[i]<<'\t'<<"w="<<weights[i]<<endl;
        //}
        //else
        //mom_vector_maker_linear_2(qvec,kmin,kmax,points1,1,firstpart_for_nonunimesh,secondpart_for_nonunimesh);
        //mom_vector_maker_1(qvec,s,kmin,kmax,a,m,eps,points1);
    
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec = -1.0*Gvec;

        //Bmat_maker_momrep(Bmat,s,qvec,qvec,a,m,eps);
        Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps);
        comp detB = Bmat.determinant();

        double relerror;

        

        //LinearSolver_2(Bmat,dsol,Gvec,relerror);
        cusolverComplex(Bmat,Gvec,dsol,size);

        relerror = abs((Bmat*dsol - Gvec).norm());
        comp result;

        interpolator_ds_integraleq_momrep(dsol,qvec,s,qval,qval,a,m,eps,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        //comp result = (s-spole)*result;
        //result = gsq*result;
        //result = rhophib(s,a,m)*result;

        cout<<"s = "<<setprecision(16)<<s<<'\t'<<"dels ="<<dels<<endl;
        cout<<"res = "<<result<<'\t'<<" run = "<<count + 1<<endl;
        cout<<"detB = "<<detB<<'\t'<<" error = "<<relerror<<endl;
        //cout<<"A = "<<someA<<'\t'<<"B = "<<someB<<endl;
        cout<<"N = "<<points1<<endl;
        fout<<setprecision(16)<<s<<'\t'
            <<real(result)<<'\t'
            <<imag(result)<<'\t'
            <<real(1.0/result)<<'\t'
            <<imag(1.0/result)<<'\t'
            <<real(detB)<<'\t'
            <<imag(detB)<<'\t'
            <<relerror<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}

void dSqqs2q2msq_belowthreshold_vs_s3_with_KDF()
{
    double a = 10000.0;
    double m = 1.0;
    comp KDF = {100.0,0.0};

    //double s = 8.65;
    double spole = 8.999975;//8.9995;//8.9246;//8.99999999916955;//8.999999440008509;//8.999983333329475;
    double delspole = real(phibthreshold(a,m)) - spole;
    
    double sinitial = 8.999999939;//8.9999996;//8.99995;//8.9995;//8.999822636602779;//8.9246;//8.99998;//8.99998;//spole - 0.0000001;//8.999000;//spole;// - delspole;//8.95;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.99;
    double spoints = 300.0;
    
    double dels = abs(sinitial-sfinal)/spoints;
    double points1 = 500;//N;
    double points2 = 500.0;
    double eps = 0.0;
    //double box = 10.0;
    //int firstpart_for_nonunimesh = someA;
    //int secondpart_for_nonunimesh = someB;

    string filename="dSqqs2q2msq_momrep_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_eps_" + to_string(eps) 
                    + "_with_KDF.dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    double eta = 20.0;


    for(int i=0;i<=spoints;++i)
    //for(double s=sinitial;s<=sfinal;s=s+dels)
    { 
        //if(i!=0) break;
        double s = sinitial + (double)i*dels;
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        cout<<"phibth = "<<phibthreshold(a,m)<<endl;
        cout<<"am = "<<a<<endl;
        //cout<<spole<<endl;
        //cout<<spole-delspole<<endl;
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp sigq = 2.0*m*m;
        comp qval = pmom(s,sigq,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        double delk = abs(kmax-kmin)/points1;
        
        //cout<<"s:"<<s<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;
        vector<comp> weights;

        //if(firstpart_for_nonunimesh==0)
        //mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
        line_maker_with_weights(qvec,weights,kmin,kmax,points1);
        //for(int i=0;i<qvec.size();++i)
        //{
            //cout<<"i="<<i<<'\t'<<"q="<<qvec[i]<<'\t'<<"w="<<weights[i]<<endl;
        //}
        //else
        //mom_vector_maker_linear_2(qvec,kmin,kmax,points1,1,firstpart_for_nonunimesh,secondpart_for_nonunimesh);
        //mom_vector_maker_1(qvec,s,kmin,kmax,a,m,eps,points1);
    
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::MatrixXcd Gmat(size,size);
        Eigen::MatrixXcd dsol(size,size);

        //Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gmat_maker_momrep(Gmat,s,qvec,qvec,m,eps);
        //Gvec = -1.0*Gvec;
        Gmat = -1.0*Gmat;

        //Bmat_maker_momrep(Bmat,s,qvec,qvec,a,m,eps);
        Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps);
        comp detB = Bmat.determinant();

        double relerror;

        

        LinearSolver_3(Bmat,dsol,Gmat,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

        //relerror = abs((Bmat*dsol - Gvec).norm());
        comp result;

        //interpolator_ds_integraleq_momrep(dsol,qvec,s,qval,qval,a,m,eps,result);

        F3infvol(dsol,qvec,weights,qvec,weights,s,a,m,eps,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        //comp result = (s-spole)*result;
        //result = gsq*result;
        //result = rhophib(s,a,m)*result;
        comp denom = result + 1.0/KDF;

        cout<<"s = "<<setprecision(16)<<s<<'\t'<<"dels ="<<dels<<endl;
        cout<<"res = "<<result<<'\t'<<" run = "<<count + 1<<endl;
        cout<<"detB = "<<detB<<'\t'<<" error = "<<relerror<<endl;
        cout<<"K^-1 + F = "<<denom<<endl;
        //cout<<"A = "<<someA<<'\t'<<"B = "<<someB<<endl;
        cout<<"N = "<<points1<<endl;
        fout<<setprecision(16)<<s<<'\t'
            <<real(result)<<'\t'
            <<imag(result)<<'\t'
            <<real(1.0/result)<<'\t'
            <<imag(1.0/result)<<'\t'
            <<real(detB)<<'\t'
            <<imag(detB)<<'\t'
            <<real(denom)<<'\t'
            <<imag(denom)<<'\t'
            <<relerror<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}



void dSqqs2q2msq_belowthreshold_vs_s3_sigk( double am, double N11, double spole, double delspole, int BSNum )
{
    double a = am;//100000.0;
    double m = 1.0;

    //double s = 8.65;
    //double spole = 8.895133263228459;//8.782848835105487;
    //double delspole = 0.0000001;//real(phibthreshold(a,m)) - spole;
    
    double delspoleBS1 = 0.001/5.0;
    double delspoleBS2 = 0.00001/5.0;
    double delspoleBS3 = 0.0000001/5.0;

    double sinitial = spole - delspole;//BS1;//8.782;//spole - 0.01;//8.75;//spole;// - delspole;//8.95;//8.72;//real(phibthreshold(a,m));//
    double sfinal = spole + delspole;//BS1;//8.788;//spole + 0.01;//real(phibthreshold(a,m));//8.99;
    
    double dels = abs(sinitial-sfinal)/10.0;
    double points1 = N11; 
    double points2 = 500.0;
    double eps = 0.0;
    //double box = 10.0;
    int sigcount = 0;
    int sigkpoints = 200.0;

    double sigkinitial = 0.001*m*m;
    double sigkfinal = 3.85*m*m;
    double delsigk = abs(sigkfinal - sigkinitial)/sigkpoints;

    double kinitial = 0.0001;
    double kfinal = 1.2;
    double kpoints = 200.0;
    double delksim = abs(kinitial - kfinal)/kpoints;

    //for(int i=0;i<sigkpoints+1;++i)
    //for(int i=0;i<kpoints+1;++i)
    {
        //double sigk = sigkinitial + i*delsigk;
        double k = 0.2;//kinitial + i*delksim;
        
        string filename="vertexfactor_momrep_a_" + to_string((int)a) 
                        + "_N_" + to_string((int)points1) 
                        + "_eps_" + to_string((int)eps)
                        + "_sigkcount_" + to_string(sigcount) 
                        + "_BS_" + to_string(BSNum)
                        + ".dat";
        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        for(double s=sinitial;s<=sfinal;s=s+dels)
        { 
            //double s = 8.78;
            //if(s>=sfinal) dels = 0.05/100.0;
            //eps = eps_above_threshold(eta,s,a,m,points);
            //cout<<"BS1 test"<<endl;
            cout<<"phibth = "<<phibthreshold(a,m)<<endl;
            //cout<<spole<<endl;
            //cout<<spole-delspole<<endl;
            double sigk = (double) real(sigma_p(s,k,m));
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = sigmab(a,m);
            comp sigq = sigk;
            comp qval = k;//pmom(s,sigq,m);
            comp qpole = k;//pmom(spole,sigq,m);

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
            double delk = abs(kmax-kmin)/points1;
        
            //cout<<"s:"<<s<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;

            mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
            //mom_vector_maker_linear_2(qvec,kmin,kmax,points1,1);
            //mom_vector_maker_1(qvec,s,kmin,kmax,a,m,eps,points1);
    
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            Gvec = -1.0*Gvec;

            Bmat_maker_momrep(Bmat,s,qvec,qvec,a,m,eps);

            double relerror;

        

            //LinearSolver_2(Bmat,dsol,Gvec,relerror);
            cusolverComplex(Bmat,Gvec,dsol,size);

            comp result;

            interpolator_ds_integraleq_momrep(dsol,qvec,s,qval,qval,a,m,eps,result);

            comp gfunc = gfuncConst(a,0.0,m);
            cout<<"gfunc = "<<gfunc<<endl;
            comp gsq = gfunc*gfunc;
            //comp result = (s-spole)*result;
            //result = gsq*result;
            //result = rhophib(s,a,m)*result;
            comp m2kres = M2kfunc(a,sigk,m,eps);
            cout<<"am = "<<a<<endl;
            cout<<"sigk = "<<sigk<<endl;
            cout<<"N = "<<qvec.size()<<endl;
            cout<<"M2kres = "<<m2kres<<endl;
            cout<<"M2kressq = "<<m2kres*m2kres<<endl;
            
            cout<<"s:"<<setprecision(16)<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            cout<<"s-pole = "<<s - spole<<'\t'<<" -(s-spole)*res = "<<-(s-spole)*result<<endl;
            fout<<setprecision(16)<<s<<'\t'
                <<sigk<<'\t'
                //<<real(qpole)<<'\t'
                <<real(k)<<'\t'
                <<imag(qpole)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(1.0/result)<<'\t'
                <<imag(1.0/result)<<'\t'
                <<real(m2kres)<<'\t'
                <<imag(m2kres)<<'\t'
                <<real(m2kres*m2kres)<<'\t'
                <<imag(m2kres*m2kres)<<endl;
            count = count + 1;
            cout<<"-------------------------"<<endl;
        }
        sigcount = sigcount + 1;

        fout.close();
    }

}


void dSqqs2q2msq_belowthreshold_vs_s3_sigk_with_weights( double k, double am, double N11, double spole, double delspole, int BSNum )
{
    double a = am;//100000.0;
    double m = 1.0;
    //double k = 0.2;

    //double s = 8.65;
    //double spole = 8.895133263228459;//8.782848835105487;
    //double delspole = 0.0000001;//real(phibthreshold(a,m)) - spole;
    
    double delspoleBS1 = 0.001/5.0;
    double delspoleBS2 = 0.00001/5.0;
    double delspoleBS3 = 0.0000001/5.0;

    double sinitial = spole - delspole;//BS1;//8.782;//spole - 0.01;//8.75;//spole;// - delspole;//8.95;//8.72;//real(phibthreshold(a,m));//
    double sfinal = spole + delspole;//BS1;//8.788;//spole + 0.01;//real(phibthreshold(a,m));//8.99;
    
    double dels = abs(sinitial-sfinal)/10.0;
    double points1 = N11; 
    double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = eps;
    //double box = 10.0;
    int sigcount = 0;
    int sigkpoints = 200.0;

    double sigkinitial = 0.001*m*m;
    double sigkfinal = 3.85*m*m;
    double delsigk = abs(sigkfinal - sigkinitial)/sigkpoints;

    double kinitial = 0.0001;
    double kfinal = 1.2;
    double kpoints = 200.0;
    double delksim = abs(kinitial - kfinal)/kpoints;

    //for(int i=0;i<sigkpoints+1;++i)
    //for(int i=0;i<kpoints+1;++i)
    {
        //double sigk = sigkinitial + i*delsigk;
        //kinitial + i*delksim;
        
        string filename="vertexfactor_momrep_a_" + to_string((int)a) 
                        + "_N_" + to_string((int)points1) 
                        + "_eps_" + to_string((int)eps)
                        + "_sigkcount_" + to_string(sigcount) 
                        + "_BS_" + to_string(BSNum)
                        + ".dat";
        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        for(double s=sinitial;s<=sfinal;s=s+dels)
        { 
            //double s = 8.78;
            //if(s>=sfinal) dels = 0.05/100.0;
            //eps = eps_above_threshold(eta,s,a,m,points);
            //cout<<"BS1 test"<<endl;
            cout<<"phibth = "<<phibthreshold(a,m)<<endl;
            //cout<<spole<<endl;
            //cout<<spole-delspole<<endl;
            double sigk = (double) real(sigma_p(s,k,m));
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = sigmab(a,m);
            comp sigq = sigk;
            comp qval = k;//pmom(s,sigq,m);
            comp qpole = k;//pmom(spole,sigq,m);

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
            double delk = abs(kmax-kmin)/points1;
        
            //cout<<"s:"<<s<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            //mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
            line_maker_with_weights(qvec,weights,kmin,kmax,points1);
            //mom_vector_maker_linear_2(qvec,kmin,kmax,points1,1);
            //mom_vector_maker_1(qvec,s,kmin,kmax,a,m,eps,points1);
    
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            Gvec = -1.0*Gvec;

            //Bmat_maker_momrep(Bmat,s,qvec,qvec,a,m,eps);
            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);


            double relerror;

        

            //LinearSolver_2(Bmat,dsol,Gvec,relerror);
            cusolverComplex(Bmat,Gvec,dsol,size);

            comp result;

            //interpolator_ds_integraleq_momrep(dsol,qvec,s,qval,qval,a,m,eps,result);
            interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            
            comp gfunc = gfuncConst(a,0.0,m);
            cout<<"gfunc = "<<gfunc<<endl;
            comp gsq = gfunc*gfunc;
            //comp result = (s-spole)*result;
            //result = gsq*result;
            //result = rhophib(s,a,m)*result;
            comp m2kres = M2kfunc(a,sigk,m,eps);
            cout<<"am = "<<a<<endl;
            cout<<"sigk = "<<sigk<<endl;
            cout<<"N = "<<qvec.size()<<endl;
            cout<<"M2kres = "<<m2kres<<endl;
            cout<<"M2kressq = "<<m2kres*m2kres<<endl;
            
            cout<<"s:"<<setprecision(16)<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            cout<<"s-pole = "<<s - spole<<'\t'<<" -(s-spole)*res = "<<-(s-spole)*result<<endl;
            fout<<setprecision(16)<<s<<'\t'
                <<sigk<<'\t'
                //<<real(qpole)<<'\t'
                <<real(k)<<'\t'
                <<imag(qpole)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(1.0/result)<<'\t'
                <<imag(1.0/result)<<'\t'
                <<real(m2kres)<<'\t'
                <<imag(m2kres)<<'\t'
                <<real(m2kres*m2kres)<<'\t'
                <<imag(m2kres*m2kres)<<endl;
            count = count + 1;
            cout<<"-------------------------"<<endl;
        }
        sigcount = sigcount + 1;

        fout.close();
    }

}

void dSqqs2q2msq_belowthreshold_vs_s3_k_vertexfactor_with_weights( double am, double N11, double spole, double delspole, int BSNum )
{
    double a = am;//100000.0;
    double m = 1.0;
    //double k = 0.2;

    //double s = 8.65;
    //double spole = 8.895133263228459;//8.782848835105487;
    //double delspole = 0.0000001;//real(phibthreshold(a,m)) - spole;
    
    double delspoleBS1 = 0.001/5.0;
    double delspoleBS2 = 0.00001/5.0;
    double delspoleBS3 = 0.0000001/5.0;

    double sinitial = spole - delspole;//BS1;//8.782;//spole - 0.01;//8.75;//spole;// - delspole;//8.95;//8.72;//real(phibthreshold(a,m));//
    double sfinal = spole + delspole;//BS1;//8.788;//spole + 0.01;//real(phibthreshold(a,m));//8.99;
    
    double dels = abs(sinitial-sfinal)/10.0;
    double points1 = N11; 
    double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = eps;
    //double box = 10.0;
    int sigcount = 0;
    int sigkpoints = 200.0;

    double sigkinitial = 0.001*m*m;
    double sigkfinal = 3.85*m*m;
    double delsigk = abs(sigkfinal - sigkinitial)/sigkpoints;

    double kinitial = 0.00001;//0.095;//0.00001;
    double kfinal = real(pmom(spole,0.0,m));//0.101;//real(pmom(spole,0.0,m));
    double kpoints = 2000.0;
    double delksim = abs(kinitial - kfinal)/kpoints;

    //for(int i=0;i<sigkpoints+1;++i)
    for(int i=0;i<kpoints+1;++i)
    {
        //double sigk = sigkinitial + i*delsigk;
        double k = kinitial + i*delksim;
        
        string filename="vertexfactor_momrep_a_" + to_string((int)a) 
                        + "_N_" + to_string((int)points1) 
                        + "_eps_" + to_string((int)eps)
                        + "_sigkcount_" + to_string(sigcount) 
                        + "_BS_" + to_string(BSNum)
                        + ".dat";
        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        for(double s=sinitial;s<=sfinal;s=s+dels)
        { 
            //double s = 8.78;
            //if(s>=sfinal) dels = 0.05/100.0;
            //eps = eps_above_threshold(eta,s,a,m,points);
            //cout<<"BS1 test"<<endl;
            cout<<"phibth = "<<phibthreshold(a,m)<<endl;
            //cout<<spole<<endl;
            //cout<<spole-delspole<<endl;
            double sigk = (double) real(sigma_p(s,k,m));
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = sigmab(a,m);
            comp sigq = sigk;
            comp qval = k;//pmom(s,sigq,m);
            comp qpole = k;//pmom(spole,sigq,m);

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
            double delk = abs(kmax-kmin)/points1;
        
            //cout<<"s:"<<s<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            //mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
            line_maker_with_weights(qvec,weights,kmin,kmax,points1);
            //mom_vector_maker_linear_2(qvec,kmin,kmax,points1,1);
            //mom_vector_maker_1(qvec,s,kmin,kmax,a,m,eps,points1);
    
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            Gvec = -1.0*Gvec;

            //Bmat_maker_momrep(Bmat,s,qvec,qvec,a,m,eps);
            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);


            double relerror;

        

            LinearSolver_2(Bmat,dsol,Gvec,relerror);
            //cusolverComplex(Bmat,Gvec,dsol,size);

            comp result;

            //interpolator_ds_integraleq_momrep(dsol,qvec,s,qval,qval,a,m,eps,result);
            interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            
            comp gfunc = gfuncConst(a,0.0,m);
            cout<<"gfunc = "<<gfunc<<endl;
            comp gsq = gfunc*gfunc;
            //comp result = (s-spole)*result;
            //result = gsq*result;
            //result = rhophib(s,a,m)*result;
            comp m2kres = M2kfunc(a,sigk,m,eps);
            cout<<"am = "<<a<<endl;
            cout<<"sigk = "<<sigk<<endl;
            cout<<"kinitial = "<<kinitial<<endl;
            cout<<"kfinal = "<<kfinal<<endl;
            cout<<"k = "<<k<<endl;
            cout<<"N = "<<qvec.size()<<endl;
            cout<<"M2kres = "<<m2kres<<endl;
            cout<<"M2kressq = "<<m2kres*m2kres<<endl;
            
            cout<<"s:"<<setprecision(16)<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            cout<<"s-pole = "<<s - spole<<'\t'<<" -(s-spole)*res = "<<-(s-spole)*result<<endl;
            fout<<setprecision(16)<<s<<'\t'
                <<sigk<<'\t'
                //<<real(qpole)<<'\t'
                <<real(k)<<'\t'
                <<imag(qpole)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(1.0/result)<<'\t'
                <<imag(1.0/result)<<'\t'
                <<real(m2kres)<<'\t'
                <<imag(m2kres)<<'\t'
                <<real(m2kres*m2kres)<<'\t'
                <<imag(m2kres*m2kres)<<endl;
            count = count + 1;
            cout<<"-------------------------"<<endl;
        }
        sigcount = sigcount + 1;

        fout.close();
    }

}

void mphib_residue_belowthreshold_vs_s3_k_vertexfactor_with_weights( double am, double N11, double spole, double delspole, int BSNum )
{
    time_t time_start, time_end;
    double a = am;//100000.0;
    double m = 1.0;
    //double k = 0.2;

    //double s = 8.65;
    //double spole = 8.895133263228459;//8.782848835105487;
    //double delspole = 0.0000001;//real(phibthreshold(a,m)) - spole;
    
    double delspoleBS1 = 0.001/5.0;
    double delspoleBS2 = 0.00001/5.0;
    double delspoleBS3 = 0.0000001/5.0;

    double sinitial = spole - delspole;//BS1;//8.782;//spole - 0.01;//8.75;//spole;// - delspole;//8.95;//8.72;//real(phibthreshold(a,m));//
    double sfinal = spole + delspole;//BS1;//8.788;//spole + 0.01;//real(phibthreshold(a,m));//8.99;
    
    double dels = abs(sinitial-sfinal)/10.0;
    double points1 = N11; 
    double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = eps;
    //double box = 10.0;
    int sigcount = 0;
    int sigkpoints = 200.0;

    double sigkinitial = 0.001*m*m;
    double sigkfinal = 3.85*m*m;
    double delsigk = abs(sigkfinal - sigkinitial)/sigkpoints;

    double kinitial = 0.00001;//0.095;//0.00001;
    double kfinal = real(pmom(spole,0.0,m));//0.101;//real(pmom(spole,0.0,m));
    double kpoints = 2000.0;
    double delksim = abs(kinitial - kfinal)/kpoints;

    //for(int i=0;i<sigkpoints+1;++i)
    //for(int i=0;i<kpoints+1;++i)
    {
        //double sigk = sigkinitial + i*delsigk;
        //double k = kinitial + i*delksim;
        
        
        string filename="vertexfactor_momrep_a_" + to_string((int)a) 
                        + "_N_" + to_string((int)points1) 
                        + "_eps_" + to_string((int)eps)
                        + "_sigkcount_" + to_string(sigcount) 
                        + "_BS_" + to_string(BSNum)
                        + ".dat";
        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        for(double s1=sinitial;s1<=sfinal;s1=s1+dels)
        { 
            double simag = 0.0;//1.0e-5;
            comp s = s1 + ii*simag;
            comp sigb = sigmab(a,m);
            comp k = pmom(s,sigb,m);
            //double s = 8.78;
            //if(s>=sfinal) dels = 0.05/100.0;
            //eps = eps_above_threshold(eta,s,a,m,points);
            //cout<<"BS1 test"<<endl;
            cout<<"phibth = "<<phibthreshold(a,m)<<endl;
            //cout<<spole<<endl;
            //cout<<spole-delspole<<endl;
            double sigk = (double) real(sigma_p(s,k,m));
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            //comp sigb = sigmab(a,m);
            comp sigq = sigk;
            comp qval = k;//pmom(s,sigq,m);
            comp qpole = k;//pmom(spole,sigq,m);

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
            double delk = abs(kmax-kmin)/points1;
        
            //cout<<"s:"<<s<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            
            int switch_for_gvec_fixer = 0; //0 means we have used contour and turn on gvec fixer
                                           //1 means we have used straight line and turn off gvec fixer
            

            //mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
            //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
            //mom_vector_maker_linear_2(qvec,kmin,kmax,points1,1);
            //mom_vector_maker_1(qvec,s,kmin,kmax,a,m,eps,points1);
            mom_vector_maker_seba_imspos_2(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            if(imag(s)<0.0)
            {
                Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            }
            else if(imag(s)>=0.0)
            {
                Gvec_maker_momrep_withtags_1(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            }

            if(switch_for_gvec_fixer==0)
                cout<<"gvec fixer switch = "<<switch_for_gvec_fixer<<",deform contour"<<endl;
            else 
                cout<<"gvec fixer switch = "<<switch_for_gvec_fixer<<",take straight line"<<endl;

            if(switch_for_gvec_fixer==0)
            {   
                //cout<<"here the problem lies 1"<<endl;
                Gvec_fixer_3(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
                //Gvec_fixer_4(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);//Gvec fixer has a contour maker inside of it
                
                //cout<<"here the problem lies 2"<<endl;
                tag_fixer_for_gvec_fixer3(Gvec,qvec,s,qval,m,eps,tag1); //this works for gvec_fixer_4 as well

                //cout<<"here the problem lies 3"<<endl;
               
            }
            Eigen::VectorXcd interpolater_Gvec = Gvec;
            Gvec = -1.0*Gvec;

            //Bmat_maker_momrep(Bmat,s,qvec,qvec,a,m,eps);
            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);


            double relerror;

        

            time(&time_start);
            LinearSolver_2(Bmat,dsol,Gvec,relerror);

            //cusolverComplex(Bmat,Gvec,dsol,size);
            time(&time_end);

            double time_taken = double(time_end - time_start);
            cout<<"Time taken to solve : "<<fixed 
                <<time_taken<<setprecision(5);
            cout<<" sec"<<endl;


            comp result;

            if(imag(s)<0.0)
            {
                //interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
                interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            }
            else if(imag(s)>=0.0)
            {
                //interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
                //interpolator_ds_integraleq_momrep_2eps_withtags_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
                interpolator_ds_integraleq_momrep_2eps_withtags_with_weights_usingGvec(dsol,qvec,weights,interpolater_Gvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
            
            }

            comp gfunc = gfuncConst(a,0.0,m);
            cout<<"gfunc = "<<gfunc<<endl;
            cout<<"gfunc^2 = "<<gfunc*gfunc<<endl;
            comp gsq = gfunc*gfunc;
            //comp result = (s-spole)*result;
            //result = gsq*result;
            //result = rhophib(s,a,m)*result;
            comp m2kres = M2kfunc(a,sigk,m,eps);
            cout<<"am = "<<a<<endl;
            cout<<"sigk = "<<sigk<<endl;
            cout<<"kinitial = "<<kinitial<<endl;
            cout<<"kfinal = "<<kfinal<<endl;
            cout<<"k = "<<k<<endl;
            cout<<"N = "<<qvec.size()<<endl;
            cout<<"M2kres = "<<m2kres<<endl;
            cout<<"M2kressq = "<<m2kres*m2kres<<endl;
            
            cout<<"s:"<<setprecision(16)<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            cout<<"s-pole = "<<s - spole<<'\t'<<" -(s-spole)*res = "<<-(s-spole)*result<<endl;
            fout<<setprecision(16)<<real(s)<<'\t'
                <<sigk<<'\t'
                //<<real(qpole)<<'\t'
                <<real(k)<<'\t'
                <<imag(qpole)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(1.0/result)<<'\t'
                <<imag(1.0/result)<<'\t'
                <<real(gfunc)<<'\t'
                <<imag(gfunc)<<'\t'
                <<real(gsq)<<'\t'
                <<imag(gsq)<<endl;
            count = count + 1;
            cout<<"-------------------------"<<endl;
        }
        sigcount = sigcount + 1;

        fout.close();
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



void sebatest_Mphib_belowthreshold_vs_s3_fixed_s3imag_contour47(double a, double pointsize)
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    //double a = 10000.0;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    double s3real = 8.70;
    double s3imag = 1.0e-5;//0.05;
    double sinitial = 0;//8.99;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;

    double delspoints = 100.0;
    double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
    double points1 = pointsize;//3000.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.01;
    //double box = 10.0;

    //string filename="fixeds3imag_Mphib_a_" + to_string((int)a) 
    //                + "_N_" + to_string((int)points1) 
    //                + "_contour47.dat";
    
    string filename;
    if(s3imag>=0.0)
    {
        //filename = "shifttest2_Mphib_s=" + to_string(s3real) + "+" + to_string(s3imag) + "i_N=" + to_string((int)points1) + ".dat";
        //filename = "Mphib_s=" + to_string(s3real) + "+" + to_string(s3imag) + "i.dat";
        //filename = "dS1_BS2_a=" + to_string((int)a) +  "_N=" + to_string((int)points1) + ".dat";
        filename = "shifttest2_kernel_4_smoothcutoff.dat";
        //filename = "tmp3.dat";
    }
    else 
    {
        //filename = "shifttest2_Mphib_s=" + to_string(s3real) + to_string(s3imag) + "i_N=" + to_string((int)points1) + ".dat";
        //filename = "Mphib_s=" + to_string(s3real) + to_string(s3imag) + "i.dat";
        //filename = "dS1_a=" + to_string((int)a) +  "_N=" + to_string((int)points1) + ".dat";
        filename = "shifttest2_kernel_4_smoothcutoff.dat";
        //filename = "tmp3.dat";
    } 
    ofstream fout;
    fout.open(filename.c_str());

    ifstream fin;
    
    

    int count = 0;

    //ofstream fout1;
    string qvecfile = "contour_s=8.6-0.05i.dat";
    //fout1.open(qvecfile.c_str());

    //double eta = 20.0;

    //for(int i=0;i<delspoints+1;++i)
    //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
    //{
    //for(int N=2000;N<=10000;N=N+500)
    //int N = 5000; 
    //for(double somepoint=-0.06;somepoint<=0.08;somepoint=somepoint+0.002)
    
    {
        int somepoint = 400;
        points1 = N;
        //double s3real = sinitial + i*dels;
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);

        comp s = s3real + ii*s3imag;
        cout<<"am = "<<a<<endl;
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp qval = pmom(s,sigb,m);

        //comp qval = pmom(s,2.0*m*m,m);
        
        //cout<<"shift = "<<somepoint<<endl;
        cout<<"qval = "<<qval<<endl;

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        cout<<"kmax = "<<kmax<<endl;
        
        //cout<<"s:"<<s<<endl;
        cout<<"run: "<<count+1<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;
        vector<comp> weights;

        int tag1=0,tag2=0;
        
        //string input_contour = "contour_N=" + to_string((int) N) + ".dat";
        string input_contour = "contour_million.dat";//"contour_" + to_string(somepoint) + ".dat";
        
        //fin.open(input_contour.c_str());
        //comp startingpoint = kmin;
        //comp endingpoint = kmax;
        //line_maker(qvec,startingpoint,endingpoint,points1);

        //mom_vector_maker_43_1(qvec,points1);

        //mom_vector_maker_43_2(somepoint, qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        
        int switch_for_gvec_fixer = 0;//0 means use the gvec_fixer_1
                                      //1 means dont use it
        //if(s3imag<=0.0)
        //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        //mom_vector_maker_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1);
                
        //else if(s3imag>0.0)
        //{
            //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            //mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
            //mom_vector_maker_47_2(somepoint, qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
            //mom_vector_maker_47_2_2(somepoint,qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
            
            //mom_vector_maker_47_2_2(somepoint,qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
            //mom_vector_maker_47_2_3(somepoint,qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
            //mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
            //mom_vector_maker_47_2_A1(somepoint,qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
            //mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
            mom_vector_maker_seba_imspos(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
            //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
            //mom_vector_maker_47_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
            //mom_vector_maker_seba_imspos_2(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                
            //mom_vector_maker_seba_imspos(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
               
            cout<<"tag1 = "<<tag1<<endl;
            cout<<"tag2 = "<<tag2<<endl;
                
                
        //}

        /*double indxbla, qvec_real,qvec_imag;
        while(fin>>indxbla>>qvec_real>>qvec_imag)
        {
            comp qvec_value = qvec_real + ii*qvec_imag;
            qvec.push_back(qvec_value);
        }
        fin.close();
        tag1 = 44;
        tag2 = 92;
        */
        for(int i=0;i<qvec.size();++i)
        {
            //cout<<" qvec found at i = "<<i<<'\t'<<" qvec = "<<qvec[i]<<endl;
        }

        cout<<"qvec created with size = "<<qvec.size()<<endl;
        int size = qvec.size();

        for(int i=0;i<qvec.size();++i)
        {
            //fout1<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
        }

        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

            
        if(s3imag<=0.0)
        {
            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        }
        else if(s3imag>0.0)
        {
            //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);

            Gvec_maker_momrep_withtags_1(Gvec,s,qvec,qval,m,eps,tag1,tag2);
        }

        if(switch_for_gvec_fixer==0)
        {
                //Gvec_fixer_1(Gvec,Gvec,qvec,s,qval,m,eps,tag1,tag2); 
                //Gvec_fixer_2(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
                //Gvec_fixer_3(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
                   
        }
        cout<<"got tags from fixer3 = "<<tag1<<'\t'<<tag2<<endl;

        //tag_fixer_for_gvec_fixer3(Gvec,qvec,s,qval,m,eps,tag1);
        //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
        Eigen::VectorXcd interpolater_Gvec = Gvec;
        //Gvec_fixer_1(Gvec,Gvec,qvec,s,qval,m,eps,tag1,tag2);
        //Gvec_fixer_2(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
        //tag_fixer(Gvec,qvec,s,qval,m,eps,tag1);
        cout<<"fixed tag1 = "<<tag1<<endl;

        ofstream fout2;//("opefile.dat");
        string opefile = "opefile.dat";
        fout2.open(opefile.c_str());
        ofstream fout3; 
        string opefile1 = "opefile1.dat";

        fout3.open(opefile1.c_str());

        for(int i=0;i<qvec.size();++i)
        {
            fout2<<i<<'\t'<<real(Gvec[i])<<'\t'<<imag(Gvec[i])<<endl;
            fout3<<i<<'\t'<<real(GS_pk(s,qvec[i],qval,m,eps))<<'\t'<<imag(GS_pk(s,qvec[i],qval,m,eps))<<'\t'<<real(GS_pk_secondsheet(s,qvec[i],qval,m,eps))<<'\t'<<imag(GS_pk_secondsheet(s,qvec[i],qval,m,eps))<<endl;
        }
        fout2.close();
        fout3.close();
        //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);

        cout<<"ran here"<<endl;
        
        Gvec = -1.0*Gvec;

        for(int i=0;i<qvec.size();++i)
        {
            //cout<<" OPE at i = "<<i<<'\t'<<" G = "<<Gvec[i]<<endl;
        }

        //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
        //Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
        Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

        Eigen::MatrixXcd onemat(qvec.size(),qvec.size());
        for(int i=0;i<qvec.size();++i)
        {
            for(int j=0;j<qvec.size();++j)
            {
                if(i==j)
                    onemat(i,j) = {1.0,0.0};
                else
                    onemat(i,j) = {0.0,0.0};
            }
        }
        //Bmat = Bmat - onemat;
        //this section will print out the bmat 
        ofstream fout4;
        
        for(int i=0;i<qvec.size();++i)
        {
            string BmatFile = "bmat_" + to_string(i) + ".dat";
            //fout4.open(BmatFile.c_str());
            for(int j=0;j<qvec.size();++j)
            {
                //cout<<"B("<<i<<","<<j<<") = "<<Bmat(i,j)<<endl;
                //fout4<<i<<'\t'<<j<<'\t'<<real(Bmat(i,j))<<'\t'<<imag(Bmat(i,j))<<'\t'<<real(qvec[j])<<'\t'<<imag(qvec[j])<<'\t'<<real(weights[j])<<'\t'<<imag(weights[j])<<endl;
            }
            //fout4.close();
        }

        for(int i=0;i<qvec.size();++i)
        {
            //comp sigk = sigma_p(s,qvec[i],m);
            //cout<<"sig k = "<<sigk<<endl;
            //cout<<"for qvec = "<<qvec[i]<<'\t'<<" M2 = "<<M2kfunc(a,sigk,m,eps_for_m2k)<<'\t'<<"omega = "<<omega_comp(qvec[i],m)<<endl;
        }
        //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);

        //cout<<"determinant of B = "<<Bmat.determinant()<<endl;
        //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
        double relerror;
            
        
        //cout<<"did this = "<<count<<" times"<<endl;
        //cout<<"i val = "<<i<<'\t';
        //cout<<"j val = "<<j<<endl;
            
        //LinearSolver_2(Bmat,dsol,Gvec,relerror);

        time(&time_start);
        //LinearSolver_2(Bmat,dsol,Gvec,relerror);
        cusolverComplex(Bmat,Gvec,dsol,size);
        time(&time_end);

        double time_taken = double(time_end - time_start);
        cout<<"Time taken to solve : "<<fixed 
            <<time_taken<<setprecision(5);
        cout<<" sec"<<endl;


        comp result;

        
        if(s3imag<=0.0)
        {
            interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            
            //interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
        }
        else if(s3imag>0.0)
        {
            interpolator_ds_integraleq_momrep_2eps_withtags_with_weights_usingGvec(dsol,qvec,weights,interpolater_Gvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
            
            
            //interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
        }
        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        comp ds = result;
        result = gsq*result;

        
        comp rhopb = rhophib(s,a,m);
        comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
        comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
        comp GSval = GS_pk(s,qval,qval,m,eps);
        //result = rhophib(s,a,m)*result;

        cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        //cout<<"somecount = "<<somepoint<<endl;
        //fout<<setprecision(16)<<real(s)<<'\t'
        //    <<imag(s)<<'\t'
        //    <<points1<<'\t'
        //    <<real(result)<<'\t'
        //    <<imag(result)<<endl;
            //<<real(ds)<<'\t'
            //<<imag(ds)<<'\t'
            //<<real(mphib2)<<'\t'
            //<<imag(mphib2)<<'\t'
            //<<real(mphib2denom)<<'\t'
            //<<imag(mphib2denom)<<'\t'
            //<<a<<'\t'
            //<<N<<endl;
        //count = count + 1;

        fout<<setprecision(16)<<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<N<<'\t'
            <<qvec.size()<<'\t'
            <<real(result)<<'\t'
            <<imag(result)<<endl;
        cout<<"-------------------------"<<endl;
        
        //comp somes = 8.6 + 0.05*ii;
        
        //comp somep = 0.5*qval;
        //comp somek = 0.5*qval;

        //cout<<"s = "<<s<<" p = "<<somep<<" k = "<<somek<<endl;
        //kernel_pk_2eps(s,somep,somek,a,m,eps,eps_for_m2k)
        //cout<<setprecision(15)<<" GS = "<<GS_pk(s,somep,somek,m,eps)<<endl;
        //cout<<setprecision(15)<<" ker = "<<kernel_pk_2eps(s,somep,somek,a,m,eps,eps_for_m2k)<<endl;
    //}
    }
    fout.close();
    //fout1.close();
    //cout<<"GS = "<<GS_pk(s,qval,qval,m,eps)<<endl;

}

void homogeneous_equation_solution_BS1( double a,
                                        int N   )
{
    //double a = 16.0;
    double m = 1.0;
    double eps = 0.0;

    double spole = 8.78285;
    

    //int N = 10000;
    poletest_bs1_vs_singleN_singlea_using_determinant(a, N,spole);
    cout<<"pole found at s = "<<spole<<endl;

    double kmin = 0.0;
    comp kmax = pmom(spole,0.0,m);//1.3;
    double kpoints = 100.0;
    //double delk = abs(kmax - kmin)/kpoints;
    double delk = abs(kmax - kmin)/N;

    ofstream fout;
    string filename = "vertexfactor_vs_k_BS1_for_a_" + to_string((int)a) + "_N_" + to_string(N) + "_withoutanalyticH.dat";
    fout.open(filename.c_str());

    //for(int i=0;i<kpoints;++i)
    //{
        //double k = kmin + i*delk;

        vector<comp> qvec;
        mom_vector_maker_linear_1(qvec,kmin,delk,N,1);
        int size = qvec.size();

        Eigen::MatrixXcd Bmat(size,size);
        Eigen::MatrixXcd eigBmat(size,size);
        
        
        Bmat_maker_momrep(Bmat,spole,qvec,qvec,a,m,eps);

        //std::cout << "Here is the matrix A:\n" << Bmat << std::endl;
        //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(Bmat);
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigensolver(Bmat);
        
        if (eigensolver.info() != Eigen::Success) abort();
        //std::cout   << "The eigenvalues of A are:\n" << eigensolver.eigenvalues() << std::endl;
        //std::cout   << "Here's a matrix whose columns are eigenvectors of A \n"
        //            << "corresponding to these eigenvalues:\n"
        //            << eigensolver.eigenvectors() << std::endl;

        Eigen::VectorXcd eigVal(size);
        eigVal = eigensolver.eigenvalues();

        eigBmat = eigensolver.eigenvectors();
        //cout<<"eigBmat = "<<eigBmat<<endl;

        //for(int i=0;i<eigVal.size();++i)
        //{
            //cout<<"eigenvalues = "<<eigVal(i)<<endl;
        //} 

        //for(int i=0;i<eigVal.size();++i)
        //{
            //for(int j=0;j<eigVal.size();++j)
            //{
                //cout<<"eigvec "<<i<<","<<j<<" = "<<eigBmat(j,i)<<endl;
            //}
        //}
        double sumvert = 0.0;
        //cout<<"0 column = "<<eigensolver.eigenvectors().col(0)<<endl;
        //cout<<"0 column printer = "<<endl;
        vector<double> m2kresvec;
        vector<double> vertexressqvec;
        vector<double> resvec;
        //for(int i=0;i<eigVal.size();++i)
        //{
            int i = 0;
            for(int j=0;j<eigVal.size();++j)
            {
                //cout<<"eigvec "<<i<<","<<j<<" = "<<eigBmat(j,i)<<endl;
                comp k = qvec[j];
                
                comp sigk = sigma_p(spole,k,m);
                comp m2kres = M2kfunc(a,sigk,m,eps);
                comp m2kressq = m2kres*m2kres;
                comp vertexres = real(eigBmat(j,i));
                comp vertexressq = vertexres*vertexres;
                comp res = m2kressq*vertexressq;
                m2kresvec.push_back(real(m2kressq));
                vertexressqvec.push_back(real(vertexressq));
                resvec.push_back(real(res));
                //fout<<real(qvec[j])<<'\t'<<real(m2kressq)<<'\t'<<real(vertexressq)<<'\t'<<real(res)<<endl;
                //cout<<j<<'\t'<<real(qvec[j])<<'\t'<<real(m2kressq)<<'\t'<<real(vertexressq)<<'\t'<<real(res)<<endl;
                if(j==0)
                {
                    sumvert = sumvert + (double)real(vertexressq)*real(qvec[j]);    
                }
                else 
                {
                    sumvert = sumvert + (double)real(vertexressq)*(real(qvec[j] - qvec[j-1]));
                }
                
            }
        //}

        for(int i=0;i<eigVal.size();++i)
        {
            fout<<real(qvec[i])<<'\t'<<m2kresvec[i]<<'\t'<<vertexressqvec[i]/sumvert<<'\t'<<resvec[i]/sumvert<<endl;
            cout<<i<<'\t'<<real(qvec[i])<<'\t'<<m2kresvec[i]<<'\t'<<vertexressqvec[i]/sumvert<<'\t'<<resvec[i]/sumvert<<endl;
            
        }

        cout<<"normalization constant = "<<setprecision(16)<<sumvert<<endl;

        comp traceval = {0.0,0.0};
        comp tracevaleig = {0.0,0.0};
        for(int i=0;i<size;++i)
        {
            traceval = traceval + Bmat(i,i);
            tracevaleig = tracevaleig + eigVal(i);
        }

        cout<<"trace from B = "<<traceval<<'\t'<<" sum of eigval = "<<tracevaleig<<endl;
        comp sigmax = sigma_p(spole,0.0,m);
        cout<<"sigmax = "<<sigmax<<endl;


    //}
}

void homogeneous_equation_solution_BS2( double a,
                                        int N)
{
    //double a = 16.0;
    double m = 1.0;
    double eps = 0.0;

    double spole = 8.78285;
    

    //int N = 10000;
    poletest_bs2_vs_singleN_singlea_using_determinant(a,N,spole);
    cout<<"pole found at s = "<<spole<<endl;

    double kmin = 0.0;
    comp kmax = pmom(spole,0.0,m);//1.3;
    double kpoints = 100.0;
    //double delk = abs(kmax - kmin)/kpoints;
    double delk = abs(kmax - kmin)/N;

    ofstream fout;
    string filename = "vertexfactor_vs_k_BS2_for_a_" + to_string((int)a) + "_N_" + to_string(N) + "_withoutanalyticH.dat";
    fout.open(filename.c_str());

    //for(int i=0;i<kpoints;++i)
    //{
        //double k = kmin + i*delk;

        vector<comp> qvec;
        mom_vector_maker_linear_1(qvec,kmin,delk,N,1);
        int size = qvec.size();

        Eigen::MatrixXcd Bmat(size,size);
        Eigen::MatrixXcd eigBmat(size,size);
        
        
        Bmat_maker_momrep(Bmat,spole,qvec,qvec,a,m,eps);

        //std::cout << "Here is the matrix A:\n" << Bmat << std::endl;
        //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(Bmat);
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigensolver(Bmat);
        
        if (eigensolver.info() != Eigen::Success) abort();
        //std::cout   << "The eigenvalues of A are:\n" << eigensolver.eigenvalues() << std::endl;
        //std::cout   << "Here's a matrix whose columns are eigenvectors of A \n"
        //            << "corresponding to these eigenvalues:\n"
        //            << eigensolver.eigenvectors() << std::endl;

        Eigen::VectorXcd eigVal(size);
        eigVal = eigensolver.eigenvalues();

        eigBmat = eigensolver.eigenvectors();
        //cout<<"eigBmat = "<<eigBmat<<endl;

        //for(int i=0;i<eigVal.size();++i)
        //{
            //cout<<"eigenvalues = "<<eigVal(i)<<endl;
        //} 

        //for(int i=0;i<eigVal.size();++i)
        //{
            //for(int j=0;j<eigVal.size();++j)
            //{
                //cout<<"eigvec "<<i<<","<<j<<" = "<<eigBmat(j,i)<<endl;
            //}
        //}
        double sumvert = 0.0;
        //cout<<"0 column = "<<eigensolver.eigenvectors().col(0)<<endl;
        //cout<<"0 column printer = "<<endl;
        vector<double> m2kresvec;
        vector<double> vertexressqvec;
        vector<double> resvec;
        //for(int i=0;i<eigVal.size();++i)
        //{
            int i = 0;
            for(int j=0;j<eigVal.size();++j)
            {
                //cout<<"eigvec "<<i<<","<<j<<" = "<<eigBmat(j,i)<<endl;
                comp k = qvec[j];
                
                comp sigk = sigma_p(spole,k,m);
                comp m2kres = M2kfunc(a,sigk,m,eps);
                comp m2kressq = m2kres*m2kres;
                comp vertexres = real(eigBmat(j,i));
                comp vertexressq = vertexres*vertexres;
                comp res = m2kressq*vertexressq;
                m2kresvec.push_back(real(m2kressq));
                vertexressqvec.push_back(real(vertexressq));
                resvec.push_back(real(res));
                //fout<<real(qvec[j])<<'\t'<<real(m2kressq)<<'\t'<<real(vertexressq)<<'\t'<<real(res)<<endl;
                //cout<<j<<'\t'<<real(qvec[j])<<'\t'<<real(m2kressq)<<'\t'<<real(vertexressq)<<'\t'<<real(res)<<endl;
                if(j==0)
                {
                    sumvert = sumvert + (double)real(vertexressq)*real(qvec[j]);    
                }
                else 
                {
                    sumvert = sumvert + (double)real(vertexressq)*(real(qvec[j] - qvec[j-1]));
                }
            }
        //}

        for(int i=0;i<eigVal.size();++i)
        {
            fout<<real(qvec[i])<<'\t'<<m2kresvec[i]<<'\t'<<vertexressqvec[i]/sumvert<<'\t'<<resvec[i]/sumvert<<endl;
            cout<<i<<'\t'<<real(qvec[i])<<'\t'<<m2kresvec[i]<<'\t'<<vertexressqvec[i]/sumvert<<'\t'<<resvec[i]/sumvert<<endl;
            
        }

        cout<<"normalization constant = "<<setprecision(16)<<sumvert<<endl;

        comp traceval = {0.0,0.0};
        comp tracevaleig = {0.0,0.0};
        for(int i=0;i<size;++i)
        {
            traceval = traceval + Bmat(i,i);
            tracevaleig = tracevaleig + eigVal(i);
        }

        cout<<"trace from B = "<<traceval<<'\t'<<" sum of eigval = "<<tracevaleig<<endl;
        comp sigmax = sigma_p(spole,0.0,m);
        cout<<"sigmax = "<<sigmax<<endl;


    //}
}

void homogeneous_equation_solution_BS3( double a,
                                        int N)
{
    //double a = 16.0;
    double m = 1.0;
    double eps = 0.0;

    double spole = 8.78285;
    

    //int N = 10000;
    poletest_bs3_vs_singleN_singlea_using_determinant(a,N,spole);
    cout<<"pole found at s = "<<spole<<endl;

    double kmin = 0.0;
    comp kmax = pmom(spole,0.0,m);//1.3;
    double kpoints = 100.0;
    //double delk = abs(kmax - kmin)/kpoints;
    double delk = abs(kmax - kmin)/N;

    ofstream fout;
    string filename = "vertexfactor_vs_k_BS3_for_a_" + to_string((int)a) + "_N_" + to_string(N) + "_withanalyticH.dat";
    fout.open(filename.c_str());

    //for(int i=0;i<kpoints;++i)
    //{
        //double k = kmin + i*delk;

        vector<comp> qvec;
        mom_vector_maker_linear_1(qvec,kmin,delk,N,1);
        int size = qvec.size();

        Eigen::MatrixXcd Bmat(size,size);
        Eigen::MatrixXcd eigBmat(size,size);
        
        
        Bmat_maker_momrep(Bmat,spole,qvec,qvec,a,m,eps);

        //std::cout << "Here is the matrix A:\n" << Bmat << std::endl;
        //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(Bmat);
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigensolver(Bmat);
        
        if (eigensolver.info() != Eigen::Success) abort();
        //std::cout   << "The eigenvalues of A are:\n" << eigensolver.eigenvalues() << std::endl;
        //std::cout   << "Here's a matrix whose columns are eigenvectors of A \n"
        //            << "corresponding to these eigenvalues:\n"
        //            << eigensolver.eigenvectors() << std::endl;

        Eigen::VectorXcd eigVal(size);
        eigVal = eigensolver.eigenvalues();

        eigBmat = eigensolver.eigenvectors();
        //cout<<"eigBmat = "<<eigBmat<<endl;

        //for(int i=0;i<eigVal.size();++i)
        //{
            //cout<<"eigenvalues = "<<eigVal(i)<<endl;
        //} 

        //for(int i=0;i<eigVal.size();++i)
        //{
            //for(int j=0;j<eigVal.size();++j)
            //{
                //cout<<"eigvec "<<i<<","<<j<<" = "<<eigBmat(j,i)<<endl;
            //}
        //}
        double sumvert = 0.0;
        //cout<<"0 column = "<<eigensolver.eigenvectors().col(0)<<endl;
        //cout<<"0 column printer = "<<endl;
        vector<double> m2kresvec;
        vector<double> vertexressqvec;
        vector<double> resvec;
        //for(int i=0;i<eigVal.size();++i)
        //{
            int i = 0;
            for(int j=0;j<eigVal.size();++j)
            {
                //cout<<"eigvec "<<i<<","<<j<<" = "<<eigBmat(j,i)<<endl;
                comp k = qvec[j];
                
                comp sigk = sigma_p(spole,k,m);
                comp m2kres = M2kfunc(a,sigk,m,eps);
                comp m2kressq = m2kres*m2kres;
                comp vertexres = real(eigBmat(j,i));
                comp vertexressq = vertexres*vertexres;
                comp res = m2kressq*vertexressq;
                m2kresvec.push_back(real(m2kressq));
                vertexressqvec.push_back(real(vertexressq));
                resvec.push_back(real(res));
                //fout<<real(qvec[j])<<'\t'<<real(m2kressq)<<'\t'<<real(vertexressq)<<'\t'<<real(res)<<endl;
                //cout<<j<<'\t'<<real(qvec[j])<<'\t'<<real(m2kressq)<<'\t'<<real(vertexressq)<<'\t'<<real(res)<<endl;
                if(j==0)
                {
                    sumvert = sumvert + (double)real(vertexressq)*real(qvec[j]);    
                }
                else 
                {
                    sumvert = sumvert + (double)real(vertexressq)*(real(qvec[j] - qvec[j-1]));
                }
                
            }
        //}

        for(int i=0;i<eigVal.size();++i)
        {
            fout<<real(qvec[i])<<'\t'<<m2kresvec[i]<<'\t'<<vertexressqvec[i]/sumvert<<'\t'<<resvec[i]/sumvert<<endl;
            cout<<i<<'\t'<<real(qvec[i])<<'\t'<<m2kresvec[i]<<'\t'<<vertexressqvec[i]/sumvert<<'\t'<<resvec[i]/sumvert<<endl;
            
        }

        cout<<"normalization constant = "<<setprecision(16)<<sumvert<<endl;

        comp traceval = {0.0,0.0};
        comp tracevaleig = {0.0,0.0};
        for(int i=0;i<size;++i)
        {
            traceval = traceval + Bmat(i,i);
            tracevaleig = tracevaleig + eigVal(i);
        }

        cout<<"trace from B = "<<traceval<<'\t'<<" sum of eigval = "<<tracevaleig<<endl;
        comp sigmax = sigma_p(spole,0.0,m);
        cout<<"sigmax = "<<sigmax<<endl;


    //}
}

void homogeneous_equation_phib_solution_BS1()
{
    ifstream fin;
    string inputfile = "bs1st_N_500.dat";
    fin.open(inputfile.c_str());
    double some_a,some_b,some_c,some_d;
    vector<double> avec;
    vector<double> spolevec;
    int N;
    while(fin>>some_a>>some_b>>some_c>>some_d)
    {
        avec.push_back(some_a);
        spolevec.push_back(some_b);
        N = some_d;
    }

    //int N=500;
    //while(fin>>some_a>>some_b>>some_c>>some_d)
    /*for(double a=16.0;a<=50.0;a=a+1.0)
    {
        double spole = 0.0;
        //poletest_bs1_vs_singleN_singlea_using_determinant(a, N,spole);
        poletest_bs1_vs_singleN_singlea_using_weights_for_vertexfactor(a, N,spole);
        avec.push_back(a);
        spolevec.push_back(spole);
        
    }*/
    for(int i=0;i<avec.size();++i)
    {
        cout<<"pole found for a = "<<setprecision(16)<<avec[i]<<" at s = "<<spolevec[i]<<endl;
    }
    //double a = 16.0;
    double m = 1.0;
    double eps = 0.0001;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    
    

    //double spole = 8.78285;
    

    //int N = 10000;
    //poletest_bs1_vs_singleN_singlea_using_determinant(a, N,spole);
    //cout<<"pole found at s = "<<spole<<endl;

    
    
    double kpoints = 100.0;
    //double delk = abs(kmax - kmin)/kpoints;
    //double delk = abs(kmax - kmin)/N;

    ofstream fout;
    string filename = "vertexfactor_k=1_vs_a_BS1_N_" + to_string(N) + "_withanalyticH_7.dat";
    fout.open(filename.c_str());

    for(int somei=0;somei<avec.size();++somei)
    {
        double a = avec[somei];
        double spole = spolevec[somei];
    //for(int i=0;i<kpoints;++i)
    //{
        //double k = kmin + i*delk;
        double s = spole;
        double kmin = 0.0;
        comp kmax = pmom(spole,0.0,m);//1.3;
        double points1 = N;
        vector<comp> qvec;
        vector<comp> weights;
        //mom_vector_maker_linear_1(qvec,kmin,delk,N,1);
        line_maker_with_weights(qvec,weights,kmin,kmax,points1);
        //mom_vector_maker_43_with_weights(qvec,weights,s,kmin,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                
        int size = qvec.size();

        Eigen::MatrixXcd Bmat(size,size);
        Eigen::MatrixXcd eigBmat(size,size);
        
        
        //Bmat_maker_momrep(Bmat,spole,qvec,qvec,a,m,eps);
        Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);


        //std::cout << "Here is the matrix A:\n" << Bmat << std::endl;
        //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(Bmat);
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigensolver(Bmat);
        
        if (eigensolver.info() != Eigen::Success) abort();
        //std::cout   << "The eigenvalues of A are:\n" << eigensolver.eigenvalues() << std::endl;
        //std::cout   << "Here's a matrix whose columns are eigenvectors of A \n"
        //            << "corresponding to these eigenvalues:\n"
        //            << eigensolver.eigenvectors() << std::endl;

        Eigen::VectorXcd eigVal(size);
        eigVal = eigensolver.eigenvalues();

        eigBmat = eigensolver.eigenvectors();
        //cout<<"eigBmat = "<<eigBmat<<endl;

        //for(int i=0;i<eigVal.size();++i)
        //{
            //cout<<"eigenvalues = "<<eigVal(i)<<endl;
        //} 

        //for(int i=0;i<eigVal.size();++i)
        //{
            //for(int j=0;j<eigVal.size();++j)
            //{
                //cout<<"eigvec "<<i<<","<<j<<" = "<<eigBmat(j,i)<<endl;
            //}
        //}
        comp sumvert = 0.0;
        comp sqrtsumvert = 0.0;
        //cout<<"0 column = "<<eigensolver.eigenvectors().col(0)<<endl;
        //cout<<"0 column printer = "<<endl;
        vector<double> m2kresvec;
        vector<double> vertexressqvec;
        vector<double> resvec;
        vector<comp> vertexresvec;
        //for(int i=0;i<eigVal.size();++i)
        //{
            int i = 0;
            for(int j=0;j<eigVal.size();++j)
            {
                //cout<<"eigvec "<<i<<","<<j<<" = "<<eigBmat(j,i)<<endl;
                comp k = qvec[j];
                cout<<"eigen value = "<<eigVal(i)<<endl;
                comp sigk = sigma_p(spole,k,m);
                comp m2kres = M2kfunc(a,sigk,m,eps);
                comp m2kressq = m2kres*m2kres;
                comp vertexres = (eigBmat(j,i));
                comp vertexressq = vertexres*vertexres;
                comp res = m2kressq*vertexressq;
                vertexresvec.push_back(vertexres);
                m2kresvec.push_back(abs(m2kressq));
                vertexressqvec.push_back(abs(vertexressq));
                resvec.push_back(real(res));
                //fout<<real(qvec[j])<<'\t'<<real(m2kressq)<<'\t'<<real(vertexressq)<<'\t'<<real(res)<<endl;
                //cout<<j<<'\t'<<real(qvec[j])<<'\t'<<real(m2kressq)<<'\t'<<real(vertexressq)<<'\t'<<real(res)<<endl;
                if(j==0)
                {
                    //sumvert = sumvert + abs(vertexressq)*(qvec[j]);
                    sumvert = sumvert + abs(vertexressq)*(weights[j]);    
                    
                }
                else 
                {
                    //sumvert = sumvert + abs(vertexressq)*((qvec[j] - qvec[j-1]));
                    sumvert = sumvert + abs(vertexressq)*(weights[j]);
                
                }
                
            }
        //}
        cout<<"normalization constant = "<<setprecision(16)<<sumvert<<endl;
        sqrtsumvert = sqrt(sumvert);
        vector<comp> normalized_vertexres;
        for(int i=0;i<vertexresvec.size();++i)
        {
            normalized_vertexres.push_back(vertexresvec[i]/sqrtsumvert);
        }

        comp vertexres_at_q = {0.0,0.0};
        comp q = 1.0;//pmom(s,sigmab(a,m),m);
        for(int i=0;i<qvec.size();++i)
        {
            //comp q = 0.000001;//pmom(s,sigmab(a,m),m);
            comp kern = kernel_pk_2eps(s,q,qvec[i],a,m,eps,eps_for_m2k);
            comp weight = weights[i];
            comp norm_vert = normalized_vertexres[i];
            vertexres_at_q = vertexres_at_q  -  weight*kern*norm_vert;

        }
        //q = qvec[0];
        //vertexres_at_q = vertexresvec[0];

        comp sigk = sigma_p(spole,q,m);//sigmab(a,m);
        //comp q = pmom(s,sigmab(a,m),m);
        comp m2kres = M2kfunc(a,sigk,m,eps);
        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        fout<<setprecision(16)<<a<<'\t'
            <<real(vertexres_at_q)<<'\t'
            <<imag(vertexres_at_q)<<'\t'
            <<abs(vertexres_at_q)<<'\t'
            <<real(m2kres*vertexres_at_q)<<'\t'
            <<imag(m2kres*vertexres_at_q)<<'\t'
            <<abs(m2kres*vertexres_at_q*m2kres*vertexres_at_q)<<'\t'
            <<real(m2kres)<<'\t'
            <<imag(m2kres)<<'\t'
            //<<real(gfunc*vertexres_at_q)<<'\t'
            //<<imag(gfunc*vertexres_at_q)<<'\t'
            //<<abs(gfunc*vertexres_at_q*gfunc*vertexres_at_q)<<'\t'
            //<<real(gfunc)<<'\t'
            //<<imag(gfunc)<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<endl;

        cout<<setprecision(16)<<"a="<<a<<'\t'
            <<"g= "<<gfunc<<'\t'
            <<"m2k= "<<m2kres<<'\t'
            <<"vert="<<vertexres_at_q<<'\t'
            <<"m2kVert="<<m2kres*vertexres_at_q<<'\t'
            <<"mVrtSq="<<m2kres*vertexres_at_q*m2kres*vertexres_at_q<<'\t'
            <<"q="<<q<<endl;

        //for(int i=0;i<eigVal.size();++i)
        //{
        //    fout<<real(qvec[i])<<'\t'<<m2kresvec[i]<<'\t'<<vertexressqvec[i]/sumvert<<'\t'<<resvec[i]/sumvert<<endl;
        //    cout<<i<<'\t'<<real(qvec[i])<<'\t'<<m2kresvec[i]<<'\t'<<vertexressqvec[i]/sumvert<<'\t'<<resvec[i]/sumvert<<endl;
            
        //}

        

        comp traceval = {0.0,0.0};
        comp tracevaleig = {0.0,0.0};
        for(int i=0;i<size;++i)
        {
            traceval = traceval + Bmat(i,i);
            tracevaleig = tracevaleig + eigVal(i);
        }

        cout<<"trace from B = "<<traceval<<'\t'<<" sum of eigval = "<<tracevaleig<<endl;
        comp sigmax = sigma_p(spole,0.0,m);
        cout<<"sigmax = "<<sigmax<<endl;


    //}
    }
}


//this is for am = 2,6,16
void homogeneous_equation_phib_residue_solution_BS1()
{
    
    double a = 16.0;
    int N = 500;

    double spole_a2_BS1_smoothcutoff = 7.252996307970249;
    double spole_a6_BS1_smoothcutoff = 8.535689549889605;
    double spole_a16_BS1_smoothcutoff = 8.782848835105487;
    double spole_a16_BS2_smoothcutoff = 8.976271055410812;
    double spole_a16_BS1_hardcutoff = 8.689974511859685;
    double spole_a16_BS2_hardcutoff = 8.975451813904034;


    double spole = spole_a16_BS2_smoothcutoff;

    double m = 1.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    
    

    //double spole = 8.78285;
    

    //int N = 10000;
    //poletest_bs1_vs_singleN_singlea_using_determinant(a, N,spole);
    //cout<<"pole found at s = "<<spole<<endl;

    
    
    double kpoints = N;
    //double delk = abs(kmax - kmin)/kpoints;
    //double delk = abs(kmax - kmin)/N;

    ofstream fout;
    string filename = "residue_phib_a_"+to_string((int)a)+"_smoothcutoff_BS2";
    fout.open(filename.c_str());

    //for(int somei=0;somei<avec.size();++somei)
    {
        
    //for(int i=0;i<kpoints;++i)
    //{
        //double k = kmin + i*delk;
        double s = spole;
        double kmin = 0.0;
        comp kmax = pmom(spole,0.0,m);//1.3;
        double points1 = N;
        vector<comp> qvec;
        vector<comp> weights;
        int tag1, tag2, switch_for_gvec_fixer;
        tag1 = 0;
        tag2 = 0;
        switch_for_gvec_fixer=0;
        //mom_vector_maker_linear_1(qvec,kmin,delk,N,1);
        //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
        //mom_vector_maker_43_with_weights(qvec,weights,s,kmin,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
        mom_vector_maker_seba_imspos_2(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
        cout<<"tag1="<<tag1<<'\t'<<"tag2="<<tag2<<endl; 
        int size = qvec.size();

        Eigen::VectorXcd Gvec(size);
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::MatrixXcd eigBmat(size,size);
        comp qval = pmom(s,sigmab(a,m),m);
        if(imag(s)<0.0)
        {
            cout<<"gvec with imag s neg selected"<<endl;
            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        }
        else if(imag(s)>=0.0)
        {
            cout<<"gvec with imag s pos selected"<<endl;
            Gvec_maker_momrep_withtags_1(Gvec,s,qvec,qval,m,eps,tag1,tag2);
        }
        
        if(switch_for_gvec_fixer==0)
            cout<<"gvec fixer switch = "<<switch_for_gvec_fixer<<",deform contour"<<endl;
        else 
            cout<<"gvec fixer switch = "<<switch_for_gvec_fixer<<",take straight line"<<endl;

        if(switch_for_gvec_fixer==0)
        {   
            //cout<<"here the problem lies 1"<<endl;
            Gvec_fixer_3(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
            //Gvec_fixer_4(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);//Gvec fixer has a contour maker inside of it
                
            //cout<<"here the problem lies 2"<<endl;
            tag_fixer_for_gvec_fixer3(Gvec,qvec,s,qval,m,eps,tag1); //this works for gvec_fixer_4 as well

            //cout<<"here the problem lies 3"<<endl;
               
        }
        //Bmat_maker_momrep(Bmat,spole,qvec,qvec,a,m,eps);
        Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);


        //std::cout << "Here is the matrix A:\n" << Bmat << std::endl;
        //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(Bmat);
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigensolver(Bmat);
        
        if (eigensolver.info() != Eigen::Success) abort();
        //std::cout   << "The eigenvalues of A are:\n" << eigensolver.eigenvalues() << std::endl;
        //std::cout   << "Here's a matrix whose columns are eigenvectors of A \n"
        //            << "corresponding to these eigenvalues:\n"
        //            << eigensolver.eigenvectors() << std::endl;

        Eigen::VectorXcd eigVal(size);
        eigVal = eigensolver.eigenvalues();

        eigBmat = eigensolver.eigenvectors();
        //cout<<"eigBmat = "<<eigBmat<<endl;

        //for(int i=0;i<eigVal.size();++i)
        //{
            //cout<<"eigenvalues = "<<eigVal(i)<<endl;
        //} 

        //for(int i=0;i<eigVal.size();++i)
        //{
            //for(int j=0;j<eigVal.size();++j)
            //{
                //cout<<"eigvec "<<i<<","<<j<<" = "<<eigBmat(j,i)<<endl;
            //}
        //}
        comp sumvert = 0.0;
        comp sqrtsumvert = 0.0;
        //cout<<"0 column = "<<eigensolver.eigenvectors().col(0)<<endl;
        //cout<<"0 column printer = "<<endl;
        vector<comp> m2kvec;
        vector<comp> zetavec;
        vector<comp> zetasqvec;
        vector<comp> gammavec;
        vector<comp> gammasqvec;
        //for(int i=0;i<eigVal.size();++i)
        //{
            int i = 0;
            cout<<"eigen value = "<<eigVal(i)<<endl;
                
            for(int j=0;j<eigVal.size();++j)
            {
                //cout<<"eigvec "<<i<<","<<j<<" = "<<eigBmat(j,i)<<endl;
                comp k = qvec[j];
                cout<<"eigen value = "<<eigVal(i)<<endl;
                comp sigk = sigma_p(spole,k,m);
                comp m2kres = M2kfunc(a,sigk,m,eps);
                comp m2kressq = m2kres*m2kres;
                comp zeta = (eigBmat(j,i));
                comp zetasq = zeta*zeta;
                comp gamma = m2kres*zeta;
                comp gammasq = m2kressq*zetasq;
                zetavec.push_back(zeta);
                m2kvec.push_back(m2kressq);
                zetasqvec.push_back(zetasq);
                gammavec.push_back(gamma);
                gammasqvec.push_back(gammasq);
                //fout<<real(qvec[j])<<'\t'<<real(m2kressq)<<'\t'<<real(vertexressq)<<'\t'<<real(res)<<endl;
                //cout<<j<<'\t'<<real(qvec[j])<<'\t'<<real(m2kressq)<<'\t'<<real(vertexressq)<<'\t'<<real(res)<<endl;
                if(j==0)
                {
                    //sumvert = sumvert + abs(vertexressq)*(qvec[j]);
                    sumvert = sumvert + zetasq*(weights[j]);    
                    
                }
                else 
                {
                    //sumvert = sumvert + abs(vertexressq)*((qvec[j] - qvec[j-1]));
                    sumvert = sumvert + zetasq*(weights[j]);
                
                }
                
            }
        //}
        cout<<"normalization constant = "<<setprecision(16)<<sumvert<<endl;
        sqrtsumvert = sqrt(sumvert);
        vector<comp> normalized_zeta;
        for(int i=0;i<zetavec.size();++i)
        {
            normalized_zeta.push_back(zetavec[i]/sqrtsumvert);
        }

        comp zeta_at_q = {0.0,0.0};
        comp q = pmom(s,sigmab(a,m),m);
        for(int i=0;i<qvec.size();++i)
        {
            //comp q = 0.000001;//pmom(s,sigmab(a,m),m);
            //comp kern = kernel_pk_2eps(s,q,qvec[i],a,m,eps,eps_for_m2k);
            comp kern = kernel_pk_2eps_using_Gvec(Gvec,i,s,q,qvec[i],a,m,eps,eps,i,i,tag1,tag2);
            
            comp weight = weights[i];
            comp norm_vert = normalized_zeta[i];
            zeta_at_q = zeta_at_q  -  weight*kern*norm_vert;

        }
        //q = qvec[0];
        //vertexres_at_q = vertexresvec[0];

        comp sigk = sigma_p(spole,q,m);//sigmab(a,m);
        //comp q = pmom(s,sigmab(a,m),m);
        comp m2kres = M2kfunc(a,sigk,m,eps);
        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        comp final_zeta = zeta_at_q;
        comp final_zetasq = zeta_at_q*zeta_at_q;
        comp final_gamma = zeta_at_q*gfunc;
        comp final_gammasq = final_gamma*final_gamma;

        fout<<setprecision(16)
            <<"ma = "<<a<<'\n'
            <<"pole = "<<spole<<'\n'
            <<"sigk = "<<sigk<<'\n'
            <<"sigb = "<<sigmab(a,m)<<'\n'
            <<"q = "<<q<<'\n'
            <<"m2k = "<<m2kres<<'\n'
            <<"|m2k|^2 = "<<abs(m2kres*m2kres)<<'\n'
            <<"g_sb = "<<gfunc<<'\n'
            <<"g_sb^2 = "<<gfunc*gfunc<<'\n'
            <<"zeta = "<<final_zeta<<'\n'
            <<"|zeta|^2 = "<<abs(final_zetasq)<<'\n'
            <<"gamma = "<<final_gamma<<'\n'
            <<"|gamma|^2 = "<<abs(final_gammasq)<<endl;

            
        cout<<setprecision(16)
            <<"ma = "<<a<<'\n'
            <<"pole = "<<spole<<'\n'
            <<"sigk = "<<sigk<<'\n'
            <<"sigb = "<<sigmab(a,m)<<'\n'
            <<"q = "<<q<<'\n'
            <<"m2k = "<<m2kres<<'\n'
            <<"|m2k|^2 = "<<abs(m2kres*m2kres)<<'\n'
            <<"g_sb = "<<gfunc<<'\n'
            <<"g_sb^2 = "<<gfunc*gfunc<<'\n'
            <<"zeta = "<<final_zeta<<'\n'
            <<"|zeta|^2 = "<<abs(final_zetasq)<<'\n'
            <<"gamma = "<<final_gamma<<'\n'
            <<"|gamma|^2 = "<<abs(final_gammasq)<<endl;

          
        //for(int i=0;i<eigVal.size();++i)
        //{
        //    fout<<real(qvec[i])<<'\t'<<m2kresvec[i]<<'\t'<<vertexressqvec[i]/sumvert<<'\t'<<resvec[i]/sumvert<<endl;
        //    cout<<i<<'\t'<<real(qvec[i])<<'\t'<<m2kresvec[i]<<'\t'<<vertexressqvec[i]/sumvert<<'\t'<<resvec[i]/sumvert<<endl;
            
        //}

        

        comp traceval = {0.0,0.0};
        comp tracevaleig = {0.0,0.0};
        for(int i=0;i<size;++i)
        {
            traceval = traceval + Bmat(i,i);
            tracevaleig = tracevaleig + eigVal(i);
        }

        cout<<"trace from B = "<<traceval<<'\t'<<" sum of eigval = "<<tracevaleig<<endl;
        comp sigmax = sigma_p(spole,0.0,m);
        cout<<"sigmax = "<<sigmax<<endl;


    //}
    }
    fout.close();
}


//this gives us the vertexfactor as a function of k
//this uses the integral equation to interpolate to 
//any momenta we want
void homogeneous_equation_integralequation_vs_k_solution_BS1()
{
    /*ifstream fin;
    string inputfile = "bs1st_N_500.dat";
    fin.open(inputfile.c_str());
    double some_a,some_b,some_c,some_d;
    vector<double> avec;
    vector<double> spolevec;
    int N;
    while(fin>>some_a>>some_b>>some_c>>some_d)
    {
        avec.push_back(some_a);
        spolevec.push_back(some_b);
        N = some_d;
    }*/

    //int N=500;
    //while(fin>>some_a>>some_b>>some_c>>some_d)
    /*for(double a=16.0;a<=50.0;a=a+1.0)
    {
        double spole = 0.0;
        //poletest_bs1_vs_singleN_singlea_using_determinant(a, N,spole);
        poletest_bs1_vs_singleN_singlea_using_weights_for_vertexfactor(a, N,spole);
        avec.push_back(a);
        spolevec.push_back(spole);
        
    }*/
    /*for(int i=0;i<avec.size();++i)
    {
        cout<<"pole found for a = "<<setprecision(16)<<avec[i]<<" at s = "<<spolevec[i]<<endl;
    }*/
    
    double a = 16.0;
    double m = 1.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    int N = 500;
    

    double spole = 8.689974511859685;//7.252996307970249;//8.535689549889605;//8.782848835105487;
    

    //int N = 10000;
    //poletest_bs1_vs_singleN_singlea_using_determinant(a, N,spole);
    //cout<<"pole found at s = "<<spole<<endl;

    
    double kmin = 0.000001;
    comp kmax = pmom(spole,0.0,m);//1.3;
    double kpoints = 2000.0;
    double delk = abs(kmax - kmin)/kpoints;
    //double delk = abs(kmax - kmin)/N;

    ofstream fout;
    //string filename = "vertexfactor_vs_a_"+to_string(a)+"BS1_N_" + to_string(N) + "_smoothcutoff.dat";
    string filename = "vertexfactor_vs_a_"+to_string(a)+"BS1_N_" + to_string(N) + "_hardcutoff.dat";
    
    fout.open(filename.c_str());

    //for(int somei=0;somei<avec.size();++somei)
    //{
        //double a = avec[somei];
        //double spole = spolevec[somei];
    for(int i=0;i<kpoints;++i)
    {
        double k = kmin + i*delk;
        double s = spole;
        
        double points1 = N;
        vector<comp> qvec;
        vector<comp> weights;
        //mom_vector_maker_linear_1(qvec,kmin,delk,N,1);
        line_maker_with_weights(qvec,weights,kmin,kmax,points1);
        //mom_vector_maker_43_with_weights(qvec,weights,s,kmin,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                
        int size = qvec.size();

        Eigen::MatrixXcd Bmat(size,size);
        Eigen::MatrixXcd eigBmat(size,size);
        
        
        //Bmat_maker_momrep(Bmat,spole,qvec,qvec,a,m,eps);
        Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);


        //std::cout << "Here is the matrix A:\n" << Bmat << std::endl;
        //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(Bmat);
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigensolver(Bmat);
        
        if (eigensolver.info() != Eigen::Success) abort();
        //std::cout   << "The eigenvalues of A are:\n" << eigensolver.eigenvalues() << std::endl;
        //std::cout   << "Here's a matrix whose columns are eigenvectors of A \n"
        //            << "corresponding to these eigenvalues:\n"
        //            << eigensolver.eigenvectors() << std::endl;

        Eigen::VectorXcd eigVal(size);
        eigVal = eigensolver.eigenvalues();

        eigBmat = eigensolver.eigenvectors();
        //cout<<"eigBmat = "<<eigBmat<<endl;

        //for(int i=0;i<eigVal.size();++i)
        //{
            //cout<<"eigenvalues = "<<eigVal(i)<<endl;
        //} 

        //for(int i=0;i<eigVal.size();++i)
        //{
            //for(int j=0;j<eigVal.size();++j)
            //{
                //cout<<"eigvec "<<i<<","<<j<<" = "<<eigBmat(j,i)<<endl;
            //}
        //}
        comp sumvert = 0.0;
        comp sqrtsumvert = 0.0;
        //cout<<"0 column = "<<eigensolver.eigenvectors().col(0)<<endl;
        //cout<<"0 column printer = "<<endl;
        vector<double> m2kresvec;
        vector<double> vertexressqvec;
        vector<double> resvec;
        vector<comp> vertexresvec;
        //for(int i=0;i<eigVal.size();++i)
        //{
            int some_i = 0;
            for(int j=0;j<eigVal.size();++j)
            {
                //cout<<"eigvec "<<i<<","<<j<<" = "<<eigBmat(j,i)<<endl;
                comp k = qvec[j];
                
                comp sigk = sigma_p(spole,k,m);
                comp m2kres = M2kfunc(a,sigk,m,eps);
                comp m2kressq = m2kres*m2kres;
                comp vertexres = (eigBmat(j,some_i));
                comp vertexressq = vertexres*vertexres;
                comp res = m2kressq*vertexressq;
                vertexresvec.push_back(vertexres);
                m2kresvec.push_back(abs(m2kressq));
                vertexressqvec.push_back(abs(vertexressq));
                resvec.push_back(real(res));
                //fout<<real(qvec[j])<<'\t'<<real(m2kressq)<<'\t'<<real(vertexressq)<<'\t'<<real(res)<<endl;
                //cout<<j<<'\t'<<real(qvec[j])<<'\t'<<real(m2kressq)<<'\t'<<real(vertexressq)<<'\t'<<real(res)<<endl;
                if(j==0)
                {
                    //sumvert = sumvert + abs(vertexressq)*(qvec[j]);
                    sumvert = sumvert + abs(vertexressq)*(weights[j]);    
                    
                }
                else 
                {
                    //sumvert = sumvert + abs(vertexressq)*((qvec[j] - qvec[j-1]));
                    sumvert = sumvert + abs(vertexressq)*(weights[j]);
                
                }
                
            }
        //}
        cout<<"normalization constant = "<<setprecision(16)<<sumvert<<endl;
        sqrtsumvert = sqrt(sumvert);
        vector<comp> normalized_vertexres;
        for(int i=0;i<vertexresvec.size();++i)
        {
            normalized_vertexres.push_back(vertexresvec[i]/sqrtsumvert);
        }

        comp vertexres_at_q = {0.0,0.0};
        comp q = k;//1.0;//pmom(s,sigmab(a,m),m);
        for(int i=0;i<qvec.size();++i)
        {
            //comp q = 0.000001;//pmom(s,sigmab(a,m),m);
            comp kern = kernel_pk_2eps(s,q,qvec[i],a,m,eps,eps_for_m2k);
            comp weight = weights[i];
            comp norm_vert = normalized_vertexres[i];
            vertexres_at_q = vertexres_at_q  -  weight*kern*norm_vert;

        }
        //q = qvec[0];
        //vertexres_at_q = vertexresvec[0];

        comp sigk = sigma_p(spole,q,m);//sigmab(a,m);
        //comp q = pmom(s,sigmab(a,m),m);
        comp m2kres = M2kfunc(a,sigk,m,eps);
        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        fout<<setprecision(16)<<a<<'\t'
            <<real(vertexres_at_q)<<'\t'
            <<imag(vertexres_at_q)<<'\t'
            <<abs(vertexres_at_q)<<'\t'
            <<real(m2kres*vertexres_at_q)<<'\t'
            <<imag(m2kres*vertexres_at_q)<<'\t'
            <<abs(m2kres*vertexres_at_q*m2kres*vertexres_at_q)<<'\t'
            <<real(m2kres)<<'\t'
            <<imag(m2kres)<<'\t'
            //<<real(gfunc*vertexres_at_q)<<'\t'
            //<<imag(gfunc*vertexres_at_q)<<'\t'
            //<<abs(gfunc*vertexres_at_q*gfunc*vertexres_at_q)<<'\t'
            //<<real(gfunc)<<'\t'
            //<<imag(gfunc)<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<endl;

        cout<<setprecision(16)<<"a="<<a<<'\t'
            <<"g= "<<gfunc<<'\t'
            <<"m2k= "<<m2kres<<'\t'
            <<"vert="<<vertexres_at_q<<'\t'
            <<"m2kVert="<<m2kres*vertexres_at_q<<'\t'
            <<"mVrtSq="<<m2kres*vertexres_at_q*m2kres*vertexres_at_q<<'\t'
            <<"q="<<q<<endl;

        //for(int i=0;i<eigVal.size();++i)
        //{
        //    fout<<real(qvec[i])<<'\t'<<m2kresvec[i]<<'\t'<<vertexressqvec[i]/sumvert<<'\t'<<resvec[i]/sumvert<<endl;
        //    cout<<i<<'\t'<<real(qvec[i])<<'\t'<<m2kresvec[i]<<'\t'<<vertexressqvec[i]/sumvert<<'\t'<<resvec[i]/sumvert<<endl;
            
        //}

        

        comp traceval = {0.0,0.0};
        comp tracevaleig = {0.0,0.0};
        for(int i=0;i<size;++i)
        {
            traceval = traceval + Bmat(i,i);
            tracevaleig = tracevaleig + eigVal(i);
        }

        cout<<"trace from B = "<<traceval<<'\t'<<" sum of eigval = "<<tracevaleig<<endl;
        comp sigmax = sigma_p(spole,0.0,m);
        cout<<"sigmax = "<<sigmax<<endl;


    }
    //}
    fout.close();
}


//this loads the referenced_ file for vertexfactor mod square that 
//we calculated using the inverse d_S method. It goes through all 
//the same momenta of that file and solves the homogeneous equation
//and gets the vertex factors. It also performs rescaling the vertex
//factors by `renormalizing' the vertex factors mod squares found from
//homogeneous equation based on the result of the inverse d_S method

void referenced_from_Inverse_dS_method_homogeneous_equation_integralequation_vs_k_solution_BS1()
{

    //input parameters:
    double a = 16.0;
    double m = 1.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    int N = 500;

    //here we note all the bound state poles for N=500
    
    double spole_a2_smoothcutoff = 7.252996307970249;
    double spole_a6_smoothcutoff = 8.535689549889605;
    double spole_a16_bs1_smoothcutoff = 8.782848835105487;
    double spole_a16_bs2_smoothcutoff = 8.976271055410812;

    double spole_a16_bs1_hardcutoff = 8.689974511859685;
    double spole_a16_bs2_hardcutoff = 8.975451813904034;

    
    
    double spole = spole_a16_bs2_smoothcutoff;

    //we load the referenced file in this section
    //this should match the file based on the input parameters 


    ifstream fin;

    //we add all the referenced file names we have used for the calculation
    //this can be modified for future use 

    string filename_a2 = "referenced_usedfor_homEqn_vertexfactor_vs_sigk_BS2_file_a=16_N=500.dat";
    
    //this is the input file name
    string inputfile = (string) filename_a2;
    
    fin.open(inputfile.c_str());
    double some_a,some_b,some_c;
    vector<double> kvec;
    vector<double> vertexsq_vec;
    vector<double> m2vertexsq_vec;

    //here we load the columns of the file into corresponding vectors
    while(fin>>some_a>>some_b>>some_c)
    {
        kvec.push_back(some_a);
        vertexsq_vec.push_back(some_b);
        m2vertexsq_vec.push_back(some_c);

    }
    fin.close();

    //here we set the output file name where the data will be saved 

    ofstream fout;
    string filename =       "vertexfactor_vs_a_"+to_string(a)
                        +   "BS2_N_" + to_string(N) 
                        +   "_smoothcutoff.dat";
    

    
    //we set the number of iterations to the length 
    //of the k vector found from the reference file]
    int kpoints = kvec.size();


    fout.open(filename.c_str());

    int somecount = 0;//for printing the normalization of the reference file 
    double ref_normalization = 0;

    for(int i=0;i<kpoints;++i)
    {
        cout<<"loaded file = "<<filename_a2<<endl;
        cout<<"cut disappears at sigma_k = "<<m*(m+sqrt(spole))<<" when sigma_max = "<<sigma_p(spole,0.0,m)<<endl;
        double k = kvec[i];
        double s = spole;
        
        double points1 = N;
        vector<comp> qvec;
        vector<comp> weights;

        //since the external momentums are set from 0 to kmax, there is 
        //no circular cut. So the integration interval can be taken as a 
        //straight line 

        double kmin = 0.0;
        comp kmax = pmom(s,0.0,m);

        line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                
        int size = qvec.size();

        Eigen::MatrixXcd Bmat(size,size);
        Eigen::MatrixXcd eigBmat(size,size);
        
        
        Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

        //we use the whole Bmat as our operator for which we are performing
        //the eigen decomposition. BV = (1+K)V = lambda V. The solution through
        //the library function below gives us all the eigen values and eigenvectors
        //we separate out the eigenvectors that correspond to the zero eigenvalue or 
        //in this case, the 'smallest' eigen value, which tends to be the first column
        //in the eigBmat matrix.

        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigensolver(Bmat);
        
        if (eigensolver.info() != Eigen::Success) abort();

        Eigen::VectorXcd eigVal(size);
        eigVal = eigensolver.eigenvalues();

        eigBmat = eigensolver.eigenvectors();


        comp sumvert_eigen = {0.0,0.0};
        comp sumverthomeqn = {0.0,0.0};
        comp sqrtsumverthomeqn = {0.0,0.0};

        vector<comp> homeqn_m2k_vec;
        vector<comp> homeqn_vertex_vec;
        vector<comp> homeqn_vertexsq_vec;
        vector<comp> homeqn_m2vertexsq_vec;
        vector<comp> vertexresvec;


        //in this section, we are separating out the eigenvector with zero 
        //eigenvalue that we need and normalize them to unity.

        int some_i = 0;
        for(int j=0;j<eigVal.size();++j)
        {
            comp q = qvec[j];
                
            comp sigk = sigma_p(spole,q,m);
            comp m2k = M2kfunc(a,sigk,m,eps);
            comp m2ksq = m2k*m2k;
            comp vertex_eigvec = (eigBmat(j,some_i));
            comp vertexsq_eigvec = vertex_eigvec*vertex_eigvec;
            comp homeqn_m2vertexsq_val = m2ksq*vertexsq_eigvec;
            homeqn_m2k_vec.push_back(abs(m2ksq));
            homeqn_vertex_vec.push_back(vertex_eigvec);
            //cout<<"solution of vertexfactor for "<<endl;
            //cout<<"k = "<<q<<'\n'
            //    <<"sol. vert = "<<vertex_eigvec<<endl;
            homeqn_vertexsq_vec.push_back(vertexsq_eigvec);
            homeqn_m2vertexsq_vec.push_back(homeqn_m2vertexsq_val);

            sumverthomeqn = sumverthomeqn + abs(vertexsq_eigvec)*(weights[j]);
            sumvert_eigen = sumvert_eigen + abs(vertexsq_eigvec);
                
        }


        //we normalize the vertexfactors found in the reference file just once and
        //print that for comparison 

        if(somecount==0)
        {
            double delk = 0.0;
            for(int i=0;i<kvec.size();++i)
            {
                if(i==0)
                {
                    delk = kvec[i];
                }
                else 
                {
                    delk = kvec[i] - kvec[i-1];
                }

                ref_normalization = ref_normalization + abs(vertexsq_vec[i])*delk;
            }

            somecount = 1;
        }

        

        cout<<"normalization of eigenvector = "<<setprecision(16)<<sumvert_eigen<<endl;
        cout<<"normalization constant, W = "<<setprecision(16)<<sumverthomeqn<<endl;
        cout<<"normalization of the ref vertexsq, Z = "<<ref_normalization<<endl;
        sqrtsumverthomeqn = sqrt(sumverthomeqn);

        //in this section, we are rescaling the vertex factor and their mod square 
        //based on the normalization 

        vector<comp> normalized_vertex;
        vector<comp> normalized_vertexsq;
        for(int i=0;i<eigVal.size();++i)
        {
            normalized_vertex.push_back(homeqn_vertex_vec[i]/sqrtsumverthomeqn);
            normalized_vertexsq.push_back(homeqn_vertexsq_vec[i]/sumverthomeqn);
        }

        //in this section, we are using the homogeneous equation to interpolate 
        //the vertex factor to our desired value of spectator momenta k 
        comp vertexres_at_q = {0.0,0.0};
        comp q = k;
        for(int i=0;i<qvec.size();++i)
        {
            comp kern = kernel_pk_2eps(s,q,qvec[i],a,m,eps,eps_for_m2k);
            comp weight = weights[i];
            comp norm_vert = normalized_vertex[i];
            vertexres_at_q = vertexres_at_q  -  weight*kern*norm_vert;

        }

        //in this section, we rescale the vertex factors squares and their m2k*vertexfactor
        //squares to the value found from the inverse d_S method

        comp sigk = sigma_p(spole,q,m);
        comp m2kres = M2kfunc(a,sigk,m,eps);

        comp rescaling_factor_of_dS_for_vertexsq = vertexsq_vec[i]/(vertexres_at_q*vertexres_at_q);
        comp rescaling_factor_of_dS_for_m2vertexsq = m2vertexsq_vec[i]/(m2kres*m2kres*vertexres_at_q*vertexres_at_q);
        
        
        cout<<"comparison between two rescaling factors = "<<rescaling_factor_of_dS_for_vertexsq<<" and "<<rescaling_factor_of_dS_for_m2vertexsq<<endl;
        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        
        
        comp final_m2vertex = sqrt(rescaling_factor_of_dS_for_m2vertexsq)*m2kres*vertexres_at_q;
        comp final_vertexsq = rescaling_factor_of_dS_for_vertexsq*vertexres_at_q*vertexres_at_q;
        comp final_m2vertexsq = rescaling_factor_of_dS_for_m2vertexsq*m2kres*m2kres*vertexres_at_q*vertexres_at_q;
        
        /*fout<<setprecision(16)<<k<<'\t'
            <<real(final_m2vertex)<<'\t'
            <<real(final_vertexsq)<<'\t'
            <<real(final_m2vertexsq)<<'\t'
            <<imag(final_m2vertex)<<'\t'
            <<imag(final_vertexsq)<<'\t'
            <<imag(final_m2vertexsq)<<'\t'
            <<abs(final_m2vertex)<<'\t'
            <<abs(final_vertexsq)<<'\t'
            <<abs(final_m2vertexsq)<<endl;
        */ // full output 

        fout<<setprecision(16)<<k<<'\t'
            <<real(final_m2vertex)<<'\t'    
            <<abs(final_m2vertex)<<'\t'
            <<abs(final_vertexsq)<<'\t'
            <<abs(final_m2vertexsq)<<endl;
        cout<<setprecision(16)<<"a="<<a<<'\n'
            <<"q="<<k<<'\n'
            <<"m2k= "<<m2kres<<'\n'
            <<"vert="<<vertexres_at_q<<'\n'
            <<"Vert HomEq = "<<real(final_vertexsq)<<'\t'
            <<"inv dS = "<<vertexsq_vec[i]<<'\t'
            <<"dif = "<<abs(real(final_vertexsq) - vertexsq_vec[i])<<'\n'
            <<"mVrtSq HomEq = "<<real(final_m2vertexsq)<<'\t'
            <<"inv dS = "<<m2vertexsq_vec[i]<<'\t'
            <<"dif = "<<abs(real(final_m2vertexsq) - m2vertexsq_vec[i])<<endl;



        

        comp traceval = {0.0,0.0};
        comp tracevaleig = {0.0,0.0};
        for(int i=0;i<size;++i)
        {
            traceval = traceval + Bmat(i,i);
            tracevaleig = tracevaleig + eigVal(i);
        }

        cout<<"trace from B = "<<traceval<<'\t'<<" sum of eigval = "<<tracevaleig<<endl;
        comp sigmax = sigma_p(spole,0.0,m);
        cout<<"sigmax = "<<sigmax<<endl;

        cout<<"----------------------------------------------"<<endl;
        cout<<endl;
    }
    fout.close();
}




void homogeneous_equation_integralequation_vs_k_solution_BS2()
{
    /*ifstream fin;
    string inputfile = "bs1st_N_500.dat";
    fin.open(inputfile.c_str());
    double some_a,some_b,some_c,some_d;
    vector<double> avec;
    vector<double> spolevec;
    int N;
    while(fin>>some_a>>some_b>>some_c>>some_d)
    {
        avec.push_back(some_a);
        spolevec.push_back(some_b);
        N = some_d;
    }*/

    //int N=500;
    //while(fin>>some_a>>some_b>>some_c>>some_d)
    /*for(double a=16.0;a<=50.0;a=a+1.0)
    {
        double spole = 0.0;
        //poletest_bs1_vs_singleN_singlea_using_determinant(a, N,spole);
        poletest_bs1_vs_singleN_singlea_using_weights_for_vertexfactor(a, N,spole);
        avec.push_back(a);
        spolevec.push_back(spole);
        
    }*/
    /*for(int i=0;i<avec.size();++i)
    {
        cout<<"pole found for a = "<<setprecision(16)<<avec[i]<<" at s = "<<spolevec[i]<<endl;
    }*/
    
    double a = 16.0;
    double m = 1.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    int N = 500;
    

    double spole = 8.976271055410812;//8.782848835105487;
    

    //int N = 10000;
    //poletest_bs1_vs_singleN_singlea_using_determinant(a, N,spole);
    //cout<<"pole found at s = "<<spole<<endl;

    
    double kmin = 0.000001;
    comp kmax = pmom(spole,0.0,m);//1.3;
    double kpoints = 2000.0;
    double delk = abs(kmax - kmin)/kpoints;
    //double delk = abs(kmax - kmin)/N;

    ofstream fout;
    string filename = "vertexfactor_vs_a_"+to_string(a)+"BS2_N_" + to_string(N) + "_smoothcutoff.dat";
    fout.open(filename.c_str());

    //for(int somei=0;somei<avec.size();++somei)
    //{
        //double a = avec[somei];
        //double spole = spolevec[somei];
    for(int i=0;i<kpoints;++i)
    {
        double k = kmin + i*delk;
        double s = spole;
        
        double points1 = N;
        vector<comp> qvec;
        vector<comp> weights;
        //mom_vector_maker_linear_1(qvec,kmin,delk,N,1);
        line_maker_with_weights(qvec,weights,kmin,kmax,points1);
        //mom_vector_maker_43_with_weights(qvec,weights,s,kmin,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                
        int size = qvec.size();

        Eigen::MatrixXcd Bmat(size,size);
        Eigen::MatrixXcd eigBmat(size,size);
        
        
        //Bmat_maker_momrep(Bmat,spole,qvec,qvec,a,m,eps);
        Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);


        //std::cout << "Here is the matrix A:\n" << Bmat << std::endl;
        //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(Bmat);
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigensolver(Bmat);
        
        if (eigensolver.info() != Eigen::Success) abort();
        //std::cout   << "The eigenvalues of A are:\n" << eigensolver.eigenvalues() << std::endl;
        //std::cout   << "Here's a matrix whose columns are eigenvectors of A \n"
        //            << "corresponding to these eigenvalues:\n"
        //            << eigensolver.eigenvectors() << std::endl;

        Eigen::VectorXcd eigVal(size);
        eigVal = eigensolver.eigenvalues();

        eigBmat = eigensolver.eigenvectors();
        //cout<<"eigBmat = "<<eigBmat<<endl;

        //for(int i=0;i<eigVal.size();++i)
        //{
            //cout<<"eigenvalues = "<<eigVal(i)<<endl;
        //} 

        //for(int i=0;i<eigVal.size();++i)
        //{
            //for(int j=0;j<eigVal.size();++j)
            //{
                //cout<<"eigvec "<<i<<","<<j<<" = "<<eigBmat(j,i)<<endl;
            //}
        //}
        comp sumvert = 0.0;
        comp sqrtsumvert = 0.0;
        //cout<<"0 column = "<<eigensolver.eigenvectors().col(0)<<endl;
        //cout<<"0 column printer = "<<endl;
        vector<double> m2kresvec;
        vector<double> vertexressqvec;
        vector<double> resvec;
        vector<comp> vertexresvec;
        //for(int i=0;i<eigVal.size();++i)
        //{
            int some_i = 0;
            for(int j=0;j<eigVal.size();++j)
            {
                //cout<<"eigvec "<<i<<","<<j<<" = "<<eigBmat(j,i)<<endl;
                comp k = qvec[j];
                
                comp sigk = sigma_p(spole,k,m);
                comp m2kres = M2kfunc(a,sigk,m,eps);
                comp m2kressq = m2kres*m2kres;
                comp vertexres = (eigBmat(j,some_i));
                comp vertexressq = vertexres*vertexres;
                comp res = m2kressq*vertexressq;
                vertexresvec.push_back(vertexres);
                m2kresvec.push_back(abs(m2kressq));
                vertexressqvec.push_back(abs(vertexressq));
                resvec.push_back(real(res));
                //fout<<real(qvec[j])<<'\t'<<real(m2kressq)<<'\t'<<real(vertexressq)<<'\t'<<real(res)<<endl;
                //cout<<j<<'\t'<<real(qvec[j])<<'\t'<<real(m2kressq)<<'\t'<<real(vertexressq)<<'\t'<<real(res)<<endl;
                if(j==0)
                {
                    //sumvert = sumvert + abs(vertexressq)*(qvec[j]);
                    sumvert = sumvert + abs(vertexressq)*(weights[j]);    
                    
                }
                else 
                {
                    //sumvert = sumvert + abs(vertexressq)*((qvec[j] - qvec[j-1]));
                    sumvert = sumvert + abs(vertexressq)*(weights[j]);
                
                }
                
            }
        //}
        cout<<"normalization constant = "<<setprecision(16)<<sumvert<<endl;
        sqrtsumvert = sqrt(sumvert);
        vector<comp> normalized_vertexres;
        for(int i=0;i<vertexresvec.size();++i)
        {
            normalized_vertexres.push_back(vertexresvec[i]/sqrtsumvert);
        }

        comp vertexres_at_q = {0.0,0.0};
        comp q = k;//1.0;//pmom(s,sigmab(a,m),m);
        for(int i=0;i<qvec.size();++i)
        {
            //comp q = 0.000001;//pmom(s,sigmab(a,m),m);
            comp kern = kernel_pk_2eps(s,q,qvec[i],a,m,eps,eps_for_m2k);
            comp weight = weights[i];
            comp norm_vert = normalized_vertexres[i];
            vertexres_at_q = vertexres_at_q  -  weight*kern*norm_vert;

        }
        //q = qvec[0];
        //vertexres_at_q = vertexresvec[0];

        comp sigk = sigma_p(spole,q,m);//sigmab(a,m);
        //comp q = pmom(s,sigmab(a,m),m);
        comp m2kres = M2kfunc(a,sigk,m,eps);
        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        fout<<setprecision(16)<<a<<'\t'
            <<real(vertexres_at_q)<<'\t'
            <<imag(vertexres_at_q)<<'\t'
            <<abs(vertexres_at_q)<<'\t'
            <<real(m2kres*vertexres_at_q)<<'\t'
            <<imag(m2kres*vertexres_at_q)<<'\t'
            <<abs(m2kres*vertexres_at_q*m2kres*vertexres_at_q)<<'\t'
            <<real(m2kres)<<'\t'
            <<imag(m2kres)<<'\t'
            //<<real(gfunc*vertexres_at_q)<<'\t'
            //<<imag(gfunc*vertexres_at_q)<<'\t'
            //<<abs(gfunc*vertexres_at_q*gfunc*vertexres_at_q)<<'\t'
            //<<real(gfunc)<<'\t'
            //<<imag(gfunc)<<'\t'
            <<real(q)<<'\t'
            <<imag(q)<<endl;

        cout<<setprecision(16)<<"a="<<a<<'\t'
            <<"g= "<<gfunc<<'\t'
            <<"m2k= "<<m2kres<<'\t'
            <<"vert="<<vertexres_at_q<<'\t'
            <<"m2kVert="<<m2kres*vertexres_at_q<<'\t'
            <<"mVrtSq="<<m2kres*vertexres_at_q*m2kres*vertexres_at_q<<'\t'
            <<"q="<<q<<endl;

        //for(int i=0;i<eigVal.size();++i)
        //{
        //    fout<<real(qvec[i])<<'\t'<<m2kresvec[i]<<'\t'<<vertexressqvec[i]/sumvert<<'\t'<<resvec[i]/sumvert<<endl;
        //    cout<<i<<'\t'<<real(qvec[i])<<'\t'<<m2kresvec[i]<<'\t'<<vertexressqvec[i]/sumvert<<'\t'<<resvec[i]/sumvert<<endl;
            
        //}

        

        comp traceval = {0.0,0.0};
        comp tracevaleig = {0.0,0.0};
        for(int i=0;i<size;++i)
        {
            traceval = traceval + Bmat(i,i);
            tracevaleig = tracevaleig + eigVal(i);
        }

        cout<<"trace from B = "<<traceval<<'\t'<<" sum of eigval = "<<tracevaleig<<endl;
        comp sigmax = sigma_p(spole,0.0,m);
        cout<<"sigmax = "<<sigmax<<endl;


    }
    //}
    fout.close();
}



void test_poles()
{
    double a = 10000.0;
    double m = 1.0;
    int N = 5000;
    double spole1, spole2, spole3;
    //poletest_bs1_vs_singleN_singlea_using_determinant(a, N, spole1);
    poletest_bs2_vs_singleN_singlea_using_determinant(a, N, spole2);
    poletest_bs3_vs_singleN_singlea_using_determinant(a, N, spole3);

    cout<<"am = "<<a<<endl;
    cout<<"N = "<<N<<endl;
    //cout<<"BS1 = "<<spole1<<endl;
    cout<<"BS2 = "<<spole2<<endl;
    cout<<"BS3 = "<<spole3<<endl;
}

void test_determinant()
{
    double a = 16.0;
    double m = 1.0;
    double sc = 8.88;
    int N = 5000;

    struct timeval stop, start;
    
    //do stuff
    
    //printf("took %lu us\n", (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec); 

    time_t time_start, time_end;
    //time(&time_start);
    gettimeofday(&start, NULL);
    cout<<dSqqs2q2msq_func_with_N_using_determinant(sc,a,N)<<endl;
    gettimeofday(&stop, NULL);
    double actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/100000.0;
    cout<<"time taken to solve = "<<actual_sol_time<<endl;
    //time(&time_end);
    //printf("took %lu us\n", (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec); 

    double time_taken = double((double)time_end - (double)time_start);
    //cout<<"Time taken to solve - cpu : " 
    //    <<setprecision(5)<<time_taken;
    //cout<<" sec"<<endl;

    //time(&time_start);
    gettimeofday(&start, NULL);
    cout<<dSqqs2q2msq_func_with_N_using_determinant_with_LUcuda(sc,a,N)<<endl;
    gettimeofday(&stop, NULL);
    actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/100000.0;
    cout<<"time taken to solve = "<<actual_sol_time<<endl;
    
    
    //printf("took %lu us\n", (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec); 

    //time(&time_end);
    //time_taken = double((double)time_end - (double)time_start);
    //cout<<"Time taken to solve - gpu : " 
    //    <<setprecision(5)<<time_taken;
    //cout<<" sec"<<endl;

}

void print_GL_energies()
{
    double m = 1.0;
    ofstream fout;
    string filename = "leftcut_energies_2.dat";
    fout.open( filename.c_str()   );

    double dela=0.01;
    for(double a=1.0001;a<=110005.0;a=a+dela)
    {
        if(a>=10) dela = 5;
        if(a>=100) dela = 20;
        if(a>=1000) dela = 100;

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));

        double phibth = real(phibthreshold(a,m));

        double en1 = sqrt(phibth) - sqrt(gvalleft);
        double en2 = sqrt(phibth) - sqrt(gvalright);
        fout<<setprecision(16)<<a<<'\t'<<en1<<'\t'<<en2<<endl;
        
    }
    fout.close();
}

//this is the function code used in forming Bmat 
//for the resonance pole search 

comp dSqqs2q2msq_resonance_func_with_N_using_weights(   comp s,
                                                        double a,
                                                        double N    )
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double m = 1.0;
    double eps_for_m2k = 0.0;
    double eps = 0.0;
    double shift = -0.1;
    double points1 = N;

    comp kmin = 0.0;
    comp kmax = pmom(s,0.0,m);

    comp sigb = 2.0*m*m;
    comp qval = pmom(s,sigb,m);

    vector<comp> qvec;
    vector<comp> weights;

    if(imag(s)<0.0)
    {
        if(real(s)>9.0)
        {
            contour_for_resonance_3(qvec,weights,s,m,kmin,kmax,eps_for_m2k,shift,points1);

        }
        else 
        line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                
    }
    else if(imag(s)>=0.0)
    {
        if(real(s)>9.0)
        {
            contour_for_resonance_3(qvec,weights,s,m,kmin,kmax,eps_for_m2k,shift,points1);

        }
        else 
        line_maker_with_weights(qvec,weights,kmin,kmax,points1);
            
    }

    cout<<"qvec created with size = "<<qvec.size()<<endl;
    int size = qvec.size();
    Eigen::MatrixXcd Bmat(size,size);
    Eigen::VectorXcd Gvec(size);
    Eigen::VectorXcd dsol(size);

    Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
    Gvec = -1.0*Gvec;

    if(real(s)>=9.0)
    {
        double pi = acos(-1.0);
        double theta = -90.0*pi/180.0;
        Bmat_maker_momrep_2eps_with_weights_rotated_m2k_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta);
        //Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

    }
    else
    Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

    struct timeval stop, start;
    gettimeofday(&start, NULL);
    comp result = Bmat.determinant();          
    gettimeofday(&stop, NULL);
    double actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/1000000.0;
    cout<<"time taken to solve = "<<actual_sol_time<<endl;


    return result;
    

}


//this is for the secant function for bound state pole
//searching 
comp dSqqs2q2msq_boundstates_func_with_N_using_weights(   comp s,
                                                        double a,
                                                        double N    )
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double m = 1.0;
    double eps_for_m2k = 0.0;
    double eps = 0.0;
    double shift = -0.1;
    double points1 = N;

    comp kmin = 0.0;
    comp kmax = pmom(s,0.0,m);

    comp sigb = 2.0*m*m;
    comp qval = pmom(s,sigb,m);

    vector<comp> qvec;
    vector<comp> weights;

    line_maker_with_weights(qvec,weights,kmin,kmax,points1);
    

    //cout<<"qvec created with size = "<<qvec.size()<<endl;
    int size = qvec.size();
    Eigen::MatrixXcd Bmat(size,size);
    Eigen::VectorXcd Gvec(size);
    Eigen::VectorXcd dsol(size);

    //Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
    //Gvec = -1.0*Gvec;

    Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

    struct timeval stop, start;
    gettimeofday(&start, NULL);
    comp result = Bmat.determinant();          
    gettimeofday(&stop, NULL);
    double actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/1000000.0;
    //cout<<"time taken to solve = "<<actual_sol_time<<endl;

    return result;
    

}


comp dSqqs2q2msq_resonance_func_with_N_using_weights_below_threebody_threshold(   comp s,
                                                        double a,
                                                        double N    )
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double m = 1.0;
    double eps_for_m2k = 0.0;
    double eps = 0.0;
    double shift = -0.1;
    double points1 = N;

    comp kmin = 0.0;
    comp kmax = pmom(s,0.0,m);

    comp sigb = 2.0*m*m;
    comp qval = pmom(s,sigb,m);

    vector<comp> qvec;
    vector<comp> weights;
    int tag_for_m2k;

    if(imag(s)<0.0)
    {
        /*if(real(s)>9.0)
        {
            contour_for_resonance_3(qvec,weights,s,m,kmin,kmax,eps_for_m2k,shift,points1);

        }
        else
        contour_for_resonance_4(qvec,weights,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k); 
        */
        //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
        contour_for_resonance_4(qvec,weights,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k); 
                
    }
    else if(imag(s)>=0.0)
    {
        /*if(real(s)>9.0)
        {
            contour_for_resonance_3(qvec,weights,s,m,kmin,kmax,eps_for_m2k,shift,points1);

        }
        else*/ 
        line_maker_with_weights(qvec,weights,kmin,kmax,points1);
            
    }

    //cout<<"qvec created with size = "<<qvec.size()<<endl;
    int size = qvec.size();
    Eigen::MatrixXcd Bmat(size,size);
    Eigen::VectorXcd Gvec(size);
    Eigen::VectorXcd dsol(size);

    //Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
    //Gvec = -1.0*Gvec;

    if(imag(s)<0.0)
    {
        double theta = 0.0;
        Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);
    
    }
    else if(imag(s)>=0.0)
    Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);


    /*if(real(s)>=9.0)
    {
        double pi = acos(-1.0);
        double theta = -90.0*pi/180.0;
        Bmat_maker_momrep_2eps_with_weights_rotated_m2k_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta);
        //Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

    }
    else if(real(s)<9.0 && imag(s)<=0.0)
    {
        double theta = 0.0;
        Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);
    }
    else 
    Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);
    */

    struct timeval stop, start;
    gettimeofday(&start, NULL);
    comp result = Bmat.determinant();          
    gettimeofday(&stop, NULL);
    double actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/1000000.0;
    //cout<<"time taken to solve = "<<actual_sol_time<<endl;


    return result;
    

}


//this is for poles lying on the upper part of the second sheet of dS
comp dSqqs2q2msq_resonance_func_with_N_using_weights_below_threebody_threshold_secondsheet_dS(   comp s,
                                                        double a,
                                                        double N    )
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double m = 1.0;
    double eps_for_m2k = 0.0;
    double eps = 0.0;
    double shift = -0.1;
    double points1 = N;

    comp kmin = 0.0;
    comp kmax = pmom(s,0.0,m);

    comp sigb = 2.0*m*m;
    comp qval = pmom(s,sigb,m);

    vector<comp> qvec;
    vector<comp> weights;
    int tag_for_m2k;
    double s3imag = imag(s);

    if(s3imag<0.0)
    {
        contour_for_resonance_4(qvec,weights,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k);

    }
    else if(s3imag>=0.0)
    {
        contour_for_resonance_6(qvec,weights,a,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k);
                
    }

    //cout<<"qvec created with size = "<<qvec.size()<<endl;
    int size = qvec.size();
    Eigen::MatrixXcd Bmat(size,size);
    Eigen::VectorXcd Gvec(size);
    Eigen::VectorXcd dsol(size);

    //Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
    //Gvec = -1.0*Gvec;

    if(imag(s)<0.0)
    {
        double pi = acos(-1.0);
        double theta = 0.0;//-90.0*pi/180.0;
        Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

    }
    else if(imag(s)>=0.0)
    {
        double pi = acos(-1.0);
        double theta = 0.0;//-90.0*pi/180.0;
        Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

    }

    /*if(real(s)>=9.0)
    {
        double pi = acos(-1.0);
        double theta = -90.0*pi/180.0;
        Bmat_maker_momrep_2eps_with_weights_rotated_m2k_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta);
        //Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

    }
    else if(real(s)<9.0 && imag(s)<=0.0)
    {
        double theta = 0.0;
        Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);
    }
    else 
    Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);
    */

    struct timeval stop, start;
    gettimeofday(&start, NULL);
    comp result = Bmat.determinant();          
    gettimeofday(&stop, NULL);
    double actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/1000000.0;
    //cout<<"time taken to solve = "<<actual_sol_time<<endl;


    return result;
    

}


//this for poles lying below real(s)<5.5
comp dSqqs2q2msq_resonance_func_with_N_using_weights_below_threebody_threshold_1(   comp s,
                                                        double a,
                                                        double N    )
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double m = 1.0;
    double eps_for_m2k = 0.0;
    double eps = 0.0;
    double shift = -0.1;
    double points1 = N;

    comp kmin = 0.0;
    comp kmax = pmom(s,0.0,m);

    comp sigb = 1.0*m*m;
    comp qval = pmom(s,sigb,m);

    vector<comp> qvec;
    vector<comp> weights;
    int tag_for_m2k;

    if(imag(s)<0.0)
    {
        /*if(real(s)>9.0)
        {
            contour_for_resonance_3(qvec,weights,s,m,kmin,kmax,eps_for_m2k,shift,points1);

        }
        else
        contour_for_resonance_4(qvec,weights,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k); 
        */
        //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
        contour_for_resonance_5(qvec,weights,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k); 
                
    }
    else if(imag(s)>=0.0)
    {
        /*if(real(s)>9.0)
        {
            contour_for_resonance_3(qvec,weights,s,m,kmin,kmax,eps_for_m2k,shift,points1);

        }
        else*/ 
        line_maker_with_weights(qvec,weights,kmin,kmax,points1);
            
    }

    cout<<"qvec created with size = "<<qvec.size()<<endl;
    int size = qvec.size();
    Eigen::MatrixXcd Bmat(size,size);
    Eigen::VectorXcd Gvec(size);
    Eigen::VectorXcd dsol(size);

    //Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
    //Gvec = -1.0*Gvec;

    if(imag(s)<0.0)
    {
        double theta = 0.0;
        Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);
    
    }
    else if(imag(s)>=0.0)
    Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);


    /*if(real(s)>=9.0)
    {
        double pi = acos(-1.0);
        double theta = -90.0*pi/180.0;
        Bmat_maker_momrep_2eps_with_weights_rotated_m2k_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta);
        //Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

    }
    else if(real(s)<9.0 && imag(s)<=0.0)
    {
        double theta = 0.0;
        Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);
    }
    else 
    Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);
    */

    struct timeval stop, start;
    gettimeofday(&start, NULL);
    comp result = Bmat.determinant();          
    gettimeofday(&stop, NULL);
    double actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/1000000.0;
    cout<<"time taken to solve = "<<actual_sol_time<<endl;


    return result;
    

}


//Here we use the secant method to find real or complex poles of dS
//We supply two initial guesses n-1 and n-2 which is used to calculate 
//the nth trial

void dSqqs2q2msq_pole_searching_secant_with_N_using_weights(    comp init_guess1,
                                                                comp init_guess2,
                                                                double a,
                                                                comp &pole, 
                                                                double N    )
{
    comp point2 = init_guess2;
    comp point1 = init_guess1; 
    int max_iteration = 100;
    double delta_s = 0.001;

    double epsilon = 1.0e-15;
    double tolerance = 1.0e-15;

    comp point0;

    for(int i=0;i<max_iteration;++i)
    {
        comp f1 = dSqqs2q2msq_resonance_func_with_N_using_weights(point1,a,N);
        comp f2 = dSqqs2q2msq_resonance_func_with_N_using_weights(point2,a,N);

        point0 = point1 - f1*(point1 - point2)/(f1 - f2);
        cout<<"iteration:"<<i<<endl;
        cout<<"point0 = "<<point0<<'\t'<<" kmax = "<<pmom(point0,0.0,1.0)<<endl;
        cout<<"point1 = "<<point1<<'\t'<<" kmax = "<<pmom(point1,0.0,1.0)<<endl;
        cout<<"point2 = "<<point2<<'\t'<<" kmax = "<<pmom(point2,0.0,1.0)<<endl;

        comp f0 = dSqqs2q2msq_resonance_func_with_N_using_weights(point0,a,N);
        cout<<"|f| = "<<abs(f0)<<endl;
        cout<<"========================================"<<endl;
        if(abs(f0)<=tolerance)
        {
            pole = point0;
            cout<<"pole found at, s = "<<pole<<endl;
            break;
        }
        else if(abs(real(point0-point1))<=epsilon && abs(imag(point0-point1))<=epsilon)
        {
            pole = point0; 
            cout<<"pole found at, s = "<<pole<<endl;
            break; 
        }

        if(i==max_iteration-1)
        {
            cout<<"did not converge to a solution, choose a diff guess or increase iteration"<<endl;
        }

        point2 = point1; 
        point1 = point0;
    }
}


//this is the secant method for finding poles of the three body 
//bound state. 
void dSqqs2q2msq_pole_searching_secant_with_N_using_weights_for_boundstates(    comp init_guess1,
                                                                comp init_guess2,
                                                                double a,
                                                                comp &pole, 
                                                                double N    )
{
    comp point2 = init_guess2;
    comp point1 = init_guess1; 
    int max_iteration = 100;
    double delta_s = 0.001;

    double epsilon = 1.0e-15;
    double tolerance = 1.0e-15;

    comp point0;

    for(int i=0;i<max_iteration;++i)
    {
        comp f1 = dSqqs2q2msq_boundstates_func_with_N_using_weights(point1,a,N);
        comp f2 = dSqqs2q2msq_boundstates_func_with_N_using_weights(point2,a,N);

        point0 = point1 - f1*(point1 - point2)/(f1 - f2);
        //cout<<"iteration:"<<i<<endl;
        //cout<<"point0 = "<<point0<<'\t'<<" kmax = "<<pmom(point0,0.0,1.0)<<endl;
        //cout<<"point1 = "<<point1<<'\t'<<" kmax = "<<pmom(point1,0.0,1.0)<<endl;
        //cout<<"point2 = "<<point2<<'\t'<<" kmax = "<<pmom(point2,0.0,1.0)<<endl;

        comp f0 = dSqqs2q2msq_boundstates_func_with_N_using_weights(point0,a,N);
        //cout<<"|f| = "<<abs(f0)<<endl;
        //cout<<"========================================"<<endl;
        if(abs(f0)<=tolerance)
        {
            pole = point0;
            //cout<<"pole found at, s = "<<pole<<endl;
            break;
        }
        else if(abs(real(point0-point1))<=epsilon && abs(imag(point0-point1))<=epsilon)
        {
            pole = point0; 
            //cout<<"pole found at, s = "<<pole<<endl;
            break; 
        }

        if(i==max_iteration-1)
        {
            cout<<"did not converge to a solution, choose a diff guess or increase iteration"<<endl;
        }

        point2 = point1; 
        point1 = point0;
    }
}


void dSqqs2q2msq_pole_searching_secant_with_N_using_weights_below_threebody_threshold(    comp init_guess1,
                                                                comp init_guess2,
                                                                double a,
                                                                comp &pole, 
                                                                double N    )
{
    comp point2 = init_guess2;
    comp point1 = init_guess1; 
    int max_iteration = 100;
    double delta_s = 0.001;

    double epsilon = 1.0e-15;
    double tolerance = 1.0e-15;

    comp point0;

    for(int i=0;i<max_iteration;++i)
    {
        comp f1 = dSqqs2q2msq_resonance_func_with_N_using_weights_below_threebody_threshold(point1,a,N);
        comp f2 = dSqqs2q2msq_resonance_func_with_N_using_weights_below_threebody_threshold(point2,a,N);

        point0 = point1 - f1*(point1 - point2)/(f1 - f2);
        //cout<<"iteration:"<<i<<endl;
        //cout<<"point0 = "<<point0<<'\t'<<" kmax = "<<pmom(point0,0.0,1.0)<<endl;
        //cout<<"point1 = "<<point1<<'\t'<<" kmax = "<<pmom(point1,0.0,1.0)<<endl;
        //cout<<"point2 = "<<point2<<'\t'<<" kmax = "<<pmom(point2,0.0,1.0)<<endl;

        comp f0 = dSqqs2q2msq_resonance_func_with_N_using_weights_below_threebody_threshold(point0,a,N);
        //cout<<"|f| = "<<abs(f0)<<endl;
        //cout<<"========================================"<<endl;
        if(abs(f0)<=tolerance)
        {
            pole = point0;
            //cout<<"pole found at, s = "<<pole<<endl;
            break;
        }
        else if(abs(real(point0-point1))<=epsilon && abs(imag(point0-point1))<=epsilon)
        {
            pole = point0; 
            //cout<<"pole found at, s = "<<pole<<endl;
            break; 
        }

        if(i==max_iteration-1)
        {
            cout<<"did not converge to a solution, choose a diff guess or increase iteration"<<endl;
        }

        point2 = point1; 
        point1 = point0;
    }
}


//this is for poles lying on the second sheet of dS
void dSqqs2q2msq_pole_searching_secant_with_N_using_weights_below_threebody_threshold_secondsheet_dS(    comp init_guess1,
                                                                comp init_guess2,
                                                                double a,
                                                                comp &pole, 
                                                                double N    )
{
    comp point2 = init_guess2;
    comp point1 = init_guess1; 
    int max_iteration = 100;
    double delta_s = 0.001;

    double epsilon = 1.0e-15;
    double tolerance = 1.0e-15;

    comp point0;

    for(int i=0;i<max_iteration;++i)
    {
        comp f1 = dSqqs2q2msq_resonance_func_with_N_using_weights_below_threebody_threshold_secondsheet_dS(point1,a,N);
        comp f2 = dSqqs2q2msq_resonance_func_with_N_using_weights_below_threebody_threshold_secondsheet_dS(point2,a,N);

        point0 = point1 - f1*(point1 - point2)/(f1 - f2);
        //cout<<"iteration:"<<i<<endl;
        //cout<<"point0 = "<<point0<<'\t'<<" kmax = "<<pmom(point0,0.0,1.0)<<endl;
        //cout<<"point1 = "<<point1<<'\t'<<" kmax = "<<pmom(point1,0.0,1.0)<<endl;
        //cout<<"point2 = "<<point2<<'\t'<<" kmax = "<<pmom(point2,0.0,1.0)<<endl;

        comp f0 = dSqqs2q2msq_resonance_func_with_N_using_weights_below_threebody_threshold_secondsheet_dS(point0,a,N);
        //cout<<"|f| = "<<abs(f0)<<endl;
        //cout<<"========================================"<<endl;
        if(abs(f0)<=tolerance)
        {
            pole = point0;
            //cout<<"pole found at, s = "<<pole<<endl;
            break;
        }
        else if(abs(real(point0-point1))<=epsilon && abs(imag(point0-point1))<=epsilon)
        {
            pole = point0; 
            //cout<<"pole found at, s = "<<pole<<endl;
            break; 
        }

        if(i==max_iteration-1)
        {
            cout<<"did not converge to a solution, choose a diff guess or increase iteration"<<endl;
        }

        point2 = point1; 
        point1 = point0;
    }
}


//this is for poles lying below real(s)<5.5
void dSqqs2q2msq_pole_searching_secant_with_N_using_weights_below_threebody_threshold_1(    comp init_guess1,
                                                                comp init_guess2,
                                                                double a,
                                                                comp &pole, 
                                                                double N    )
{
    comp point2 = init_guess2;
    comp point1 = init_guess1; 
    int max_iteration = 100;
    double delta_s = 0.001;

    double epsilon = 1.0e-15;
    double tolerance = 1.0e-15;

    comp point0;

    for(int i=0;i<max_iteration;++i)
    {
        comp f1 = dSqqs2q2msq_resonance_func_with_N_using_weights_below_threebody_threshold_1(point1,a,N);
        comp f2 = dSqqs2q2msq_resonance_func_with_N_using_weights_below_threebody_threshold_1(point2,a,N);

        point0 = point1 - f1*(point1 - point2)/(f1 - f2);
        cout<<"iteration:"<<i<<endl;
        cout<<"point0 = "<<point0<<'\t'<<" kmax = "<<pmom(point0,0.0,1.0)<<endl;
        cout<<"point1 = "<<point1<<'\t'<<" kmax = "<<pmom(point1,0.0,1.0)<<endl;
        cout<<"point2 = "<<point2<<'\t'<<" kmax = "<<pmom(point2,0.0,1.0)<<endl;

        comp f0 = dSqqs2q2msq_resonance_func_with_N_using_weights_below_threebody_threshold_1(point0,a,N);
        cout<<"|f| = "<<abs(f0)<<endl;
        cout<<"========================================"<<endl;
        if(abs(f0)<=tolerance)
        {
            pole = point0;
            cout<<"pole found at, s = "<<pole<<endl;
            break;
        }
        else if(abs(real(point0-point1))<=epsilon && abs(imag(point0-point1))<=epsilon)
        {
            pole = point0; 
            cout<<"pole found at, s = "<<pole<<endl;
            break; 
        }

        if(i==max_iteration-1)
        {
            cout<<"did not converge to a solution, choose a diff guess or increase iteration"<<endl;
        }

        point2 = point1; 
        point1 = point0;
    }
}


//here we test the pole searching secant code for resonances 
void poletest_res1_vs_N_using_weights()
{
    double m = 1.0;
    comp pole = 0.0;

    comp ii = {0.0,1.0};

    comp init_guess1 = 9.00000004 - 0.0000012*ii;
    comp init_guess2 = 9.00000004 - 0.0000016*ii; 

    double delta_s = 0.000001;
    //double a = -2.0; 
    double N = 750.0;

    ofstream fout;
    string filename1 = "res_3_pole_N_" + to_string((int)N) + "_1.dat";

    double ainitial = -1337;
    double afinal = -201;
    double apoints = (abs(ainitial-afinal));//150.0;
    double dela = abs(afinal - ainitial)/apoints;

    fout.open(filename1.c_str());
    for(int i=0;i<apoints+1;++i)
    {
        double a = ainitial + i*dela;
        dSqqs2q2msq_pole_searching_secant_with_N_using_weights(init_guess1, init_guess2, a, pole, N);

        fout<<setprecision(16)<<a<<'\t'<<real(pole)<<'\t'<<imag(pole)<<'\t'<<N<<endl;
        init_guess1 = pole;
        init_guess2 = pole - delta_s;

        cout<<"***********************************"<<endl;
        cout<<"a="<<a<<'\t'<<"pole = "<<real(pole)<<'\t'<<imag(pole)<<"i"<<'\t'<<"N="<<N<<endl;
        cout<<"***********************************"<<endl;
        
    }
    fout.close();


}


void poletest_res1_vs_N_using_weights_below_threebody_threshold()
{
    double m = 1.0;
    comp pole = 0.0;

    comp ii = {0.0,1.0};

    comp init_guess1 = 5.91	- 3.5920928*ii;
    comp init_guess2 = 5.91 - 3.5930928*ii; 

    double delta_s = 0.000001;
    //double a = -2.0; 
    double N = 500.0;

    ofstream fout;
    string filename1 = "res_1_pole_N_" + to_string((int)N) + "_5.dat";

    double ainitial = -0.552918;//-2.15;
    double afinal = -0.005;
    double apoints = 500.0;//(abs(ainitial-afinal));//150.0;
    double dela = abs(afinal - ainitial)/apoints;

    fout.open(filename1.c_str());
    for(int i=0;i<apoints+1;++i)
    {
        double a = ainitial + i*dela;
        dSqqs2q2msq_pole_searching_secant_with_N_using_weights_below_threebody_threshold(init_guess1, init_guess2, a, pole, N);

        fout<<setprecision(16)<<a<<'\t'<<real(pole)<<'\t'<<imag(pole)<<'\t'<<N<<endl;
        init_guess1 = pole;
        init_guess2 = pole - delta_s;

        cout<<"***********************************"<<endl;
        cout<<"a="<<a<<'\t'<<"pole = "<<real(pole)<<'\t'<<imag(pole)<<"i"<<'\t'<<"N="<<N<<endl;
        cout<<"***********************************"<<endl;
        
    }
    fout.close();


}


//this one locates the poles for the resonances entirely on the upper 
//second riemann sheet for dS. 

void poletest_res2_vs_N_using_weights_below_threebody_threshold()
{
    double m = 1.0;
    comp pole = 0.0;

    comp ii = {0.0,1.0};

    comp init_guess1 = 8.994376134878165	+ 0.00873595*ii;//9.0 + (5.75164e-10)*ii;
    comp init_guess2 = 8.994376134878165	+ 0.00873596*ii;//9.0 + (5.7517e-10)*ii; 

    double delta_s = 0.000001;
    //double a = -2.0; 
    double N = 500.0;

    ofstream fout;
    string filename1 = "res_2_pole_N_" + to_string((int)N) + "_sheet2_4.dat";

    double ainitial = -8.1;//-2.15;
    double afinal = -11.0;
    double apoints = 100.0;//(abs(ainitial-afinal));//150.0;
    double dela = abs(afinal - ainitial)/apoints;

    fout.open(filename1.c_str());
    for(int i=0;i<apoints+1;++i)
    {
        double a = ainitial - i*dela;
        dSqqs2q2msq_pole_searching_secant_with_N_using_weights_below_threebody_threshold_secondsheet_dS(init_guess1, init_guess2, a, pole, N);

        fout<<setprecision(16)<<a<<'\t'<<real(pole)<<'\t'<<imag(pole)<<'\t'<<N<<endl;
        init_guess1 = pole;
        init_guess2 = pole - delta_s;

        cout<<"***********************************"<<endl;
        cout<<"a="<<a<<'\t'<<"pole = "<<real(pole)<<'\t'<<imag(pole)<<"i"<<'\t'<<"N="<<N<<endl;
        cout<<"***********************************"<<endl;
        
    }
    fout.close();


}


//this is for 1st resonances whose pole lies below real(s)<5.5
void poletest_res1_vs_N_using_weights_below_threebody_threshold_1()
{
    double m = 1.0;
    comp pole = 0.0;

    comp ii = {0.0,1.0};

    comp init_guess1 = 9.000000040971532 -1.054537716476989e-06*ii;
    comp init_guess2 = 9.000000050971532 -1.054557716476989e-06*ii; 

    double delta_s = 0.000001;
    //double a = -2.0; 
    double N = 500.0;

    ofstream fout;
    string filename1 = "res_3_pole_N_" + to_string((int)N) + "_8.dat";

    double ainitial = -1337;		//-2.15;
    double afinal = -5000;
    double apoints = 50.0;//(abs(ainitial-afinal));//150.0;
    double dela = abs(afinal - ainitial)/apoints;

    fout.open(filename1.c_str());
    for(int i=0;i<apoints+1;++i)
    {
        double a = ainitial - i*dela;
        dSqqs2q2msq_pole_searching_secant_with_N_using_weights_below_threebody_threshold_1(init_guess1, init_guess2, a, pole, N);

        fout<<setprecision(16)<<a<<'\t'<<real(pole)<<'\t'<<imag(pole)<<'\t'<<N<<endl;
        init_guess1 = pole;
        init_guess2 = pole - delta_s;

        cout<<"***********************************"<<endl;
        cout<<"a="<<a<<'\t'<<"pole = "<<real(pole)<<'\t'<<imag(pole)<<"i"<<'\t'<<"N="<<N<<endl;
        cout<<"***********************************"<<endl;
        
    }
    fout.close();


}


//this is pole finder using secant method for specific a
//you change the guesses and get the 1st, 2nd or 3rd 
//bound state poles. 
void poletest_BS_vs_N_using_weights_secant_specific_a()
{
    double m = 1.0;
    comp pole = 0.0;

    comp ii = {0.0,1.0};

    comp init_guess1 = 9 - 5e-10;
    comp init_guess2 = 9 - 6e-10; 

    double delta_s = 0.000001;
    //double a = -2.0; 
    double N = 500.0;
    double a = 1000000;

    ofstream fout;
    string filename1 = "bs_4_pole_a_" + to_string(a) + ".dat";


    fout.open(filename1.c_str());
    dSqqs2q2msq_pole_searching_secant_with_N_using_weights_for_boundstates(init_guess1, init_guess2, a, pole, N);

    //a, pole, binding energy, N
    fout<<setprecision(16)<<a<<'\t'<<real(pole)<<'\t'<<abs(sqrt(pole) - 3.0)<<'\t'<<N<<endl;

    cout<<"***********************************"<<endl;
    cout<<setprecision(20)<<"a="<<a<<'\t'<<"pole = "<<real(pole)<<'\t'<<imag(pole)<<"i"<<'\t'<<"N="<<N<<endl;
    cout<<"***********************************"<<endl;
        
    fout.close();


}


//void D_second_sheet_test()
//{

//}


void F3_testing()
{
    double a = 1000000.0;
    double m = 1.0;
    comp KDF = {100.0,0.0};

    //double s = 8.65;
    double spole = 8.999975;//8.9995;//8.9246;//8.99999999916955;//8.999999440008509;//8.999983333329475;
    double delspole = real(phibthreshold(a,m)) - spole;
    
    double sinitial = 8.0;//8.999999939;//8.9999996;//8.99995;//8.9995;//8.999822636602779;//8.9246;//8.99998;//8.99998;//spole - 0.0000001;//8.999000;//spole;// - delspole;//8.95;//8.72;//real(phibthreshold(a,m));//
    double sfinal = 8.9;//real(phibthreshold(a,m));//8.99;
    double spoints = 1000.0;
    double sub_spoints = 100.0;
    
    double dels = abs(sinitial-sfinal)/spoints;
    double points1 = 500;//N;
    double points2 = 500.0;
    double eps = 0.0;
    //double box = 10.0;
    //int firstpart_for_nonunimesh = someA;
    //int secondpart_for_nonunimesh = someB;

    string filename="F3test.dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    double eta = 20.0;

    double dels_start = abs(sinitial-8.9)/sub_spoints;//1.0e-1;
    dels = dels_start; 
    int sub_spoint_counter = 0;
    comp s;
    double pow9 = 1;

    for(int i=0;i<=spoints;++i)
    //for(double s=sinitial;s<=sfinal;s=s+dels)
    { 
        //if(i!=0) break;
        sub_spoint_counter += 1;
        if(sub_spoint_counter==100)
        {
            pow9 = pow9 + 1;
            sinitial = real(s); 
            sfinal = sfinal + 9*pow(10,-pow9);
            dels = abs(sinitial-sfinal)/sub_spoints;
            sub_spoint_counter = 1;
        }
        //s = sinitial + (double)i*dels + ii*1.0e-15;
        s = sinitial + (double)sub_spoint_counter*dels + ii*1.0e-15;

        
        //double s = 8.99999985402895;
        //comp s = 8.99999985402895+1.0e-11*ii;//9.006-0.002*ii;//8.99999985402895 - ii*1.0e-16;

        //comp some_p = -0.8827566756850663+0.00027971206371478176*ii;
        //comp some_sigp = sigma_p(s,some_p,m);
        //cout<<"the sigk = "<<some_sigp<<endl;
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        //cout<<"phibth = "<<phibthreshold(a,m)<<endl;
        cout<<"am = "<<a<<endl;
        //cout<<spole<<endl;
        //cout<<spole-delspole<<endl;
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp sigq = 2.0*m*m;
        comp qval = pmom(s,sigq,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        double delk = abs(kmax-kmin)/points1;
        
        //cout<<"s:"<<s<<endl;
        //cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        //cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;
        vector<comp> weights;

        //if(firstpart_for_nonunimesh==0)
        //mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
        line_maker_with_weights(qvec,weights,kmin,kmax,points1);
        //for(int i=0;i<qvec.size();++i)
        //{
            //cout<<"i="<<i<<'\t'<<"q="<<qvec[i]<<'\t'<<"w="<<weights[i]<<endl;
        //}
        //else
        //mom_vector_maker_linear_2(qvec,kmin,kmax,points1,1,firstpart_for_nonunimesh,secondpart_for_nonunimesh);
        //mom_vector_maker_1(qvec,s,kmin,kmax,a,m,eps,points1);
    
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::MatrixXcd Gmat(size,size);
        Eigen::MatrixXcd dsol(size,size);

        //Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gmat_maker_momrep(Gmat,s,qvec,qvec,m,eps);
        //Gvec = -1.0*Gvec;
        Gmat = -1.0*Gmat;

        //Bmat_maker_momrep(Bmat,s,qvec,qvec,a,m,eps);
        Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps);
        comp detB = Bmat.determinant();

        double relerror;

        

        LinearSolver_3(Bmat,dsol,Gmat,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

        //relerror = abs((Bmat*dsol - Gvec).norm());
        comp result;

        //interpolator_ds_integraleq_momrep(dsol,qvec,s,qval,qval,a,m,eps,result);

        F3infvol(dsol,qvec,weights,qvec,weights,s,a,m,eps,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        //comp result = (s-spole)*result;
        //result = gsq*result;
        //result = rhophib(s,a,m)*result;
        comp denom = 1.0/result;// + 1.0/KDF;
        comp F3res = result; 
        comp minusF3invres = -1.0/F3res; 

        cout<<"s = "<<setprecision(16)<<s<<'\t'<<"dels ="<<dels<<endl;
        cout<<"res = "<<result<<'\t'<<" run = "<<count + 1<<endl;
        cout<<"detB = "<<detB<<'\t'<<" error = "<<relerror<<endl;
        cout<<setprecision(20)<<"F3 = "<<result<<endl;
        cout<<setprecision(20)<<"F3inv = "<<1.0/result<<endl;
        //cout<<"A = "<<someA<<'\t'<<"B = "<<someB<<endl;
        cout<<"N = "<<points1<<endl;
        fout<<setprecision(16)<<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<abs(real(sqrt(s) - 3.0))<<'\t'
            <<real(F3res)<<'\t'
            <<imag(F3res)<<'\t'
            <<real(minusF3invres)<<'\t'
            <<imag(minusF3invres)<<endl;//'\t'
            //<<real(detB)<<'\t'
            //<<imag(detB)<<'\t'
            //<<real(denom)<<'\t'
            //<<imag(denom)<<'\t'
            //<<relerror<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

}

void resonance_vs_real_s_sheet1()
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    int count = 0;

    double a = -8.0;
    double m = 1.0;
    double eps_for_m2k = 0.0;
    double eps = 0.0;
    double shift = -0.1;
    int N = 500.0;
    double points1 = N;

    double sinitial = 9.00000000000001;
    double sfinal = 9.03;
    double spoints = 500.0;
    double dels = abs(sinitial - sfinal)/spoints;

    string filename = "res_vs_s_for_a_" + to_string(a) + "_N_" + to_string(N) + "_sheet1.dat";
    ofstream fout;
    fout.open(filename.c_str());

    for(int i=0;i<=spoints;++i)
    {
        comp s = sinitial + i*dels - ii*0.002;

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);

        comp sigb = 2.0*m*m;
        comp qval = pmom(s,sigb,m);

        vector<comp> qvec;
        vector<comp> weights;

        line_maker_with_weights(qvec,weights,kmin,kmax,points1);

        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec = -1.0*Gvec;

       
        Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

        double relerror; 
        struct timeval stop, start;
        gettimeofday(&start, NULL);
        LinearSolver_2(Bmat,dsol,Gvec,relerror);
        
        gettimeofday(&stop, NULL);
        double actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/1000000.0;
        cout<<"time taken to solve = "<<actual_sol_time<<endl;

        comp result;

        interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);


           
        comp m2k = M2kfunc(a,sigb,m, eps_for_m2k);
        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        comp ds = result;
        result = m2k*result*m2k;
        comp Ds = result;
        cout<<setprecision(20)<<"s:"<<s<<'\t'<<"dS:"<<ds<<'\t'<<"run:"<<count + 1<<endl;
        cout<<setprecision(20)<<"Ds:"<<Ds<<endl;
        fout<<setprecision(20)<<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<real(Ds)<<'\t'
            <<imag(Ds)<<'\t'
            <<real(ds)<<'\t'
            <<imag(ds)<<'\t'
            <<a<<'\t'
            <<N<<'\t'
            <<eps<<endl;
        count = count + 1;
            
        cout<<"-------------------------"<<endl;
    }
    fout.close();
}



void resonance_vs_real_s_sheet2()
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    int count = 0;

    double a = -8.0;
    double m = 1.0;
    double eps_for_m2k = 0.0;
    double eps = 0.0;
    double shift = -0.1;
    int N = 500.0;
    double points1 = N;

    double sinitial = 9.00000000000001;
    double sfinal = 9.03;
    double spoints = 500.0;
    double dels = abs(sinitial - sfinal)/spoints;

    string filename = "res_vs_s_for_a_" + to_string(a) + "_N_" + to_string(N) + "_sheet2.dat";
    ofstream fout;
    fout.open(filename.c_str());

    int tag_for_m2k;

    for(int i=0;i<=spoints;++i)
    {
        comp s = sinitial + i*dels - ii*0.002;

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);

        comp sigb = 2.0*m*m;
        comp qval = pmom(s,sigb,m);

        vector<comp> qvec;
        vector<comp> weights;

        if(imag(s)<0.0)
        {
            if(real(s)>9.0)
            {
                //contour_for_resonance_3(qvec,weights,s,m,kmin,kmax,eps_for_m2k,shift,points1);
                contour_for_resonance_4(qvec,weights,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k);

            }
        else 
            line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                
        }
        else if(imag(s)>=0.0)
        {
            if(real(s)>9.0)
            {
                contour_for_resonance_3(qvec,weights,s,m,kmin,kmax,eps_for_m2k,shift,points1);

            }
            else 
            line_maker_with_weights(qvec,weights,kmin,kmax,points1);
            
        }
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec = -1.0*Gvec;
        
        if(real(s)>=9.0)
        {
            double pi = acos(-1.0);
            double theta = 0.0;//-90.0*pi/180.0;
            
            //Bmat_maker_momrep_2eps_with_weights_rotated_m2k_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta);
            Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

        }
        else
        Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

        double relerror; 
        struct timeval stop, start;
        gettimeofday(&start, NULL);
        LinearSolver_2(Bmat,dsol,Gvec,relerror);
        
        gettimeofday(&stop, NULL);
        double actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/1000000.0;
        cout<<"time taken to solve = "<<actual_sol_time<<endl;

        comp result;

            //if(real(s)>=9.0 && imag(s)<=0.0)
        if(real(s)>=9.0)
        {
            double pi = acos(-1.0);
            double theta = 0.0;//-90.0*pi/180.0;
            //m2k_cut_tagger(qvec,s,m,tag_for_m2k);
            //interpolator_ds_integraleq_momrep_2eps_with_weights_rotated_m2k_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,result);
            interpolator_ds_integraleq_momrep_2eps_with_weights_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,tag_for_m2k,result);
            
        
        }
        else 
            interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
           
        comp m2k = M2kfunc(a,sigb,m, eps_for_m2k);
        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        comp ds = result;
        result = m2k*result*m2k;
        comp Ds = result;
        cout<<setprecision(20)<<"s:"<<s<<'\t'<<"dS:"<<ds<<'\t'<<"run:"<<count + 1<<endl;
        cout<<setprecision(20)<<"Ds:"<<Ds<<endl;
        fout<<setprecision(20)<<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<real(Ds)<<'\t'
            <<imag(Ds)<<'\t'
            <<real(ds)<<'\t'
            <<imag(ds)<<'\t'
            <<a<<'\t'
            <<N<<'\t'
            <<eps<<endl;
        count = count + 1;
            
        cout<<"-------------------------"<<endl;
    }
    fout.close();
}

void Mphib_analytic_continuation_results_with_weights()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = 8.1;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    double s3imag = -1.0e-5;
    

    double delspoints = 500.0;
    double points1 = 500.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    //double box = 10.0;
    int acount = 0;

    double dela = 0.051;
    //for(double a=2.401;a<=16.03;a=a+dela)
    vector<double> epsvec(5);
    epsvec[0] = 0.000010;
    epsvec[1] = 0.000007;
    epsvec[2] = 0.000005;
    epsvec[3] = 0.000003;
    epsvec[4] = 0.000001;

    //for(int epsi=0;epsi<epsvec.size();++epsi)
    //{

    //    eps = epsvec[epsi];

    //int N = 500;
    //for(int N=900;N<=1000;N=N+20)
    {
        //if(a>=3.7) dela = 0.101;
        //else if(a>=17.0) dela = 0.51;
        //else if(a>=1000.0) dela = 300.01;

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        double sinitial = 8.0;//6.5;//8.3;//5.0;//8.974;//8.6;//8.975;//8.96815;//7.5;//8.7;//7.5;//8.75;//1.001;//gvalright - 0.5;//(gvalleft + (gvalleft+gvalright)/2.0)/2.0;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
        double sfinal = 9.0;//real(phibthreshold(a,m));//8.82;
        double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
        points1 = N;
        //string filename="Mphib_secondsheet_acount_" + to_string(acount) + ".dat";
        string filename="Mphib_a_" + to_string(a) + "_N_" + to_string((int)N) + "_1.dat";
        
        //string filename="testGL_Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_2.dat";
        //string filename="tmp.dat";
        

        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        //double eta = 20.0;

        for(int i=0;i<delspoints+1;++i) 
        //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
        { 
            double s3real = sinitial + i*dels;
            //comp s = 8.6 + ii*0.00001;
            //if(s>=sfinal) dels = 0.05/100.0;
            //eps = eps_above_threshold(eta,s,a,m,points);
            //s3real = 8.75; //this was added
            //s3imag = 0.00001; //this was added 

            comp s = s3real + ii*s3imag;
            cout<<"am = "<<a<<endl;
            //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = sigmab(a,m);
            comp qval = pmom(s,sigb,m);
            cout<<"q = "<<qval<<endl;

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        

            //cout<<"s:"<<s<<endl;
            cout<<"run: "<<count+1<<endl;
            cout<<"a = "<<a<<'\t'<<" dela = "<<dela<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            
            int switch_for_gvec_fixer = 1; //0 means we have used contour and turn on gvec fixer
                                           //1 means we have used straight line and turn off gvec fixer
            
            if(s3imag<0.0)
            {
                //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                cout<<"contour_43 chosen"<<endl;
                //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                mom_vector_maker_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1);
                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
            }
            else if(s3imag>=0.0)
            {
                //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                cout<<"contour_47 chosen"<<endl;
                //mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
                int changed_points = points1 + 5;
                vector<comp> tmp_qvec;
                vector<comp> tmp_weights;
                //mom_vector_maker_seba_imspos(tmp_qvec,tmp_weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)changed_points,tag1,tag2,switch_for_gvec_fixer);
                //mom_vector_maker_seba_imspos(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                //mom_vector_maker_47_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
                //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                
                //mom_vector_maker_seba_imspos_2(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer); //change this
            
                //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                
                /*
                mom_vector_maker_seba_imsneg(tmp_qvec,tmp_weights,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1);
                for(int i=0;i<tmp_qvec.size();++i)
                {
                    comp qv = real(tmp_qvec[i]) - ii*imag(tmp_qvec[i]);
                    comp w = real(tmp_weights[i]) - ii*imag(tmp_weights[i]);
                    qvec.push_back(qv);
                    weights.push_back(w);
                }
                cout<<"tag1 = "<<tag1<<endl;
                cout<<"tag2 = "<<tag2<<endl;
                */

                //mom_vector_maker_seba_imspos_5(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                mom_vector_maker_seba_imspos_5_ext1(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)points1,tag1,tag2,switch_for_gvec_fixer);
        
                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                
                //this portion is for loading ready made contours:
                /*ifstream fin;
                string contour_file = "for_digonto.txt";
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
                qvec = contour;*/
                
                
            }

            cout<<"qvec created with size = "<<qvec.size()<<endl;
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            //Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            if(s3imag<0.0)
            {
                Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            }
            else if(s3imag>=0.0)
            {
                Gvec_maker_momrep_withtags_1(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            }

            //switch_for_gvec_fixer=1;
            if(switch_for_gvec_fixer==0)
            {
                //Gvec_fixer_1(Gvec,Gvec,qvec,s,qval,m,eps,tag1,tag2); 
                //Gvec_fixer_2(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
                
                //this two were used for the paper calculations
                //Gvec_fixer_3(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
                //tag_fixer_for_gvec_fixer3(Gvec,qvec,s,qval,m,eps,tag1);
                //--------------------------------------------//


                Gvec_fixer_5(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
                //using this will make the tags redundant if we pass the 
                //fixed Gvec to the interpolator
                tag_fixer_for_gvec_fixer5(Gvec,qvec,s,qval,m,eps,tag1);
                

                //this did not work for im(s)>0 with hard cutoff

                //Gvec_fixer_6(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
               
                
            }

            cout<<"tag fixer = "<<tag1<<'\t'<<tag2<<endl;
            
            
            //sebatest part ---------
            ofstream fout2;//("opefile.dat");
            string opefile = "opefile.dat";
            fout2.open(opefile.c_str());
            ofstream fout3; 
            string opefile1 = "opefile1.dat";

            fout3.open(opefile1.c_str());

            for(int i=0;i<qvec.size();++i)
            {
                fout2<<i<<'\t'<<real(Gvec[i])<<'\t'<<imag(Gvec[i])<<endl;
                fout3<<i<<'\t'<<real(GS_pk(s,qvec[i],qval,m,eps))<<'\t'<<imag(GS_pk(s,qvec[i],qval,m,eps))<<'\t'<<real(GS_pk_secondsheet(s,qvec[i],qval,m,eps))<<'\t'<<imag(GS_pk_secondsheet(s,qvec[i],qval,m,eps))<<endl;
            }
            fout2.close();
            fout3.close();
            
            //-----------------------
            
            //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            Eigen::VectorXcd interpolater_Gvec = Gvec;
            Gvec = -1.0*Gvec;

            //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
            //Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
            //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);
            //Bmat_maker_momrep_2eps_withtags_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,tag1,tag2);
            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

            //cout<<"determinant of B = "<<Bmat.determinant()<<endl;
            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        
            //cout<<"did this = "<<count<<" times"<<endl;
            //cout<<"i val = "<<i<<'\t';
            //cout<<"j val = "<<j<<endl;
            
            //LinearSolver_2(Bmat,dsol,Gvec,relerror);

            time(&time_start);
            LinearSolver_2(Bmat,dsol,Gvec,relerror);

            //cusolverComplex(Bmat,Gvec,dsol,size);
            time(&time_end);

            double time_taken = double(time_end - time_start);
            cout<<"Time taken to solve : "<<fixed 
                <<time_taken<<setprecision(5);
            cout<<" sec"<<endl;


            comp result;

            //interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            
            if(s3imag<0.0)
            {
                //interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
                interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            }
            else if(s3imag>=0.0)
            {
                //interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
                //interpolator_ds_integraleq_momrep_2eps_withtags_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
                
                //the tags become redundant here//
                interpolator_ds_integraleq_momrep_2eps_withtags_with_weights_usingGvec(dsol,qvec,weights,interpolater_Gvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
            
            }
            //result = dsol[10];
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            comp ds = result;
            result = gsq*result;

        
            comp rhopb = rhophib(s,a,m);
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;

            cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            fout<<setprecision(16)<<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(ds)<<'\t'
                <<imag(ds)<<'\t'
                <<real(mphib2)<<'\t'
                <<imag(mphib2)<<'\t'
                <<real(mphib2denom)<<'\t'
                <<imag(mphib2denom)<<'\t'
                <<a<<'\t'
                <<N<<'\t'
                <<gvalleft<<'\t'
                <<gvalright<<'\t'
                <<real(phibthreshold(a,m))<<'\t'
                <<eps<<endl;
            count = count + 1;
            
            cout<<"-------------------------"<<endl;
        }
        acount = acount + 1;
        fout.close();

    }

    //}

}

void ds_qq2msq_analytic_continuation_results_with_weights()
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = 8.1;
    double m = 1.0;

    cout<<"am = "<<a<<endl;
    //double s = 8.65;

    double s3imag = 0.0;//1.0e-5;
    

    double delspoints = 500.0;
    double points1 = 500.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    //double box = 10.0;
    int acount = 0;

    double dela = 0.051;
    //for(double a=2.401;a<=16.03;a=a+dela)
    vector<double> epsvec(5);
    epsvec[0] = 0.000010;
    epsvec[1] = 0.000007;
    epsvec[2] = 0.000005;
    epsvec[3] = 0.000003;
    epsvec[4] = 0.000001;

    //for(int epsi=0;epsi<epsvec.size();++epsi)
    //{

    //    eps = epsvec[epsi];

    //int N = 500;
    //for(int N=900;N<=1000;N=N+20)
    {
        //if(a>=3.7) dela = 0.101;
        //else if(a>=17.0) dela = 0.51;
        //else if(a>=1000.0) dela = 300.01;

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        double sinitial = 8.0;//6.5;//8.3;//5.0;//8.974;//8.6;//8.975;//8.96815;//7.5;//8.7;//7.5;//8.75;//1.001;//gvalright - 0.5;//(gvalleft + (gvalleft+gvalright)/2.0)/2.0;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
        double sfinal = 9.0;//real(phibthreshold(a,m));//8.82;
        double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
        points1 = N;
        //string filename="Mphib_secondsheet_acount_" + to_string(acount) + ".dat";
        string filename="dsqq_a_" + to_string(a) + "_N_" + to_string((int)N) + ".dat";
        
        //string filename="testGL_Mphib_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_2.dat";
        //string filename="tmp.dat";
        

        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

        //double eta = 20.0;

        for(int i=0;i<delspoints+1;++i) 
        //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
        { 
            double s3real = sinitial + i*dels;
            //comp s = 8.6 + ii*0.00001;
            //if(s>=sfinal) dels = 0.05/100.0;
            //eps = eps_above_threshold(eta,s,a,m,points);
            //s3real = 8.75; //this was added
            //s3imag = 0.00001; //this was added 

            comp s = s3real + ii*s3imag;
            cout<<"am = "<<a<<endl;
            //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = 2.0;//sigmab(a,m);
            comp qval = pmom(s,sigb,m);
            cout<<"q = "<<qval<<endl;

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        

            //cout<<"s:"<<s<<endl;
            cout<<"run: "<<count+1<<endl;
            cout<<"a = "<<a<<'\t'<<" dela = "<<dela<<endl;
            cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            
            int switch_for_gvec_fixer = 0; //0 means we have used contour and turn on gvec fixer
                                           //1 means we have used straight line and turn off gvec fixer
            
            if(s3imag<0.0)
            {
                //mom_vector_maker_41(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                cout<<"contour_43 chosen"<<endl;
                //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                mom_vector_maker_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1);
                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
            }
            else if(s3imag>=0.0)
            {
                //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                cout<<"contour_47 chosen"<<endl;
                //mom_vector_maker_47(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
                int changed_points = points1 + 5;
                vector<comp> tmp_qvec;
                vector<comp> tmp_weights;
                //mom_vector_maker_seba_imspos(tmp_qvec,tmp_weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)changed_points,tag1,tag2,switch_for_gvec_fixer);
                //mom_vector_maker_seba_imspos(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                //mom_vector_maker_47_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r,tag1,tag2);
                //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                
                //mom_vector_maker_seba_imspos_2(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer); //change this
            
                //mom_vector_maker_seba_imspos_1(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                
                /*
                mom_vector_maker_seba_imsneg(tmp_qvec,tmp_weights,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1);
                for(int i=0;i<tmp_qvec.size();++i)
                {
                    comp qv = real(tmp_qvec[i]) - ii*imag(tmp_qvec[i]);
                    comp w = real(tmp_weights[i]) - ii*imag(tmp_weights[i]);
                    qvec.push_back(qv);
                    weights.push_back(w);
                }
                cout<<"tag1 = "<<tag1<<endl;
                cout<<"tag2 = "<<tag2<<endl;
                */

                //mom_vector_maker_seba_imspos_5(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)points1,tag1,tag2,switch_for_gvec_fixer);
        
                line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                
                //this portion is for loading ready made contours:
                /*ifstream fin;
                string contour_file = "for_digonto.txt";
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
                qvec = contour;*/
                
                
            }

            cout<<"qvec created with size = "<<qvec.size()<<endl;
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            /*if(s3imag<0.0)
            {
                Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            }
            else if(s3imag>=0.0)
            {
                Gvec_maker_momrep_withtags_1(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            }*/

            switch_for_gvec_fixer=1;
            if(switch_for_gvec_fixer==0)
            {
                //Gvec_fixer_1(Gvec,Gvec,qvec,s,qval,m,eps,tag1,tag2); 
                //Gvec_fixer_2(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
                
                //this two were used for the paper calculations
                //Gvec_fixer_3(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
                //tag_fixer_for_gvec_fixer3(Gvec,qvec,s,qval,m,eps,tag1);
                //--------------------------------------------//


                Gvec_fixer_5(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
                //using this will make the tags redundant if we pass the 
                //fixed Gvec to the interpolator
                tag_fixer_for_gvec_fixer5(Gvec,qvec,s,qval,m,eps,tag1);
                

                //this did not work for im(s)>0 with hard cutoff

                //Gvec_fixer_6(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
               
                
            }

            cout<<"tag fixer = "<<tag1<<'\t'<<tag2<<endl;
            
            
            //sebatest part ---------
            /*ofstream fout2;//("opefile.dat");
            string opefile = "opefile.dat";
            fout2.open(opefile.c_str());
            ofstream fout3; 
            string opefile1 = "opefile1.dat";

            fout3.open(opefile1.c_str());

            for(int i=0;i<qvec.size();++i)
            {
                fout2<<i<<'\t'<<real(Gvec[i])<<'\t'<<imag(Gvec[i])<<endl;
                fout3<<i<<'\t'<<real(GS_pk(s,qvec[i],qval,m,eps))<<'\t'<<imag(GS_pk(s,qvec[i],qval,m,eps))<<'\t'<<real(GS_pk_secondsheet(s,qvec[i],qval,m,eps))<<'\t'<<imag(GS_pk_secondsheet(s,qvec[i],qval,m,eps))<<endl;
            }
            fout2.close();
            fout3.close();
            */
            //-----------------------
            
            //Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            Eigen::VectorXcd interpolater_Gvec = Gvec;
            Gvec = -1.0*Gvec;

            //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
            //Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
            //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);
            //Bmat_maker_momrep_2eps_withtags_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,tag1,tag2);
            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

            //cout<<"determinant of B = "<<Bmat.determinant()<<endl;
            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        
            //cout<<"did this = "<<count<<" times"<<endl;
            //cout<<"i val = "<<i<<'\t';
            //cout<<"j val = "<<j<<endl;
            
            //LinearSolver_2(Bmat,dsol,Gvec,relerror);

            time(&time_start);
            LinearSolver_2(Bmat,dsol,Gvec,relerror);

            //cusolverComplex(Bmat,Gvec,dsol,size);
            time(&time_end);

            double time_taken = double(time_end - time_start);
            cout<<"Time taken to solve : "<<fixed 
                <<time_taken<<setprecision(5);
            cout<<" sec"<<endl;


            comp result;

            interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            
            /*if(s3imag<0.0)
            {
                //interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
                interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            }
            else if(s3imag>=0.0)
            {
                //interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
                //interpolator_ds_integraleq_momrep_2eps_withtags_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
                
                //the tags become redundant here//
                interpolator_ds_integraleq_momrep_2eps_withtags_with_weights_usingGvec(dsol,qvec,weights,interpolater_Gvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
            
            }*/
            //result = dsol[10];
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            comp ds = result;
            result = result;

        
            comp rhopb = rhophib(s,a,m);
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;

            cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            fout<<setprecision(16)<<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(ds)<<'\t'
                <<imag(ds)<<'\t'
                <<real(mphib2)<<'\t'
                <<imag(mphib2)<<'\t'
                <<real(mphib2denom)<<'\t'
                <<imag(mphib2denom)<<'\t'
                <<a<<'\t'
                <<N<<'\t'
                <<gvalleft<<'\t'
                <<gvalright<<'\t'
                <<real(phibthreshold(a,m))<<'\t'
                <<eps<<endl;
            count = count + 1;
            
            cout<<"-------------------------"<<endl;
        }
        acount = acount + 1;
        fout.close();

    }

    //}

}


void Mphib_around_pole_with_weights(    double a,
                                        comp sinitial,
                                        comp sfinal,
                                        double delspoints,
                                        comp polelocation,
                                        string filename     )
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    //double a = 16;
    double m = 1.0;
    double simag = -1.0e-5;

    double points1 = 500.0;
    double N = points1;
    //double points2 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    //double box = 10.0;
    int acount = 0;

    {

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        
        comp dels = abs(sinitial-sfinal)/delspoints;
        points1 = N;
          

        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;


        for(int i=0;i<delspoints+1;++i) 
        { 
            comp s = sinitial + ((comp)i)*dels + ii*simag;
            //cout<<"am = "<<a<<endl;
            //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = sigmab(a,m);
            comp qval = pmom(s,sigb,m);
            //cout<<"q = "<<qval<<endl;

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        

            //cout<<"s:"<<s<<endl;
            //cout<<"run: "<<count+1<<endl;
            //cout<<"a = "<<a<<'\t'<<" dela = "<<dela<<endl;
            //cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            //cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            
            int switch_for_gvec_fixer = 1; //0 means we have used contour and turn on gvec fixer
                                           //1 means we have used straight line and turn off gvec fixer
            
            if(imag(s)<0.0)
            {
                mom_vector_maker_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1);
            }
            else if(imag(s)>=0.0)
            {
                mom_vector_maker_seba_imspos_5_ext1(qvec,weights,s,0.0,kmax,a,m,eps,eps,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                
            }

            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            
            if(imag(s)<0.0)
            {
                Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            }
            else if(imag(s)>=0.0)
            {
                Gvec_maker_momrep_withtags_1(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            }

            if(switch_for_gvec_fixer==0)
            {
                Gvec_fixer_5(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
                //using this will make the tags redundant if we pass the 
                //fixed Gvec to the interpolator
                tag_fixer_for_gvec_fixer5(Gvec,qvec,s,qval,m,eps,tag1);
                

                //this did not work for im(s)>0 with hard cutoff

                //Gvec_fixer_6(Gvec,Gvec,qvec,s,qval,a,m,eps,0.0,tag1,tag2);
               
                
            }

            //-----------------------
            
            Eigen::VectorXcd interpolater_Gvec = Gvec;
            Gvec = -1.0*Gvec;

            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

            double relerror;
            
        

            time(&time_start);
            LinearSolver_2(Bmat,dsol,Gvec,relerror);

            //cusolverComplex(Bmat,Gvec,dsol,size);
            time(&time_end);

            double time_taken = double(time_end - time_start);
            //cout<<"Time taken to solve : "<<fixed 
            //    <<time_taken<<setprecision(5);
            //cout<<" sec"<<endl;


            comp result;

        
            if(imag(s)<0.0)
            {
                //interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
                interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            }
            else if(imag(s)>=0.0)
            {
                //interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
                //interpolator_ds_integraleq_momrep_2eps_withtags_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
                
                //the tags become redundant here//
                interpolator_ds_integraleq_momrep_2eps_withtags_with_weights_usingGvec(dsol,qvec,weights,interpolater_Gvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
            
            }
            //result = dsol[10];
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            comp ds = result;
            result = gsq*result;

        
            comp rhopb = rhophib(s,a,m);
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;

            //cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
            fout<<setprecision(16)<<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(1.0/result)<<'\t'
                <<imag(1.0/result)<<'\t'
                <<real(gsq)<<'\t'
                <<imag(gsq)<<'\t'
                <<real(polelocation)<<'\t'
                <<imag(polelocation)<<'\t'
                <<a<<'\t'
                <<N<<'\t'
                <<gvalleft<<'\t'
                <<gvalright<<'\t'
                <<real(phibthreshold(a,m))<<'\t'
                <<eps<<endl;
            count = count + 1;
            
            //cout<<"-------------------------"<<endl;
        }
        acount = acount + 1;
        fout.close();

    }

    //}

}

void check_gamma_phib_vs_a()
{
    //first find the poles 
    double init_guess1 = 8.96;//2.87;
    double init_guess2 = 8.99999999;//2.9;
    double m = 1.0;
    double N = 500;
    comp pole;
    double ainitial = 13.03;
    double afinal = 1000000;
    double apoints = 2000;
    double dela = (afinal - ainitial)/apoints;//(1.0 - exp(apoints-1.0));
    dela = 0.0;//1.01;
    double a = ainitial;
    //for(int i=0;i<apoints;++i)
    for(int i=0; ;++i)
    //int i = 2;
    {
        //double a = ainitial + (1.0-exp(i))*dela;
        //double a = ainitial + i*dela;
        a = a + dela;

        if(a>afinal) break; 
        if(a>=ainitial && a<10.0) dela = 0.2501;
        else if(a>10.0 && a<100.0) dela = 2.51;
        else if(a>100.0 && a<1000.0) dela = 25.101;
        else if(a>1000.0 && a<10000.0) dela = 251.01;
        else if(a>10000.0 && a<100000.0) dela = 1257.01;
        else if(a>100000.0 && a<1000000.0) dela = 2501.01;
        cout<<"next dela = "<<dela<<endl;

        if(i==0)
        {
            init_guess1 = 8.96;
            init_guess2 = real(phibthreshold(a,m));
            
        }
        else if(i>0)
        {
            init_guess1 = real(pole) - 0.0001;
            init_guess2 = real(phibthreshold(a,m));//real(pole);// + 0.00005;
        }

        cout<<"running for am = "<<a<<endl;
        dSqqs2q2msq_pole_searching_secant_with_N_using_weights_for_boundstates(init_guess1,init_guess2,a,pole,N);

        comp sinitial = pole - 0.000001;
        comp sfinal = pole + 0.000001;
        comp polelocation = pole;
        double delspoints = 10.0;
        string filename = "BS2_res_acount_" + to_string(i) + ".dat";
        
        /*string file_for_d = "dtest.dat";
        ofstream fout;
        fout.open(file_for_d.c_str());
        comp dels = (sfinal - sinitial)/delspoints;
        for(int i=0;i<delspoints;++i)
        {
            comp s = sinitial + ((comp)i)*dels; 
            comp dres = dSqqs2q2msq_boundstates_func_with_N_using_weights(s,a,N);
            fout<<real(s)<<'\t'<<imag(s)<<'\t'<<real(dres)<<'\t'<<imag(dres)<<endl;

        }
        fout.close();
        */
        Mphib_around_pole_with_weights(a,sinitial,sfinal,delspoints,polelocation,filename);
        cout<<"found pole at s = "<<setprecision(16)<<pole<<endl;
        cout<<"phib threshold = "<<phibthreshold(a,m)<<endl;
        cout<<"file = "<<filename<<" created"<<endl;
        cout<<"======================================"<<endl;
    }
    
}

//residue matching portion for BS1, res1, and res2
//this can also serve as a datagenerator for residue 
//or vertex factors, works the same way

void ds_qq2msq_BS_dataprinter_for_residue_matching(    string filename,
                                                    double a,
                                                    comp sinitial, 
                                                    comp sfinal,
                                                    double delspoints,
                                                    int N    )
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double m = 1.0;


    double s3imag = 0.0;//1.0e-5;
    

    double points1 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    {

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
        points1 = N;

        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;


        for(int i=0;i<delspoints+1;++i) 
        { 
            double s3real = real(sinitial) + i*dels;

            comp s = s3real + ii*s3imag;
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = 2.0;//sigmab(a,m);
            comp qval = pmom(s,sigb,m);

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        


            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            
            int switch_for_gvec_fixer = 0; //0 means we have used contour and turn on gvec fixer
                                           //1 means we have used straight line and turn off gvec fixer
            
            if(s3imag<0.0)
            {
                line_maker_with_weights(qvec,weights,kmin,kmax,points1);
            }
            else if(s3imag>=0.0)
            {
                line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                
            }

            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            
            Eigen::VectorXcd interpolater_Gvec = Gvec;
            Gvec = -1.0*Gvec;

            Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

            double relerror;

            time(&time_start);
            LinearSolver_2(Bmat,dsol,Gvec,relerror);
            time(&time_end);

            double time_taken = double(time_end - time_start);

            comp result;

            interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
        
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            comp m2k = M2kfunc(a,sigb,m,eps);
            comp ds = result;
            result = m2k*result*m2k;

        
            comp rhopb = rhophib(s,a,m);
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;

            fout<<setprecision(16)<<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(ds)<<'\t'
                <<imag(ds)<<'\t'
                <<real(1.0/result)<<'\t'
                <<imag(1.0/result)<<'\t'
                <<real(1.0/ds)<<'\t'
                <<imag(1.0/ds)<<'\t'
                <<a<<'\t'
                <<N<<'\t'
                <<gvalleft<<'\t'
                <<gvalright<<'\t'
                <<real(phibthreshold(a,m))<<'\t'
                <<eps<<endl;
            count = count + 1;
            
        }
        fout.close();

    }

}

void ds_qq2msq_RES_dataprinter_for_residue_matching(    string filename,
                                                    double a,
                                                    comp sinitial, 
                                                    comp sfinal,
                                                    double delspoints,
                                                    int N    )
{
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double m = 1.0;


    double points1 = 500.0;
    double eps = 0.0;
    double eps_for_m2k = 0.0;
    double qvec_r = 0.0101;
    int acount = 0;

    double dela = 0.051;

    {

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        comp dels = (sfinal - sinitial)/delspoints;//abs(8.70275246887-8.70115128532);//
        points1 = N;
        

        ofstream fout;
        fout.open(filename.c_str());

        int count = 0;

       
        for(int i=0;i<delspoints+1;++i)
        //for(double s3=sinitial;s3<=sfinal;s3=s3+dels)
        { 

            comp s =  sinitial + ((comp) i)*dels; 
            comp sigmamin = {0.0,0.0};
            comp sigmamax = sigmax(s,m);
            comp sigb = 2.0*m*m;//sigmab(a,m);
            comp qval = pmom(s,sigb,m);

            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        

            vector<comp> qvec;
            vector<comp> weights;

            int tag1=0,tag2=0;
            
            int switch_for_gvec_fixer = 0; //0 means we have used contour and turn on gvec fixer
                                           //1 means we have used straight line and turn off gvec fixer
            
            
            int tag_for_m2k = 0;

            if(imag(s)<0.0)
            {
                contour_for_resonance_4(qvec,weights,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k);

            }
            else if(imag(s)>=0.0)
            {
                contour_for_resonance_6(qvec,weights,a,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k);
                
            }

            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
          
            Eigen::VectorXcd interpolater_Gvec = Gvec;
            Gvec = -1.0*Gvec;

            

            if(imag(s)<0.0)
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }
            else if(imag(s)>=0.0)
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                Bmat_maker_momrep_2eps_with_weights_for_resonances(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k);

            }

            double relerror;
            
        

            time(&time_start);
            LinearSolver_2(Bmat,dsol,Gvec,relerror);

            //cusolverComplex(Bmat,Gvec,dsol,size);
            time(&time_end);

            double time_taken = double(time_end - time_start);


            comp result;


            if(imag(s)<0.0)
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                interpolator_ds_integraleq_momrep_2eps_with_weights_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,tag_for_m2k,result);
            }
            else 
            {
                double pi = acos(-1.0);
                double theta = 0.0;//-90.0*pi/180.0;
                interpolator_ds_integraleq_momrep_2eps_with_weights_for_resonances(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,tag_for_m2k,result);
            }

            
            comp m2k = M2kfunc(a,sigb,m, eps_for_m2k);
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            comp ds = result;
            result = m2k*result*m2k;

            comp rhopb = rhophib(s,a,m);
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);

            fout<<setprecision(16)<<real(s)<<'\t'
                <<imag(s)<<'\t'
                <<real(result)<<'\t'
                <<imag(result)<<'\t'
                <<real(ds)<<'\t'
                <<imag(ds)<<'\t'
                <<real(1.0/result)<<'\t'
                <<imag(1.0/result)<<'\t'
                <<real(1.0/ds)<<'\t'
                <<imag(1.0/ds)<<'\t'
                <<a<<'\t'
                <<N<<'\t'
                <<gvalleft<<'\t'
                <<gvalright<<'\t'
                <<real(phibthreshold(a,m))<<'\t'
                <<eps<<endl;
            count = count + 1;
            
        }
        fout.close();

    }

    //}

}



void residue_matching_between_1st_bs_and_1st_2nd_res_datagenerator()
{
    double m = 1.0;
    comp pole = 0.0;

    comp ii = {0.0,1.0};

    comp init_guess1 = 8.999;//9 - 5e-5;
    comp init_guess2 = 8.99999999;//8.9999999;//9 - 6e-6; 

    double delta_s = 0.000001;
    //double a = -2.0; 
    double N = 500.0;
    double a = -8.8;

    double ainitial = -8.8;
    double afinal = -8.7;
    double apoints = 112;

    double dela = abs(ainitial - afinal)/apoints; 

    comp prevpole = 0.0;

    //generate data for the BS1 
    /*
    for(int i=0;i<apoints;++i)
    {
        comp newpole = 0;
        double a = ainitial + i*dela;
        dSqqs2q2msq_pole_searching_secant_with_N_using_weights_for_boundstates(init_guess1, init_guess2, a, pole, N);
        newpole = pole; 
        cout<<setprecision(16)<<"a = "<<a<<'\t'<<"pole = "<<pole<<endl;
        string filename = "BS1_residue_data_" + to_string(i) + ".dat";
        comp gap = 1.0e-6;
        comp sinitial = pole - gap; 
        comp sfinal = pole + gap; 
        ds_qq2msq_BS_dataprinter_for_residue_matching(filename,a,sinitial,sfinal,10,500);
        cout<<"file = "<<filename<<" created"<<endl; 
        init_guess1 = real(pole);
         
        if(abs(prevpole) == abs(newpole))
        {
            cout<<"broken!!"<<endl;
            break; 
        }
        prevpole = pole;

    }
    */
    afinal = -8.6;
    dela = abs(ainitial - afinal)/apoints; 
    //generate data for res1 and 2 
    init_guess1 = 9.0000001 - ii*5e-5;
    init_guess2 = 9.0000001 - ii*6e-5;//8.9999999;//9 - 6e-6; 

    for(int i=0;i<apoints;++i)
    {
        comp newpole = 0.0;
        double a = ainitial + i*dela;
        
        //dSqqs2q2msq_pole_searching_secant_with_N_using_weights_below_threebody_threshold(init_guess1, init_guess2, a, pole, N);
        dSqqs2q2msq_pole_searching_secant_with_N_using_weights_below_threebody_threshold_secondsheet_dS(init_guess1, init_guess2, a, pole, N);
        cout<<setprecision(16)<<"i="<<i<<'\t'<<"a = "<<a<<'\t'<<"pole = "<<pole<<endl;
        newpole = pole;
        init_guess1 = pole - ii*5e-5;
        init_guess2 = pole - ii*6e-6; 
        //fout<<setprecision(16)<<a<<'\t'<<real(pole)<<'\t'<<imag(pole)<<endl;
        string filename = "Res_1_2_residue_data_" + to_string(i) + ".dat";
        
        comp gap = 1.0e-6;
        comp sinitial = pole - gap; 
        comp sfinal = pole + gap; 
        ds_qq2msq_RES_dataprinter_for_residue_matching(filename, a, sinitial, sfinal, 10, 500);
        
        cout<<"file = "<<filename<<" created"<<endl; 

         
        if(abs(prevpole) == abs(newpole))
        {
            cout<<"broken!!"<<endl;
            //break; 
        }
        prevpole = pole;
    }


    
}

void secondsheet_poletesting_new_datagenerator()
{
    double m = 1.0;
    comp pole = 0.0;

    comp ii = {0.0,1.0};

    comp init_guess1 = 8.9998210944585;//9 - 5e-5;
    comp init_guess2 = 8.9999210944585;//8.9999999;//9 - 6e-6; 

    double delta_s = 0.000001;
    //double a = -2.0; 
    double N = 1000.0;
    double a = -8.8;

    double ainitial = -586.70;
    double afinal = -211.0;
    double apoints = 1870;

    double dela = abs(ainitial - afinal)/apoints; 

    comp prevpole = 0.0;

    //generate data for the BS1 

    string filename = "BS_2_extradata_3.dat";
    ofstream fout; 
    fout.open(filename.c_str());
    comp newpole = 0;
    
    for(int i=0;i<apoints+1;++i)
    {
        double a = ainitial + i*dela;
        dSqqs2q2msq_pole_searching_secant_with_N_using_weights_for_boundstates(init_guess1, init_guess2, a, pole, N);
        cout<<setprecision(16)<<"i = "<<i<<'\t'<<"a = "<<a<<'\t'<<"pole = "<<pole<<endl;
        init_guess1 = real(pole);
        init_guess2 = 8.99999999999;//pole - ii*6e-6; 
        if(real(newpole)!=real(pole))
        {   
            fout<<setprecision(16)<<a<<'\t'<<real(pole)<<'\t'<<imag(pole)<<endl;
        }
        else 
        {

        }    
        newpole = pole;

    }
    fout.close();
    /*
    ainitial = -238.6;
    afinal = -206.2;
    dela = abs(ainitial - afinal)/apoints; 
    //generate data for res1 and 2 
    init_guess1 = 8.999992125321745 - ii*6.165799659773153e-06;
    init_guess2 = 8.999992125321745 - ii*6.166799659773153e-06;//8.9999999;//9 - 6e-6; 

    string filename = "Res_3_extradata_1.dat";
    ofstream fout; 
    fout.open(filename.c_str());

    comp newpole = 0.0;

    for(int i=0;i<apoints;++i)
    {
        double a = ainitial + i*dela;
        
        //dSqqs2q2msq_pole_searching_secant_with_N_using_weights_below_threebody_threshold(init_guess1, init_guess2, a, pole, N);
        dSqqs2q2msq_pole_searching_secant_with_N_using_weights_below_threebody_threshold_secondsheet_dS(init_guess1, init_guess2, a, pole, N);
        cout<<setprecision(16)<<"i="<<i<<'\t'<<"a = "<<a<<'\t'<<"pole = "<<pole<<endl;
        init_guess1 = pole - ii*5e-7;
        init_guess2 = pole - ii*6e-6; 
        if(real(newpole)!=real(pole))
        {   
            fout<<setprecision(16)<<a<<'\t'<<real(pole)<<'\t'<<imag(pole)<<endl;
        }
        else 
        {

        }    
        newpole = pole;
        
    }
    fout.close();

    */
    
}


int main(int argc, char *argv[])
{
    //test_determinant();

    //--------------reference vertex factor ---------------------//

    //using BF method
    //double a = atof(argv[1]);
    
    //double N = atof(argv[2]);
    //double pole = 0;//atof(argv[2]);
    //double delspoleBS1 = 0.001/5.0;
    //double delspoleBS2 = 0.00001/5.0;
    //double delspoleBS3 = 0.0000001/5.0;
    //double delspole = 0.0002;//atof(argv[3]);
    //poletest_bs1_vs_singleN_singlea(N,pole);
    //poletest_bs1_vs_singleN_singlea_using_determinant(a,N,pole);
    //dSqqs2q2msq_belowthreshold_vs_s3_sigk(a, N, pole, delspoleBS1,1);
    //poletest_bs2_vs_singleN_singlea_using_determinant(a,N,pole);
    //dSqqs2q2msq_belowthreshold_vs_s3_sigk(a, N, pole, delspoleBS2,2);
    //poletest_bs3_vs_singleN_singlea_using_determinant(a,N,pole);
    //dSqqs2q2msq_belowthreshold_vs_s3_sigk(a, N, pole, delspoleBS3,3);
    
    //using GL method
    //this portion is for making vertex functions 
    //using the inverse dS method 

    double am = 16.0;
    double N11 = 500.0;
    double spole_a2_smooth = 7.252996307970249;//8.535689549889605;//8.976271055410812;//8.782848835105487;
    double spole_a2_hard = 6.849677772527381;
    double spole_a6_smooth = 8.535689549889605;
    double spole_a6_hard = 8.386039338306581;
    double spole_a16_1_smooth = 8.782848835105487;
    double spole_a16_2_smooth = 8.976271055410812;
    double spole_a16_1_hard = 8.689974511859685;
    double spole_a16_2_hard = 8.975451813904034;
    double spole = spole_a16_2_smooth; 
    //double k = 0.2002582157465981;//0.2005422667075143;//0.200336586980402;//1.000566527914533;//0.1004514762276402;//0.200245413394969;
    int BSNum = 2;
    double delspoleBS1 = 0.001/5.0;
    double delspoleBS2 = 0.00001/5.0;
    double delspoleBS3 = 0.0000001/5.0;
    double delspoleX = 1.0e-9;
    double delspole = delspoleX;//delspoleBS2;//delspoleBS2/100.0;

    //cout<<phibthreshold(am,1)<<endl;
    //dSqqs2q2msq_belowthreshold_vs_s3_sigk_with_weights(k,am,N11,spole,delspole,BSNum );
    //dSqqs2q2msq_belowthreshold_vs_s3_k_vertexfactor_with_weights( am,N11,spole,delspole,BSNum );
    
    //this is for calculating the residue factors from the inverse dS method for mphib a=2,6,16
    //---------------------------------------------------------------//
    //mphib_residue_belowthreshold_vs_s3_k_vertexfactor_with_weights( am,N11,spole,delspole,BSNum );
    //------------------------------------------------------------//

    //----------pole location smooth cutoff : ------------// 
    //poletest_bs1_vs_N_using_determinant();

    //poletest_bs2_vs_N_using_determinant();

    //poletest_bs3_vs_N_using_determinant();
    //---------------------------------------------------//

    //--Vertex Factor Creators using Homogeneous Eqn Solving Method:----//

    //double a = atof(argv[1]);
    
    //double N = atof(argv[2]);
    
    //homogeneous_equation_solution_BS1(a,N);
    //homogeneous_equation_solution_BS2(a,N);
    //homogeneous_equation_solution_BS3(a,N);
    //-----------------------------------------------------------------//

    //test_poles();

    //Mphib_belowthreshold_vs_s3_3d_omp();
    //Mphib_belowthreshold_vs_s3_3d_omp_qvecsecondsheet();

    //Mphib_belowthreshold_vs_s3_for_different_am();

    //Mphib_belowthreshold_vs_s3_fixed_s3imag_contour45();

    //Mphib_belowthreshold_vs_s3_fixed_s3imag_contour46();

    //Mphib_belowthreshold_vs_s3_fixed_s3imag_contour47();

    //Mphib_belowthreshold_vs_s3_fixed_s3imag_contourlinear();

    //Mphib_belowthreshold_vs_s3_3d_omp_qvecsecondsheet();

    //Mphib_belowthreshold_vs_s3_3d_omp_qvec_undersheet();

    //Mphib_belowthreshold_vs_s3_fixed_s3imag_contour43();

    double a = 16;
    int N = 5000;
    //sebatest_Mphib_belowthreshold_vs_s3_fixed_s3imag_contour47(a,N);

    //M_belowthreshold_vs_s3_3d_omp_sigkloop();
    
    //dSqqs2q2msq_belowthreshold_vs_s3();
    
    //Mphib_belowthreshold_vs_s3();
    
    //Mphib_belowthreshold_vs_eps();
    
    //Mphib_belowthreshold_vs_s3();
    
    //Mphib_belowthreshold_vs_s3_N_eps();
    
    //Mphib_belowthreshold_vs_s3_opposite();
    
    //poletest_bs1_vs_N();

    //poletest_bs1_vs_singleN();
    //poletest_bs2_vs_singleN();

    
    
    /*dSqqs2q2msq_belowthreshold_vs_s3( 0, 0, 500);
    dSqqs2q2msq_belowthreshold_vs_s3( 4, 1, 3000);*/
    
    //poletest_bs2_vs_N();
    
    //poletest_bs3_vs_N();
    
    //poletest_bs3_vs_N();
    //poletest_bs3_vs_N_fix();
    //poletest_bs3_vs_singleN_fix();
    //poletest_bs3_vs_singleN();
    //Mphib_belowthreshold_vs_s3();
    
    //Mphib_belowthreshold_vs_eps();
    
    //poletest_with_Mphib_bs1_vs_N();
    
    //poletest_with_Mphib_bs2_vs_N();
    
    //testing_Mphib_belowthreshold_vs_s3();

    //------------------GL Method Solutions------------------------ //

    //Mphib_secondsheet_belowthreshold_vs_s3_fixed_s3imag_contour43_with_weights();
    //Mphib_secondsheet_belowthreshold_vs_s3_a_fixed_s3imag_with_weights();
    //Mphib_secondsheet_belowthreshold_vs_s3_N_fixed_s3imag_with_weights(); //this is for firstsheet as well
    //Mphib_secondsheet_belowthreshold_vs_s3_N_3d_with_weights(); //this is for the 3d plots of Mphib
    
    //Mphib_analytic_continuation_results_with_weights();
    //ds_qq2msq_analytic_continuation_results_with_weights();

    //check_gamma_phib_vs_a();
    
    //-------------------3 body resonances---------------------------//
    
    //ds_belowthreshold_vs_s3_N_3d_with_weights();
    //ds_belowthreshold_vs_s3_N_3d_with_weights_for_resonances();
    //ds_belowthreshold_vs_s3_N_3d_with_weights_for_resonances_1();

    //this is the final one, extends dS to the second sheet from top
    //towards the direction of s = 0
    //ds_belowthreshold_vs_s3_N_3d_with_weights_for_resonances_secondsheet_bottom_half();

    //ds_belowthreshold_vs_s3_N_3d_with_weights_for_resonances_firstsheet();
    //ds_belowthreshold_vs_s3_N_3d_with_weights_for_resonances_secondsheet_full(); //this correctly portrays solutions on the second sheet
    //ds_belowthreshold_vs_s3_N_3d_with_weights_for_resonances_secondsheet_full_a_loop();
    
    //3d amplitude for sheet -1 
    //ds_belowthreshold_vs_s3_N_3d_with_weights_for_resonances_sheet_minus1_full();

    //check_m2k_cut_transition();

    //poletest_res1_vs_N_using_weights();
    //poletest_res1_vs_N_using_weights_below_threebody_threshold();
    //poletest_res1_vs_N_using_weights_below_threebody_threshold();
    
    // for poles on the second riemann sheet on the upper portion
    //poletest_res2_vs_N_using_weights_below_threebody_threshold(); //this is the last one we were checking for Efimov
    
    //this is for poles lying below real(s)<5.5
    //poletest_res1_vs_N_using_weights_below_threebody_threshold_1();

    //new data generator for poles of resonances vs a 
    secondsheet_poletesting_new_datagenerator();
    //here we test the residues of 1st bs, 1st res, 2nd res
    //the data for inverse amplitude near the pole is generated here
    //the actual residue data will generated after fitting using a python script
    //residue_matching_between_1st_bs_and_1st_2nd_res_datagenerator();

    //-----------resonance testing on the second sheet of D---------//
    //F3_testing();
    //poletest_BS_vs_N_using_weights_secant_specific_a();
    //resonance_vs_real_s_sheet1();
    //resonance_vs_real_s_sheet2();

    //--------------------------------------------------------------//


    //Mphib_firstsheet_abovethreshold_vs_s3_N_fixed_s3imag_with_weights();
    //Mphib_firstsheet_abovethreshold_SAmethod_vs_s3_N_fixed_s3imag_with_weights();
    //dSqqsigma2msq_belowthreshold_vs_s3_N_fixed_s3imag_with_weights();
    //dSqqs2q2msq_belowthreshold_vs_s3_with_KDF();
    //poletest_virtual_vs_N_using_weights();
    //poletest_BS1_vs_N_using_weights_KDF();
    //poletest_BS2_vs_N_using_weights_KDF();
    //poletest_BS3_vs_N_using_weights_KDF();
    //poletest_BS1_vs_N_singlea_using_weights_KDF();
    //poletest_BS2_vs_N_singlea_using_weights_KDF();
    //poletest_BS3_vs_N_singlea_using_weights_KDF();

    //poletest_bs1_vs_N_using_weights();
    //poletest_bs1_vs_N_using_weights1();
    //poletest_bs2_vs_N_using_weights();
    //poletest_bs2_vs_N_using_weights1();
    poletest_bs3_vs_N_using_weights();
    //poletest_bs3_vs_N_using_weights1();
    //poletest_bs4_vs_N_using_weights();

    //-----vertex_factor maker using GL method ------------// 

    //homogeneous_equation_phib_solution_BS1();
    //homogeneous_equation_integralequation_vs_k_solution_BS1();
    //homogeneous_equation_integralequation_vs_k_solution_BS2();

    //------vertex_factor maker using GL method and using the reference file 
    //      from inverse d_S method for rescaling the vertex factors ----------//
    //referenced_from_Inverse_dS_method_homogeneous_equation_integralequation_vs_k_solution_BS1();

    //Mphib Residue//
    //homogeneous_equation_phib_residue_solution_BS1();

    //----------------------//

    //print_GL_energies();
    return 0;
}
