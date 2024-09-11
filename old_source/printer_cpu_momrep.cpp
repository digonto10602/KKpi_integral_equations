//macos:g++-11 printer_cpu_momrep.cpp -o printer -O3 -std=c++14 -I /usr/local/opt/eigen/include/eigen3
//wahab:g++ printer.cpp -o printer -O3 -std=c++14 -I /cm/shared/applications/eigen/3.3.7/include/eigen3/
//wahab_gpu:nvcc -I /cm/shared/applications/cuda-toolkit/10.0.130/include -I /cm/shared/applications/eigen/3.3.7/include/eigen3/ printer_cpu_momrep.cpp  -O3 -std=c++14 -L /cm/shared/applications/cuda-toolkit/10.0.130/lib64 -lcublas -lcusolver -lcudart -o printer_gpu
//wahab_cpu:g++ -I /cm/shared/applications/eigen/3.3.7/include/eigen3/ printer_cpu_momrep.cpp  -O3 -std=c++14 -fopenmp -o printer_cpu

#include<bits/stdc++.h>
#include "functions_momrep_based.h"
#include "momentum_vector_maker.h"
#include "integralequation_momrep_based.h"
#include "interpolator_momrep.h"
#include "solvers.h"
#include<omp.h> 

using namespace std;

typedef complex<double> comp;

void Mphib_belowthreshold_vs_s3()
{
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.70;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/100.0;
    double points1 = 7000.0;
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

        mom_vector_maker_4(qvec,s,kmin,kmax,a,m,eps,eps_for_m2k,points1,qvec_r);
    
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec = -1.0*Gvec;

        Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);

        double relerror;

        

        LinearSolver(Bmat,dsol,Gvec,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

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
    }

    fout.close();

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
                    + ".dat";
    ofstream fout;
    

    int count = 0;

    double eta = 20.0;

    double s3realinitial = 8.72;
    double s3realfinal = real(phibthreshold(a,m));
    double dels3real = abs(s3realinitial-s3realfinal)/s3realpoints;

    double s3imaginitial = -0.1;
    double s3imagfinal = 0.1;
    double dels3imag = abs(s3imaginitial - s3imagfinal)/s3imagpoints;
    cout<<"dels3real = "<<dels3real<<'\t'<<"dels3imag = "<<dels3imag<<endl;

    //vector<vector<comp> > resultvec(points1,points1);
    comp **resultvec = new comp*[(int)s3imagpoints+1];
    for(int i=0; i<(int)s3imagpoints+1; ++i)
    {
        resultvec[i] = new comp[(int)s3realpoints+1];
    }

    //comp resultvec[(int)s3imagpoints+1][(int)s3realpoints+1];

    int i;
    #pragma omp parallel for shared(resultvec)
    for(i=0;i<(int)s3imagpoints+1;++i)
    {
        //for(double s=sinitial;s<=sfinal;s=s+dels)
        
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

            //#pragma omp critical
            {
                //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
                //fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<'\t'<<real(1.0/result)<<'\t'<<imag(1.0/result)<<endl;
                resultvec[i][j] = result;
                
                count = count + 1;
                //cout<<"run = "<<count<<endl;
                //cout<<"-------------------------"<<endl;
            }
        }
    }

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

            cout<<s3real<<'\t'<<s3imag<<'\t'<<real(resultvec[k][l])<<'\t'<<imag(resultvec[k][l])<<endl;
            fout<<s3real<<'\t'<<s3imag<<'\t'<<real(resultvec[k][l])<<'\t'<<imag(resultvec[k][l])<<endl;
        }
    }

    fout.close();

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

    string filename="Mphib_momrep_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_3d"
                    //+ "_epspow_" + to_string((int)abs(log10(eps))) 
                    + "_withoutomp_qss_cpu.dat";//_followq_poss.dat";
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
    omp_set_num_threads(20);
    //#pragma omp distribute parallel for shared(resultvec,count)

    
    //cudaStream_t nstreams[(int)s3imagpoints+1];
    
    
    #pragma omp target teams distribute shared(resultvec)
    for(i=0;i<(int)s3imagpoints+1;++i)
    {
        //for(double s=sinitial;s<=sfinal;s=s+dels)
        //std::cout<<"from Thread = "<<omp_get_thread_num()<<endl;
        //std::cout<<"total Thread = "<<omp_get_num_threads()<<endl;
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
            //if(s3imag<=0.0) qval = -qval;


            comp kmin = 0.0;
            comp kmax = pmom(s,0.0,m);
        
            //cout<<"s:"<<s<<endl;
            //cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
            //cout<<"eps:"<<eps<<endl;

            vector<comp> qvec;

            comp delk = abs(kmin - kmax)/points1;
            int tag1,tag2;
            
            mom_vector_maker_6(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,points1,qvec_r,tag1,tag2);
            //mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
            //mom_vector_maker_4(qvec,s,kmin,kmax,a,m,eps,eps_for_m2k,points1,qvec_r);
    
            int size = qvec.size();
            Eigen::MatrixXcd Bmat(size,size);
            Eigen::VectorXcd Gvec(size);
            Eigen::VectorXcd dsol(size);

            
            //Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
            Gvec_maker_momrep_withtags(Gvec,s,qvec,qval,m,eps,tag1,tag2);
            Gvec = -1.0*Gvec;

            //std::cout<<"i = "<<i<<'\t'<<" j = "<<j<<endl;
            Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
            //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);


            //std::cout<<"Bmat made from i = "<<i<<'\t'<<" j = "<<j<<endl;
            double relerror;
            
        
            //cout<<"did this = "<<count<<" times"<<endl;
            //cout<<"i val = "<<i<<'\t';
            //cout<<"j val = "<<j<<endl;
            
            #pragma omp critical
            {
                //std::cout<<"running iteration i = "<<i<<" from thread = "<<omp_get_thread_num()<<endl;
                LinearSolver_2(Bmat,dsol,Gvec,relerror);
            }
            

            
            //cudaStream_t somethingstream;
            //cudaStreamCreate(&somethingstream);

            //#pragma omp critical
            //{
                //cusolverComplex(Bmat,Gvec,dsol,size);
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

            //interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
            interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
            comp gfunc = gfuncConst(a,0.0,m);
            comp gsq = gfunc*gfunc;
            result = gsq*result;
            //result = rhophib(s,a,m)*result;
            
            resultvec[i][j] = result;

            //#pragma omp critical
            {
                //cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
                //fout<<s<<'\t'<<real(result)<<'\t'<<imag(result)<<'\t'<<real(1.0/result)<<'\t'<<imag(1.0/result)<<endl;
                
                
                //count = count + 1;
                //std::cout<<count<<" out of "<<((s3realpoints+1)*(s3imagpoints+1))<<" runs finished"<<endl;
                //cout<<"run = "<<count<<endl;
                //cout<<"-------------------------"<<endl;
            }

            //if(count%10==0 || count == ((s3realpoints+1)*(s3imagpoints+1)))
            //std::cout<<count<<" out of "<<((s3realpoints+1)*(s3imagpoints+1))<<" runs finished"<<endl;

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

            std::cout<<s3real<<'\t'<<s3imag<<'\t'<<real(resultvec[k][l])<<'\t'<<imag(resultvec[k][l])<<endl;
            fout<<s3real<<'\t'<<s3imag<<'\t'<<real(resultvec[k][l])<<'\t'<<imag(resultvec[k][l])<<endl;
        }
    }

    fout.close();

    for(int i=0; i<s3imagpoints+1; ++i)
    {
        delete [] resultvec[i];
    }

    delete [] resultvec;
}



void dSqqs2q2msq_belowthreshold_vs_s3()
{
    double a = 1000.0;
    double m = 1.0;

    //double s = 8.65;
    double spole = 8.9985;
    double delspole = real(phibthreshold(a,m)) - spole;
    
    double sinitial = 8.99998;//spole;// - delspole;//8.95;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.99;
    
    double dels = abs(sinitial-sfinal)/100.0;
    double points1 = 500.0;
    double points2 = 500.0;
    double eps = 0.0;
    //double box = 10.0;

    string filename="dSqqs2q2msq_momrep_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_eps_" + to_string(eps) 
                    + "_cpu.dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    double eta = 20.0;


    for(double s=sinitial;s<=sfinal;s=s+dels)
    { 
        //double s = 8.78;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        cout<<phibthreshold(a,m)<<endl;
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
        //result = gsq*result;
        //result = rhophib(s,a,m)*result;

        cout<<"s:"<<setprecision(16)<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        fout<<setprecision(16)<<s<<'\t'
            <<real(result)<<'\t'
            <<imag(result)<<'\t'
            <<real(1.0/result)<<'\t'
            <<imag(1.0/result)<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    }

    fout.close();

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
    
    double epsilon = 1.0e-10;
    double tolerance = 1.0e-15;
    double sa = apoint;
    double sb = bpoint;

    for(int i=0;i<max_iteration;++i)
    {
        
        double sc = (sa + sb)/2.0;
        double fc = (double) real(dSqqs2q2msq_func_with_N(sc,a,N));
        double fa = (double) real(dSqqs2q2msq_func_with_N(sa,a,N));
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



void dSqqs2q2msq_belowthreshold_vs_N()
{
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.78;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/100.0;
    //double points1 = 1000.0;
    double points2 = 500.0;
    double eps = 0.01;
    //double box = 10.0;
    double s = 8.78;

    string filename="dSqqs2q2msq_momrep_a_" + to_string((int)a) 
                    + "_s_" + to_string(s) 
                    + "_eps_" + to_string(eps) 
                    + ".dat";
    ofstream fout;
    fout.open(filename.c_str());

    int count = 0;

    double eta = 20.0;

    

    //for(double s=sinitial;s<=sfinal;s=s+dels)
    for(double points1=500;points1<=3000;points1=points1 + 500 )
    { 
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
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;

        mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
        //mom_vector_maker_1(qvec,s,kmin,kmax,a,m,eps,points1);
        cout<<"qvecsize:"<<qvec.size()<<endl;
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
        //result = gsq*result;
        //result = rhophib(s,a,m)*result;

        cout<<"s:"<<s<<'\t'<<"dels:"<<dels<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        fout<<points1<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
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

        

        LinearSolver(Bmat,dsol,Gvec,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

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


void gamma_test()
{
    double a = 16.0;
    double m = 1.0;

    //double s = 8.65;
    double sinitial = 8.72;//8.72;//real(phibthreshold(a,m));//
    double sfinal = real(phibthreshold(a,m));//8.82;
    double dels = abs(sinitial-sfinal)/100.0;
    double points1 = 10.0;
    double points2 = 500.0;
    double eps = 0.0;
    double sigvalinitial = 1.0*m*m;
    double sigvalfinal = 2.0*m*m;
    double delsigval = abs(sigvalinitial - sigvalfinal)/10.0;
    //double box = 10.0;

    string filename="gammatest_momrep_a_" + to_string((int)a) 
                    + "_N_" + to_string((int)points1) 
                    + "_eps_" + to_string(eps) 
                    + ".dat";
    //ofstream fout;
    //fout.open(filename.c_str());

    int count = 0;

    double eta = 20.0;


    //for(double sigval=sigvalinitial;sigval<=sigvalfinal;sigval=sigval+delsigval)
    //{ 
        double sigval = sigvalfinal;
        double s = 8.78287;
        //if(s>=sfinal) dels = 0.05/100.0;
        //eps = eps_above_threshold(eta,s,a,m,points);
        comp sigmamin = {0.0,0.0};
        comp sigmamax = sigmax(s,m);
        comp sigb = sigmab(a,m);
        comp sigq = sigval;
        comp qval = pmom(s,sigq,m);

        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        cout<<"kmax = "<<kmax<<endl;
        double delk = abs(kmax-kmin)/points1;
        
        //cout<<"s:"<<s<<endl;
        cout<<"sigmin:"<<sigmamin<<'\t'<<"sigmax:"<<sigmamax<<'\t'<<"sigb:"<<sigb<<endl;
        cout<<"eps:"<<eps<<endl;

        vector<comp> qvec;

        mom_vector_maker_linear_1(qvec,kmin,delk,points1,1);
        //mom_vector_maker_1(qvec,s,kmin,kmax,a,m,eps,points1);
    
        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        //Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        //Gvec = -1.0*Gvec;
        for(int i=0;i<size;++i) Gvec(i) = 0;

        Bmat_maker_momrep(Bmat,s,qvec,qvec,a,m,eps);

        double relerror;

        

        LinearSolver_1(Bmat,dsol,Gvec,relerror);
        //cusolverComplex(Bmat,Gvec,dsol,size);

        comp result;

        interpolator_gamma_integraleq_momrep(dsol,qvec,s,qval,qval,a,m,eps,result);

        comp gfunc = gfuncConst(a,0.0,m);
        comp gsq = gfunc*gfunc;
        //result = gsq*result;
        //result = rhophib(s,a,m)*result;

        cout<<"sigval:"<<sigval<<'\t'<<"dels:"<<delsigval<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
        //fout<<sigval<<'\t'<<real(qval)<<'\t'<<imag(qval)<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
        count = count + 1;
        cout<<"-------------------------"<<endl;
    //}

    //fout.close();

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

                    cout<<"s:"<<s<<'\t'<<"dels:"<<size<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
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
    double init_guess1 = 8.70;
    double init_guess2 = 8.99;//(double) real(phibthreshold(a,m));
    
    for(int N=500; N<=4000; N=N+500)
    {
        vector<double> polevec;
        string filename = "bs1st_N_" + to_string(N) + ".dat";
        
        for(double a = 50.0; a<3000.0; a = a + 50.0)
        {
            //double dels = 0.00001;
            //double tempval = 0.0;
            //double set_s = (double) real(phibthreshold(a,m));
            //init_guess2 = (double) real(phibthreshold(a,m));
    
            //(double) real(phibthreshold(a,m));
            
        
        
            dSqqs2k2msq_pole_searching_bissection_with_N(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
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



void poletest_bs2_vs_N()
{
    ofstream fout;
    double m = 1.0;
    double pole = 0.0;
    
    //vector<double> init_guessvec(20);
    //init_guessvec[0] = 0.0;
    double init_guess1 = 8.99;
    double init_guess2 = 0.0;//(double) real(phibthreshold(a,m));
    
    for(int N=500; N<=4000; N=N+500)
    {
        vector<double> polevec;
        string filename = "bs2nd_N" + to_string(N) + ".dat";
        
        for(double a = 50.0; a<3000.0; a = a + 50.0)
        {
            double dels = 0.00001;
            double tempval = 0.0;
            double set_s = (double) real(phibthreshold(a,m));
            init_guess2 = (double) real(phibthreshold(a,m));
    
            //(double) real(phibthreshold(a,m));
            for(int i=0;i<100;++i)
            {
                double s = set_s - i*dels;
                tempval = real(dSqqs2q2msq_func_with_N(s,a,N));
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
            }
        
        
            dSqqs2k2msq_pole_searching_bissection_with_N(init_guess1,init_guess2,a,pole,N);
            polevec.push_back(pole);
            //init_guess2 = pole;
        }

        fout.open(filename.c_str());
        for(int i=0;i<polevec.size();++i)
        {
            double a1 = 50.0 + i*50.0;
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



void test_bissection_endpoint()
{
    double a = 50.0;
    double m = 1.0;

    double sinitial = 8.9;
    double sfinal = (double) real(phibthreshold(a,m));
    double dels = abs(sinitial-sfinal)/200.0;

    for(int i=0;i<200;++i)
    {
        double s = sinitial + i*dels;

        cout<<setprecision(16)<<s<<'\t'<<real(dSqqs2q2msq_func(s,a))<<endl;
    }
}

int main()
{
    //Mphib_belowthreshold_vs_s3();
    //Mphib_belowthreshold_vs_eps();
    //Mphib_belowthreshold_vs_s3();
    //dSqqs2q2msq_belowthreshold_vs_s3();
    //dSqqs2q2msq_belowthreshold_vs_N();
    //gamma_test();
    //poletest();
    //poletest_bs1_vs_N();
    //poletest_bs2_vs_N();
    //Mphib_belowthreshold_vs_s3_3d_omp();
    Mphib_belowthreshold_vs_s3_3d_omp_qvecsecondsheet();
    //dSqqs2q2msq_belowthreshold_vs_s3();
    //test_bissection_endpoint();
    
    //Mphib_belowthreshold_vs_s3_opposite();

    return 0;
}
