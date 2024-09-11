//macos:g++-11 printer.cpp -o printer -O3 -std=c++14 -I ~/homebrew/include/eigen3/ -L ~/homebrew/lib/
//wahab:g++ printer.cpp -o printer -O3 -std=c++14 -I /cm/shared/applications/eigen/3.3.7/include/eigen3/
//wahab_gpu:nvcc -I /cm/shared/applications/cuda-toolkit/10.0.130/include -I /cm/shared/applications/eigen/3.3.7/include/eigen3/ printer_gpu_momrep.cpp  -O3 -std=c++14 -L /cm/shared/applications/cuda-toolkit/10.0.130/lib64 -lcublas -lcusolver -lcudart -o printer_gpu
//wahab_gpu_omp:nvcc -I /cm/shared/applications/cuda-toolkit/10.0.130/include -I /cm/shared/applications/eigen/3.3.7/include/eigen3/ printer_gpu_momrep.cpp  -O3 -std=c++14 -L /cm/shared/applications/cuda-toolkit/10.0.130/lib64 -lcublas -lcusolver -lcudart -lgomp -Xcompiler -fopenmp -o printer_gpu
//ubuntupc_omp:nvcc117 -I/usr/include/eigen3/ -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/include -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/lib64 printer_gpu_momrep.cpp  -O3 -std=c++14 -lcublas -lcusolver -lcudart -lgomp -Xcompiler -fopenmp -o printer_gpu
//ubuntupc_omp1:nvcc -I/usr/include/eigen3/ -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/include -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/lib64 printer_gpu_momrep.cpp  -O3 -std=c++14 -lcublas -lcusolver -lcudart -lgomp -Xcompiler -fopenmp -o printer_gpu

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

void raul_test()
{
    int count = 0;
    //comp ii = {0.0,1.0};
    double relerror = 0.0;
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = 16;
    double m = 1.0;
    

    double points1 = 500.0;
    double N = points1;
    double qvec_r = 0.001;
    double sreal = 8.6;
    double simag = 1.0e-3;

    double eps = 0.0;
    double eps_for_m2k = 0.0;

    comp s = sreal + ii*simag;

    comp sigmamin = {0.0,0.0};
    comp sigmamax = sigmax(s,m);
    comp sigb = sigmab(a,m);
    comp qval = pmom(s,sigb,m);
    cout<<"q = "<<qval<<endl;

    comp kmin = 0.0;
    comp kmax = pmom(s,0.0,m);

    vector<comp> qvec;
    vector<comp> weights;

    int tag1=0,tag2=0;

    vector<comp> tmp_qvec;
    vector<comp> tmp_weights;
    //mom_vector_maker_seba_imsneg(tmp_qvec,tmp_weights,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1);
    mom_vector_maker_43_with_weights_with_seba_imsneg(tmp_qvec,tmp_weights,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps,(double)points1,qvec_r);
        
    
    for(int i=0;i<tmp_qvec.size();++i)
    {
        comp qv = real(tmp_qvec[i]) - ii*imag(tmp_qvec[i]);
        comp w = real(tmp_weights[i]) - ii*imag(tmp_weights[i]);
        qvec.push_back(qv);
        weights.push_back(w);
    }
    cout<<"tag1 = "<<tag1<<endl;
    cout<<"tag2 = "<<tag2<<endl;
    //this portion is for loading ready made contours:
    ifstream fin;
    string contour_file = "for_digonto_sr8.8_si0.001_new.txt";
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
    //qvec = contour;

    

    cout<<"qvec created with size = "<<qvec.size()<<endl;
    int size = qvec.size();
    Eigen::MatrixXcd Bmat(size,size);
    Eigen::VectorXcd Gvec(size);
    Eigen::VectorXcd dsol(size);

    Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
    Eigen::VectorXcd interpolater_Gvec = Gvec;
    Gvec = -1.0*Gvec;

    for(int i=0;i<qvec.size();++i)
    {
        cout<<i<<'\t'<<setprecision(16)<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
    }

    Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
    //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);
    //Bmat_maker_momrep_2eps_withtags_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,tag1,tag2);
    //Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);
    cout<<"Bmat = "<<endl;
    //cout<<setprecision(4)<<Bmat<<endl;
    cout<<"================================="<<endl;

    Eigen::MatrixXcd Gmat(size,size);

    for(int i=0;i<qvec.size();++i)
    {
        for(int j=0;j<qvec.size();++j)
        {
            comp p = qvec[i];
            comp k = qvec[j];
            Gmat(i,j) = GS_pk(s,p,k,m,eps);
        }
    }

    cout<<"Gmat"<<endl;
    
    cout<<"------------------------"<<endl;

    time(&time_start);
    LinearSolver_2(Bmat,dsol,Gvec,relerror);

    //cusolverComplex(Bmat,Gvec,dsol,size);
    time(&time_end);

    double time_taken = double(time_end - time_start);
    cout<<"Time taken to solve : "<<fixed 
        <<time_taken<<setprecision(5);
    cout<<" sec"<<endl;


    comp result;
    interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
    //interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
    //interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
    //interpolator_ds_integraleq_momrep_2eps_withtags_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
    //interpolator_ds_integraleq_momrep_2eps_withtags_with_weights_usingGvec(dsol,qvec,weights,interpolater_Gvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
    comp gfunc = gfuncConst(a,0.0,m);
    comp gsq = gfunc*gfunc;
    comp ds = result;
    result = gsq*result;

        
    comp rhopb = rhophib(s,a,m);
    comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
    comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
    comp GSval = GS_pk(s,qval,qval,m,eps);
    

    /*for(int i=0;i<qvec.size();++i)
    {
        cout<<i<<'\t'<<setprecision(16)<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
    }*/

    for(int i=0;i<qvec.size();++i)
    {
        comp sigk = sigma_p(s,qvec[i],m);
        comp m2k = M2kfunc(a,sigk,m,eps);
        cout<<setprecision(16)<<"qreal = "<<real(qvec[i])<<'\t'
            <<"qimag = "<<imag(qvec[i])<<'\t'
            <<"|m2k| = "<<abs(m2k)<<'\t'
            <<"|G_S| = "<<abs(Gvec[i])<<endl;
    }

    cout<<"M2k = "<<endl;
    for(int i=0;i<qvec.size();++i)
    {
        comp sigk = sigma_p(s,qvec[i],m);
        comp m2k = M2kfunc(a,sigk,m,eps);
        cout<<qvec[i]<<'\t'<<sigk<<'\t'<<m2k<<endl;
    }

    cout<<"-------------------------"<<endl;
    cout<<"s = "<<s<<'\t'<<"Mphib = "<<result<<'\t'<<"run:"<<count + 1<<endl;
    cout<<"-------------------------"<<endl;

    for(int i=0;i<qvec.size();++i)
    {
        comp sigp = sigma_p(s,qvec[i],m);
        comp sigk = sigp;
        cout<<qvec[i]<<'\t'<<Hfunc_comp(sigp,sigk,m)<<endl;
    }
}

void raul_test2()
{
    int count = 0;
    //comp ii = {0.0,1.0};
    double relerror = 0.0;
    time_t time_start, time_end;
    comp ii = {0.0,1.0};
    double a = 16;
    double m = 1.0;
    

    double points1 = 500.0;
    double N = points1;
    double qvec_r = 0.001;
    double sreal = 8.6;
    double simag = 1.0e-3;

    double eps = 0.0;
    double eps_for_m2k = 0.0;

    comp s = sreal + ii*simag;

    ofstream fout;
    string filename = "for_raul.dat";
    fout.open(filename.c_str());

    for(int N=300;N<=2000;N=N+100)
    {
    points1 = N;
    comp sigmamin = {0.0,0.0};
    comp sigmamax = sigmax(s,m);
    comp sigb = sigmab(a,m);
    comp qval = pmom(s,sigb,m);
    cout<<"q = "<<qval<<endl;

    comp kmin = 0.0;
    comp kmax = pmom(s,0.0,m);

    vector<comp> qvec;
    vector<comp> weights;

    int tag1=0,tag2=0;

    vector<comp> tmp_qvec;
    vector<comp> tmp_weights;
    //mom_vector_maker_seba_imsneg(tmp_qvec,tmp_weights,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1);
    mom_vector_maker_43_with_weights_with_seba_imsneg(tmp_qvec,tmp_weights,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps,(double)points1,qvec_r);
        
    
    for(int i=0;i<tmp_qvec.size();++i)
    {
        comp qv = real(tmp_qvec[i]) - ii*imag(tmp_qvec[i]);
        comp w = real(tmp_weights[i]) - ii*imag(tmp_weights[i]);
        qvec.push_back(qv);
        weights.push_back(w);
    }
    cout<<"tag1 = "<<tag1<<endl;
    cout<<"tag2 = "<<tag2<<endl;
    //this portion is for loading ready made contours:
    ifstream fin;
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
    //qvec = contour;

    

    cout<<"qvec created with size = "<<qvec.size()<<endl;
    int size = qvec.size();
    Eigen::MatrixXcd Bmat(size,size);
    Eigen::VectorXcd Gvec(size);
    Eigen::VectorXcd dsol(size);

    Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
    Eigen::VectorXcd interpolater_Gvec = Gvec;
    Gvec = -1.0*Gvec;

    Bmat_maker_momrep_2eps(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k);
    //Bmat_maker_momrep_2eps_withtags(Bmat,s,qvec,qvec,a,m,eps,eps_for_m2k,tag1,tag2);
    //Bmat_maker_momrep_2eps_withtags_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,tag1,tag2);
    //Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

    time(&time_start);
    LinearSolver_2(Bmat,dsol,Gvec,relerror);

    //cusolverComplex(Bmat,Gvec,dsol,size);
    time(&time_end);

    double time_taken = double(time_end - time_start);
    cout<<"Time taken to solve : "<<fixed 
        <<time_taken<<setprecision(5);
    cout<<" sec"<<endl;


    comp result;
    interpolator_ds_integraleq_momrep_2eps(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result);
    //interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
    //interpolator_ds_integraleq_momrep_2eps_withtags(dsol,qvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
    //interpolator_ds_integraleq_momrep_2eps_withtags_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
    //interpolator_ds_integraleq_momrep_2eps_withtags_with_weights_usingGvec(dsol,qvec,weights,interpolater_Gvec,s,qval,qval,a,m,eps,eps_for_m2k,result,tag1,tag2);
    comp gfunc = gfuncConst(a,0.0,m);
    comp gsq = gfunc*gfunc;
    comp ds = result;
    result = gsq*result;

        
    comp rhopb = rhophib(s,a,m);
    comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
    comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
    comp GSval = GS_pk(s,qval,qval,m,eps);
    cout<<"s:"<<s<<'\t'<<"res:"<<result<<'\t'<<"run:"<<count + 1<<endl;
    cout<<"-------------------------"<<endl;
    fout<<N<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;

    /*for(int i=0;i<qvec.size();++i)
    {
        cout<<i<<'\t'<<setprecision(16)<<real(qvec[i])<<'\t'<<imag(qvec[i])<<endl;
    }*/

    /*for(int i=0;i<qvec.size();++i)
    {
        comp sigk = sigma_p(s,qvec[i],m);
        comp m2k = M2kfunc(a,sigk,m,eps);
        cout<<setprecision(16)<<"qreal = "<<real(qvec[i])<<'\t'
            <<"qimag = "<<imag(qvec[i])<<'\t'
            <<"|m2k| = "<<abs(m2k)<<'\t'
            <<"|G_S| = "<<abs(Gvec[i])<<endl;
    }*/

    
    }
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
    

    double spole = 8.782848835105487;//8.689974511859685;//7.252996307970249;//8.535689549889605;//8.782848835105487;
    

    //int N = 10000;
    //poletest_bs1_vs_singleN_singlea_using_determinant(a, N,spole);
    //cout<<"pole found at s = "<<spole<<endl;

    
    double kmin = 0.0;
    comp kmax = pmom(spole,0.0,m);//1.3;
    double kpoints = 2000.0;
    double delk = abs(kmax - kmin)/kpoints;
    //double delk = abs(kmax - kmin)/N;

    ofstream fout;
    //string filename = "vertexfactor_vs_a_"+to_string(a)+"BS1_N_" + to_string(N) + "_smoothcutoff.dat";
    string filename = "temp.txt";
    
    fout.open(filename.c_str());

    //for(int somei=0;somei<avec.size();++somei)
    //{
        //double a = avec[somei];
        //double spole = spolevec[somei];
    //for(int i=0;i<kpoints;++i)
    {
        double k = 0.1;//kmin + i*delk;
        double s = spole;
        
        double points1 = N;
        vector<comp> qvec;
        vector<comp> weights;
        //mom_vector_maker_linear_1(qvec,kmin,delk,N,1);
        //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
        /*1.00000000e-08 1.45897688e-01 2.91795366e-01 4.37693044e-01
 5.83590721e-01 7.29488399e-01 8.75386077e-01 1.02128375e+00
 1.16718143e+00*/
        qvec.push_back(1.00000000e-08);
        qvec.push_back(1.45897688e-01);
        qvec.push_back(2.91795366e-01);
        qvec.push_back(4.37693044e-01);
        qvec.push_back(5.83590721e-01);
        qvec.push_back(7.29488399e-01);
        qvec.push_back(8.75386077e-01);
        qvec.push_back(1.02128375e+00);
        qvec.push_back(1.16718143e+00);

        
        //mom_vector_maker_43_with_weights(qvec,weights,s,kmin,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                
        int size = qvec.size();

        Eigen::MatrixXcd Bmat(size,size);
        Eigen::MatrixXcd eigBmat(size,size);
        
        
        Bmat_maker_momrep(Bmat,s,qvec,qvec,a,m,eps);
        //Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);


        cout<<"Bmat Determinant = "<<Bmat.determinant()<<endl;
        cout<<"=========================="<<endl;
        cout<<"Bmat = "<<endl;
        cout<<Bmat<<endl;


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

        for(int i=0;i<eigVal.size();++i)
        {
            cout<<"eigenvalues = "<<eigVal(i)<<endl;
        } 

        for(int i=0;i<eigVal.size();++i)
        {
            for(int j=0;j<eigVal.size();++j)
            {
                cout<<"eigvec "<<i<<","<<j<<" = "<<eigBmat(j,i)<<endl;
            }
        }


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
            //for(int j=0;j<eigVal.size();++j)
            for(int j=0;j<qvec.size();++j)
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
                cout<<"k = "<<qvec[j]<<'\t'<<"m2k = "<<m2kres<<'\t'<<"Eigvec = "<<vertexres<<'\t'<<"Gamma(k) = "<<vertexres*m2kres<<endl;
                //fout<<real(qvec[j])<<'\t'<<real(m2kressq)<<'\t'<<real(vertexressq)<<'\t'<<real(res)<<endl;
                //cout<<j<<'\t'<<real(qvec[j])<<'\t'<<real(m2kressq)<<'\t'<<real(vertexressq)<<'\t'<<real(res)<<endl;
                if(j==0)
                {
                    //sumvert = sumvert + abs(vertexressq)*(qvec[j]);
                    //last
                    //sumvert = sumvert + abs(vertexressq)*(weights[j]);    
                    
                }
                else 
                {
                    //sumvert = sumvert + abs(vertexressq)*((qvec[j] - qvec[j-1]));
                    //last
                    //sumvert = sumvert + abs(vertexressq)*(weights[j]);
                
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
        comp detB = {1.0,0.0};
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



int main()
{
    raul_test();
    //homogeneous_equation_integralequation_vs_k_solution_BS1();
    return 0;
}