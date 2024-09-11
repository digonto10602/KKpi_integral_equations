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
#include "d_secondsheet_func.h"
#include <sys/time.h>

using namespace std;

typedef complex<double> comp;

void test_functions()
{
    double a = -8;
    double m = 1.0;

    comp s = {8.8,0.0001};

    int N = 500;

    comp sig_q = 2.0*m*m;
    comp q = pmom(s,sig_q, m);

    double eps_for_m2k = 0.0;
    double eps = 0.0;

    //int N = 500;

    //solution for Ds
    comp kmin = 0.0;
    comp kmax = pmom(s,0.0,m);
    comp sigb = 2.0*m*m;
    comp qval = pmom(s,sigb,m);

    vector<comp> qvec;
    vector<comp> weights;

    line_maker_with_weights(qvec,weights,kmin,kmax,N);

    int size = qvec.size();
    Eigen::MatrixXcd Bmat(size,size);
    Eigen::MatrixXcd Gmat(size,size);
    Eigen::VectorXcd Gvec(size);
    Eigen::VectorXcd dsol(size);
    Eigen::MatrixXcd dmat(size,size);

    Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
    Gvec = -1.0*Gvec;
    Gmat_maker_momrep(Gmat,s,qvec,qvec,m,eps);
    Gmat = -1.0*Gmat;

        
    Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

    double relerror; 
    struct timeval stop, start;
    gettimeofday(&start, NULL);
    LinearSolver_2(Bmat,dsol,Gvec,relerror);
    LinearSolver_3(Bmat,dmat,Gmat,relerror);
        
        
    gettimeofday(&stop, NULL);
    double actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/1000000.0;
    cout<<"time taken to solve = "<<actual_sol_time<<endl;

    Eigen::VectorXcd Ds(size);
    Eigen::MatrixXcd Dmats(size,size);
    for(int i=0;i<size;++i)
    {
        comp sig_i = sigma_p(s,qvec[i],m);
        comp sig_j = sigb;

        comp m2k_i = M2kfunc(a,sig_i,m, eps_for_m2k);
        comp m2k_j = M2kfunc(a,sigb,m, eps_for_m2k);
        Ds(i) = m2k_i*dsol(i)*m2k_j;
        for(int j=0;j<size;++j)
        {
            sig_j = sigma_p(s,qvec[j],m);
            m2k_j = M2kfunc(a,sig_j,m, eps_for_m2k);
        
            //Dmats(i,j) = m2k_i*dmat(i,j)*m2k_j;
        }
    }

    
    //solution for D_dagger
    Eigen::MatrixXcd Cmat(size,size);
    Eigen::VectorXcd Bvec(size);
    Eigen::VectorXcd Ddagsol(size);

    Bvec_Dsecsht_maker_with_weights(Bvec,Ds,qvec,weights,s,qval,m,a,eps);
    cmat_maker_with_weights(Cmat,Dmats,qvec,weights,s,m,a,eps);

    double relerror1; 
    //struct timeval stop, start;
    gettimeofday(&start, NULL);
    LinearSolver_2(Cmat,Ddagsol,Bvec,relerror1);
    
        
    gettimeofday(&stop, NULL);
    actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/1000000.0;
    cout<<"time taken to solve = "<<actual_sol_time<<endl;

    for(int i=0;i<size;++i)
    {
        cout<<"q = "<<qvec[i]<<'\t'<<"Ds = "<<Ds(i)<<'\t'<<"D_dag_s = "<<Ddagsol(i)<<endl;
    }
       
}

void test_functions_2()
{
    comp ii = {0.0,1.0};
    double a = -8;
    double m = 1.0;

    //comp s = {8.8,0.0001};

    double simag = 0.0075;

    double sreal_initial = 8.99;
    double sreal_final = 9.01;
    double spoints = 500.0;
    double dels = abs(sreal_initial - sreal_final)/spoints;

    int N = 500;

    string filename = "D_dag_test.dat";
    ofstream fout;
    fout.open(filename.c_str());

    for(int i=0;i<spoints+1;++i)
    {
        comp s = sreal_initial + i*dels + ii*simag;
        comp sig_q = 2.0*m*m;
        comp q = pmom(s,sig_q, m);

        double eps_for_m2k = 0.0;
        double eps = 0.0;

        //int N = 500;

        //solution for Ds
        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        comp sigb = 2.0*m*m;
        comp qval = pmom(s,sigb,m);

        vector<comp> qvec;
        vector<comp> weights;

        line_maker_with_weights(qvec,weights,kmin,kmax,N);

        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::MatrixXcd Gmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);
        Eigen::MatrixXcd dmat(size,size);

        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec = -1.0*Gvec;
        Gmat_maker_momrep(Gmat,s,qvec,qvec,m,eps);
        Gmat = -1.0*Gmat;

        
        Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

        double relerror; 
        struct timeval stop, start;
        gettimeofday(&start, NULL);
        LinearSolver_2(Bmat,dsol,Gvec,relerror);
        LinearSolver_3(Bmat,dmat,Gmat,relerror);
        
        
        gettimeofday(&stop, NULL);
        double actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/1000000.0;
        cout<<"time taken to solve = "<<actual_sol_time<<endl;

        Eigen::VectorXcd Ds(size);
        Eigen::MatrixXcd Dmats(size,size);
        for(int i=0;i<size;++i)
        {
            comp sig_i = sigma_p(s,qvec[i],m);
            comp sig_j = sigb;

            comp m2k_i = M2kfunc(a,sig_i,m, eps_for_m2k);
            comp m2k_j = M2kfunc(a,sigb,m, eps_for_m2k);
            Ds(i) = m2k_i*dsol(i)*m2k_j;
            for(int j=0;j<size;++j)
            {
                sig_j = sigma_p(s,qvec[j],m);
                m2k_j = M2kfunc(a,sig_j,m, eps_for_m2k);
        
                //Dmats(i,j) = m2k_i*dmat(i,j)*m2k_j;
            }
        }

    
        //solution for D_dagger
        Eigen::MatrixXcd Cmat(size,size);
        Eigen::VectorXcd Bvec(size);
        Eigen::VectorXcd Ddagsol(size);

        
        double relerror1; 
        //struct timeval stop, start;
        gettimeofday(&start, NULL);

        Bvec_Dsecsht_maker_with_weights(Bvec,Ds,qvec,weights,s,qval,m,a,eps);
        cmat_maker_with_weights(Cmat,Dmats,qvec,weights,s,m,a,eps);

        LinearSolver_2(Cmat,Ddagsol,Bvec,relerror1);
    
        
        gettimeofday(&stop, NULL);
        actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/1000000.0;
        cout<<"time taken to solve = "<<actual_sol_time<<endl;

        double temp_val=0.0;
        int ind_q = 0;
        comp selected_val = 0.0;

        for(int j=0;j<size;++j)
        {
            comp some_p = qvec[j];
            //cout<<"q = "<<qvec[i]<<'\t'<<"Ds = "<<Ds(i)<<'\t'<<"D_dag_s = "<<Ddagsol(i)<<endl;
            if(j==0)
            {
                temp_val = abs(some_p - qval);
                ind_q = j;
                continue; 
            }
            
            if(abs(some_p - qval)<temp_val)
            {
                temp_val = abs(some_p - qval);
                ind_q = j;
                selected_val = some_p;
            }
        
        
        }

        cout<<"run = "<<i+1<<endl;
        cout<<"qval = "<<qval<<'\t'<<"selected_val = "<<selected_val<<endl;
        cout<<"diff = "<<abs(qval - selected_val)<<endl;
        cout<<"s = "<<s<<'\t'<<"Ds = "<<Ds(ind_q)<<'\t'<<"D_dag_s = "<<Ddagsol(ind_q)<<endl;
        fout<<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<real(Ds(ind_q))<<'\t'
            <<imag(Ds(ind_q))<<'\t'
            <<real(Ddagsol(ind_q))<<'\t'
            <<imag(Ddagsol(ind_q))<<endl;
        cout<<"======================================"<<endl;

    }
    fout.close();  
}

void test_functions_3()
{
    comp ii = {0.0,1.0};
    double a = -8;
    double m = 1.0;

    //comp s = {8.8,0.0001};

    double sreal = 9.01;
    double simag = 0.0075;

    double sreal_initial = 8.99;
    double sreal_final = 9.01;
    double spoints = 100.0;
    double dels = abs(sreal_initial - sreal_final)/spoints;

    double simag_initial = -0.01;
    double simag_final = 0.01;
    double dels_imag = abs(simag_initial - simag_final)/spoints;
    int N = 500;

    string filename = "D_dag_test_simag.dat";
    ofstream fout;
    fout.open(filename.c_str());

    for(int i=0;i<spoints+1;++i)
    {
        //comp s = sreal_initial + i*dels + ii*simag;

        comp s = sreal + ii*simag_initial +  (comp (i))*ii*dels_imag;
        
        comp sig_q = 2.0*m*m;
        comp q = pmom(s,sig_q, m);

        double eps_for_m2k = 0.0;
        double eps = 0.0;

        //int N = 500;

        //solution for Ds
        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        comp sigb = 2.0*m*m;
        comp qval = pmom(s,sigb,m);

        vector<comp> qvec;
        vector<comp> weights;

        line_maker_with_weights(qvec,weights,kmin,kmax,N);

        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::MatrixXcd Gmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);
        Eigen::MatrixXcd dmat(size,size);

        Gvec_maker_momrep(Gvec,s,qvec,qval,m,eps);
        Gvec = -1.0*Gvec;
        Gmat_maker_momrep(Gmat,s,qvec,qvec,m,eps);
        Gmat = -1.0*Gmat;

        
        Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

        double relerror; 
        struct timeval stop, start;
        gettimeofday(&start, NULL);
        LinearSolver_2(Bmat,dsol,Gvec,relerror);
        LinearSolver_3(Bmat,dmat,Gmat,relerror);
        
        
        gettimeofday(&stop, NULL);
        double actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/1000000.0;
        cout<<"time taken to solve = "<<actual_sol_time<<endl;

        Eigen::VectorXcd Ds(size);
        Eigen::MatrixXcd Dmats(size,size);
        for(int i=0;i<size;++i)
        {
            comp sig_i = sigma_p(s,qvec[i],m);
            comp sig_j = sigb;

            comp m2k_i = M2kfunc(a,sig_i,m, eps_for_m2k);
            comp m2k_j = M2kfunc(a,sigb,m, eps_for_m2k);
            Ds(i) = m2k_i*dsol(i)*m2k_j;
            for(int j=0;j<size;++j)
            {
                sig_j = sigma_p(s,qvec[j],m);
                m2k_j = M2kfunc(a,sig_j,m, eps_for_m2k);
        
                //Dmats(i,j) = m2k_i*dmat(i,j)*m2k_j;
            }
        }

    
        //solution for D_dagger
        Eigen::MatrixXcd Cmat(size,size);
        Eigen::VectorXcd Bvec(size);
        Eigen::VectorXcd Ddagsol(size);

        
        double relerror1; 
        //struct timeval stop, start;
        gettimeofday(&start, NULL);

        Bvec_Dsecsht_maker_with_weights(Bvec,Ds,qvec,weights,s,qval,m,a,eps);
        cmat_maker_with_weights(Cmat,Dmats,qvec,weights,s,m,a,eps);

        LinearSolver_2(Cmat,Ddagsol,Bvec,relerror1);
    
        
        gettimeofday(&stop, NULL);
        actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/1000000.0;
        cout<<"time taken to solve = "<<actual_sol_time<<endl;

        double temp_val=0.0;
        int ind_q = 0;
        comp selected_val = 0.0;

        for(int j=0;j<size;++j)
        {
            comp some_p = qvec[j];
            //cout<<"q = "<<qvec[i]<<'\t'<<"Ds = "<<Ds(i)<<'\t'<<"D_dag_s = "<<Ddagsol(i)<<endl;
            if(j==0)
            {
                temp_val = abs(some_p - qval);
                ind_q = j;
                continue; 
            }
            
            if(abs(some_p - qval)<temp_val)
            {
                temp_val = abs(some_p - qval);
                ind_q = j;
                selected_val = some_p;
            }
        
        
        }

        cout<<"run = "<<i+1<<endl;
        cout<<"qval = "<<qval<<'\t'<<"selected_val = "<<selected_val<<endl;
        cout<<"diff = "<<abs(qval - selected_val)<<endl;
        cout<<"s = "<<s<<'\t'<<"Ds = "<<Ds(ind_q)<<'\t'<<"D_dag_s = "<<Ddagsol(ind_q)<<endl;
        fout<<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<real(Ds(ind_q))<<'\t'
            <<imag(Ds(ind_q))<<'\t'
            <<real(Ddagsol(ind_q))<<'\t'
            <<imag(Ddagsol(ind_q))<<endl;
        cout<<"======================================"<<endl;

    }
    fout.close();  
}

void test_functions_4()
{
    comp ii = {0.0,1.0};
    double a = -8;
    double m = 1.0;

    //comp s = {8.8,0.0001};

    double simag = 0.0;

    double sreal_initial = 8.99;
    double sreal_final = 9.01;
    double spoints = 500.0;
    double dels = abs(sreal_initial - sreal_final)/spoints;

    int N = 500;

    string filename = "dS_discontinuity_test.dat";
    ofstream fout;
    fout.open(filename.c_str());

    for(int i=0;i<spoints+1;++i)
    {
        comp s = sreal_initial + i*dels + ii*simag;
        comp sig_q = 2.0*m*m;
        comp q = pmom(s,sig_q, m);

        
        double eps_for_m2k = 0.0;
        double eps = 0.0;
        double eps_for_s = 0.000001;

        comp s_plus = s + ii*eps_for_s;
        comp s_minus = s - ii*eps_for_s;
        //int N = 500;

        //solution for Ds
        comp kmin = 0.0;
        comp kmax = pmom(s,0.0,m);
        comp kmax_plus = pmom(s_plus, 0.0, m);
        comp kmax_minus = pmom(s_minus, 0.0, m);

        comp sigb = 2.0*m*m;
        comp qval = pmom(s,sigb,m);
        comp qval_plus = pmom(s_plus, sigb, m);
        comp qval_minus = pmom(s_minus, sigb, m);



        vector<comp> qvec_plus;
        vector<comp> weights_plus;
        vector<comp> qvec_minus;
        vector<comp> weights_minus;

        line_maker_with_weights(qvec_plus,weights_plus,kmin,kmax_plus,N);
        line_maker_with_weights(qvec_minus,weights_minus,kmin,kmax_minus,N);


        int size = qvec.size();
        Eigen::MatrixXcd Bmat(size,size);
        Eigen::VectorXcd Gvec(size);
        Eigen::VectorXcd dsol(size);

        Gvec_maker_momrep(Gvec,s_plus,qvec_plus,qval_plus,m,eps);
        Gvec = -1.0*Gvec;
        
        
        Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec_plus,qvec_plus,weights_plus,a,m,eps,eps_for_m2k);

        double relerror; 
        struct timeval stop, start;
        gettimeofday(&start, NULL);
        LinearSolver_2(Bmat,dsol,Gvec,relerror);
        
        
        gettimeofday(&stop, NULL);
        double actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/1000000.0;
        cout<<"time taken to solve = "<<actual_sol_time<<endl;

        

    
        
        double relerror1; 
        //struct timeval stop, start;
        gettimeofday(&start, NULL);

        
        
        gettimeofday(&stop, NULL);
        actual_sol_time = ((stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec)/1000000.0;
        cout<<"time taken to solve = "<<actual_sol_time<<endl;

        

        cout<<"run = "<<i+1<<endl;
        cout<<"qval = "<<qval<<'\t'<<"selected_val = "<<selected_val<<endl;
        cout<<"diff = "<<abs(qval - selected_val)<<endl;
        cout<<"s = "<<s<<'\t'<<"Ds = "<<Ds(ind_q)<<'\t'<<"D_dag_s = "<<Ddagsol(ind_q)<<endl;
        fout<<real(s)<<'\t'
            <<imag(s)<<'\t'
            <<real(Ds(ind_q))<<'\t'
            <<imag(Ds(ind_q))<<'\t'
            <<real(Ddagsol(ind_q))<<'\t'
            <<imag(Ddagsol(ind_q))<<endl;
        cout<<"======================================"<<endl;

    }
    fout.close();  
}



int main()
{
    //this is where first tested D_dag for different k values
    //test_functions();

    //this is where we first tested D_dag for different s on the real axis
    //test_functions_2();

    //this is where we test D_dag accross the imaginary s axis to see
    //if there is a discontinuity present or not
    //test_functions_3();
    
    
    test_functions_4();
    
    return 0;
}