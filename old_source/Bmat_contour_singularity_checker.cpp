//macos:g++-12 ../testing_on_mac.cpp -o printer -O3 -std=c++14 -I /usr/local/Cellar/eigen/3.4.0_1/include/eigen3
//wahab:g++ printer.cpp -o printer -O3 -std=c++14 -I /cm/shared/applications/eigen/3.3.7/include/eigen3/
//wahab_gpu:nvcc -I /cm/shared/applications/cuda-toolkit/10.0.130/include -I /cm/shared/applications/eigen/3.3.7/include/eigen3/ printer_gpu_momrep.cpp  -O3 -std=c++14 -L /cm/shared/applications/cuda-toolkit/10.0.130/lib64 -lcublas -lcusolver -lcudart -o printer_gpu
//wahab_gpu_omp:nvcc -I /cm/shared/applications/cuda-toolkit/10.0.130/include -I /cm/shared/applications/eigen/3.3.7/include/eigen3/ printer_gpu_momrep.cpp  -O3 -std=c++14 -L /cm/shared/applications/cuda-toolkit/10.0.130/lib64 -lcublas -lcusolver -lcudart -lgomp -Xcompiler -fopenmp -o printer_gpu
//ubuntupc_omp:nvcc117 -I/usr/include/eigen3/ -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/include -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/lib64 printer_gpu_momrep.cpp  -O3 -std=c++14 -lcublas -lcusolver -lcudart -lgomp -Xcompiler -fopenmp -o printer_gpu
//ubuntupc_omp1:nvcc -I/usr/include/eigen3/ -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/include -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/lib64 printer_gpu_momrep.cpp  -O3 -std=c++14 -lcublas -lcusolver -lcudart -lgomp -Xcompiler -fopenmp -o printer_gpu
//ubuntu cpu: g++ ../Bmat_contour_singularity_checker.cpp -o contoursing -O3 -std=c++14 -I /usr/include/eigen3

#include<bits/stdc++.h>
#include "functions_momrep_based.h"
#include "momentum_vector_maker.h"
#include "integralequation_momrep_based.h"
#include "interpolator_momrep.h"
#include "solvers.h"
#include<omp.h>
#include <Eigen/Eigenvalues> //for vertex factor calculations
//#include "LUSolver_modified.h"
#include "fastGL_based_functions.h"
#include "SA_method_functions.h"
#include <sys/time.h>

using namespace std;

typedef complex<double> comp;

//the 4 roots of k in the argument of the logarithmic function 
void Bmat_singularity_for_contour_11(   vector<comp> qvec,
                                        vector<comp> &singularity_contour,
                                        comp s, 
                                        double m )
{
    for(int i=0;i<qvec.size();++i)
    {
        comp p = qvec[i];
        comp C = s - 2.0*sqrt(s)*omega_comp(qvec[i],m) + m*m; 
        comp gamma = 2.0*p*(omega_comp(p,m) - sqrt(s));
        comp beta = m*m*pow(omega_comp(p,m) - sqrt(s),2.0) - C*C/4.0;
        comp alpha = p*p + pow(omega_comp(p,m) - sqrt(s),2.0);

        comp num =  - (m*m*gamma*gamma  - 2.0*alpha*beta) 
                    + sqrt(pow(m*m*gamma*gamma - 2.0*alpha*beta,2.0)
                    - 4.0*(gamma*gamma - alpha*alpha)*beta*beta);
        comp denom = 2.0*(gamma*gamma - alpha*alpha);

        singularity_contour.push_back(+sqrt(num/denom));
    }
}

void Bmat_singularity_for_contour_12(   vector<comp> qvec,
                                        vector<comp> &singularity_contour,
                                        comp s, 
                                        double m )
{
    for(int i=0;i<qvec.size();++i)
    {
        comp p = qvec[i];
        comp C = s - 2.0*sqrt(s)*omega_comp(qvec[i],m) + m*m; 
        comp gamma = 2.0*p*(omega_comp(p,m) - sqrt(s));
        comp beta = m*m*pow(omega_comp(p,m) - sqrt(s),2.0) - C*C/4.0;
        comp alpha = p*p + pow(omega_comp(p,m) - sqrt(s),2.0);

        comp num =  - (m*m*gamma*gamma  - 2.0*alpha*beta) 
                    - sqrt(pow(m*m*gamma*gamma - 2.0*alpha*beta,2.0)
                    - 4.0*(gamma*gamma - alpha*alpha)*beta*beta);
        comp denom = 2.0*(gamma*gamma - alpha*alpha);

        singularity_contour.push_back(+sqrt(num/denom));
    }
}

void Bmat_singularity_for_contour_21(   vector<comp> qvec,
                                        vector<comp> &singularity_contour,
                                        comp s, 
                                        double m )
{
    for(int i=0;i<qvec.size();++i)
    {
        comp p = qvec[i];
        comp C = s - 2.0*sqrt(s)*omega_comp(qvec[i],m) + m*m; 
        comp gamma = 2.0*p*(omega_comp(p,m) - sqrt(s));
        comp beta = m*m*pow(omega_comp(p,m) - sqrt(s),2.0) - C*C/4.0;
        comp alpha = p*p + pow(omega_comp(p,m) - sqrt(s),2.0);

        comp num =  - (m*m*gamma*gamma  - 2.0*alpha*beta) 
                    + sqrt(pow(m*m*gamma*gamma - 2.0*alpha*beta,2.0)
                    - 4.0*(gamma*gamma - alpha*alpha)*beta*beta);
        comp denom = 2.0*(gamma*gamma - alpha*alpha);

        singularity_contour.push_back(-sqrt(num/denom));
    }
}

void Bmat_singularity_for_contour_22(   vector<comp> qvec,
                                        vector<comp> &singularity_contour,
                                        comp s, 
                                        double m )
{
    for(int i=0;i<qvec.size();++i)
    {
        comp p = qvec[i];
        comp C = s - 2.0*sqrt(s)*omega_comp(qvec[i],m) + m*m; 
        comp gamma = 2.0*p*(omega_comp(p,m) - sqrt(s));
        comp beta = m*m*pow(omega_comp(p,m) - sqrt(s),2.0) - C*C/4.0;
        comp alpha = p*p + pow(omega_comp(p,m) - sqrt(s),2.0);

        comp num =  - (m*m*gamma*gamma  - 2.0*alpha*beta) 
                    - sqrt(pow(m*m*gamma*gamma - 2.0*alpha*beta,2.0)
                    - 4.0*(gamma*gamma - alpha*alpha)*beta*beta);
        comp denom = 2.0*(gamma*gamma - alpha*alpha);

        singularity_contour.push_back(-sqrt(num/denom));
    }
}

void singularity_check()
{
    double a = 16.0;//-8.1;
    double m = 1.0;

    double sreal = 8.605;
    double simag = -0.005;//1.0e-5;

    comp s = sreal + ii*simag;

    double eps = 0.0;
    double eps_for_m2k = 0.0;

    int N = 500;

    comp kmin = 0.0;
    comp sigk = 0.0;//2.0*m*m;
    comp kmax = pmom(s,sigk,m);

    vector<comp> qvec;
    vector<comp> weights;
    int tag_for_m2k_1;
    int tag_for_m2k_2; 

    //previous analytic continuation paper
    int tag1, tag2, switch_for_gvec_fixer;
    //for im(s)>0
    //mom_vector_maker_seba_imspos_5(qvec,weights,s,kmin,kmax,a,m,eps,eps,N,tag1,tag2,switch_for_gvec_fixer);
    //for im(s)<0
    mom_vector_maker_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,N);
                
    //contour for third sheet of dS     
    //contour_for_resonance_7(qvec,weights,a,s,m,kmin,kmax,eps_for_m2k,N,tag_for_m2k_1,tag_for_m2k_2);

    cout<<"tag_for_m2k 1 = "<<tag_for_m2k_1<<endl;
    cout<<"tag_for_m2k 2 = "<<tag_for_m2k_2<<endl;
    
    int size = qvec.size();

    vector<comp> singularity_qvec_11;
    vector<comp> singularity_qvec_12;
    vector<comp> singularity_qvec_21;
    vector<comp> singularity_qvec_22;

    Bmat_singularity_for_contour_11(qvec,singularity_qvec_11,s,m);
    Bmat_singularity_for_contour_12(qvec,singularity_qvec_12,s,m);
    Bmat_singularity_for_contour_21(qvec,singularity_qvec_21,s,m);
    Bmat_singularity_for_contour_22(qvec,singularity_qvec_22,s,m);

    ofstream fout;
    string filename = "Bmat_contour_singularity_check_prevcheck.dat";
    fout.open(filename.c_str());

    for(int i=0;i<size;++i)
    {
        fout<<real(qvec[i])<<'\t'
            <<imag(qvec[i])<<'\t'
            <<real(singularity_qvec_11[i])<<'\t'
            <<imag(singularity_qvec_11[i])<<'\t'
            <<real(singularity_qvec_12[i])<<'\t'
            <<imag(singularity_qvec_12[i])<<'\t'
            <<real(singularity_qvec_21[i])<<'\t'
            <<imag(singularity_qvec_21[i])<<'\t'
            <<real(singularity_qvec_22[i])<<'\t'
            <<imag(singularity_qvec_22[i])<<endl;
    }
    fout.close();


}

int main()
{
    singularity_check();
    return 0;
}

