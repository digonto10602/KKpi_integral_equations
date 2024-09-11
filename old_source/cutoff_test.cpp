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
#include "solvers.h"
#include<omp.h>
#include <Eigen/Eigenvalues> //for vertex factor calculations
//#include "LUSolver_modified.h"
#include "fastGL_based_functions.h"
#include "SA_method_functions.h"
#include <sys/time.h>

using namespace std;

typedef complex<double> comp;


comp s0_based_Jfunc_comp_smooth(    double s0, 
                                    comp x  )
{
    if(real(x)<=0.0) return 0.0;
    else if(real(x)>0.0 && real(x)<s0) return exp((-1.0/x)*exp(-1.0/(1.0-x)));
    else return 1.0;
}

double s0_based_Jfunc_comp_smooth_1(    double s0, 
                                        comp x  )
{
    if(real(x)<=0.0) return 0.0;
    else if(real(x)>0.0 && real(x)<s0) return exp((-1.0/real(x))*exp(-1.0/(1.0-real(x))));
    else return 1.0;
}

comp s0_based_Jfunc_comp_hard(      double s0, 
                                    comp x  )
{
    if(real(x)<=s0) return 0.0;
    //else if(real(x)>0.0 && real(x)<s0) return exp((-1.0/x)*exp(-1.0/(1.0-x)));
    else return 1.0;
}

comp s0_based_Hfunc_comp_smooth(    comp sigmap,
                                    double s0,
                                    double m    )
{
    return s0_based_Jfunc_comp_smooth(s0, sigmap/(4.0*m*m));
}

double s0_based_Hfunc_comp_smooth_1(    comp sigmap,
                                    double s0,
                                    double m    )
{
    return s0_based_Jfunc_comp_smooth_1(s0, sigmap/(4.0*m*m));
}

comp s0_based_Hfunc_comp_hard(      comp sigmap,
                                    double s0,
                                    double m    )
{
    return s0_based_Jfunc_comp_hard(s0, sigmap/(4.0*m*m));
}


comp K2_scattlength(    comp s,
                        double scattering_length)
{
    comp pi = acos(-1.0);

    return -8.0*pi*sqrt(s)*scattering_length;
}

comp G_00_smooth(   comp s,
                    double px,
                    double py,
                    double pz,
                    double kx,
                    double ky,
                    double kz,
                    double s0,
                    double m    )

{
    comp p = sqrt(px*px + py*py + pz*pz);
    comp k = sqrt(kx*kx + ky*ky + kz*kz);
    comp p_plus_k = sqrt(       px*px + py*py + pz*pz
                            +   kx*kx + ky*ky + kz*kz 
                            +  2.0*px*kx + 2.0*py*ky + 2.0*pz*kz    );
    comp s2p = sigma_p(s,p,m);
    comp s2k = sigma_p(s,k,m);

    comp cutoff1 = s0_based_Hfunc_comp_smooth(s2p,s0,m);
    comp cutoff2 = s0_based_Hfunc_comp_smooth(s2k,s0,m);

    comp omegap = omega_comp(p,m);
    comp omegak = omega_comp(k,m);
    comp omegakp = omega_comp(p_plus_k,m);

    comp res = cutoff1*cutoff2/(2.0*omegakp*(sqrt(s) - omegak - omegap - omegakp));

    return res;

}

comp G_00_hard(   comp s,
                    double px,
                    double py,
                    double pz,
                    double kx,
                    double ky,
                    double kz,
                    double s0,
                    double m    )

{
    comp p = sqrt(px*px + py*py + pz*pz);
    comp k = sqrt(kx*kx + ky*ky + kz*kz);
    comp p_plus_k = sqrt(       px*px + py*py + pz*pz
                            +   kx*kx + ky*ky + kz*kz 
                            +  2.0*px*kx + 2.0*py*ky + 2.0*pz*kz    );
    comp s2p = sigma_p(s,p,m);
    comp s2k = sigma_p(s,k,m);

    comp cutoff1 = s0_based_Hfunc_comp_hard(s2p,s0,m);
    comp cutoff2 = s0_based_Hfunc_comp_hard(s2k,s0,m);

    comp omegap = omega_comp(p,m);
    comp omegak = omega_comp(k,m);
    comp omegakp = omega_comp(p_plus_k,m);

    comp res = cutoff1*cutoff2/(2.0*omegakp*(sqrt(s) - omegak - omegap - omegakp));

    return res;

}

comp F_ieps_00_smooth(  comp s,
                        double kx,
                        double ky,
                        double kz, 
                        double s0,
                        double m,
                        double epsilon,
                        double L,
                        double kpoints,
                        double n_max    )
{
    comp ii = {0.0,1.0};
    comp res_I = {0.0,0.0};
    double pi = acos(-1.0);

    double kmin = 0.0;
    comp kmax = pmom(s,0.0,m);
    double real_kmax = real(kmax);
    //double kpoints = 500.0;
    double delk = abs(real_kmax - kmin)/kpoints;
    for(int nx=0;nx<n_max;++nx)
    {
        for(int ny=0;ny<n_max;++ny)
        {
            for(int nz=0;nz<n_max;++nz)
            {
                double n = sqrt(nx*nx + ny*ny + nz*nz);
                if(n<=0) continue;
                if(n>=n_max) continue;
                
                for(double ax=kmin;ax<=real_kmax;ax=ax+delk)
                {
                    for(double ay=kmin;ay<=real_kmax;ay=ay+delk)
                    {
                        for(double az=kmin;az<=real_kmax;az=az+delk)
                        {
                            double amom = sqrt(ax*ax + ay*ay + az*az);

                            double phasespace = delk*delk*delk/pow(2.0*pi,3.0);

                            double amom_n = ax*nx + ay*ny + az*nz;
                            comp expon = exp(ii*amom_n*L);

                            double k = sqrt(kx*kx + ky*ky + kz*kz);
                            comp sigk = sigma_p(s,k,m);
                            comp siga = sigma_p(s,amom,m);

                            double k_plus_a = sqrt(     kx*kx + ky*ky + kz*kz 
                                                    +   ax*ax + ay*ay + az*az   
                                                    +  2.0*kx*ax + 2.0*ky*ay + 2.0*kz*az    );
                            comp sigka = sigma_p(s,k_plus_a,m);

                            comp cutoffs = s0_based_Hfunc_comp_smooth(sigk,s0,m)
                                            *s0_based_Hfunc_comp_smooth(siga,s0,m)
                                            *s0_based_Hfunc_comp_smooth(sigka,s0,m);

                            comp omegak = omega_comp(k,m);
                            comp omegaa = omega_comp(amom, m);
                            comp omegaka = omega_comp(k_plus_a, m);

                            comp denom = 2.0*omegaa*2.0*omegaka*(sqrt(s) - omegak - omegaa - omegaka + ii*epsilon);
                            
                            /*
                            cout<<"nx = "<<nx<<'\t'<<"ny = "<<ny<<'\t'<<"nz = "<<nz<<endl;
                            cout<<"n = "<<n<<endl;
                            
                            cout<<"ax = "<<ax<<'\t'<<"ay = "<<ay<<'\t'<<"az = "<<az<<endl;
                            cout<<"a = "<<amom<<endl;

                            cout<<"exponent = "<<expon<<'\t'<<" cutoffs = "<<cutoffs<<endl;
                            */
                            comp tot = phasespace * expon * cutoffs / denom; 
                            
                            //cout<<"result = "<<tot<<endl;

                            //cout<<"====================================="<<endl;

                            res_I = res_I + tot;
                        }
                    }
                }
            }
        }
    }

    return res_I*0.5;
}


comp F_ieps_00_hard(    comp s,
                        double kx,
                        double ky,
                        double kz, 
                        double s0,
                        double m,
                        double epsilon,
                        double L,
                        double kpoints,
                        double n_max    )
{
    comp ii = {0.0,1.0};
    comp res_I = {0.0,0.0};
    double pi = acos(-1.0);

    double kmin = 0.0;
    comp kmax = pmom(s,0.0,m);
    double real_kmax = real(kmax);
    //double kpoints = 500.0;
    double delk = abs(real_kmax - kmin)/kpoints;
    for(int nx=0;nx<n_max;++nx)
    {
        for(int ny=0;ny<n_max;++ny)
        {
            for(int nz=0;nz<n_max;++nz)
            {
                double n = sqrt(nx*nx + ny*ny + nz*nz);
                if(n<=0) continue;
                if(n>=n_max) continue;
                
                for(double ax=kmin;ax<=real_kmax;ax=ax+delk)
                {
                    for(double ay=kmin;ay<=real_kmax;ay=ay+delk)
                    {
                        for(double az=kmin;az<=real_kmax;az=az+delk)
                        {
                            double amom = sqrt(ax*ax + ay*ay + az*az);

                            double phasespace = delk*delk*delk/pow(2.0*pi,3.0);

                            double amom_n = ax*nx + ay*ny + az*nz;
                            comp expon = exp(ii*amom_n*L);

                            double k = sqrt(kx*kx + ky*ky + kz*kz);
                            comp sigk = sigma_p(s,k,m);
                            comp siga = sigma_p(s,amom,m);

                            double k_plus_a = sqrt(     kx*kx + ky*ky + kz*kz 
                                                    +   ax*ax + ay*ay + az*az   
                                                    +  2.0*kx*ax + 2.0*ky*ay + 2.0*kz*az    );
                            comp sigka = sigma_p(s,k_plus_a,m);

                            comp cutoffs = s0_based_Hfunc_comp_hard(sigk,s0,m)
                                            *s0_based_Hfunc_comp_hard(siga,s0,m)
                                            *s0_based_Hfunc_comp_hard(sigka,s0,m);

                            comp omegak = omega_comp(k,m);
                            comp omegaa = omega_comp(amom, m);
                            comp omegaka = omega_comp(k_plus_a, m);

                            comp denom = 2.0*omegaa*2.0*omegaka*(sqrt(s) - omegak - omegaa - omegaka + ii*epsilon);
                            
                            /*
                            cout<<"nx = "<<nx<<'\t'<<"ny = "<<ny<<'\t'<<"nz = "<<nz<<endl;
                            cout<<"n = "<<n<<endl;
                            
                            cout<<"ax = "<<ax<<'\t'<<"ay = "<<ay<<'\t'<<"az = "<<az<<endl;
                            cout<<"a = "<<amom<<endl;

                            cout<<"exponent = "<<expon<<'\t'<<" cutoffs = "<<cutoffs<<endl;
                            */
                            comp tot = phasespace * expon * cutoffs / denom; 
                            
                            //cout<<"result = "<<tot<<endl;

                            //cout<<"====================================="<<endl;

                            res_I = res_I + tot;
                        }
                    }
                }
            }
        }
    }

    return res_I*0.5;
}

comp rho_bar(   comp s2k, 
                double m    )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    if(real(s2k)>4.0*m*m)
    {
        comp num = -ii*sqrt(s2k/4.0 - m*m);
        comp denom = 16.0*pi*sqrt(s2k);

        return num/denom;
    }
    else 
    {
        comp num = abs(sqrt(s2k/4.0 - m*m));
        comp denom = 16.0*pi*sqrt(s2k);

        return num/denom; 
    }
}

comp rho00_smooth(  comp s,
                    double kx, 
                    double ky, 
                    double kz,
                    double s0, 
                    double m    )
{
    comp k = sqrt(kx*kx + ky*ky + kz*kz);

    comp s2k = sigma_p(s,k,m);

    comp Hk = s0_based_Hfunc_comp_smooth(s2k, s0, m);

    return Hk*rho_bar(s2k,m);
}

comp rho00_hard(    comp s,
                    double kx, 
                    double ky, 
                    double kz,
                    double s0, 
                    double m    )
{
    comp k = sqrt(kx*kx + ky*ky + kz*kz);

    comp s2k = sigma_p(s,k,m);

    comp Hk = s0_based_Hfunc_comp_hard(s2k, s0, m);

    return Hk*rho_bar(s2k,m);
}

comp F00_smooth(    comp s,
                    double kx,
                    double ky,
                    double kz, 
                    double s0,
                    double m,
                    double epsilon,
                    double L,
                    double kpoints,
                    double n_max  )
{
    return F_ieps_00_smooth(s, kx, ky, kz, s0, m, epsilon, L, kpoints, n_max) + rho00_smooth(s, kx, ky, kz, s0, m);
}

comp F00_hard(      comp s,
                    double kx,
                    double ky,
                    double kz, 
                    double s0,
                    double m,
                    double epsilon,
                    double L,
                    double kpoints,
                    double n_max  )
{
    return F_ieps_00_hard(s, kx, ky, kz, s0, m, epsilon, L, kpoints, n_max) + rho00_hard(s, kx, ky, kz, s0, m);
}

comp F3_smooth( comp s,
                double scattering_length, 
                double kx,
                double ky,
                double kz, 
                double s0,
                double m,
                double epsilon,
                double L,
                double kpoints,
                double n_max    )
{
    comp F = F00_smooth(s, kx, ky, kz, s0, m, epsilon, L, kpoints, n_max);
    comp k = sqrt(kx*kx + ky*ky + kz*kz);
    comp omegak = omega_comp(k, m);

    comp G = G_00_smooth(s, kx, ky, kz, kx, ky, kz, s0, m);

    comp K2 = K2_scattlength(s, scattering_length);

    comp firstterm = F/(2.0*omegak*L);
    double secondterm = -2.0/3.0;
    comp thirdterm_inv = 1.0/(1.0 + K2*G);

    comp forthterm = 1.0/(1.0 + thirdterm_inv*K2*F);

    return firstterm*(secondterm + forthterm);

}

comp F3_hard(   comp s,
                double scattering_length, 
                double kx,
                double ky,
                double kz, 
                double s0,
                double m,
                double epsilon,
                double L,
                double kpoints,
                double n_max    )
{
    comp F = F00_hard(s, kx, ky, kz, s0, m, epsilon, L, kpoints, n_max);
    comp k = sqrt(kx*kx + ky*ky + kz*kz);
    comp omegak = omega_comp(k, m);

    comp G = G_00_hard(s, kx, ky, kz, kx, ky, kz, s0, m);

    comp K2 = K2_scattlength(s, scattering_length);

    comp firstterm = F/(2.0*omegak*L);
    double secondterm = -2.0/3.0;
    comp thirdterm_inv = 1.0/(1.0 + K2*G);

    comp forthterm = 1.0/(1.0 + thirdterm_inv*K2*F);

    return firstterm*(secondterm + forthterm);

}


//int main()
void cutoff_test()
{

    // Inputs
    double pi = acos(-1.0);
    double a = 2.0;
    double s = 9.2;
    double kx = 0.0;
    double ky = 0.0;
    double kz = 0.0;
    comp k = sqrt(kx*kx + ky*ky + kz*kz);
    //double s0 = 1.0;
    double m = 1.0;
    double epsilon = 0.001;
    //double L = 4.0;
    double kpoints = 100.0;
    double n_max = 10.0;

    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb,m);
    cout<<"pole at q = "<<q<<endl;

    ofstream fout;
    int s0count = 0;
    for(double s0=0.90;s0>0.25;s0=s0-0.25)
    {
        string filename = "F3_check_s0_" + to_string(s0count) + ".dat";
        fout.open(filename.c_str());
        for(int L=2;L<20;++L)
        {
            //comp F = F_ieps_00_smooth(  s, kx, ky, kz, s0, m, epsilon, L, kpoints, n_max);


            comp F3smooth = F3_smooth(s, a, kx, ky, kz, s0, m, epsilon, L, kpoints, n_max);
            comp F3hard   = F3_hard  (s, a, kx, ky, kz, s0, m, epsilon, L, kpoints, n_max);
        
            comp delF3 = F3smooth - F3hard;
            cout<<"s0="<<s0<<'\t'<<"L="<<L<<'\t'<<"F3smooth="<<F3smooth<<'\t'<<"F3hard="<<F3hard<<'\t'<<"delF3="<<delF3<<endl;
            fout<<a<<'\t'<<s0<<'\t'<<real(k)<<'\t'<<imag(k)<<'\t'<<L<<'\t'
                <<real(F3smooth)<<'\t'<<imag(F3smooth)<<'\t'
                <<real(F3hard)<<'\t'<<imag(F3hard)<<'\t'
                <<real(delF3)<<'\t'<<imag(delF3)<<endl;
        }
        fout.close();
        s0count = s0count + 1;
    }
    //return 0;
}

void testing_hardcutoff()
{
    // Inputs
    double pi = acos(-1.0);
    double a = 2.0;
    double s = 8.5;
    double kx = 0.0;
    double ky = 0.0;
    double kz = 0.0;
    comp k = sqrt(kx*kx + ky*ky + kz*kz);

    
    //double s0 = 1.0;
    double m = 1.0;
    comp sigk = sigma_p(s,k,m);


    double epsilon = 0.001;
    //double L = 4.0;
    double kpoints = 100.0;
    double n_max = 10.0;

    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb,m);
    cout<<"pole at q = "<<q<<endl;

    comp cutoff = s0_based_Hfunc_comp_hard(sigk, 0.9, m);

    cout<<"cutoff = "<<cutoff<<endl;

    //return 0;
}

void F3_vs_s()
{
    //Inputs
    double pi = acos(-1.0);
    double a = -10.0;
    double kx = 0.0;
    double ky = 0.0;
    double kz = 0.0;
    comp k = sqrt(kx*kx + ky*ky + kz*kz);

    double m = 1.0;
    double epsilon = 0.001;
    double L = 6.0;
    double kpoints = 100.0;
    double n_max = 10.0;
    double s0 = 0.90;
    double sinitial = 6.25;
    double sfinal = 20.25;
    double spoints = 200.0;
    double dels = abs(sfinal - sinitial)/spoints; 

    ofstream fout;
    string filename = "F3_vs_s_for_a_"+to_string(a)+"_L_"+to_string((int)L)+"_s0_"+to_string(s0)+".dat";
    fout.open(filename.c_str());

    int scount = 0;
    for(int i=0;i<spoints+1;++i)
    {
        double s = sinitial + i*dels;

        comp F3smooth = F3_smooth(s, a, kx, ky, kz, s0, m, epsilon, L, kpoints, n_max);
        //comp F3hard   = F3_hard  (s, a, kx, ky, kz, s0, m, epsilon, L, kpoints, n_max);
        
        //comp delF3 = F3smooth - F3hard;

        cout<<"s = "<<s<<endl;
        cout<<"F3smooth = "<<F3smooth<<endl;
        //cout<<"F3hard = "<<F3hard<<endl;
        cout<<"cutoff_s0 = "<<s0<<endl;
        cout<<"========================="<<endl;
        cout<<endl;

        fout<<setprecision(16)
            <<s<<'\t'<<a<<'\t'<<s0<<'\t'
            <<real(k)<<'\t'<<imag(k)<<'\t'<<L<<'\t'
            <<real(F3smooth)<<'\t'<<imag(F3smooth)<<endl;
            //<<real(F3hard)<<'\t'<<imag(F3hard)<<'\t'
            //<<real(delF3)<<'\t'<<imag(delF3)<<endl;
            
        cout<<"run = "<<scount+1<<endl;
        scount = scount + 1;
    }
    fout.close();
}


// This portion of the code is made based 
// on https://arxiv.org/pdf/1803.04169.pdf


//in these functions we assume p=k
comp onebyomega_Kmat_smooth(    double scat_length,
                                comp s,
                                comp k,
                                double s0,
                                double m  )
{
    double pi = acos(-1.0);
    comp sigk = sigma_p(s,k,m);
    comp omegak = omega_comp(k,m);

    comp Hk = s0_based_Hfunc_comp_smooth(sigk,s0,m);

    comp q2k = sqrt(sigk/4.0 - m*m);

    return (-1.0/scat_length + abs(q2k)*(1.0 - Hk))/(32.0*pi*omegak*sqrt(sigk));

}

comp Gs_00_smooth(  comp s,
                    comp p,
                    comp k,
                    double m,
                    double s0,
                    double L    )
{
    comp sigp = sigma_p(s,p,m);
    comp sigk = sigma_p(s,k,m);
    comp Hk = s0_based_Hfunc_comp_smooth(sigk,s0,m);
    comp Hp = s0_based_Hfunc_comp_smooth(sigp,s0,m);

    comp omegak = omega_comp(k,m);
    comp omegap = omega_comp(p,m);
    comp omegakp = omega_comp(p+k,m);

    return (Hk*Hp)/(8.0*L*L*L*omegak*omegap*omegakp*(sqrt(s) - omegak - omegap - omegakp));
}

comp Fs_smooth_finite(  comp s,
                        comp k,
                        double m,
                        double s0,
                        double n_max,
                        double L)
{
    double pi = acos(-1.0);
    comp sigk = sigma_p(s,k,m);
    comp Hk = s0_based_Hfunc_comp_smooth(sigk,s0,m);
    comp omegak = omega_comp(k,m);

    comp summ = {0.0,0.0};
    for(int nx=0;nx<n_max;++nx)
    {
        for(int ny=0;ny<n_max;++ny)
        {
            for(int nz=0;nz<n_max;++nz)
            {
                double n = sqrt(nx*nx + ny*ny + nz*nz);

                comp az = (2.0*pi/L)*nz; 
                comp a = (2.0*pi/L)*n;

                comp omegaa = omega_comp(a,m);

                comp kplusa = k*k + a*a + 2.0*k*az;

                comp omegaka = omega_comp(kplusa,m);

                comp denom = 4.0*omegaa*omegaka*(sqrt(s)- omegak - omegaa - omegaka);

                comp res = (Hk/(4.0*omegak*L*L*L))*(1.0/denom);
                summ = summ + res;
            }
        }
    }

    return summ;
}

comp Fs_smooth_infinite(  comp s,
                        comp k,
                        double m,
                        double s0,
                        double points)
{
    double pi = acos(-1.0);
    comp sigk = sigma_p(s,k,m);
    comp Hk = s0_based_Hfunc_comp_smooth_1(sigk,s0,m);
    comp omegak = omega_comp(k,m);

    comp summ = {0.0,0.0};
    comp a_min = 0.0;
    comp a_max = pmom(s,0.0,m);
    vector<comp> qvec;
    vector<comp> weights;
    line_maker_with_weights(qvec,weights,a_min,a_max,points);

    for(int nx=0;nx<qvec.size();++nx)
    {
        for(int ny=0;ny<qvec.size();++ny)
        {
            for(int nz=0;nz<qvec.size();++nz)
            {
                comp ax = qvec[nx];
                comp ay = qvec[ny];
                comp az = qvec[nz];

                comp dax = weights[nx];
                comp day = weights[ny];
                comp daz = weights[nz];

                comp a = sqrt(ax*ax + ay*ay + az*az);

                comp omegaa = omega_comp(a,m);

                comp kplusa = k*k + a*a + 2.0*k*az;

                comp omegaka = omega_comp(kplusa,m);

                comp denom = 4.0*omegaa*omegaka*(sqrt(s)- omegak - omegaa - omegaka);
                //comp denom = 4.0*omegaa*omegaka;//*(sqrt(s)- omegak - omegaa - omegaka);

                //comp res = (Hk/(4.0*omegak*2.0*pi*2.0*pi*2.0*pi))*dax*day*daz*(1.0/denom);
                comp res = (1.0/(4.0*omegak*2.0*pi*2.0*pi*2.0*pi))*dax*day*daz*(1.0/denom);
                
                double check_val = abs(sqrt(s)- omegak - omegaa - omegaka);
                double tolerance = 1.0e-6;
                if(check_val>tolerance)
                {
                    summ = summ + res;
                }
                else 
                {
                    summ = summ;
                }
            }
        }
    }

    return summ;
}

comp Fs_smooth( comp s,
                comp k,
                double m,
                double s0,
                double n_max,
                double L, 
                double points   )

{
    comp Fsfinite = Fs_smooth_finite(s,k,m,s0,n_max,L);
    comp Fsinfinite = Fs_smooth_infinite(s,k,m,s0,points);

    return Fsfinite - Fsinfinite;
}

comp F3s_smooth(    double scat_length,
                    comp s,
                    comp p,
                    comp k,
                    double m,
                    double s0,
                    double n_max,
                    double L, 
                    double points   )
{
    comp Fs00 = Fs_smooth(s,k,m,s0,n_max,L,points);

    comp Gs00 = Gs_00_smooth(s,p,k,m,s0,L);

    comp onebKmat = onebyomega_Kmat_smooth(scat_length,s,k,s0,m);

    comp denom = onebKmat + Fs00 + Gs00;

    return (1.0/L*L*L)*(Fs00/3.0 - Fs00*(1.0/denom)*Fs00);
}

void F3s_vs_s()
{
    //Inputs
    double pi = acos(-1.0);
    double a = -10.0;
    double kx = 0.0;
    double ky = 0.0;
    double kz = 0.0;
    comp k = sqrt(kx*kx + ky*ky + kz*kz);

    double m = 1.0;
    double epsilon = 0.001;
    double L = 6.0;
    double kpoints = 100.0;
    double n_max = 15.0;
    double s0 = 1.0;
    double sinitial = 6.25;
    double sfinal = 20.25;
    double spoints = 200.0;
    double dels = abs(sfinal - sinitial)/spoints; 

    ofstream fout;
    string filename = "F3_vs_s_for_a_"+to_string(a)+"_L_"+to_string((int)L)+"_s0_"+to_string(s0)+".dat";
    fout.open(filename.c_str());

    int scount = 0;
    for(int i=0;i<spoints+1;++i)
    {
        double s = sinitial + i*dels;

        comp k = sqrt(kx*kx + ky*ky + kz*kz);

        comp F3smooth = F3s_smooth(a, s, k, k, m, s0, n_max, L, kpoints);
        //comp F3hard   = F3_hard  (s, a, kx, ky, kz, s0, m, epsilon, L, kpoints, n_max);
        
        //comp delF3 = F3smooth - F3hard;

        cout<<"s = "<<s<<endl;
        cout<<"F3smooth = "<<F3smooth<<endl;
        //cout<<"F3hard = "<<F3hard<<endl;
        cout<<"cutoff_s0 = "<<s0<<endl;
        cout<<"========================="<<endl;
        cout<<endl;

        fout<<setprecision(16)
            <<s<<'\t'<<a<<'\t'<<s0<<'\t'
            <<real(k)<<'\t'<<imag(k)<<'\t'<<L<<'\t'
            <<real(F3smooth)<<'\t'<<imag(F3smooth)<<endl;
            //<<real(F3hard)<<'\t'<<imag(F3hard)<<'\t'
            //<<real(delF3)<<'\t'<<imag(delF3)<<endl;
            
        cout<<"run = "<<scount+1<<endl;
        scount = scount + 1;
    }
    fout.close();
}

void number_test()
{
    double a = -2.0;
    double m = 1.0;

    double s = 9.01;
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb, m);

    
}

void raul_test()
{
    comp s = 7.0;
    double m = 1.0;

    cout<<pmom(s,0.0,m)<<endl;
}

void Fs_KSS(    comp s,
                comp k,
                double m    )
{
    
}

int main()
{
    //cutoff_test();
    //F3_vs_s();
    F3s_vs_s();
    return 0;
}