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

//here in these packages above, we specifically set H = 1 everywhere first

comp s_based_cutoff(    comp s,
                        double threshold  )
{
    if(real(s)<threshold)
    return 0.0;
    else return 1.0;
}

comp s2k_based_cutoff(  comp s2k,
                        double threshold    )
{
    if(real(s2k)>threshold) return 1.0;
    else return 0.0;
}

void Gvec_newcutoff(    Eigen::VectorXcd &Gvec,
                        comp s,
                        double threshold,
                        vector<comp> qvec,
                        comp k,
                        double m,
                        double epsilon
                        )
{
    for(int i=0;i<qvec.size();++i)
    {
        //Gvec(i) = s_based_cutoff(s,threshold)*GS_pk(s,qvec[i],p,m,epsilon);
        comp s2p = sigma_p(s,qvec[i],m);
        comp s2k = sigma_p(s,k,m);
        Gvec(i) = s2k_based_cutoff(s2p,threshold)*s2k_based_cutoff(s2k,threshold)*GS_pk(s,qvec[i],k,m,epsilon);
        
    }
}

void Gmat_newcutoff(    Eigen::MatrixXcd &Gmat,
                        comp s,
                        double threshold,
                        vector<comp> qvec,
                        vector<comp> qvec1,
                        double m,
                        double epsilon
                        )
{
    for(int i=0;i<qvec.size();++i)
    {
        for(int j=0;j<qvec1.size();++j)
        {
            Gmat(i,j) = s_based_cutoff(s,threshold)*GS_pk(s,qvec[i],qvec1[j],m,epsilon);
        }
        
    }
}

comp kernel_pk_2eps_newcutoff(    comp s,
                        comp p,
                        comp k,
                        double a,
                        double m,
                        double epsilon,
                        double eps_for_m2k,
                        double threshold )
{
    double pi = acos(-1.0);
    comp picomp = (comp) pi;
    //cout<<"Gs:"<<GS(s,sigp,sigk,m,epsilon)<<endl;
    //cout<<"tau:"<<tau(s,sigk,m)<<endl;
    //cout<<"M2:"<<M2kfunc(a,sigk,m,epsilon)<<endl;
    comp sigk = sigma_p(s,k,m);
    comp sigp = sigma_p(s,p,m);
    
    //return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*GS_pk(s,p,k,m,epsilon)*M2kfunc(a,sigk,m,epsilon);
    return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*s2k_based_cutoff(sigp,threshold)*s2k_based_cutoff(sigk,threshold)*GS_pk(s,p,k,m,epsilon)*M2kfunc(a,sigk,m,eps_for_m2k);
    
    //return (k*k/(pow(2.0*pi,2.0)*omega_comp(k,m)))*M2kfunc(a,sigk,m,eps_for_m2k);
    
    
    //return ((comp)1.0/(2.0*picomp))*GS(s,sigp,sigk,m,epsilon)*tau1(s,sigk,m)*M2kfunc(a,sigk,m,epsilon);
    //return GS(s,sigp,sigk,m,epsilon);
}



void Bmat_newcutoff(       Eigen::MatrixXcd &Bmat,
                                                comp s,
                                                vector<comp> &p,
                                                vector<comp> &k,
                                                vector<comp> &weights,
                                                double a,
                                                double m,
                                                double epsilon,
                                                double eps_for_m2k,
                                                double threshold  )
{
    comp ii = {0.0,1.0};
    comp pi = acos(-1.0);
    comp delk = {0.0,0.0};

    for(int i=0;i<p.size();++i)
    {
        for(int j=0;j<k.size();++j)
        {
            if(j==0)
            {
                delk = k[1] - k[0];
            }
            else
            {
                delk = k[j] - k[j-1];
            }

            comp weight = weights[j];
            //make sure you supply both the invariant vectors
            //as the same, otherwise this condition will fail
            if(j==i)
            {
                comp one = {1.0,0.0};
                Bmat(i,j) = one + weight*kernel_pk_2eps_newcutoff(s,p[i],k[j],a,m,epsilon,eps_for_m2k,threshold);
                //Bmat(i,j) =  weight*kernel_pk_2eps(s,p[i],k[j],a,m,epsilon,eps_for_m2k);
            
            }
            else 
            {
                Bmat(i,j) = weight*kernel_pk_2eps_newcutoff(s,p[i],k[j],a,m,epsilon,eps_for_m2k,threshold);
            }
        }
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

    double s3imag = 1.0e-5;
    

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
    int s2kcutoffcount = 0;
    double s2kinitial = 0.0;
    double s2kfinal = 3.5;
    double dels2k = 0.5;

    for(int s2kcount=0;s2kcount<10;++s2kcount)
    {
        double s2kcutoff = s2kinitial + s2kcount*dels2k;
        //if(a>=3.7) dela = 0.101;
        //else if(a>=17.0) dela = 0.51;
        //else if(a>=1000.0) dela = 300.01;

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        double sinitial = 8.5;//3.5;//6.5;//8.3;//5.0;//8.974;//8.6;//8.975;//8.96815;//7.5;//8.7;//7.5;//8.75;//1.001;//gvalright - 0.5;//(gvalleft + (gvalleft+gvalright)/2.0)/2.0;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
        double sfinal = real(phibthreshold(a,m));//9.0;//real(phibthreshold(a,m));//8.82;
        
        double dels = abs(sinitial-sfinal)/delspoints;//abs(8.70275246887-8.70115128532);//
        points1 = N;
        //string filename="Mphib_secondsheet_acount_" + to_string(acount) + ".dat";
        string filename="Mphib_a_" + to_string(a) + "_N_" + to_string((int)N) + "_s2kcutoff_" + to_string(s2kcount) + ".dat";
        
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
            if(real(s)>=real(phibthreshold(a,m)))
            {
                double eta = 15.0;
                eps_for_m2k = eps_above_threshold(eta,s,a,m,points1);
            }
            else 
            {
                eps_for_m2k = 0.0;
            }
            cout<<"am = "<<a<<endl;
            //cout<<"gLeft = "<<setprecision(16)<<gvalleft<<'\t'<<"gRight = "<<gvalright<<endl;
            comp sigmamin = s2kcutoff;//{0.0,0.0};
            comp sigmamax = sigmax(s,m);
            //if(((double)abs(real(sigmamax)))<s2kcutoff) continue;
            comp sigb = sigmab(a,m);
            comp qval = pmom(s,sigb,m);
            cout<<"q = "<<qval<<endl;

            comp kmin = 0.0;
            comp kmax = pmom(s,sigmamin,m);
        

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
                //cout<<"contour_43 chosen"<<endl;
                //mom_vector_maker_43(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43_with_weights(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1);
                line_maker_with_weights(qvec,weights,kmin,kmax,points1);
            }
            else if(s3imag>=0.0)
            {
                //mom_vector_maker_41(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_43(qvec,real(s) - ii*imag(s),0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //mom_vector_maker_44(qvec,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,qvec_r);
                //cout<<"contour_47 chosen"<<endl;
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

                if(real(s)<=real(phibthreshold(a,m)))
                mom_vector_maker_seba_imspos_5(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                else
                {
                    switch_for_gvec_fixer = 1;
                    line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                }
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

            //switch_for_gvec_fixer=0;
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
            /*
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
            comp sigk = sigma_p(s,qval,m);
            comp m2k = M2kfunc(a,sigk,m,eps_for_m2k);
            //result = m2k*result*m2k;
            comp rhopb = rhophib(s,a,m);
            comp mphib2 = result/(1.0 + 2.0*ii*rhopb*result);
            comp mphib2denom = (1.0 + 2.0*ii*rhopb*result);
            comp GSval = GS_pk(s,qval,qval,m,eps);
            //result = rhophib(s,a,m)*result;

            cout<<"s2kcount = "<<s2kcount<<endl;
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
                <<eps<<'\t'
                <<s2kcutoff<<'\t'
                <<pow(sqrt(s2kcutoff) + m,2.0)<<endl;
            count = count + 1;
            
            cout<<"-------------------------"<<endl;
        }
        acount = acount + 1;
        fout.close();

    }

    //}

}



int main()
{
    Mphib_secondsheet_belowthreshold_vs_s3_N_fixed_s3imag_with_weights();
    return 0;
}