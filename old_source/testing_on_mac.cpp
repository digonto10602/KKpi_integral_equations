//macos:g++-12 ../testing_on_mac.cpp -o printer -O3 -std=c++14 -I /usr/local/Cellar/eigen/3.4.0_1/include/eigen3
//wahab:g++ printer.cpp -o printer -O3 -std=c++14 -I /cm/shared/applications/eigen/3.3.7/include/eigen3/
//wahab_gpu:nvcc -I /cm/shared/applications/cuda-toolkit/10.0.130/include -I /cm/shared/applications/eigen/3.3.7/include/eigen3/ printer_gpu_momrep.cpp  -O3 -std=c++14 -L /cm/shared/applications/cuda-toolkit/10.0.130/lib64 -lcublas -lcusolver -lcudart -o printer_gpu
//wahab_gpu_omp:nvcc -I /cm/shared/applications/cuda-toolkit/10.0.130/include -I /cm/shared/applications/eigen/3.3.7/include/eigen3/ printer_gpu_momrep.cpp  -O3 -std=c++14 -L /cm/shared/applications/cuda-toolkit/10.0.130/lib64 -lcublas -lcusolver -lcudart -lgomp -Xcompiler -fopenmp -o printer_gpu
//ubuntupc_omp:nvcc117 -I/usr/include/eigen3/ -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/include -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/lib64 printer_gpu_momrep.cpp  -O3 -std=c++14 -lcublas -lcusolver -lcudart -lgomp -Xcompiler -fopenmp -o printer_gpu
//ubuntupc_omp1:nvcc -I/usr/include/eigen3/ -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/include -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/lib64 printer_gpu_momrep.cpp  -O3 -std=c++14 -lcublas -lcusolver -lcudart -lgomp -Xcompiler -fopenmp -o printer_gpu
//ubuntu cpu: g++ ../testing_on_mac.cpp -o printer -O3 -std=c++14 -I /usr/include/eigen3/


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

void discontinuity_Gs_fixed_q()
{
    comp ii = {0.0,1.0};
    double a = -1600000.0;
    double m = 1.0;

    double epsilon = 1.0e-10;

    double sinitial = 8.5;
    double sfinal = 9.9;
    double spoints = 500.0;
    double dels = abs(sinitial - sfinal)/spoints;

    

    for(int i=0;i<spoints;++i)
    {
        comp s = sinitial + i*dels;
        comp sigb = 2.0;//sigmab(a,m);
        comp q = pmom(s,sigb,m);
        comp Gsplus = GS_pk(s+ii*epsilon,q,q,m,0.0);
        comp Gsminus = GS_pk(s-ii*epsilon,q,q,m,0.0);

        cout<<setprecision(20);
        cout<<real(s)<<'\t'
            <<real(Gsplus)<<'\t'
            <<imag(Gsplus)<<'\t'
            <<real(Gsminus)<<'\t'
            <<imag(Gsminus)<<'\t'
            <<real(Gsplus-Gsminus)<<'\t'
            <<imag(Gsplus-Gsminus)<<'\t'
            <<abs(Gsplus-Gsminus)<<endl;

    }

}

void discontinuity_Gs_fixed_s()
{
    double a = 16.0;
    double m = 1.0;

    double epsilon = 1.0e-10;

    double sinitial = 8.5;
    double sfinal = 9.9;
    double spoints = 500.0;
    double dels = abs(sinitial - sfinal)/spoints;


    double sreal = 9.01;
    comp splus = sreal + ii*epsilon;
    comp sminus = sreal - ii*epsilon; 
    comp qinitial = 0.0;
    comp qfinalplus = pmom(splus,0.0,m);
    comp qfinalminus = pmom(sminus,0.0,m);

    vector<comp> qvecplus;
    vector<comp> qvecminus; 
    vector<comp> weights; 
    double qvecpoints = 500;
    line_maker_with_weights(qvecplus,weights,qinitial,qfinalplus,qvecpoints);
    line_maker_with_weights(qvecminus,weights,qinitial,qfinalminus,qvecpoints);


    for(int i=0;i<qvecpoints;++i)
    {
        comp sigb = sigmab(a,m);
        comp pplus = pmom(splus,sigb,m);
        comp pminus = pmom(sminus,sigb,m);
        comp qplus = qvecplus[i];
        comp qminus = qvecminus[i];
        comp Gsplus = GS_pk(splus,pplus,qplus,m,0.0);
        comp Gsminus = GS_pk(sminus,pminus,qminus,m,0.0);

        cout<<setprecision(20);
        cout<<real(qplus)<<'\t'
            <<imag(qplus)<<'\t'
            <<real(qminus)<<'\t'
            <<imag(qminus)<<'\t'
            <<real(Gsplus)<<'\t'
            <<imag(Gsplus)<<'\t'
            <<real(Gsminus)<<'\t'
            <<imag(Gsminus)<<'\t'
            <<real(Gsplus-Gsminus)<<'\t'
            <<imag(Gsplus-Gsminus)<<'\t'
            <<abs(Gsplus-Gsminus)<<endl;

    }

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
                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
                vector<comp> temp_qvec;
                vector<comp> temp_weights;
                /*contour_for_resonance_4(temp_qvec,temp_weights,real(s) - ii*imag(s),m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k);

                for(int qi=0;qi<temp_qvec.size();++qi)
                {
                    comp temp_q1 = real(temp_qvec[qi]) - ii*abs(imag(temp_qvec[qi]));
                    qvec.push_back(temp_q1);
                    weights.push_back(temp_weights[qi]);
                }*/
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
            //Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);


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
            //interpolator_ds_integraleq_momrep_2eps_with_weights(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,result);
            

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

void ds_belowthreshold_vs_s3_N_3d_with_weights_for_resonances_thirdsheet_bottom_half()
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
        string filename="dS_sheet3_bottomhalf_a_" + to_string((int)a) + "_N_" + to_string((int)N) + "_3d.dat";
        
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
            int tag_for_m2k_1 = 0;
            int tag_for_m2k_2 = 0;

            if(s3imag<0.0)
            {
                contour_for_resonance_7(qvec,weights,a,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k_1,tag_for_m2k_2);

            }
            else if(s3imag>=0.0)
            {
                contour_for_resonance_6(qvec,weights,a,s,m,kmin,kmax,eps_for_m2k,points1,tag_for_m2k);

                //line_maker_with_weights(qvec,weights,kmin,kmax,points1);
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
                Bmat_maker_momrep_2eps_with_weights_for_resonances_sheet3(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k_1,tag_for_m2k_2);
                

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
                interpolator_ds_integraleq_momrep_2eps_with_weights_for_resonances_sheet3(dsol,qvec,weights,s,qval,qval,a,m,eps,eps_for_m2k,theta,tag_for_m2k_1,tag_for_m2k_2,result);
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

void testing_Bmat_for_sheet3()
{
    double m = 1.0;
    double a = -8.1;

    double sreal = 9.005;
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
    //mom_vector_maker_seba_imsneg(qvec,weights,s,0.0,kmax,a,m,eps,eps_for_m2k,N);
                
    //contour for third sheet of dS     
    contour_for_resonance_7(qvec,weights,a,s,m,kmin,kmax,eps_for_m2k,N,tag_for_m2k_1,tag_for_m2k_2);

    cout<<"tag_for_m2k 1 = "<<tag_for_m2k_1<<endl;
    cout<<"tag_for_m2k 2 = "<<tag_for_m2k_2<<endl;
    
    int size = qvec.size();

    Eigen::MatrixXcd Bmat(size,size);
    double theta = 0;
    //Bmat_maker_momrep_2eps_with_weights(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k);

    Bmat_maker_momrep_2eps_with_weights_for_resonances_sheet3(Bmat,s,qvec,qvec,weights,a,m,eps,eps_for_m2k,theta,tag_for_m2k_1,tag_for_m2k_2);

    string Bmat_file = "test_Bmat.dat";
    ofstream fout;
    fout.open(Bmat_file.c_str());
    for(int i=0;i<size;++i)
    {
        for(int j=0;j<size;++j)
        {
            comp q1 = qvec[i];
            comp q2 = qvec[j];
            fout<<i<<'\t'
                <<j<<'\t'
                <<real(q1)<<'\t'
                <<imag(q1)<<'\t'
                <<real(q2)<<'\t'
                <<imag(q2)<<'\t'
                <<real(Bmat(i,j))<<'\t'
                <<imag(Bmat(i,j))<<'\t'
                <<size<<endl;
        }
        cout<<"i = "<<i<<" done"<<endl;
    }
    fout.close();

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
    {
        //if(a>=3.7) dela = 0.101;
        //else if(a>=17.0) dela = 0.51;
        //else if(a>=1000.0) dela = 300.01;

        double gvalleft = real(Gleftbranchpoint(a,m));
        double gvalright = real(Grightbranchpoint(a,m));
    
        double sinitial = 8.5;//6.5;//8.3;//5.0;//8.974;//8.6;//8.975;//8.96815;//7.5;//8.7;//7.5;//8.75;//1.001;//gvalright - 0.5;//(gvalleft + (gvalleft+gvalright)/2.0)/2.0;//8.965;//8.70115128532;//8.72;//real(phibthreshold(a,m));//
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

            //switch_for_gvec_fixer=0;
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



int main(int argc, char *argv[])
{
    //resonance_vs_real_s_sheet2();

    //discontinuity_Gs_fixed_q();

    //ds_belowthreshold_vs_s3_N_3d_with_weights_for_resonances_thirdsheet_bottom_half();
    //testing_Bmat_for_sheet3();
    //ds_belowthreshold_vs_s3_N_3d_with_weights_for_resonances_secondsheet_full();
    //discontinuity_Gs_fixed_s();


    //testing previous results
    Mphib_secondsheet_belowthreshold_vs_s3_N_fixed_s3imag_with_weights();
    return 0;
}