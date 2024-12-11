/* Here we write down all the functions needed to calculate d(q,q)
which is needed to test for Mphib amplitude. This will serve as a test 
for the integral equation which will then be used to solve for KKpi amplitude */

#ifndef MPHIB_TESTING_H
#define MPHIB_TESTING_H

#include<bits/stdc++.h>
#include "momentum_vector.h"
#include "integral_equations_2plus1_system.h"
#include "linear_solvers_cpu.h"


//This is to calculate the GS(p,qb) matrix where qb is the 
//spectator momentum for which the two-body system will have 
//a bound-state. Depending on which qb we pass, only one of 
//the two-body system will form a bound-state. 

void negative_GSpqb_mat(    Eigen::MatrixXcd &negGpqb_mat,
                            comp En,
                            std::vector<comp> &pvec_for_m1m2, 
                            std::vector<comp> &kvec_for_m1m1,
                            std::vector<comp> &weights_for_pvec_for_m1m2, 
                            std::vector<comp> &weights_for_kvec_for_m1m1, 
                            comp qb, //passed spectator momentum for two-body bound-state
                            double m1, 
                            double m2, 
                            double eps_for_ope,
                            double eps_for_cutoff,
                            comp total_P )
{
    char debug = 'n';
    double mi,mj,mk; 
    comp ii = {0.0,1.0};
    comp pi = std::acos(-1.0); 

    int size1 = pvec_for_m1m2.size();
    int size2 = kvec_for_m1m1.size(); 

    Eigen::MatrixXcd Filler0_22(size2, 1);
    Filler0_22 = Eigen::MatrixXcd::Zero(size2, 1);

    //for (i,j) = 1 1 and k = 2
    
    mi = m1; 
    mj = m1;
    mk = m2; 

    Eigen::MatrixXcd G_11(size1, 1);

    for(int i=0; i<size1; ++i)
    {
        comp p = pvec_for_m1m2[i];
        comp G = GS_pk(En, p, qb, mi, mj, mk, eps_for_ope, eps_for_cutoff);
            
        G_11(i,0) = G;
    } 

    //for (i,j) = 1 2 and k = 1
    
    mi = m1; 
    mj = m2;
    mk = m1; 

    Eigen::MatrixXcd G_12(size1, 1);

    for(int i=0; i<size1; ++i)
    {
        comp p = pvec_for_m1m2[i];
        comp G = GS_pk(En, p, qb, mi, mj, mk, eps_for_ope, eps_for_cutoff);
            
        G_12(i,0) = G;
    } 
    //for (i,j) = 2 1 and k = 1
    
    mi = m2; 
    mj = m1;
    mk = m1; 

    Eigen::MatrixXcd G_21(size2, 1);

    for(int i=0; i<size2; ++i)
    {
        comp p = kvec_for_m1m1[i];
        comp G = GS_pk(En, p, qb, mi, mj, mk, eps_for_ope, eps_for_cutoff);
            
        G_21(i,0) = G;
    } 

    G_12 = std::sqrt(2.0)*G_12; 
    G_21 = std::sqrt(2.0)*G_21; 

    Eigen::MatrixXcd Gmat(size1 + size2, 2); 

    Gmat << G_11, G_12, 
            G_21, Filler0_22;

    negGpqb_mat = -Gmat; 

    //std::cout<<"here"<<std::endl; 


    if(debug=='y')
    {
        std::cout<<"======== negative Gmat(p,qb) ========="<<std::endl; 
        std::cout<<negGpqb_mat<<std::endl; 
        std::cout<<"======================================"<<std::endl; 
    }

}

void vec_printer(std::vector<comp> vec)
{
    for(int i=0; i<vec.size(); ++i)
    std::cout<<"vec["<<i<<"] = "<<vec[i]<<std::endl; 
}



void dpqb_solver_ERE(   Eigen::MatrixXcd &dmatpqb,
                        comp En,
                        double m1, 
                        double m2, 
                        std::vector<comp> &pvec_m1m2,
                        std::vector<comp> &weights_pvec,
                        std::vector<comp> &kvec_m1m1, 
                        std::vector<comp> &weights_kvec, 
                        comp &qb_val, //relative spectator momentum for two-body bound-state
                        double eps_for_m2k, 
                        double eps_for_ope, 
                        double eps_for_cutoff, 
                        comp total_P, 
                        double a0_m1,
                        double r0_m1,
                        double a0_m2,
                        double r0_m2, 
                        int number_of_points    )
{
    char debug = 'n';
    double mi = m1;
    double mj = m1; 
    double mk = m2; 


    comp kmax_for_m1 = pmom(En, 0.0, m1);
    comp kmax_for_m2 = pmom(En, 0.0, m2); 

    comp epsilon_for_kvec = 1.0e-5; 

    if(debug=='y')
    {
        std::cout<<"kmax_for_m1 = "<<kmax_for_m1<<std::endl; 
        std::cout<<"kmax_for_m2 = "<<kmax_for_m2<<std::endl; 

    }

    std::vector<comp> pvec_for_m1m2;
    std::vector<comp> weights_for_pvec_for_m1m2; 
    std::vector<comp> kvec_for_m1m1; 
    std::vector<comp> weights_for_kvec_for_m1m1; 

    flavor_based_momentum_vector(pvec_for_m1m2, weights_for_pvec_for_m1m2, En, m1, number_of_points);


    flavor_based_momentum_vector(kvec_for_m1m1, weights_for_kvec_for_m1m1, En, m2, number_of_points);

    pvec_m1m2 = pvec_for_m1m2; 
    weights_pvec = weights_for_pvec_for_m1m2; 
    kvec_m1m1 = kvec_for_m1m1; 
    weights_kvec = weights_for_kvec_for_m1m1; 

    if(debug=='y')
    {
        std::cout<<"pvec for m1m2 = "<<std::endl; 
        vec_printer(pvec_for_m1m2); 
        std::cout<<"kvec for m1m1 = "<<std::endl; 
        vec_printer(kvec_for_m1m1); 
        std::cout<<"=============================="<<std::endl; 
    }


    Eigen::MatrixXcd B_mat; 
    Eigen::MatrixXcd negG_mat; 
    Eigen::MatrixXcd d_mat; 

    comp sigb = sigma_b_plus(a0_m1, mj, mk);
    comp qb = qb_i(En, sigb, mi);
    qb_val = qb; 
    double relerr; 

    Bmat_2plus1_system_ERE( B_mat, En, pvec_for_m1m2, kvec_for_m1m1, weights_for_pvec_for_m1m2, weights_for_kvec_for_m1m1, m1, m2, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_m1, r0_m1, a0_m2, r0_m2 );
    //std::cout<<"here"<<std::endl; 
    
    negative_GSpqb_mat(negG_mat, En, pvec_for_m1m2, kvec_for_m1m1, weights_for_pvec_for_m1m2, weights_for_kvec_for_m1m1, qb, m1, m2, eps_for_ope, eps_for_cutoff, total_P);
    //std::cout<<"here"<<std::endl; 
    
    LinearSolver_2(B_mat, d_mat, negG_mat, relerr); 

    if(debug=='y')
    {
        std::cout<<"Bmat = "<<std::endl; 
        std::cout<<B_mat<<std::endl; 
        std::cout<<"==========================="<<std::endl; 
        std::cout<<"negG_mat = "<<std::endl; 
        std::cout<<negG_mat<<std::endl; 
        std::cout<<"==========================="<<std::endl;
        std::cout<<"d_mat = "<<std::endl; 
        std::cout<<d_mat<<std::endl; 
        std::cout<<"==========================="<<std::endl; 
        std::cout<<"relative error of inversion = "<<relerr<<std::endl; 
    
    }

    dmatpqb = d_mat; 

}

void dqq_interpolator(  Eigen::MatrixXcd &result_dqbqb_mat,
                        Eigen::MatrixXcd &dpqb_mat, 
                        comp En,
                        double m1, 
                        double m2, 
                        std::vector<comp> &pvec_for_m1m2, 
                        std::vector<comp> &weights_for_pvec_for_m1m2, 
                        std::vector<comp> &kvec_for_m1m1,
                        std::vector<comp> &weights_for_kvec_for_m1m1, 
                        comp &qb, //relative spectator momentum for the two-body bound-state
                        double eps_for_m2k, 
                        double eps_for_ope, 
                        double eps_for_cutoff, 
                        comp total_P, 
                        double a0_m1,
                        double r0_m1,
                        double a0_m2,
                        double r0_m2, 
                        int number_of_points    )
{
    char debug = 'n';
    comp ii = {0.0,1.0};
    comp pi = std::acos(-1.0); 

    double mi = m1; 
    double mj = m1; 
    double mk = m2; 

    int size1 = pvec_for_m1m2.size();
    int size2 = kvec_for_m1m1.size(); 

    Eigen::MatrixXcd Filler0_22(1, 1);
    Filler0_22 = Eigen::MatrixXcd::Zero(1, 1);

    //for (i,j) = 1 1 and k = 2
    
    mi = m1; 
    mj = m1;
    mk = m2; 

    Eigen::MatrixXcd G_11(1, 1);
    
    comp G = GS_pk(En, qb, qb, mi, mj, mk, eps_for_ope, eps_for_cutoff);
            
    G_11(0,0) = G;

    //for (i,j) = 1 2 and k = 1
    
    mi = m1; 
    mj = m2;
    mk = m1; 

    Eigen::MatrixXcd G_12(1, 1);
    
    G = GS_pk(En, qb, qb, mi, mj, mk, eps_for_ope, eps_for_cutoff);
            
    G_12(0,0) = G;

    //for (i,j) = 2 1 and k = 1
    
    mi = m2; 
    mj = m1;
    mk = m1; 

    Eigen::MatrixXcd G_21(1, 1);
    
    G = GS_pk(En, qb, qb, mi, mj, mk, eps_for_ope, eps_for_cutoff);
            
    G_21(0,0) = G;

    G_12 = std::sqrt(2.0)*G_12; 
    G_21 = std::sqrt(2.0)*G_21; 

    Eigen::MatrixXcd Gmat(2, 2); 

    Gmat << G_11, G_12, 
            G_21, Filler0_22;

    if(debug=='y')
    {
        std::cout<<"======== Gqq ==========="<<std::endl; 
        std::cout<<Gmat<<std::endl;
        std::cout<<"========================"<<std::endl; 
    }

    Eigen::MatrixXcd negGqbqb_mat;

    negGqbqb_mat = -Gmat;

    //First we build the M2 matrix 
    //int size1 = pvec_for_m1m2.size(); 
    //int size2 = kvec_for_m1m1.size(); 

    Eigen::MatrixXcd M2k_m1m2(size1, size1);
    Eigen::MatrixXcd M2k_m1m1(size2, size2);

    M2k_m1m2 = Eigen::MatrixXcd::Zero(size1,size1);
    M2k_m1m1 = Eigen::MatrixXcd::Zero(size2,size2); 

    Eigen::MatrixXcd M2k_filler0_12(size1, size2); 
    Eigen::MatrixXcd M2k_filler0_21(size2, size1); 

    M2k_filler0_12 = Eigen::MatrixXcd::Zero(size1, size2); 
    M2k_filler0_21 = Eigen::MatrixXcd::Zero(size2, size1); 

    //When m1 is in the spectator
    mi = m1; 
    mj = m1; 
    mk = m2; 

    for(int i=0; i<size1; ++i)
    {
        comp p = pvec_for_m1m2[i];
        comp M2k = M2k_ERE(En, p, total_P, a0_m1, r0_m1, mi, mj, mk, eps_for_m2k);
        M2k_m1m2(i,i) = M2k; 
    }

    //When m2 is in the spectator
    mi = m2; 
    mj = m1; 
    mk = m1; 

    for(int i=0; i<size2; ++i)
    {
        comp k = kvec_for_m1m1[i]; 
        comp M2k = M2k_ERE(En, k, total_P, a0_m2, r0_m2, mi, mj, mk, eps_for_m2k);
        M2k_m1m1(i,i) = M2k; 
    }

    M2k_m1m1 = 0.5*M2k_m1m1; //Multiplied with 0.5 because of identical particle symmetry

    Eigen::MatrixXcd M2k_mat(size1 + size2, size1 + size2); 

    M2k_mat <<  M2k_m1m2, M2k_filler0_12,
                M2k_filler0_21, M2k_m1m1; 
    
    if(debug=='y')
    {
        std::cout<<"======== M(k)mat ==========="<<std::endl; 
        std::cout<<M2k_mat<<std::endl;
        std::cout<<"========================"<<std::endl; 
    }
    /*--------------------------------------------------*/

    //Now we build the G mat inside the integral with the weights 
    
    //for (i,j) = 1 1 and k = 2
    mi = m1;
    mj = m1; 
    mk = m2; 
    Eigen::MatrixXcd WG_11(1, size1);

    for(int i=0; i<size1; ++i)
    {
        comp k = pvec_for_m1m2[i]; 
        comp weights = weights_for_pvec_for_m1m2[i]; 
        comp omgk = omega_func(k, mj); 

        comp G = GS_pk(En, qb, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
        comp W = weights*k*k/(pow(2.0*pi,2.0)*omgk);
        WG_11(0,i) = W*G; 
    }

    //for (i,j) = 1 2 and k = 1
    mi = m1;
    mj = m2; 
    mk = m1; 
    Eigen::MatrixXcd WG_12(1, size2);

    for(int i=0; i<size2; ++i)
    {
        comp k = kvec_for_m1m1[i]; 
        comp weights = weights_for_kvec_for_m1m1[i]; 
        comp omgk = omega_func(k, mj); 

        comp G = GS_pk(En, qb, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
        comp W = weights*k*k/(pow(2.0*pi,2.0)*omgk);
        WG_12(0,i) = W*G; 
    }

    //for (i,j) = 2 1 and k = 1
    mi = m2;
    mj = m1; 
    mk = m1; 
    Eigen::MatrixXcd WG_21(1, size1);

    for(int i=0; i<size1; ++i)
    {
        comp k = pvec_for_m1m2[i]; 
        comp weights = weights_for_pvec_for_m1m2[i]; 
        comp omgk = omega_func(k, mj); 

        comp G = GS_pk(En, qb, k, mi, mj, mk, eps_for_ope, eps_for_cutoff);
        comp W = weights*k*k/(pow(2.0*pi,2.0)*omgk);
        WG_21(0,i) = W*G; 
    }

    WG_12 = std::sqrt(2.0)*WG_12; 
    WG_21 = std::sqrt(2.0)*WG_21; 

    Eigen::MatrixXcd WG_mat(2, size1 + size2); 
    Eigen::MatrixXcd WG_Filler0_22(1, size2);
    WG_Filler0_22 = Eigen::MatrixXcd::Zero(1, size2);
    
    WG_mat <<   WG_11, WG_12,
                WG_21, WG_Filler0_22; 
    
    if(debug=='y')
    {
        std::cout<<"======== WGqk ==========="<<std::endl; 
        std::cout<<WG_mat<<std::endl;
        std::cout<<"========================="<<std::endl; 
    }


    if(debug=='y')
    {
        Eigen::MatrixXcd first_op = WG_mat*M2k_mat;

        std::cout<<"======== First Op WGM2 ==========="<<std::endl; 
        std::cout<<first_op<<std::endl;
        std::cout<<"=================================="<<std::endl; 
    
        Eigen::MatrixXcd second_op = first_op*dpqb_mat; 

        std::cout<<"======== Second Op WGM2dpqb ==========="<<std::endl; 
        std::cout<<second_op<<std::endl;
        std::cout<<"======================================="<<std::endl; 
    }


    result_dqbqb_mat = negGqbqb_mat - WG_mat*M2k_mat*dpqb_mat; 

    if(debug=='y')
    {
        std::cout<<"======== result dqbqb mat ========="<<std::endl; 
        std::cout<<result_dqbqb_mat<<std::endl; 
        std::cout<<"==================================="<<std::endl; 
    }

}

/*  SA Method */
comp delta_epsilon( comp sigk, 
                    comp sigb, 
                    double eps_for_m2k  )
{
    double pi = std::acos(-1.0); 

    comp res = eps_for_m2k/(pi*(std::pow(sigk - sigb,2.0) + std::pow(eps_for_m2k,2.0)));

    return res; 
}

//comp dela_M2k_ERE(  )



#endif