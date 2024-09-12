#ifndef INTEGRALEQN_MOM_BASED_H
#define INTEGRALEQN_MOM_BASED_H

#include "functions_mom_based.h"


/* Here the convention is the 2+1 system is defined
with 2 particles of mass m1 and one particle of mass
m2, */
/* For the analysis of KKpi system, mK is m1 
and mpi is m2 */
void kernel_2plus1_system_ERE(  comp E, 
                                std::vector<comp> &pvec_for_m1m2,
                                std::vector<comp> &kvec_for_m1m1, 
                                double m1, 
                                double m2,
                                double eps_for_m2k,
                                double eps_for_ope,
                                double eps_for_cutoff,
                                comp total_P,
                                double a,
                                double r    )
{
    //Boiler Plate
    int size_for_M2_1 = pvec_for_m1m2.size(); 
    int size_for_M2_2 = kvec_for_m1m1.size(); 

    //when m1 and m2 are in twobody subchannel
    Eigen::MatrixXcd M2k_for_12(size_for_M2_1, size_for_M2_1);
    //when m1 and m1 are in twobody subchannel 
    Eigen::MatrixXcd M2k_for_11(size_for_M2_2, size_for_M2_2);

    Eigen::MatrixXcd Filler0_12(size_for_M2_1, size_for_M2_2);
    Eigen::MatrixXcd Filler0_21(size_for_M2_2, size_for_M2_1);
    Eigen::MatrixXcd Filler0_22(size_for_M2_2, size_for_M2_2);

    Filler0_12 = Eigen::MatrixXcd::Zero(size1,size2);
    Filler0_21 = Eigen::MatrixXcd::Zero(size2,size1);
    Filler0_22 = Eigen::MatrixXcd::Zero(size2,size2);


    //First we setup the matrix for M2
    //if m1 is spectator then the 2body sub-channel
    //have m1 and m2 
    double mi = m1;
    double mj = m1;
    double mk = m2; 
    
    

    for(int i=0; i<size_for_M2_1; ++i)
    {
        comp mom_p = pvec_for_m1m2[i];

        comp M2k = M2k_ERE(E, mom_p, total_P, a, r, mi, mj, mk, eps_for_m2k);

        M2k_for_12(i,i) = M2k; 
    }

    //Now m2 is the spectator, so the two body subchannel 
    //have m1 and m1 
    mi = m2;
    mj = m1; 
    mk = m1; 
    for(int i=0; i<size_for_M2_2; ++i)
    {
        comp mom_k = kvec_for_m1m1[i];

        comp M2k = M2k_ERE(E, mom_k, total_P, a, r, mi, mj, mk, eps_for_m2k);
    }



}

#endif 