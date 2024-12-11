#include <bits/stdc++.h>
#include "integral_equations_2plus1_system.h"
#include "mphib_testing.h"

void test_M2()
{
    double pmin = 1e-10;
    double pmax = 1.0; 
    double ppoints = 5; 

    double delp = std::abs(pmax - pmin)/ppoints; 

    std::vector<comp> pvec; 
    for(int i=0; i<ppoints; ++i)
    {
        double p = pmin + i*delp; 

        pvec.push_back(p); 
    }

    double E = 3.1; 
    double a = -2.0; 
    double r = 0.0; 
    double m1 = 1.0;
    double m2 = 0.5; 

    double eps = 0.0; 
    double total_P = 0.0; 

    std::vector<comp> weights; 
    Eigen::MatrixXcd kernmat;

    kernel_2plus1_system_ERE(kernmat, E, pvec, pvec, weights, weights, m1, m2, eps, eps, eps, total_P, a, r, a, r); 

}

void test_q2k_sigk()
{
    double En = 3.1; 
    double pmin = -1.5;
    double pmax = 1.5;
    double ppoints = 20;
    double mi = 1.0;
    double total_P = 0.0; 

    double delp = std::abs(pmax - pmin)/ppoints; 

    std::vector<comp> repvec;
    std::vector<comp> impvec; 

    for(int i=0;i<ppoints; ++i)
    {
        double p = pmin + i*delp; 

        repvec.push_back(p);
        impvec.push_back(ii*p); 
    }

    for(int i=0; i<repvec.size(); ++i)
    {
        comp sig1 = sigma(En, repvec[i], mi, total_P);
        comp sig2 = sigma(En, impvec[i], mi, total_P); 

        std::cout<<repvec[i]<<'\t'<<sig1<<'\t'<<impvec[i]<<'\t'<<sig2<<std::endl; 
    }
}

void test_sigk_q2k()
{
    double En = 3.1;
    double mi = 1.0; 

    double re_sigmin = 2.1;
    double im_sigmin = -1.0;
    double im_sigmax = 1.0;
    double sig_points = 20;
    double del_sig = std::abs(im_sigmin - im_sigmax)/sig_points; 

    for(int i=0; i<sig_points; ++i)
    {
        double im_sig = im_sigmin + i*del_sig; 
        comp sigk = re_sigmin + ii*im_sig; 

        comp p = pmom(En, sigk, mi);

        std::cout<<"sigk:"<<sigk<<'\t'<<"p="<<p<<std::endl;  
    }
}

void M2k_Plotter_KK()
{
    // atmpi = 0.06906
    // atmK = 0.09698
    // scattering_length_1_piK = 4.04;// - 0.2; //total uncertainty 0.05 stat 0.15 systematic 
    // scattering_length_2_KK = 4.07;// - 0.07; //total uncertainty 0.07 stat 
    // atinv = 5.666 GeV 
    // anisotropy xi = 3.444

    //double 
}

void test_kinematic_variables()
{
    double m1 = 1.0;
    double m2 = 1.0;//0.999;

    double mi = m1;
    double mj = m1;
    double mk = m2; 

    double a01 = 2.0; 
    double a02 = 2.0; 

    double En = 2.6;

    comp kmax1 = pmom(En, 0.0, m1);
    comp kmax2 = pmom(En, 0.0, m2);
    comp sigb1 = sigma_b_plus(a01, mj, mk); 
    comp sigb2 = sigma_b_minus(a01, mj, mk); 

    comp qb = qb_i(En, sigb1, mi);

    std::cout<<"kmax1 = "<<kmax1<<std::endl; 
    std::cout<<"kmax2 = "<<kmax2<<std::endl; 
    std::cout<<"sigb1 = "<<sigb1<<std::endl; 
    std::cout<<"sigb2 = "<<sigb2<<std::endl; 
    std::cout<<"qb = "<<qb<<std::endl; 


}

void test_dpqb_building()
{
    double m1 = 1.0;
    double m2 = 0.99999999990;//0.999; 
    double a0_1 = 2.0; 
    double a0_2 = 2.0; 
    double r0_1 = 0.0; 
    double r0_2 = 0.0; 

    //double En = 1.95;
    double total_P = 0.0; 
    double r = 0; 
    int number_of_points = 100;//200;//250;  

    double eps_for_m2k = 0.001;
    double eps_for_ope = 0.001; 
    double eps_for_cutoff = 0.0; 

    /*-----------------------------------------*/

    comp sigb1plus = sigma_b_plus(a0_1, m1, m2);
    comp sigb1minus = sigma_b_minus(a0_1, m1, m2); 
    comp sigb2plus = sigma_b_plus(a0_2, m1, m1); 
    comp sigb2minus = sigma_b_minus(a0_2, m1, m1); 
    comp phib1plus = std::sqrt(sigb1plus) + m1;
    comp phib1minus = std::sqrt(sigb1minus) + m1;
    comp phib2plus = std::sqrt(sigb2plus) + m2;
    comp phib2minus = std::sqrt(sigb2minus) + m2;

    std::cout<<"sigb1+ = "<<sigb1plus<<std::endl; 
    std::cout<<"sigb1- = "<<sigb1minus<<std::endl; 
    
    std::cout<<"sigb2+ = "<<sigb2plus<<std::endl; 
    std::cout<<"sigb2- = "<<sigb2minus<<std::endl; 

    std::cout<<"phib+ threshold 1 = "<<phib1plus<<std::endl; 
    std::cout<<"phib- threshold 1 = "<<phib1minus<<std::endl; 
    std::cout<<"phib+ threshold 2 = "<<phib2plus<<std::endl; 
    std::cout<<"phib- threshold 2 = "<<phib2minus<<std::endl; 

    comp threeparticle_threshold = (m1 + m1 + m2); 

    //comp En = (phib1 + threeparticle_threshold)/2.0;

    std::cout<<"threeparticle threshold = "<<threeparticle_threshold<<std::endl; 

    
    //abort();

    
    /* Generating a file here to check */
    std::ofstream fout;
    std::string filename =  "check_dqq_test1.dat";

    fout.open(filename.c_str()); 

    comp En_initial = phib1plus; 
    comp En_final = threeparticle_threshold; 

    double En_points = 3;//100.0; 
    comp delEn = std::abs(En_initial - En_final)/En_points; 

    //for(int i=0; i<(int)En_points; ++i)
    {
        int i = 1; 
        comp En = En_initial + ((comp)i)*delEn; 
        //comp qb; 
        Eigen::MatrixXcd dpqbmat; 
        Eigen::MatrixXcd dqqmat; 
        std::vector<comp> pvec_for_m1m2; 
        std::vector<comp> weights_for_pvec_for_m1m2; 
        std::vector<comp> kvec_for_m1m1; 
        std::vector<comp> weights_for_kvec_for_m1m1; 

        double eta = 25; 

        comp qb = qb_i(En, sigb1plus, m1);
        comp kmax = pmom(En, 0.0, m1); 

        eps_for_m2k = energy_dependent_epsilon(eta, En, qb, sigb1plus, kmax, m1, number_of_points);
        eps_for_ope = 0.01; //eps_for_m2k; 
        dpqb_solver_ERE(dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points); 
        dqq_interpolator(dqqmat, dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points);

        //std::cout<<"pvec size = "<<pvec_for_m1m2.size()<<std::endl; 
        /*for(int i=0; i<pvec_for_m1m2.size(); ++i)
        {
            std::cout<<"pvec i="<<i<<'\t'<<"val = "<<pvec_for_m1m2[i]<<std::endl;
        }*/

        comp dqq11 = dqqmat(0,0); 
        
        comp gval = gfunc(sigb1plus,a0_1);

        comp mphib = gval*gval*dqq11; 
        comp rhophib_mphib = rhophib(qb, En)*mphib; 

        std::cout<<En.real()<<'\t'<<(En*En).real()<<'\t'<<rhophib_mphib.real()<<'\t'<<rhophib_mphib.imag()<<std::endl; 
        fout<<std::real(En)<<'\t'<<std::real(En*En)<<'\t'<<rhophib_mphib.real()<<'\t'<<rhophib_mphib.imag()<<std::endl; 
        
    
    }
    fout.close(); 
}



void test_dpqb_building_1()
{
    double m1 = 1.0;
    double m2 = 1.0;//0.99999999990;//0.999; 
    double a0_1 = 16.0; 
    double a0_2 = 2.0; 
    double r0_1 = 0.0; 
    double r0_2 = 0.0; 

    //double En = 1.95;
    double total_P = 0.0; 
    double r = 0; 
    int number_of_points = 250;  

    double eps_for_m2k = 0.001;
    double eps_for_ope = 0.001; 
    double eps_for_cutoff = 0.0; 

    /*-----------------------------------------*/

    comp sigb1plus = sigma_b_plus(a0_1, m1, m2);
    comp sigb1minus = sigma_b_minus(a0_1, m1, m2); 
    comp sigb2plus = sigma_b_plus(a0_2, m1, m1); 
    comp sigb2minus = sigma_b_minus(a0_2, m1, m1); 
    comp phib1plus = std::sqrt(sigb1plus) + m1;
    comp phib1minus = std::sqrt(sigb1minus) + m1;
    comp phib2plus = std::sqrt(sigb2plus) + m2;
    comp phib2minus = std::sqrt(sigb2minus) + m2;

    std::cout<<"sigb1+ = "<<sigb1plus<<std::endl; 
    std::cout<<"sigb1- = "<<sigb1minus<<std::endl; 
    
    std::cout<<"sigb2+ = "<<sigb2plus<<std::endl; 
    std::cout<<"sigb2- = "<<sigb2minus<<std::endl; 

    std::cout<<"phib+ threshold 1 = "<<phib1plus<<std::endl; 
    std::cout<<"phib- threshold 1 = "<<phib1minus<<std::endl; 
    std::cout<<"phib+ threshold 2 = "<<phib2plus<<std::endl; 
    std::cout<<"phib- threshold 2 = "<<phib2minus<<std::endl; 

    comp threeparticle_threshold = (m1 + m1 + m2); 

    //comp En = (phib1 + threeparticle_threshold)/2.0;

    std::cout<<"threeparticle threshold = "<<threeparticle_threshold<<std::endl; 

    
    //abort();

    
    /* Generating a file here to check */
    std::ofstream fout;
    std::string filename =  "mphib_above_phibthreshold_a" + std::to_string((int)a0_1) + ".dat";

    fout.open(filename.c_str()); 

    comp En_initial = phib1plus; 
    comp En_final = threeparticle_threshold; 

    double En_points = 3000.0; 
    comp delEn = std::abs(En_initial - En_final)/En_points; 

    for(int i=0; i<(int)En_points; ++i)
    {
        comp En = En_initial + ((comp)i)*delEn; 
        //comp qb; 
        Eigen::MatrixXcd dpqbmat; 
        Eigen::MatrixXcd dqqmat; 
        std::vector<comp> pvec_for_m1m2; 
        std::vector<comp> weights_for_pvec_for_m1m2; 
        std::vector<comp> kvec_for_m1m1; 
        std::vector<comp> weights_for_kvec_for_m1m1; 

        double eta = 15; 

        comp qb = qb_i(En, sigb1plus, m1);
        comp kmax = pmom(En, 0.0, m1); 

        eps_for_m2k = energy_dependent_epsilon(eta, En, qb, sigb1plus, kmax, m1, number_of_points);
        eps_for_ope = eps_for_m2k; 
        dpqb_solver_ERE(dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points); 
        dqq_interpolator(dqqmat, dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points);

        /*for(int i=0; i<pvec_for_m1m2.size(); ++i)
        {
            std::cout<<"pvec i="<<i<<'\t'<<"val = "<<pvec_for_m1m2[i]<<std::endl;
        }*/

        comp dqq11 = dqqmat(0,0); 
        
        comp gval = gfunc(sigb1plus,a0_1);

        comp mphib = gval*gval*dqq11; 
        comp rhophib_mphib = rhophib(qb, En)*mphib; 
        comp rhophib_val = rhophib(qb,En); 

        double diff = (double)std::abs((std::imag(1.0/mphib) + std::real(rhophib_val))/(std::real(rhophib_val)))*100;
        
        std::cout<<std::setprecision(20)<<En.real()<<'\t'<<(En*En).real()<<'\t'<<rhophib_mphib.real()<<'\t'<<rhophib_mphib.imag()<<'\t'<<diff<<std::endl; 
        fout<<std::setprecision(20)<<std::real(En)<<'\t'<<std::real(En*En)<<'\t'<<rhophib_mphib.real()<<'\t'<<rhophib_mphib.imag()<<'\t'<<diff<<std::endl; 
        
    
    }
    fout.close(); 
}


//We change m2 value by 0.1 starting from 0.9
//and check rhoMphib value 
void test_d11_building_1()
{
    comp ii = {0.0,1.0};
    double m1 = 1.0;
    double m2 = 1.0;//0.99999999990;//0.999; 
    double a0_1 = 2.0; 
    double a0_2 = 2.0; 
    double r0_1 = 0.0; 
    double r0_2 = 0.0; 

    double m2_initial = 0.7; 
    double m2_final = 1.0; 
    double del_m2 = 0.1; 
    double m2_points = std::abs(m2_initial - m2_final)/del_m2; 
    //double En = 1.95;
    double total_P = 0.0; 
    double r = 0; 
    int number_of_points = 250;  

    double eps_for_m2k = 0.001;
    double eps_for_ope = 0.001; 
    double eps_for_cutoff = 0.0; 

    /*-----------------------------------------*/

    
    
    //abort();

    
    /* Generating a file here to check */
    std::ofstream fout;
    std::ofstream fout1; 
    

    /* file index naming*/
    fout<<"En"<<'\t'
        <<"s"<<'\t'
        <<"Re(rhoMphib)"<<'\t'
        <<"Im(rhoMphib)"<<'\t'
        <<"perc_diff"
        <<std::endl; 

    fout1<<"En"<<'\t'
         <<"s"<<'\t'
         <<"Redqq11"<<'\t'
         <<"Imdqq11"<<'\t'
         <<"Redqq12"<<'\t'
         <<"Imdqq12"<<'\t'
         <<"Redqq21"<<'\t'
         <<"Imdqq21"<<'\t'
         <<"Redqq22"<<'\t'
         <<"Imdqq22"
         <<std::endl; 

    for(int j=0; j<(int)m2_points; ++j)
    {
        m2 = m2_initial + j*del_m2; 

        comp sigb1plus = sigma_b_plus(a0_1, m1, m2);
        comp sigb1minus = sigma_b_minus(a0_1, m1, m2); 
        comp sigb2plus = sigma_b_plus(a0_2, m1, m1); 
        comp sigb2minus = sigma_b_minus(a0_2, m1, m1); 
        comp phib1plus = std::sqrt(sigb1plus) + m1;
        comp phib1minus = std::sqrt(sigb1minus) + m1;
        comp phib2plus = std::sqrt(sigb2plus) + m2;
        comp phib2minus = std::sqrt(sigb2minus) + m2;

        std::cout<<"sigb1+ = "<<sigb1plus<<std::endl; 
        std::cout<<"sigb1- = "<<sigb1minus<<std::endl; 
    
        std::cout<<"sigb2+ = "<<sigb2plus<<std::endl; 
        std::cout<<"sigb2- = "<<sigb2minus<<std::endl; 

        std::cout<<"phib+ threshold 1 = "<<phib1plus<<std::endl; 
        std::cout<<"phib- threshold 1 = "<<phib1minus<<std::endl; 
        std::cout<<"phib+ threshold 2 = "<<phib2plus<<std::endl; 
        std::cout<<"phib- threshold 2 = "<<phib2minus<<std::endl; 

        comp threeparticle_threshold = (m1 + m1 + m2); 

        //comp En = (phib1 + threeparticle_threshold)/2.0;

        std::cout<<"threeparticle threshold = "<<threeparticle_threshold<<std::endl; 

        std::string filename =    "mphib_above_phibthreshold_m2_" + std::to_string(m2) 
                                + "_a1_" + std::to_string(a0_1)
                                + ".dat";
        std::string filename1 =   "dqq_m2_" + std::to_string(m2)
                                + "_a1_" + std::to_string(a0_1) 
                                + ".dat"; 

        fout.open(filename.c_str()); 
        fout1.open(filename1.c_str()); 

        comp En_initial = phib1plus; 
        comp En_final = threeparticle_threshold; 

        double En_points = 2000.0; 
        comp delEn = std::abs(En_initial - En_final)/En_points; 
        //abort();
    
    for(int i=0; i<(int)En_points; ++i)
    {
        double img_En = 0.0; 
        comp En = std::real(En_initial + ((comp)i)*delEn) + ii*img_En; 
        //comp qb; 
        Eigen::MatrixXcd dpqbmat; 
        Eigen::MatrixXcd dqqmat; 
        std::vector<comp> pvec_for_m1m2; 
        std::vector<comp> weights_for_pvec_for_m1m2; 
        std::vector<comp> kvec_for_m1m1; 
        std::vector<comp> weights_for_kvec_for_m1m1; 

        double eta = 55; 

        comp qb = qb_i(En, sigb1plus, m1);
        comp kmax = pmom(En, 0.0, m1); 

        eps_for_m2k = energy_dependent_epsilon(eta, En, qb, sigb1plus, kmax, m1, number_of_points);
        eps_for_ope = eps_for_m2k; 
        dpqb_solver_ERE(dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points); 
        dqq_interpolator(dqqmat, dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points);

        /*for(int i=0; i<pvec_for_m1m2.size(); ++i)
        {
            std::cout<<"pvec i="<<i<<'\t'<<"val = "<<pvec_for_m1m2[i]<<std::endl;
        }*/

        comp dqq11 = dqqmat(0,0); 
        comp dqq12 = dqqmat(0,1); 
        comp dqq21 = dqqmat(1,0); 
        comp dqq22 = dqqmat(1,1); 

        
        comp gval = gfunc(sigb1plus,a0_1);

        comp mphib = gval*gval*dqq11; 
        comp rhophib_mphib = rhophib(qb, En)*mphib; 
        comp rhophib_val = rhophib(qb,En); 

        double diff = (double)std::abs((std::imag(1.0/mphib) + std::real(rhophib_val))/(std::real(rhophib_val)))*100;
        
        std::cout<<std::setprecision(6)
                 <<"m2 = "<<m2<<'\t'
                 <<"sigb = "<<sigb1plus<<'\t'
                 <<"En = "<<En<<'\t'
                 <<"s = "<<(En*En)<<'\t'
                 <<"rhoM = "<<rhophib_mphib<<'\t'
                 <<"diff = "<<diff<<'\t'
                 <<"d12 = "<<dqq12<<'\t'
                 <<"d21 = "<<dqq21<<'\t'
                 <<"d22 = "<<dqq22
                 <<std::endl; 
        
        
        fout<<std::setprecision(20)
            <<std::real(En)<<'\t'
            <<std::real(En*En)<<'\t'
            <<rhophib_mphib.real()<<'\t'
            <<rhophib_mphib.imag()<<'\t'
            <<diff
            <<std::endl; 
        
        fout1<<std::setprecision(20)
             <<std::real(En)<<'\t'
             <<std::real(En*En)<<'\t'
             <<std::real(dqq11)<<'\t'
             <<std::imag(dqq11)<<'\t'
             <<std::real(dqq12)<<'\t'
             <<std::imag(dqq12)<<'\t'
             <<std::real(dqq21)<<'\t'
             <<std::imag(dqq21)<<'\t'
             <<std::real(dqq22)<<'\t'
             <<std::imag(dqq22)
             <<std::endl; 
    
    }
    fout.close(); 
    fout1.close(); 
    }
}


int main()
{
    //test_M2(); 
    //test_q2k_sigk();
    //test_sigk_q2k();
    //test_kinematic_variables();
    //test_dpqb_building();
    //test_dpqb_building_1();
    test_d11_building_1();
    
    return 0;
}