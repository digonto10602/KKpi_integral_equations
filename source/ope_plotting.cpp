/* We use this script to generate ope plots 
   for different kinematic points */



#include<bits/stdc++.h>
#include "functions_mom_based.h"
#include "integral_equations_2plus1_system.h"
#include<omp.h> 

typedef std::complex<double> comp; 

//This is to print the function name and value 
#define Print_Func(x) std::cout << #x << " = " << x << std::endl;

/* This plots the GS11, GS12, GS21 for a set spectator 
momentum in the outgoing channel for different energies 
on the real axis, this is a full code and does not take 
input from outside */
/* We make sure that the outgoing spectator is such 
that the outgoing two particle sub-channel always forms a 
bound-state */
void Gs_plot_1( char debug  )
{
    comp ii = {0.0,1.0}; 
    double m1 = 1.0; 
    double m2 = 0.9; 

    double a0_1 = 2.0; //scattering length between m1 and m2
    double a0_2 = 2.0; //scattering length between m1 and m1

    double eps_for_ope = 0.0; 
    double eps_for_cutoff = 0.0; 

    comp sigb1 = sigma_b_plus(a0_1, m1, m2); //bound-state mass sq of m1 and m2
    comp sigb2 = sigma_b_plus(a0_2, m1, m1); //bound-state mass sq of m1 and m1

    double threeparticle_threshold_in_E = m1 + m1 + m2; 
    comp phib_threshold_1_in_E = std::sqrt(sigb1) + m1;
    comp phib_threshold_2_in_E = std::sqrt(sigb2) + m2; 

    if(debug=='y')
    {
        Print_Func(phib_threshold_1_in_E);
        Print_Func(phib_threshold_2_in_E); 
    }

    double En_initial = 0.1; 
    double En_final = threeparticle_threshold_in_E + m2; 
    double En_points = 100; 

    double delEn = std::abs(En_initial - En_final)/En_points; 

    int omp_ind = 0; 
    int countdown = 1; 
    int prev_count_val = 0;

    #pragma omp parallel for 

    for(omp_ind=0; omp_ind<(int)En_points; ++omp_ind)
    {
        double En = En_initial + omp_ind*delEn; 
        int mom_points = 250; 
        std::string filename = "ope_file_" + std::to_string(omp_ind) + ".dat"; 

        std::ofstream fout; 
        comp qb_1 = qb_i(En, sigb1, m1); 
        comp qb_2 = qb_i(En, sigb2, m2); 

        comp preal_val = 0.5; 
        comp preal_initial = -std::abs(preal_val);//-1.0; 
        comp preal_final = std::abs(preal_val); 
        comp del_preal = std::abs(preal_initial - preal_final)/mom_points; 

        comp pimag_val = 0.5; 
        comp pimag_initial = -std::abs(pimag_val); 
        comp pimag_final = std::abs(pimag_val); 
        comp del_pimag = std::abs(pimag_initial - pimag_final)/mom_points; 
        double mi, mj, mk; 

        fout.open(filename.c_str()); 
        
            
        fout<<"Re_En"<<'\t'
            <<"Im_En"<<'\t'
            <<"Re_qb_1"<<'\t'
            <<"Im_qb_1"<<'\t'
            <<"Re_qb_2"<<'\t'
            <<"Im_qb_2"<<'\t'
            <<"Re_p"<<'\t'
            <<"Im_p"<<'\t'
            <<"Re_G11"<<'\t'
            <<"Im_G11"<<'\t'
            <<"Re_G12"<<'\t'
            <<"Im_G12"<<'\t'
            <<"Re_G21"<<'\t'
            <<"Im_G21"<<std::endl; 

        for(int j=0; j<mom_points; ++j)
        {
            for(int k=0; k<mom_points; ++k)
            {
                comp preal = preal_initial + ((comp)j)*del_preal; 
                comp pimag = pimag_initial + ((comp)k)*del_pimag; 

                comp p = preal + ii*pimag; 

                //(i,j) = 1 1 and k = 2
                mi = m1; 
                mj = m1; 
                mk = m2; 
                comp GS11 = GS_pk(En, p, qb_1, mi, mj, mk, eps_for_ope, eps_for_cutoff);

                //(i,j) = 1 2 and k = 1 
                mi = m1; 
                mj = m2; 
                mk = m1; 
                comp GS12 = GS_pk(En, p, qb_2, mi, mj, mk, eps_for_ope, eps_for_cutoff); 

                //(i,j) = 2 1 and k = 1
                mi = m2; 
                mj = m1; 
                mk = m1; 
                comp GS21 = GS_pk(En, p, qb_1, mi, mj, mk, eps_for_ope, eps_for_cutoff); 

                fout<<std::setprecision(20);
                fout<<std::real(En)<<'\t'
                    <<std::imag(En)<<'\t'
                    <<std::real(qb_1)<<'\t'
                    <<std::imag(qb_1)<<'\t'
                    <<std::real(qb_2)<<'\t'
                    <<std::imag(qb_2)<<'\t'
                    <<std::real(p)<<'\t'
                    <<std::imag(p)<<'\t'
                    <<std::real(GS11)<<'\t'
                    <<std::imag(GS11)<<'\t'
                    <<std::real(GS12)<<'\t'
                    <<std::imag(GS12)<<'\t'
                    <<std::real(GS21)<<'\t'
                    <<std::imag(GS21)<<std::endl; 

            }
        }
        fout.close(); 

        #pragma omp critical
        {
            countdown = countdown + 1; 
        }

        int count_val = ((int)(countdown/En_points * 100.0))%10;
        int prog_val = (int)(countdown/En_points * 100.0);

        if(count_val==0)
        {
            if(prev_count_val==0)
            {
                prev_count_val = prog_val; 
                std::cout<<"Progress = "<<(countdown/En_points * 100.0)<<std::endl; 
            }
            else if(prev_count_val!=prog_val)
            {
                prev_count_val = prog_val; 
                std::cout<<"Progress = "<<(countdown/En_points * 100.0)<<std::endl; 
            }

        }


    }

    

}

/* This one takes input for mass m1 and m2, 
scattering length between m1&m2 called a0_1 and m1&m1 called a0_2
eps for ope and eps for cutoff 
initial energy and how much extra from the three particle threshold it is going to calculate,
along with the energy points, 
imaginary part of Energy, 
how many momentum points will be used for surface plots, 
preal value for setting the x axis 
pimag value for setting the y axis 
file initial
*/
void Gs_plot_2( char debug,
                double m1, 
                double m2, 
                double a0_1,
                double a0_2, 
                double eps_for_ope, 
                double eps_for_cutoff, 
                comp En_initial, 
                comp del_ext_En, 
                double En_points, 
                comp En_imag, 
                int mom_points, 
                comp preal_val, 
                comp pimag_val,
                std::string file_initial   )
{
    comp ii = {0.0,1.0}; 
    //double m1 = 1.0; 
    //double m2 = 0.9; 

    //double a0_1 = 2.0; //scattering length between m1 and m2
    //double a0_2 = 2.0; //scattering length between m1 and m1

    //double eps_for_ope = 0.0; 
    //double eps_for_cutoff = 0.0; 

    comp sigb1 = sigma_b_plus(a0_1, m1, m2); //bound-state mass sq of m1 and m2
    comp sigb2 = sigma_b_plus(a0_2, m1, m1); //bound-state mass sq of m1 and m1

    comp threeparticle_threshold_in_E = m1 + m1 + m2; 
    comp phib_threshold_1_in_E = std::sqrt(sigb1) + m1;
    comp phib_threshold_2_in_E = std::sqrt(sigb2) + m2; 

    if(debug=='y')
    {
        Print_Func(phib_threshold_1_in_E);
        Print_Func(phib_threshold_2_in_E); 
    }

    //double En_initial = 0.1; 
    comp En_final = threeparticle_threshold_in_E + del_ext_En; 
    //double En_points = 100; 

    comp delEn = std::abs(En_initial - En_final)/En_points; 

    int omp_ind = 0; 
    int countdown = 1; 
    int prev_count_val = 0;

    #pragma omp parallel for 

    for(omp_ind=0; omp_ind<(int)En_points; ++omp_ind)
    {
        comp En = En_initial + ((comp)omp_ind)*delEn + ii*En_imag; 
        //int mom_points = 250; 
        std::string filename = file_initial + std::to_string(omp_ind) + ".dat"; 

        std::ofstream fout; 
        comp qb_1 = qb_i(En, sigb1, m1); 
        comp qb_2 = qb_i(En, sigb2, m2); 

        //comp preal_val = 0.5; 
        comp preal_initial = -std::abs(preal_val);//-1.0; 
        comp preal_final = std::abs(preal_val); 
        comp del_preal = std::abs(preal_initial - preal_final)/mom_points; 

        //comp pimag_val = 0.5; 
        comp pimag_initial = -std::abs(pimag_val); 
        comp pimag_final = std::abs(pimag_val); 
        comp del_pimag = std::abs(pimag_initial - pimag_final)/mom_points; 
        double mi, mj, mk; 

        fout.open(filename.c_str()); 
        
            
        fout<<"Re_En"<<'\t'
            <<"Im_En"<<'\t'
            <<"Re_qb_1"<<'\t'
            <<"Im_qb_1"<<'\t'
            <<"Re_qb_2"<<'\t'
            <<"Im_qb_2"<<'\t'
            <<"Re_p"<<'\t'
            <<"Im_p"<<'\t'
            <<"Re_G11"<<'\t'
            <<"Im_G11"<<'\t'
            <<"Re_G12"<<'\t'
            <<"Im_G12"<<'\t'
            <<"Re_G21"<<'\t'
            <<"Im_G21"<<std::endl; 

        for(int j=0; j<mom_points; ++j)
        {
            for(int k=0; k<mom_points; ++k)
            {
                comp preal = preal_initial + ((comp)j)*del_preal; 
                comp pimag = pimag_initial + ((comp)k)*del_pimag; 

                comp p = preal + ii*pimag; 

                //(i,j) = 1 1 and k = 2
                mi = m1; 
                mj = m1; 
                mk = m2; 
                comp GS11 = GS_pk(En, p, qb_1, mi, mj, mk, eps_for_ope, eps_for_cutoff);

                //(i,j) = 1 2 and k = 1 
                mi = m1; 
                mj = m2; 
                mk = m1; 
                comp GS12 = GS_pk(En, p, qb_2, mi, mj, mk, eps_for_ope, eps_for_cutoff); 

                //(i,j) = 2 1 and k = 1
                mi = m2; 
                mj = m1; 
                mk = m1; 
                comp GS21 = GS_pk(En, p, qb_1, mi, mj, mk, eps_for_ope, eps_for_cutoff); 

                fout<<std::setprecision(20);
                fout<<std::real(En)<<'\t'
                    <<std::imag(En)<<'\t'
                    <<std::real(qb_1)<<'\t'
                    <<std::imag(qb_1)<<'\t'
                    <<std::real(qb_2)<<'\t'
                    <<std::imag(qb_2)<<'\t'
                    <<std::real(p)<<'\t'
                    <<std::imag(p)<<'\t'
                    <<std::real(GS11)<<'\t'
                    <<std::imag(GS11)<<'\t'
                    <<std::real(GS12)<<'\t'
                    <<std::imag(GS12)<<'\t'
                    <<std::real(GS21)<<'\t'
                    <<std::imag(GS21)<<std::endl; 

            }
        }
        fout.close(); 

        #pragma omp critical
        {
            countdown = countdown + 1; 
        }

        int count_val = ((int)(countdown/En_points * 100.0))%10;
        int prog_val = (int)(countdown/En_points * 100.0);

        if(count_val==0)
        {
            if(prev_count_val==0)
            {
                prev_count_val = prog_val; 
                std::cout<<"Progress = "<<(countdown/En_points * 100.0)<<"%"<<std::endl; 
            }
            else if(prev_count_val!=prog_val)
            {
                prev_count_val = prog_val; 
                std::cout<<"Progress = "<<(countdown/En_points * 100.0)<<"%"<<std::endl; 
            }

        }


    }
    std::cout<<"OPE File Generation Finished"<<std::endl; 
    std::cout<<"============================"<<std::endl; 

    

}

void GsM2_plot( char debug,
                double m1, 
                double m2, 
                double a0_1,
                double a0_2, 
                double eps_for_m2k, 
                double eps_for_ope, 
                double eps_for_cutoff,
                comp total_P, 
                comp En_initial, 
                comp del_ext_En, 
                double En_points, 
                comp En_imag, 
                int mom_points, 
                comp preal_val, 
                comp pimag_val,
                std::string file_initial   )
{
    comp ii = {0.0,1.0}; 
    //double total_P = 0.0; 
    //double m1 = 1.0; 
    //double m2 = 0.9; 

    //double a0_1 = 2.0; //scattering length between m1 and m2
    //double a0_2 = 2.0; //scattering length between m1 and m1

    //double eps_for_ope = 0.0; 
    //double eps_for_cutoff = 0.0; 

    comp sigb1 = sigma_b_plus(a0_1, m1, m2); //bound-state mass sq of m1 and m2
    comp sigb2 = sigma_b_plus(a0_2, m1, m1); //bound-state mass sq of m1 and m1

    comp threeparticle_threshold_in_E = m1 + m1 + m2; 
    comp phib_threshold_1_in_E = std::sqrt(sigb1) + m1;
    comp phib_threshold_2_in_E = std::sqrt(sigb2) + m2; 

    if(debug=='y')
    {
        Print_Func(phib_threshold_1_in_E);
        Print_Func(phib_threshold_2_in_E); 
    }

    //double En_initial = 0.1; 
    comp En_final = threeparticle_threshold_in_E + del_ext_En; 
    //double En_points = 100; 

    comp delEn = std::abs(En_initial - En_final)/En_points; 

    int omp_ind = 0; 
    int countdown = 1; 
    int prev_count_val = 0;

    #pragma omp parallel for 

    for(omp_ind=0; omp_ind<(int)En_points; ++omp_ind)
    {
        comp En = En_initial + ((comp)omp_ind)*delEn + ii*En_imag; 
        //int mom_points = 250; 
        std::string filename = file_initial + std::to_string(omp_ind) + ".dat"; 

        std::ofstream fout; 
        comp qb_1 = qb_i(En, sigb1, m1); 
        comp qb_2 = qb_i(En, sigb2, m2); 

        //comp preal_val = 0.5; 
        comp preal_initial = -std::abs(preal_val);//-1.0; 
        comp preal_final = std::abs(preal_val); 
        comp del_preal = std::abs(preal_initial - preal_final)/mom_points; 

        //comp pimag_val = 0.5; 
        comp pimag_initial = -std::abs(pimag_val); 
        comp pimag_final = std::abs(pimag_val); 
        comp del_pimag = std::abs(pimag_initial - pimag_final)/mom_points; 
        double mi, mj, mk; 
        double eta_i; 

        fout.open(filename.c_str()); 
        
            
        fout<<"Re_En"<<'\t'
            <<"Im_En"<<'\t'
            <<"Re_qb_1"<<'\t'
            <<"Im_qb_1"<<'\t'
            <<"Re_qb_2"<<'\t'
            <<"Im_qb_2"<<'\t'
            <<"Re_p"<<'\t'
            <<"Im_p"<<'\t'
            <<"Re_G11"<<'\t'
            <<"Im_G11"<<'\t'
            <<"Re_G12"<<'\t'
            <<"Im_G12"<<'\t'
            <<"Re_G21"<<'\t'
            <<"Im_G21"<<std::endl; 

        for(int j=0; j<mom_points; ++j)
        {
            for(int k=0; k<mom_points; ++k)
            {
                comp preal = preal_initial + ((comp)j)*del_preal; 
                comp pimag = pimag_initial + ((comp)k)*del_pimag; 

                comp p = preal + ii*pimag; 

                //(i,j) = 1 1 and k = 2
                mi = m1; 
                mj = m1; 
                mk = m2; 
                eta_i = 1.0; 
                //M211 = M2_1
                comp M211 = M2k_ERE(eta_i, En, p, total_P, a0_1, 0.0, mi, mj, mk, eps_for_m2k);
                comp GS11 = GS_pk(En, p, qb_1, mi, mj, mk, eps_for_ope, eps_for_cutoff);
                comp kern11 = GS11*M211; 
                //(i,j) = 1 2 and k = 1 
                mi = m1; 
                mj = m2; 
                mk = m1; 
                eta_i = 0.5; 
                //M212 = M2_2 
                comp M212 = M2k_ERE(eta_i, En, p, total_P, a0_2, 0.0, mi, mj, mk, eps_for_m2k);
                comp GS12 = GS_pk(En, p, qb_2, mi, mj, mk, eps_for_ope, eps_for_cutoff); 
                comp kern12 = GS12*M212; 
                //(i,j) = 2 1 and k = 1
                mi = m2; 
                mj = m1; 
                mk = m1; 
                eta_i = 1.0;
                //M221 = M2_1
                comp M221 = M2k_ERE(eta_i, En, p, total_P, a0_1, 0.0, mi, mj, mk, eps_for_m2k);
                comp GS21 = GS_pk(En, p, qb_1, mi, mj, mk, eps_for_ope, eps_for_cutoff); 
                comp kern21 = GS21*M221; 

                fout<<std::setprecision(20);
                fout<<std::real(En)<<'\t'
                    <<std::imag(En)<<'\t'
                    <<std::real(qb_1)<<'\t'
                    <<std::imag(qb_1)<<'\t'
                    <<std::real(qb_2)<<'\t'
                    <<std::imag(qb_2)<<'\t'
                    <<std::real(p)<<'\t'
                    <<std::imag(p)<<'\t'
                    <<std::real(kern11)<<'\t'
                    <<std::imag(kern11)<<'\t'
                    <<std::real(kern12)<<'\t'
                    <<std::imag(kern12)<<'\t'
                    <<std::real(kern21)<<'\t'
                    <<std::imag(kern21)<<std::endl; 

            }
        }
        fout.close(); 

        #pragma omp critical
        {
            countdown = countdown + 1; 
        }

        int count_val = ((int)(countdown/En_points * 100.0))%10;
        int prog_val = (int)(countdown/En_points * 100.0);

        if(count_val==0)
        {
            if(prev_count_val==0)
            {
                prev_count_val = prog_val; 
                std::cout<<"Progress = "<<(countdown/En_points * 100.0)<<"%"<<std::endl; 
            }
            else if(prev_count_val!=prog_val)
            {
                prev_count_val = prog_val; 
                std::cout<<"Progress = "<<(countdown/En_points * 100.0)<<"%"<<std::endl; 
            }

        }


    }
    std::cout<<"OPE File Generation Finished"<<std::endl; 
    std::cout<<"============================"<<std::endl; 

    

}


void parallel_ope_plotting( int number_of_files, 
                            int mom_points, 
                            std::string file_initial, 
                            std::string python_script )
{
    char* endptr;  
    //int number_of_files = std::strtol(argv[1], &endptr, 10);
    //std::string python_script = argv[2];
    int omp_ind = 0;
    //std::cout<<"running from "<<omp_ind<<std::endl; 

    int countdown = 1; 
    int prev_count_val = 0;
    std::cout<<"======================"<<std::endl; 
    std::cout<<"Plotting Files Started"<<std::endl; 

    #pragma omp parallel for 
    for(omp_ind=0; omp_ind<number_of_files; ++omp_ind)
    {
        std::string runcommand_str = "python  " + python_script + " " + file_initial + " " + std::to_string(omp_ind) + " " + std::to_string(mom_points);
        
        char const *runcommand = runcommand_str.c_str();
        //std::cout<<"running from "<<omp_ind<<std::endl; 
        
        int runfile = std::system(runcommand);

        #pragma omp critical
        {
            countdown = countdown + 1; 
        }

        int count_val = ((int)(countdown/((double)number_of_files) * 100.0))%10;
        int prog_val = (int)(countdown/((double)number_of_files) * 100.0);

        if(count_val==0)
        {
            if(prev_count_val==0)
            {
                prev_count_val = prog_val; 
                std::cout<<"Progress = "<<(countdown/((double)number_of_files) * 100.0)<<"%"<<std::endl; 
            }
            else if(prev_count_val!=prog_val)
            {
                prev_count_val = prog_val; 
                std::cout<<"Progress = "<<(countdown/((double)number_of_files) * 100.0)<<"%"<<std::endl; 
            }

        }
    }

    std::cout<<"Plot Files Generated"<<std::endl;
    std::cout<<"===================="<<std::endl; 
}



int main(int argc, char* argv[])
{
    char* endptr; 
    
    if(argc==1)
    {
        std::cout<<"input needed: <debug_y/n> <m1> <m2> "
                 <<"<a0_1> <a0_2> <eps_for_m2> <eps_for_ope> <eps_for_cutoff> "
                 <<"<En_initial> <total_P> <del_ext_En> <En_points> <En_imag> "
                 <<"<mom_points> <preal_val> <pimag_val> "
                 <<"<file_initial> <dest to python plot script>"<<std::endl; 
        std::exit(121); 
    }
    else 
    {
        std::string debug1 = argv[1];
        int debug1_length = debug1.length(); 
        char debug;
        if(debug1_length>1)
        {
            std::cout<<"debug is one letter, either y or n"<<std::endl; 
            std::exit(11);
        }
        else 
        {
            debug = debug1[0];
        } 
        //std::strcpy(debug, debug1.c_str());  

        double m1 = std::strtod(argv[2], &endptr); 
        double m2 = std::strtod(argv[3], &endptr); 
        double a0_1 = std::strtod(argv[4], &endptr); 
        double a0_2 = std::strtod(argv[5], &endptr); 
        double eps_for_m2k = std::strtod(argv[6], &endptr); 
        double eps_for_ope = std::strtod(argv[7], &endptr); 
        double eps_for_cutoff = std::strtod(argv[8], &endptr); 
        double En_initial = std::strtod(argv[9], &endptr); 
        double total_P = std::strtod(argv[10], &endptr); 
        double del_ext_En = std::strtod(argv[11], &endptr); 
        double En_points = std::strtod(argv[12], &endptr); 
        double En_imag = std::strtod(argv[13], &endptr); 
        double mom_points = std::strtod(argv[14], &endptr); 
        double preal_val = std::strtod(argv[15], &endptr); 
        double pimag_val = std::strtod(argv[16], &endptr); 
        std::string file_initial = argv[17]; 
        std::string dest_to_py_script = argv[18]; 

        /*Gs_plot_2(  debug, m1, m2, a0_1, a0_2, eps_for_ope, eps_for_cutoff, 
                    En_initial, del_ext_En, En_points, En_imag, 
                    mom_points, preal_val, pimag_val, file_initial );
        */
        
        GsM2_plot(  debug, m1, m2, a0_1, a0_2, eps_for_m2k, eps_for_ope, eps_for_cutoff, 
                    En_initial, total_P, del_ext_En, En_points, En_imag, 
                    mom_points, preal_val, pimag_val, file_initial );

        
        parallel_ope_plotting(En_points, mom_points, file_initial, dest_to_py_script);
        
    }
    
    //char debug = 'y';

    //Gs_plot_1(debug);

    return 0; 
}