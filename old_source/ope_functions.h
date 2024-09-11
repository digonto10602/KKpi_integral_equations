#ifndef OPE_FUNCTIONS_H
#define OPE_FUNCTIONS_H
#include<bits/stdc++.h>
#include<Eigen/Dense>
//#include "functions_momrep_based.h"



//This GS uses 2 indices and 2 tags, if the indices are between the 
//tags then the OPE is evaluated on the unphysical sheet by adding the 
//discontinuity with the original OPE 
comp GS_pk_withtags(    comp s,
                        comp p,
                        comp k,
                        double m,
                        double eps,
                        int index1,
                        int index2,
                        int tag1,
                        int tag2  )
{
    comp ii = {0.0,1.0};
    comp sigp = sigma_p(s,p,m);
    comp sigk = sigma_p(s,k,m);
    comp alphapk = (sqrt(s)-omega_comp(p,m)-omega_comp(k,m))*(sqrt(s)-omega_comp(p,m)-omega_comp(k,m)) - p*p - k*k - m*m;
    comp num = alphapk - 2.0*p*k + ii*eps;
    comp denom = alphapk + 2.0*p*k + ii*eps;
    double pi = acos(-1.0);

    //if(abs(p)<1.0e-10 || abs(k)<1.0e-10)
    //return Hfunc_comp(sigp,sigk,m)/(alphapk + ii*eps);
    //return -1.0/(alphapk + ii*eps);
    
    //else
    {
        if((index1>=tag1 && index1<=tag2) || (index2>=tag1 && index2<=tag2))
        //smooth cutoff
        return -(Hfunc_comp(sigp,sigk,m)/(4.0*p*k))*(log(num/denom)) -(Hfunc_comp(sigp,sigk,m)/(4.0*p*k))*(2.0*pi*ii)  ;
        //hard cutoff
        //return -(1.0/(4.0*p*k))*(log(num/denom)) -(1.0/(4.0*p*k))*(2.0*pi*ii)  ;
        
        else 
        //smooth cutoff
        return -(Hfunc_comp(sigp,sigk,m)/(4.0*p*k))*(log(num/denom));
        
        //hard cutoff 
        //return -(1.0/(4.0*p*k))*(log(num/denom));
    
    }

    
}

comp GS_pk_secondsheet(    comp s,
                        comp p,
                        comp k,
                        double m,
                        double eps )
{
    comp ii = {0.0,1.0};
    comp sigp = sigma_p(s,p,m);
    comp sigk = sigma_p(s,k,m);
    comp alphapk = (sqrt(s)-omega_comp(p,m)-omega_comp(k,m))*(sqrt(s)-omega_comp(p,m)-omega_comp(k,m)) - p*p - k*k - m*m;
    comp num = alphapk - 2.0*p*k + ii*eps;
    comp denom = alphapk + 2.0*p*k + ii*eps;
    double pi = acos(-1.0);

    //if(abs(p)<1.0e-10 || abs(k)<1.0e-10)
    //return Hfunc_comp(sigp,sigk,m)/(alphapk + ii*eps);
    //return -1.0/(alphapk + ii*eps);
    
    //smooth cutoff
    return -(Hfunc_comp(sigp,sigk,m)/(4.0*p*k))*(log(num/denom)) -(Hfunc_comp(sigp,sigk,m)/(4.0*p*k))*(2.0*pi*ii)  ;
    //hard cutoff 
    //return -(1.0/(4.0*p*k))*(log(num/denom)) -(1.0/(4.0*p*k))*(2.0*pi*ii)  ;
        
    
}


//------------------------------------------------------------//

comp sigmab(    double a,
                double m    )
{
    return 4.0*(m*m - 1.0/(a*a));
}


void Gvec_maker_momrep(     Eigen::VectorXcd &Gvec,
                            comp s,
                            vector<comp> &qvec,
                            comp p,
                            double m,
                            double epsilon  )
{
    for(int i=0;i<qvec.size();++i)
    {
        Gvec(i) = GS_pk(s,p,qvec[i],m,epsilon);
    }
}

void Gmat_maker_momrep(     Eigen::MatrixXcd &Gmat,
                            comp s,
                            vector<comp> &qvec,
                            vector<comp> &qvec1,
                            double m,
                            double epsilon  )
{
    for(int i=0;i<qvec.size();++i)
    {
        for(int j=0;j<qvec1.size();++j)
        {
            Gmat(i,j) = GS_pk(s,qvec[i],qvec1[j],m,epsilon);
        }
        
    }
}

void Gvec_maker_momrep_withtags(        Eigen::VectorXcd &Gvec,
                                        comp s,
                                        vector<comp> &qvec,
                                        comp p,
                                        double m,
                                        double epsilon,
                                        int tag1,
                                        int tag2   )
{
    for(int i=0;i<qvec.size();++i)
    {
        int index = i;
        if(tag1==0 && tag2==0)
        {
            Gvec(i) = GS_pk(s,p,qvec[i],m,epsilon);
        }
        else 
        {
            Gvec(i) = GS_pk_withtags(s,p,qvec[i],m,epsilon,index,index,tag1,tag2);
        }
    }
}

//this is a much simpler definition
void Gvec_maker_momrep_withtags_1(        Eigen::VectorXcd &Gvec,
                                        comp s,
                                        vector<comp> &qvec,
                                        comp p,
                                        double m,
                                        double epsilon,
                                        int tag1,
                                        int tag2   )
{
    for(int i=0;i<qvec.size();++i)
    {
        int index = i;
        if(tag1==0 && tag2==0)
        {
            Gvec(i) = GS_pk(s,p,qvec[i],m,epsilon);
        }
        else 
        {
            //Gvec(i) = GS_pk_withtags(s,p,qvec[i],m,epsilon,index,index,tag1,tag2);
            if(i>=tag1 && i<=tag2)
            {
                Gvec(i) = GS_pk_secondsheet(s,p,qvec[i],m,epsilon);
            }
            else
            {
                Gvec(i) = GS_pk(s,p,qvec[i],m,epsilon);   
            }
        }
    }
}


//this gvec_fixer is based on contour_imspos_2
//where we have one discontinuity to left 
//and one similarity to right
void Gvec_fixer_3(    Eigen::VectorXcd &Gvec,
                    Eigen::VectorXcd &fixed_Gvec,
                    vector<comp> qvec,
                    comp s,
                    comp k,
                    double a,
                    double m,
                    double eps,
                    double eps_for_m2k,
                    int &tag1,
                    int &tag2    )
{

    comp kmax = pmom(s,0.0,m);
    int switch_for_gvec_fixer = 0;// this is not used here
    vector<comp> temp_qvec;
    vector<comp> temp_weights;
    int points1 = 10000;
    int temp_tag1,temp_tag2;
    //mom_vector_maker_seba_imspos(temp_qvec,temp_weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,temp_tag1,temp_tag2,switch_for_gvec_fixer);
    //mom_vector_maker_47_with_weights(temp_qvec,temp_weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,0.00101,temp_tag1,temp_tag2);
    mom_vector_maker_seba_imspos_2(temp_qvec,temp_weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,temp_tag1,temp_tag2,switch_for_gvec_fixer);
                                
    Eigen::VectorXcd temp_Gvec(temp_qvec.size());
    temp_Gvec[0] = Gvec[0];

    double past_difference = 0.0;
    vector<int> close_indices;
    vector<double> gap;
    for(int i=0;i<temp_qvec.size()-1;++i)
    {
        comp p = temp_qvec[i];
        comp future_p = temp_qvec[i+1];
        comp GS_firstsheet = GS_pk(s,p,k,m,eps);
        //comp GS_secondsheet = GS_pk_secondsheet(s,p,k,m,eps);
        comp future_GS_firstsheet = GS_pk(s,future_p,k,m,eps);
        //comp future_GS_secondsheet = GS_pk_secondsheet(s,future_p,k,m,eps);
        double realGS_firstsheet = real(GS_firstsheet);
        //double realGS_secondsheet = real(GS_secondsheet);
        double future_realGS_firstsheet = real(future_GS_firstsheet);
        //double future_realGS_secondsheet = real(future_GS_secondsheet);
        
        double present_difference = abs(realGS_firstsheet - future_realGS_firstsheet);
        //double future_difference = abs(future_realGS_firstsheet - future_realGS_secondsheet);

        close_indices.push_back(i);
        gap.push_back(present_difference);
        temp_Gvec[i] = GS_firstsheet;
        
    }

    //cout<<"ran here"<<endl;

    double temp_max1 = 0.0;
    int temp_ind1 = 0;
    vector<int> selected_indices;
    vector<double> selected_gap;
    for(int i=0;i<gap.size();++i)
    {
        if(gap[i]>temp_max1)
        {
            temp_max1 = gap[i];
            temp_ind1 = i;
        }
    }

    double temp_max2 = INT_MAX;
    int temp_ind2 = 0;

    for(int i=temp_tag1;i<=temp_tag2;++i)
    {
        comp p = temp_qvec[i];
        comp future_p = temp_qvec[i+1];
        comp GS_firstsheet = GS_pk(s,p,k,m,eps);
        comp GS_secondsheet = GS_pk_secondsheet(s,p,k,m,eps);
        //comp future_GS_firstsheet = GS_pk(s,future_p,k,m,eps);
        //comp future_GS_secondsheet = GS_pk_secondsheet(s,future_p,k,m,eps);
        double realGS_firstsheet = real(GS_firstsheet);
        double realGS_secondsheet = real(GS_secondsheet);
        //double future_realGS_firstsheet = real(future_GS_firstsheet);
        //double future_realGS_secondsheet = real(future_GS_secondsheet);
        
        //double present_difference = abs(realGS_firstsheet - future_realGS_firstsheet);
        double present_difference = abs(realGS_firstsheet - realGS_secondsheet);
        
        //double future_difference = abs(future_realGS_firstsheet - future_realGS_secondsheet);

        //close_indices.push_back(i);
        //gap.push_back(present_difference);
        //temp_Gvec[i] = GS_firstsheet;
        if(present_difference<=temp_max2)
        {
            temp_ind2 = i;
            temp_max2 = present_difference;
        }
        
    }

    /*for(int i=0;i<gap.size();++i)
    {
        if(i==temp_ind1)
        {
            continue;
        }
        else
        {
            if(gap[i]>temp_max2)
            {
                temp_max2 = gap[i];
                temp_ind2 = i;
            }
        }
    }*/

    if(temp_ind1>temp_ind2)
    {
        int sometemp = temp_ind2;
        temp_ind2 = temp_ind1;
        temp_ind1 = sometemp;
    }

    cout<<"tags before rescaling = "<<temp_ind1<<'\t'<<temp_ind2<<endl;

    comp temp_gvec_tag1_val = temp_Gvec[temp_ind1];
    comp temp_gvec_tag2_val = temp_Gvec[temp_ind2];
    
    int rescaled_tag1 = ceil(temp_ind1*qvec.size()/temp_qvec.size()) + 1 ;
    int rescaled_tag2 = ceil(temp_ind2*qvec.size()/temp_qvec.size()) - 1;
    //int rescaled_tag1 = floor(temp_ind1*qvec.size()/temp_qvec.size());
    //int rescaled_tag2 = floor(temp_ind2*qvec.size()/temp_qvec.size());
    
    
    cout<<"min ind = "<<rescaled_tag1<<'\t'<<" max ind = "<<rescaled_tag2<<endl;
    cout<<"gap 1 = "<<temp_max1<<'\t'<<" gap 2 = "<<temp_max2<<endl;
    for(int i=0;i<qvec.size();++i)
    {
        comp p = qvec[i];
        comp GS_firstsheet = GS_pk(s,p,k,m,eps);
        comp GS_secondsheet = GS_pk_secondsheet(s,p,k,m,eps);

        if(i>=rescaled_tag1 && i<=rescaled_tag2)
            fixed_Gvec[i] = GS_secondsheet;
        else    
            fixed_Gvec[i] = GS_firstsheet;
        
       //if(real(GS_firstsheet)>)
    }

    tag1 = rescaled_tag1;//temp_ind1;
    tag2 = rescaled_tag2;//temp_ind2;

    /*for(int i=0;i<fixed_Gvec.size()-1;++i)
    {
        double real_gvec = real(fixed_Gvec[i]);
        double real_gvec_next = real(fixed_Gvec[i+1]);
        comp p = qvec[i+1];
        double real_GS_one = real(GS_pk(s,p,k,m,eps));

        double difference1 = abs(real_gvec - real_gvec_next);
        double difference2 = abs(real_gvec - real_GS_one);
        if(difference1>difference2)
        {
            fixed_Gvec[i+1] = GS_pk(s,p,k,m,eps);
        }
    }

    for(int i=0;i<fixed_Gvec.size()-1;++i)
    {
        double real_gvec = real(fixed_Gvec[i]);
        double real_gvec_next = real(fixed_Gvec[i+1]);
        comp p = qvec[i+1];
        double real_GS_one = real(GS_pk_secondsheet(s,p,k,m,eps));

        double difference1 = abs(real_gvec - real_gvec_next);
        double difference2 = abs(real_gvec - real_GS_one);
        if(difference1>difference2)
        {
            fixed_Gvec[i+1] = GS_pk_secondsheet(s,p,k,m,eps);
        }
    }*/
}


void tag_fixer_for_gvec_fixer3( Eigen::VectorXcd &Gvec,
                vector<comp> &qvec,
                comp s, 
                comp k,
                double m,
                double eps,
                int &problem_tag     )
{
    int tag_initial = problem_tag - 5;
    int tag_final = problem_tag + 5;
    int temp_tag = 0;
    int somecount = 0;

    if(tag_initial<0)
        tag_initial=0;
    if(tag_final<0)
        tag_final=0;
    if(tag_final>qvec.size()-1)
        tag_final=qvec.size()-1;

    for(int i=tag_initial;i<tag_final;++i)
    {
        comp p = qvec[i];
        
        comp Gvecval = Gvec[i];
        comp Gvecvalnext = Gvec[i+1]; 
        comp Gs1stsheet = GS_pk(s,p,k,m,eps);
        comp Gs2ndsheet = GS_pk_secondsheet(s,p,k,m,eps);

        //comp GS_secondsheet = GS_pk_secondsheet(s,p,k,m,eps);
        //comp future_GS_firstsheet = GS_pk(s,future_p,k,m,eps);
        //comp future_GS_secondsheet = GS_pk_secondsheet(s,future_p,k,m,eps);
        //double realGS_firstsheet = real(GS_firstsheet);
        //double realGS_secondsheet = real(GS_secondsheet);

        comp pnext = qvec[i+1];
        comp Gs1stsheetnext = GS_pk(s,pnext,k,m,eps);
        comp Gs2ndsheetnext = GS_pk_secondsheet(s,pnext,k,m,eps);
        
        double first_difference = abs(real(Gvecval) - real(Gvecvalnext));
        double second_difference = abs(real(Gvecval) - real(Gs2ndsheetnext));
        
        if(first_difference<=second_difference)
        {
            continue;
        }
        else
        {
            Gvec[i+1] = GS_pk_secondsheet(s,pnext,k,m,eps);
            problem_tag = i+1;

        }
        
        
    }

    
}


#endif