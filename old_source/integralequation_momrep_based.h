#ifndef INTEGRALEQUATION_MOMREP_H
#define INTEGRALEQUATION_MOMREP_H

#include<bits/stdc++.h>
#include<Eigen/Dense>
using namespace std;

typedef complex<double> comp;

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



void Bmat_maker_momrep(     Eigen::MatrixXcd &Bmat,
                            comp s,
                            vector<comp> &p,
                            vector<comp> &k,
                            double a,
                            double m,
                            double epsilon  )
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
                delk = k[1]- k[0];
            }
            else
            {
                delk = k[j] - k[j-1];
            }

            //make sure you supply both the invariant vectors
            //as the same, otherwise this condition will fail
            if(j==i)
            {
                comp one = {1.0,0.0};
                Bmat(i,j) = one + delk*kernel_pk(s,p[i],k[j],a,m,epsilon);
                //Bmat(i,j) = one + kernel_pk(s,p[i],k[j],a,m,epsilon);
            }
            else 
            {
                Bmat(i,j) = delk*kernel_pk(s,p[i],k[j],a,m,epsilon);
                //Bmat(i,j) = kernel_pk(s,p[i],k[j],a,m,epsilon);
            }
        }
    }

}

void Bmat_maker_momrep_2eps(        Eigen::MatrixXcd &Bmat,
                                    comp s,
                                    vector<comp> &p,
                                    vector<comp> &k,
                                    double a,
                                    double m,
                                    double epsilon,
                                    double eps_for_m2k  )
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

            //make sure you supply both the invariant vectors
            //as the same, otherwise this condition will fail
            if(j==i)
            {
                comp one = {1.0,0.0};
                Bmat(i,j) = one + delk*kernel_pk_2eps(s,p[i],k[j],a,m,epsilon,eps_for_m2k);
                //Bmat(i,j) = one + kernel_pk_2eps(s,p[i],k[j],a,m,epsilon,eps_for_m2k);
            }
            else 
            {
                Bmat(i,j) = delk*kernel_pk_2eps(s,p[i],k[j],a,m,epsilon,eps_for_m2k);
                //Bmat(i,j) = kernel_pk_2eps(s,p[i],k[j],a,m,epsilon,eps_for_m2k);
            }
        }
    }

}

void Bmat_maker_momrep_2eps_withtags(   Eigen::MatrixXcd &Bmat,
                                        comp s,
                                        vector<comp> &p,
                                        vector<comp> &k,
                                        double a,
                                        double m,
                                        double epsilon,
                                        double eps_for_m2k,
                                        int tag1,
                                        int tag2  )
{
    comp ii = {0.0,1.0};
    comp pi = acos(-1.0);
    comp delk = {0.0,0.0};

    for(int i=0;i<p.size();++i)
    {
        for(int j=0;j<k.size();++j)
        {
            int index1 = i;
            int index2 = j;
            if(j==0)
            {
                delk = k[1] - k[0];
            }
            else
            {
                delk = k[j] - k[j-1];
            }

            //make sure you supply both the invariant vectors
            //as the same, otherwise this condition will fail
            if(j==i)
            {
                comp one = {1.0,0.0};
                if(tag1==0 && tag2==0)
                {
                    Bmat(i,j) = one + delk*kernel_pk_2eps(s,p[i],k[j],a,m,epsilon,eps_for_m2k);
                }
                else 
                {
                    Bmat(i,j) = one + delk*kernel_pk_2eps_withtags(s,p[i],k[j],a,m,epsilon,eps_for_m2k,index1,index2,tag1,tag2);
                }
            }
            else 
            {
                if(tag1==0 && tag2==0)
                {
                    Bmat(i,j) = delk*kernel_pk_2eps(s,p[i],k[j],a,m,epsilon,eps_for_m2k);
                }
                else 
                {
                    Bmat(i,j) = delk*kernel_pk_2eps_withtags(s,p[i],k[j],a,m,epsilon,eps_for_m2k,index1,index2,tag1,tag2);
                }
            }
        }
    }

}

void Bmat_maker_momrep_2eps_withtags_with_weights(   Eigen::MatrixXcd &Bmat,
                                        comp s,
                                        vector<comp> &p,
                                        vector<comp> &k,
                                        vector<comp> &weights,
                                        double a,
                                        double m,
                                        double epsilon,
                                        double eps_for_m2k,
                                        int tag1,
                                        int tag2  )
{
    comp ii = {0.0,1.0};
    comp pi = acos(-1.0);
    comp delk = {0.0,0.0};

    for(int i=0;i<p.size();++i)
    {
        for(int j=0;j<k.size();++j)
        {
            int index1 = i;
            int index2 = j;
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
                if(tag1==0 && tag2==0)
                {
                    Bmat(i,j) = one + weight*kernel_pk_2eps(s,p[i],k[j],a,m,epsilon,eps_for_m2k);
                }
                else 
                {
                    Bmat(i,j) = one + weight*kernel_pk_2eps_withtags(s,p[i],k[j],a,m,epsilon,eps_for_m2k,index1,index2,tag1,tag2);
                }
            }
            else 
            {
                if(tag1==0 && tag2==0)
                {
                    Bmat(i,j) = weight*kernel_pk_2eps(s,p[i],k[j],a,m,epsilon,eps_for_m2k);
                }
                else 
                {
                    Bmat(i,j) = weight*kernel_pk_2eps_withtags(s,p[i],k[j],a,m,epsilon,eps_for_m2k,index1,index2,tag1,tag2);
                }
            }
        }
    }

}


void Bmat_maker_momrep_2eps_with_weights(       Eigen::MatrixXcd &Bmat,
                                                comp s,
                                                vector<comp> &p,
                                                vector<comp> &k,
                                                vector<comp> &weights,
                                                double a,
                                                double m,
                                                double epsilon,
                                                double eps_for_m2k  )
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
                Bmat(i,j) = one + weight*kernel_pk_2eps(s,p[i],k[j],a,m,epsilon,eps_for_m2k);
                //Bmat(i,j) = one + weight*kernel_pk_2eps_1(s,p[i],k[j],a,m,epsilon,eps_for_m2k);
                //Bmat(i,j) =  weight*kernel_pk_2eps(s,p[i],k[j],a,m,epsilon,eps_for_m2k);
            
            }
            else 
            {
                Bmat(i,j) = weight*kernel_pk_2eps(s,p[i],k[j],a,m,epsilon,eps_for_m2k);
                //Bmat(i,j) = weight*kernel_pk_2eps_1(s,p[i],k[j],a,m,epsilon,eps_for_m2k);
            }
        }
    }

}


void Bmat_maker_momrep_2eps_with_weights_rotated_m2k_for_resonances(       Eigen::MatrixXcd &Bmat,
                                                comp s,
                                                vector<comp> &p,
                                                vector<comp> &k,
                                                vector<comp> &weights,
                                                double a,
                                                double m,
                                                double epsilon,
                                                double eps_for_m2k,
                                                double theta )
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
                //sheet_tag=0 for first sheet
                //sheet_tag=1 for second sheet
                comp one = {1.0,0.0};
                Bmat(i,j) = one + weight*kernel_pk_2eps_rotated_m2k_for_resonances(s,p[i],k[j],a,m,epsilon,eps_for_m2k,theta);

            
            }
            else 
            {
                Bmat(i,j) = weight*kernel_pk_2eps_rotated_m2k_for_resonances(s,p[i],k[j],a,m,epsilon,eps_for_m2k,theta);
                
            }
        }
    }

}


//this handles sheet2 of dS
void Bmat_maker_momrep_2eps_with_weights_for_resonances(       Eigen::MatrixXcd &Bmat,
                                                comp s,
                                                vector<comp> &p,
                                                vector<comp> &k,
                                                vector<comp> &weights,
                                                double a,
                                                double m,
                                                double epsilon,
                                                double eps_for_m2k,
                                                double theta,
                                                int tag_for_m2k  )
{
    comp ii = {0.0,1.0};
    comp pi = acos(-1.0);
    comp delk = {0.0,0.0};

    for(int i=0;i<p.size();++i)
    {
        for(int j=0;j<k.size();++j)
        {
            

            comp weight = weights[j];
            //make sure you supply both the invariant vectors
            //as the same, otherwise this condition will fail
            if(j==i)
            {
                //sheet_tag=0 for first sheet
                //sheet_tag=1 for second sheet
                comp one = {1.0,0.0};
                if(j<=tag_for_m2k)
                Bmat(i,j) = one + weight*kernel_pk_2eps_for_resonances(s,p[i],k[j],a,m,epsilon,eps_for_m2k,theta,1);
                else 
                Bmat(i,j) = one + weight*kernel_pk_2eps_for_resonances(s,p[i],k[j],a,m,epsilon,eps_for_m2k,theta,0);

            
            }
            else 
            {
                if(j<=tag_for_m2k)
                Bmat(i,j) = weight*kernel_pk_2eps_for_resonances(s,p[i],k[j],a,m,epsilon,eps_for_m2k,theta,1);
                else 
                Bmat(i,j) = weight*kernel_pk_2eps_for_resonances(s,p[i],k[j],a,m,epsilon,eps_for_m2k,theta,0);//0 before
                
            }
        }
    }

}

void Bmat_maker_momrep_2eps_with_weights_for_resonances_sheet3(       Eigen::MatrixXcd &Bmat,
                                                comp s,
                                                vector<comp> &p,
                                                vector<comp> &k,
                                                vector<comp> &weights,
                                                double a,
                                                double m,
                                                double epsilon,
                                                double eps_for_m2k,
                                                double theta,
                                                int tag_for_m2k_1,
                                                int tag_for_m2k_2  )
{
    comp ii = {0.0,1.0};
    comp pi = acos(-1.0);
    comp delk = {0.0,0.0};

    for(int i=0;i<p.size();++i)
    {
        for(int j=0;j<k.size();++j)
        {
            

            comp weight = weights[j];
            //make sure you supply both the invariant vectors
            //as the same, otherwise this condition will fail
            if(j==i)
            {
                //sheet_tag=0 for first sheet
                //sheet_tag=1 for second sheet
                comp one = {1.0,0.0};
                if(j>=tag_for_m2k_1 && j<=tag_for_m2k_2)
                Bmat(i,j) = one + weight*kernel_pk_2eps_for_resonances(s,p[i],k[j],a,m,epsilon,eps_for_m2k,theta,1);
                else 
                Bmat(i,j) = one + weight*kernel_pk_2eps_for_resonances(s,p[i],k[j],a,m,epsilon,eps_for_m2k,theta,0);

            
            }
            else 
            {
                if(j>=tag_for_m2k_1 && j<=tag_for_m2k_2)
                Bmat(i,j) = weight*kernel_pk_2eps_for_resonances(s,p[i],k[j],a,m,epsilon,eps_for_m2k,theta,1);
                else 
                Bmat(i,j) = weight*kernel_pk_2eps_for_resonances(s,p[i],k[j],a,m,epsilon,eps_for_m2k,theta,0);//0 before
                
            }
        }
    }

}


comp rho3_func( comp s,
                comp k,
                double m    )
{
    comp s2k = sigma_p(s,k,m);
    comp ii = {0.0,1.0};
    if(real(s2k)>4.0*m*m)
    {
        return (1.0/(16.0*pi*sqrt(s2k)))*(-ii*sqrt(s2k/4.0 - m*m));
    }
    else 
    {
        return (1.0/(16.0*pi*sqrt(s2k)))*(abs(sqrt(s2k/4.0 - m*m)));
    }
}

//based on our efimov paper
comp rho3_func_1( comp s,
                comp k,
                double m    )
{
    comp s2k = sigma_p(s,k,m);
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    return (1.0/(16.0*pi*sqrt(s2k)))*(-ii*mysqrt(s2k/4.0 - m*m));
    
}

void F3infvol(  Eigen::MatrixXcd &dsol,
                vector<comp> &qvec1,
                vector<comp> &weights1,
                vector<comp> &qvec2,
                vector<comp> &weights2,
                comp s,
                double a,
                double m,
                double eps,
                comp &result  )
{
    comp firstterm = {0.0,0.0};
    comp secondterm = {0.0,0.0};
    double pi = acos(-1.0);
    
    for(int i=0; i<qvec1.size(); ++i)
    {
        comp rho3_i = rho3_func(s,qvec1[i],m);
        comp sig_i = sigma_p(s,qvec1[i],m);
        comp m2k_i = M2kfunc(a,sig_i,m,eps);
        
        comp temp_firstterm = weights1[i]*qvec1[i]*qvec1[i]/(pow(2.0*pi,2.0)*omega_comp(qvec1[i],m))*(rho3_i/3.0 - rho3_i*m2k_i*rho3_i);
        firstterm = firstterm + temp_firstterm;
        for(int j=0; j<qvec2.size(); ++j)
        {
            comp rho3_j = rho3_func(s,qvec2[j],m);
            comp sig_j = sigma_p(s,qvec2[j],m);
            comp m2k_j = M2kfunc(a,sig_j,m,eps);
            //cout<<rho3_i<<'\t'<<rho3_j<<endl;

            comp temp_phasespace1 = weights1[i]*qvec1[i]*qvec1[i]/(pow(2.0*pi,2.0)*omega_comp(qvec1[i],m));
            comp temp_phasespace2 = weights2[j]*qvec2[j]*qvec2[j]/(pow(2.0*pi,2.0)*omega_comp(qvec2[j],m));

            comp temp_secondterm = temp_phasespace1*temp_phasespace2*rho3_i*m2k_i*dsol(i,j)*m2k_j*rho3_j;

            secondterm = secondterm + temp_secondterm;
        }
    }


    result = firstterm - secondterm;


}


//the loops are separated
void F3infvol_1(  Eigen::MatrixXcd &dsol,
                vector<comp> &qvec1,
                vector<comp> &weights1,
                vector<comp> &qvec2,
                vector<comp> &weights2,
                comp s,
                double a,
                double m,
                double eps,
                comp &result  )
{
    comp firstterm = {0.0,0.0};
    comp secondterm = {0.0,0.0};
    double pi = acos(-1.0);
    
    for(int i=0; i<qvec1.size(); ++i)
    {
        comp rho3_i = rho3_func(s,qvec1[i],m);
        comp sig_i = sigma_p(s,qvec1[i],m);
        comp m2k_i = M2kfunc(a,sig_i,m,eps);
        
        comp temp_firstterm = weights1[i]*qvec1[i]*qvec1[i]/(pow(2.0*pi,2.0)*omega_comp(qvec1[i],m))*(rho3_i/3.0 - rho3_i*m2k_i*rho3_i);
        firstterm = firstterm + temp_firstterm;
    }

    for(int i=0; i<qvec1.size(); ++i)
    {   
        comp rho3_i = rho3_func(s,qvec1[i],m);
        comp sig_i = sigma_p(s,qvec1[i],m);
        comp m2k_i = M2kfunc(a,sig_i,m,eps);
        for(int j=0; j<qvec2.size(); ++j)
        {
            comp rho3_j = rho3_func(s,qvec2[j],m);
            comp sig_j = sigma_p(s,qvec2[j],m);
            comp m2k_j = M2kfunc(a,sig_j,m,eps);
            //cout<<rho3_i<<'\t'<<rho3_j<<endl;

            comp temp_phasespace1 = weights1[i]*qvec1[i]*qvec1[i]/(pow(2.0*pi,2.0)*omega_comp(qvec1[i],m));
            comp temp_phasespace2 = weights2[j]*qvec2[j]*qvec2[j]/(pow(2.0*pi,2.0)*omega_comp(qvec2[j],m));

            comp temp_secondterm = temp_phasespace1*temp_phasespace2*rho3_i*m2k_i*dsol(i,j)*m2k_j*rho3_j;

            secondterm = secondterm + temp_secondterm;
        }
    }


    result = firstterm - secondterm;


}




#endif