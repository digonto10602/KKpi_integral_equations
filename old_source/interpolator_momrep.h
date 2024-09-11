#ifndef INTERPOLATOR_MOMREP_H
#define INTERPOLATOR_MOMREP_H

#include<bits/stdc++.h>
#include<Eigen/Dense>
using namespace std;

typedef complex<double> comp;

void interpolator_ds_integraleq_momrep(     Eigen::VectorXcd &dsol,
                                            vector<comp> &qvec,
                                            comp s,
                                            comp p,
                                            comp k,
                                            double a,
                                            double m,
                                            double eps,
                                            comp &result  )
{
    comp ii = {0.0,1.0};

    comp Gs = GS_pk(s,p,k,m,eps);
    comp delq;
    comp integralsum = {0.0,0.0};

    for(int i=0;i<qvec.size();++i)
    {
        if(i==0)
        {
            delq = qvec[1] - qvec[0];
        }
        else
        {
            delq = qvec[i] - qvec[i-1];
        }

        integralsum = integralsum + delq*kernel_pk(s,p,qvec[i],a,m,eps)*dsol(i);


    }

    comp one = {1.0,0.0};
    result = -one*Gs - one*integralsum;
}

void interpolator_ds_integraleq_momrep_2eps(     Eigen::VectorXcd &dsol,
                                            vector<comp> &qvec,
                                            comp s,
                                            comp p,
                                            comp k,
                                            double a,
                                            double m,
                                            double eps,
                                            double eps_for_m2k,
                                            comp &result  )
{
    comp ii = {0.0,1.0};

    comp Gs = GS_pk(s,p,k,m,eps);
    comp delq;
    comp integralsum = {0.0,0.0};

    for(int i=0;i<qvec.size();++i)
    {
        if(i==0)
        {
            delq = qvec[1] - qvec[0];
        }
        else
        {
            delq = qvec[i] - qvec[i-1];
        }

        integralsum = integralsum + delq*kernel_pk_2eps(s,p,qvec[i],a,m,eps,eps_for_m2k)*dsol(i);


    }

    comp one = {1.0,0.0};
    result = -one*Gs - one*integralsum;
}

void interpolator_ds_integraleq_momrep_2eps_with_weights(       Eigen::VectorXcd &dsol,
                                                                vector<comp> &qvec,
                                                                vector<comp> &weights,
                                                                comp s,
                                                                comp p,
                                                                comp k,
                                                                double a,
                                                                double m,
                                                                double eps,
                                                                double eps_for_m2k,
                                                                comp &result  )
{
    comp ii = {0.0,1.0};

    comp Gs = GS_pk(s,p,k,m,eps);
    comp delq;
    comp integralsum = {0.0,0.0};

    for(int i=0;i<qvec.size();++i)
    {
        if(i==0)
        {
            delq = qvec[1] - qvec[0];
        }
        else
        {
            delq = qvec[i] - qvec[i-1];
        }

        integralsum = integralsum + weights[i]*kernel_pk_2eps(s,p,qvec[i],a,m,eps,eps_for_m2k)*dsol(i);


    }

    comp one = {1.0,0.0};
    result = -one*Gs - one*integralsum;
}


void interpolator_ds_integraleq_momrep_2eps_with_weights_for_resonances(       Eigen::VectorXcd &dsol,
                                                                vector<comp> &qvec,
                                                                vector<comp> &weights,
                                                                comp s,
                                                                comp p,
                                                                comp k,
                                                                double a,
                                                                double m,
                                                                double eps,
                                                                double eps_for_m2k,
                                                                double theta, 
                                                                int tag_for_m2k,
                                                                comp &result  )
{
    comp ii = {0.0,1.0};

    comp Gs = GS_pk(s,p,k,m,eps);
    comp delq;
    comp integralsum = {0.0,0.0};

    for(int i=0;i<qvec.size();++i)
    {
        if(i==0)
        {
            delq = qvec[1] - qvec[0];
        }
        else
        {
            delq = qvec[i] - qvec[i-1];
        }


        if(i<=tag_for_m2k)
        integralsum = integralsum + weights[i]*kernel_pk_2eps_for_resonances(s,p,qvec[i],a,m,eps,eps_for_m2k,theta,1)*dsol(i);
        else 
        integralsum = integralsum + weights[i]*kernel_pk_2eps_for_resonances(s,p,qvec[i],a,m,eps,eps_for_m2k,theta,0)*dsol(i); //0 before
        

    }

    comp one = {1.0,0.0};
    result = -one*Gs - one*integralsum;
}

void interpolator_ds_integraleq_momrep_2eps_with_weights_for_resonances_sheet3(       Eigen::VectorXcd &dsol,
                                                                vector<comp> &qvec,
                                                                vector<comp> &weights,
                                                                comp s,
                                                                comp p,
                                                                comp k,
                                                                double a,
                                                                double m,
                                                                double eps,
                                                                double eps_for_m2k,
                                                                double theta, 
                                                                int tag_for_m2k_1,
                                                                int tag_for_m2k_2,
                                                                comp &result  )
{
    comp ii = {0.0,1.0};

    comp Gs = GS_pk(s,p,k,m,eps);
    comp delq;
    comp integralsum = {0.0,0.0};

    for(int i=0;i<qvec.size();++i)
    {
        if(i==0)
        {
            delq = qvec[1] - qvec[0];
        }
        else
        {
            delq = qvec[i] - qvec[i-1];
        }


        if(i>=tag_for_m2k_1 && i<=tag_for_m2k_2)
        integralsum = integralsum + weights[i]*kernel_pk_2eps_for_resonances(s,p,qvec[i],a,m,eps,eps_for_m2k,theta,1)*dsol(i);
        else 
        integralsum = integralsum + weights[i]*kernel_pk_2eps_for_resonances(s,p,qvec[i],a,m,eps,eps_for_m2k,theta,0)*dsol(i); //0 before
        

    }

    comp one = {1.0,0.0};
    result = -one*Gs - one*integralsum;
}



void interpolator_ds_integraleq_momrep_2eps_with_weights_rotated_m2k_for_resonances(       Eigen::VectorXcd &dsol,
                                                                vector<comp> &qvec,
                                                                vector<comp> &weights,
                                                                comp s,
                                                                comp p,
                                                                comp k,
                                                                double a,
                                                                double m,
                                                                double eps,
                                                                double eps_for_m2k,
                                                                double theta, 
                                                                comp &result  )
{
    comp ii = {0.0,1.0};

    comp Gs = GS_pk(s,p,k,m,eps);
    comp delq;
    comp integralsum = {0.0,0.0};

    for(int i=0;i<qvec.size();++i)
    {
        if(i==0)
        {
            delq = qvec[1] - qvec[0];
        }
        else
        {
            delq = qvec[i] - qvec[i-1];
        }


        integralsum = integralsum + weights[i]*kernel_pk_2eps_rotated_m2k_for_resonances(s,p,qvec[i],a,m,eps,eps_for_m2k,theta)*dsol(i);
        

    }

    comp one = {1.0,0.0};
    result = -one*Gs - one*integralsum;
}


void interpolator_ds_integraleq_momrep_2eps_with_weights_2mysqrt(       Eigen::VectorXcd &dsol,
                                                                vector<comp> &qvec,
                                                                vector<comp> &weights,
                                                                comp s,
                                                                comp p,
                                                                comp k,
                                                                double a,
                                                                double m,
                                                                double eps,
                                                                double eps_for_m2k,
                                                                comp &result  )
{
    comp ii = {0.0,1.0};

    comp Gs = GS_pk(s,p,k,m,eps);
    comp delq;
    comp integralsum = {0.0,0.0};

    for(int i=0;i<qvec.size();++i)
    {
        if(i==0)
        {
            delq = qvec[1] - qvec[0];
        }
        else
        {
            delq = qvec[i] - qvec[i-1];
        }

        if(real(s)<9.0*m*m)
        integralsum = integralsum + weights[i]*kernel_pk_2eps(s,p,qvec[i],a,m,eps,eps_for_m2k)*dsol(i);
        else 
        integralsum = integralsum + weights[i]*kernel_pk_2eps_1(s,p,qvec[i],a,m,eps,eps_for_m2k)*dsol(i);

    }

    comp one = {1.0,0.0};
    result = -one*Gs - one*integralsum;
}



void interpolator_ds_integraleq_momrep_2eps_withtags(     Eigen::VectorXcd &dsol,
                                            vector<comp> &qvec,
                                            comp s,
                                            comp p,
                                            comp k,
                                            double a,
                                            double m,
                                            double eps,
                                            double eps_for_m2k,
                                            comp &result,
                                            int tag1,
                                            int tag2  )
{
    comp ii = {0.0,1.0};

    comp Gs = GS_pk(s,p,k,m,eps);
    comp delq;
    comp integralsum = {0.0,0.0};

    for(int i=0;i<qvec.size();++i)
    {
        if(i==0)
        {
            delq = qvec[1] - qvec[0];
        }
        else
        {
            delq = qvec[i] - qvec[i-1];
        }

        int index1 = i;
        int index2 = i;

        if(tag1==0 && tag2==0)
        {
            integralsum = integralsum + delq*kernel_pk_2eps(s,p,qvec[i],a,m,eps,eps_for_m2k)*dsol(i);
        }
        else 
        {
            integralsum = integralsum + delq*kernel_pk_2eps_withtags(s,p,qvec[i],a,m,eps,eps_for_m2k,index1,index2,tag1,tag2)*dsol(i);
        }

    }

    comp one = {1.0,0.0};
    result = -one*Gs - one*integralsum;
}

void interpolator_ds_integraleq_momrep_2eps_withtags_with_weights(     Eigen::VectorXcd &dsol,
                                            vector<comp> &qvec,
                                            vector<comp> &weights,
                                            comp s,
                                            comp p,
                                            comp k,
                                            double a,
                                            double m,
                                            double eps,
                                            double eps_for_m2k,
                                            comp &result,
                                            int tag1,
                                            int tag2  )
{
    comp ii = {0.0,1.0};

    comp Gs = GS_pk(s,p,k,m,eps);
    comp delq;
    comp integralsum = {0.0,0.0};

    for(int i=0;i<qvec.size();++i)
    {
        if(i==0)
        {
            delq = qvec[1] - qvec[0];
        }
        else
        {
            delq = qvec[i] - qvec[i-1];
        }

        int index1 = i;
        int index2 = i;

        if(tag1==0 && tag2==0)
        {
            integralsum = integralsum + weights[i]*kernel_pk_2eps(s,p,qvec[i],a,m,eps,eps_for_m2k)*dsol(i);
        }
        else 
        {
            integralsum = integralsum + weights[i]*kernel_pk_2eps_withtags(s,p,qvec[i],a,m,eps,eps_for_m2k,index1,index2,tag1,tag2)*dsol(i);
        }

    }

    comp one = {1.0,0.0};
    result = -one*Gs - one*integralsum;
}

void interpolator_ds_integraleq_momrep_2eps_withtags_with_weights_usingGvec(     Eigen::VectorXcd &dsol,
                                            vector<comp> &qvec,
                                            vector<comp> &weights,
                                            Eigen::VectorXcd &Gvec,
                                            comp s,
                                            comp p,
                                            comp k,
                                            double a,
                                            double m,
                                            double eps,
                                            double eps_for_m2k,
                                            comp &result,
                                            int tag1,
                                            int tag2  )
{
    comp ii = {0.0,1.0};

    comp Gs = GS_pk(s,p,k,m,eps);
    comp delq;
    comp integralsum = {0.0,0.0};

    for(int i=0;i<qvec.size();++i)
    {
        if(i==0)
        {
            delq = qvec[1] - qvec[0];
        }
        else
        {
            delq = qvec[i] - qvec[i-1];
        }

        int index1 = i;
        int index2 = i;

        /*if(tag1==0 && tag2==0)
        {
            integralsum = integralsum + weights[i]*kernel_pk_2eps(s,p,qvec[i],a,m,eps,eps_for_m2k)*dsol(i);
        }
        else 
        {
            integralsum = integralsum + weights[i]*kernel_pk_2eps_withtags(s,p,qvec[i],a,m,eps,eps_for_m2k,index1,index2,tag1,tag2)*dsol(i);
        }*/

        integralsum = integralsum + weights[i]*kernel_pk_2eps_using_Gvec(Gvec,i,s,p,qvec[i],a,m,eps,eps_for_m2k,index1,index2,tag1,tag2)*dsol(i);
        

    }

    comp one = {1.0,0.0};
    result = -one*Gs - one*integralsum;
}


void interpolator_gamma_integraleq_momrep(     Eigen::VectorXcd &dsol,
                                            vector<comp> &qvec,
                                            comp s,
                                            comp p,
                                            comp k,
                                            double a,
                                            double m,
                                            double eps,
                                            comp &result  )
{
    comp ii = {0.0,1.0};

    comp Gs = GS_pk(s,p,k,m,eps);
    comp delq;
    comp integralsum = {0.0,0.0};

    for(int i=0;i<qvec.size();++i)
    {
        if(i==0)
        {
            delq = qvec[0];
        }
        else
        {
            delq = qvec[i] - qvec[i-1];
        }

        integralsum = integralsum + delq*kernel_pk(s,p,qvec[i],a,m,eps)*dsol(i);


    }

    comp one = {1.0,0.0};
    result = - one*integralsum;
}

#endif