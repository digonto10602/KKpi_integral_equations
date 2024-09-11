#ifndef INVARIANT_VECTOR_MAKER_H
#define INVARIANT_VECTOR_MAKER_H

#include<bits/stdc++.h>
#include<Eigen/Dense>
#include "fastGL_based_functions.h"


using namespace std;

typedef complex<double> comp;


void mom_vector_maker_linear_1(   vector<comp> &qvec,
                                    comp startingpoint,
                                    comp delsig,
                                    double points,
                                    int a               )
{
    comp ii = {0.0,1.0};
    if(a==1)
    {
        for(int i=1;i<=(int)points;++i)
        {
            comp tempsig = startingpoint + ((comp)i)*delsig;
            qvec.push_back(tempsig);
        }
    }
    else if(a==2)
    {
        for(int i=1;i<=(int)points;++i)
        {
            comp tempsig = startingpoint + ii*((comp)i)*delsig;
            qvec.push_back(tempsig);
        }
    }

}

//non-uniform linear momentum vector maker 
void mom_vector_maker_linear_2(   vector<comp> &qvec,
                                    comp kmin,
                                    comp kmax, 
                                    double points,
                                    int a,
                                    int firstpart,
                                    int secondpart               )
{
    comp ii = {0.0,1.0};
    if(a==1)
    {
        int newpoints1 = firstpart*points/5;
        int newpoints2 = secondpart*points/5;

        double halfpoint = 0.002;//abs(kmax + kmin)/8.0;
        double delk1 = (double) abs(halfpoint - kmin)/newpoints1;
        cout<<"delk1 = "<<delk1<<endl;
        double delk2 = (double) abs(kmax - halfpoint)/newpoints2;
        cout<<"delk2 = "<<delk2<<endl;
        for(int i=1;i<(int)newpoints1+1;++i)
        {
            comp tempsig = kmin + ((comp)i)*delk1;
            qvec.push_back(tempsig);
        }
        comp startingpoint = qvec[qvec.size() - 1];
        for(int i=1;i<(int)newpoints2+1;++i)
        {
            comp tempsig = startingpoint + ((comp)i)*delk2;
            qvec.push_back(tempsig);
        }
    }
    else if(a==2)
    {
        int newpoints1 = 3*points/5;
        int newpoints2 = 2*points/5;

        double halfpoint = abs(kmax + kmin)/10.0;
        double delk1 = (double) abs(halfpoint - kmin)/newpoints1;
        double delk2 = (double) abs(kmax - halfpoint)/newpoints2;
        for(int i=1;i<(int)newpoints1+1;++i)
        {
            comp tempsig = kmin + ii*((comp)i)*delk1;
            qvec.push_back(tempsig);
        }
        comp startingpoint = qvec[qvec.size() - 1];
        for(int i=1;i<(int)newpoints2+1;++i)
        {
            comp tempsig = startingpoint + ii*((comp)i)*delk2;
            qvec.push_back(tempsig);
        }
    }

    cout<<"qvec generated with size = "<<qvec.size()<<endl;

}

//The linear_vector_maker_3 rescale the kmax - kmin between
//[0 : 1], then for all the chosen points selects the points
//following f(x) = x^p, where p is some power, p=4 for BS1,BS2
//and p=7 for BS3. We then scale back to the original kmin->kmax 

void mom_vector_maker_linear_3( vector<comp> &qvec,
                                comp kmax,
                                double N,
                                double p   )
{
    //comp delk = abs(kmax - kmin)/N;
    double deluni = 1.0/N;
    for(int i=1;i<N+1;++i)
    {
        //comp k = kmin + (comp)i*delk;
        double uni = 0.0 + i*deluni;
        double rescale = pow(uni,p);
        comp rerescale = (comp)rescale*kmax;
        qvec.push_back(rerescale);
        //cout<<k<<'\t'<<k/abs(kmax)<<'\t'<<uni<<'\t'<<uni*kmax<<endl;
    }
}

void mom_vector_maker_1(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double points1    )
{
    comp ii = {0.0,1.0};
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    comp delq = {0.0,0.0};

    if(real(rad)>0.0)
    {
        /*  this is the initial straight line along the negative imaginary axis */
        comp startingpoint = min;
        delq = abs((abs(imag(qpls)) + abs(diffqplsqmns) + eps)-min)/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,-delq,points1,2);
        
        /*  this is the contour going right in the imaginary axis   */
        startingpoint = qvec[qvec.size()-1];
        delq = abs((abs(real(qc1)) + eps) - real(startingpoint))/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,delq,points1,1);

        /*  This contour now goes up in the imaginary axis                       */
        startingpoint = qvec[qvec.size()-1];
        delq = abs(imag(startingpoint))/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,delq,points1,2);

        /*  This contour now goes right in the real axis                       */
        startingpoint = qvec[qvec.size()-1];
        delq = abs(real(max) - real(startingpoint))/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,delq,points1,1);

        
        /* print the size of the created sigma_vector */
        cout<<"sigvec created with size = "<<qvec.size()<<endl;

    }
    else
    {
        double length = real(max - min);
        comp delsig = length/(4.0*points1) ;

        int totpoints = 4.0*points1;
        for(int i=1;i<=totpoints;++i)
        {
            comp qmom = min + ii*eps + ((comp)i)*delsig;

            qvec.push_back(qmom);
        }
    }

}

void mom_vector_maker_opposite_1(   vector<comp> &qvec,   
                                    comp s,
                                    comp min,
                                    comp max,
                                    double a,
                                    double m,
                                    double eps,
                                    double points1    )
{
    comp ii = {0.0,1.0};
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qmns = q_minus(s,a,m,eps);
    comp qc1 = q_c1(s,a,m,eps);
    comp delq = {0.0,0.0};

    if(real(rad)>0.0)
    {
        /*  this is the initial straight line along the negative imaginary axis */
        comp startingpoint = min;
        delq = abs((abs(imag(qmns)) + 2.0*eps)-min)/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,delq,points1,2);
        
        /*  this is the contour going right in the imaginary axis   */
        startingpoint = qvec[qvec.size()-1];
        delq = abs((abs(real(qc1)) + eps) - real(startingpoint))/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,delq,points1,1);

        /*  This contour now goes up in the imaginary axis                       */
        startingpoint = qvec[qvec.size()-1];
        delq = abs(imag(startingpoint))/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,-delq,points1,2);

        /*  This contour now goes right in the real axis                       */
        startingpoint = qvec[qvec.size()-1];
        delq = abs(real(max) - real(startingpoint))/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,delq,points1,1);

        
        /* print the size of the created sigma_vector */
        cout<<"sigvec created with size = "<<qvec.size()<<endl;

    }
    else
    {
        double length = real(max - min);
        comp delsig = length/(4.0*points1) ;

        int totpoints = 4.0*points1;
        for(int i=1;i<=totpoints;++i)
        {
            comp qmom = min + ii*eps + ((comp)i)*delsig;

            qvec.push_back(qmom);
        }
    }

}


void xy_coordinates_m(    double x1,
                        double y1,
                        double x2,
                        double y2,
                        double r,
                        double &x4,
                        double &y4   )
{
    double x3,y3;
    x3 = (x1 + x2)/2.0;
    y3 = (y1 + y2)/2.0;

    double m1 = (y2 - y1)/(x2 - x1);
    double m2 = -1.0/m1;

    if(m1==0.0||abs(m1)<1.0e-16)
    {
        m2 = INT_MAX;
    }
    
    
    x4 = x3 - r/sqrt(1.0 + m2*m2);
    y4 = y3 - m2*(x3 - x4);

    //cout<<"for m: m1:"<<m1<<'\t'<<" m2:"<<m2<<endl;
    //cout<<"       x3:"<<x3<<'\t'<<" y3:"<<y3<<endl;

}

void xy_coordinates_p(    double x1,
                        double y1,
                        double x2,
                        double y2,
                        double r,
                        double &x4,
                        double &y4   )
{
    double x3,y3;
    x3 = (x1 + x2)/2.0;
    y3 = (y1 + y2)/2.0;

    double m1 = (y2 - y1)/(x2 - x1);
    double m2 = -1.0/m1;

    if(m1==0.0||abs(m1)<1.0e-8)
    {
        m2 = INT_MAX;
    }
    
    
    x4 = x3 + r/sqrt(1.0 + m2*m2);
    y4 = y3 - m2*(x3 - x4);

    //cout<<"for p: m1:"<<m1<<'\t'<<" m2:"<<m2<<endl;
    //cout<<"       x3:"<<x3<<'\t'<<" y3:"<<y3<<endl;

}

int m1_checker( double x1,
                double y1,
                double x2,
                double y2   )
{
    double m1 = (y2 - y1)/(x2 - x1);

    if(m1>=0.0) return 1;
    else return -1;
}

void mom_vector_maker_along_pminus_xmtoxp(  vector<comp> &qvec,
                                            comp s,
                                            double a,
                                            double m,
                                            double eps,
                                            double eps_for_m2k,
                                            double r,
                                            double points   )
{
    comp ii = {0.0,1.0};
    double xinitial = -1.0;
    double xfinal = 1.0;
    comp qc1 = q_c1(s,a,m,eps);
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp x1 = x_for_pk(s,qc1,q,m,eps);
    xinitial = real(x1);
    double delx1 = abs(xinitial-xfinal)/(points + 1.0);
    double delx2 = abs(xinitial-xfinal)/(points);

    vector<comp> tempvec;
    
    for(int i=0;i<(points+1);++i)
    //for(double x=xinitial;x<=xfinal;x=x+delx1)
    {
        double x = xinitial + i*delx1;
        comp pcut = pcut_minus_comp(x,s,q,m,eps);

        tempvec.push_back(pcut);
    }

    for(int i=1;i<tempvec.size();++i)
    {
        double x4,y4;
        double x1 = (double)real(tempvec[i-1]);
        double y1 = (double)imag(tempvec[i-1]);
        double x2 = (double)real(tempvec[i]);
        double y2 = (double)imag(tempvec[i]);

        int m1val = m1_checker(x1,y1,x2,y2);
        if(m1val==1)
        {
            xy_coordinates_m(x1,y1,x2,y2,r,x4,y4);
        }
        else
        {
            xy_coordinates_p(x1,y1,x2,y2,r,x4,y4);
        }
        comp qpoint = x4 + ii*y4;
        //cout<<"from up to down: x:"<<x4<<'\t'<<" y4:"<<y4<<endl;
        qvec.push_back(qpoint);

    }
}

void mom_vector_maker_along_pminus_xptoxm(  vector<comp> &qvec,
                                            comp s,
                                            double a,
                                            double m,
                                            double eps,
                                            double eps_for_m2k,
                                            double r,
                                            double points   )
{
    comp ii = {0.0,1.0};
    double xinitial = -1.0;
    double xfinal = 1.0;
    comp qc1 = q_c1(s,a,m,eps);
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp x1 = x_for_pk(s,qc1,q,m,eps);
    xinitial = real(x1);
    double delx1 = abs(xinitial-xfinal)/(points + 1.0);
    double delx2 = abs(xinitial-xfinal)/(points);

    vector<comp> tempvec;
    
    for(int i=0;i<points+1;++i)
    //for(double x=xfinal;x>xinitial;x=x-delx1)
    {
        double x = xfinal - i*delx1;
        comp pcut = pcut_minus_comp(x,s,q,m,eps);

        tempvec.push_back(pcut);
    }

    for(int i=1;i<tempvec.size();++i)
    {
        double x4,y4;
        double x1 = (double)real(tempvec[i-1]);
        double y1 = (double)imag(tempvec[i-1]);
        double x2 = (double)real(tempvec[i]);
        double y2 = (double)imag(tempvec[i]);

        int m1val = m1_checker(x1,y1,x2,y2);
        if(m1val==1)
        {
            xy_coordinates_p(x1,y1,x2,y2,r,x4,y4);
        }
        else
        {
            xy_coordinates_m(x1,y1,x2,y2,r,x4,y4);
        }

        comp qpoint = x4 + ii*y4;

        //cout<<"from down to up: x4:"<<x4<<'\t'<<" y4:"<<y4<<endl;
        //cout<<qpoint<<endl;
        qvec.push_back(qpoint);

    }
}

void half_circle_drawer_momrep(             vector<comp> &qvec,
                                            comp s,
                                            double a,
                                            double m,
                                            double eps,
                                            double eps_for_m2k,
                                            double r,
                                            double points     )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    double thetainitial = 0.0+0.0000001;
    double thetafinal = pi-0.0000001;
    double deltheta = abs(thetainitial - thetafinal)/(points);

    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp pcutm = pcut_minus(1.0,s,q,m,eps);

    

    for(int i=1;i<points+1;++i)
    //for(double theta=thetainitial+deltheta;theta<thetafinal;theta=theta+deltheta)
    {
        double theta = thetainitial + i*deltheta;
        comp qpoint = pcutm + r*exp(ii*theta);
        qvec.push_back(qpoint);
    }
}


void half_circle_drawer_momrep_1(           vector<comp> &qvec,
                                            comp s,
                                            double a,
                                            double m,
                                            double eps,
                                            double eps_for_m2k,
                                            double r,
                                            double points     )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    

    comp sigb = sigmab(a,m);
    comp qc1 = q_c1(s,a,m,eps);    
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp pcutm = pcut_minus(1.0,s,q,m,eps);


    double xinitial = -1.0;
    double xfinal = 1.0;
    comp x1angle = x_for_pk(s,qc1,q,m,eps);
    xinitial = real(x1angle);
    double delx1 = abs(xinitial-xfinal)/(points + 1.0);
    double delx2 = abs(xinitial-xfinal)/(points);

    comp firstval = qvec[qvec.size()-1];
    double firstangle = atan((imag(firstval) - imag(pcutm))/(real(firstval)-real(pcutm)));

    double thetainitial = abs(firstangle) + 0.00001;
    double thetafinal = pi + abs(firstangle) - 0.00001;
    double deltheta = abs(thetainitial - thetafinal)/(points);

    for(int i=1;i<points+1;++i)
    //for(double theta=thetainitial+deltheta;theta<thetafinal;theta=theta+deltheta)
    {
        double theta = thetainitial + i*deltheta;
        comp qpoint = pcutm + r*exp(ii*theta);
        qvec.push_back(qpoint);
    }
}



void mom_vector_maker_2(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points1,
                            double r    )
{
    comp ii = {0.0,1.0};
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    comp delq = {0.0,0.0};

    if(real(rad)>0.05)
    {
        /*  this is the initial straight line along the positive axis */
        comp startingpoint = min;
        delq = abs((abs(real(qc1)) - 1.5*r)-min)/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,+delq,points1,1);
        /*  this is the contour going along pcut down   */
        startingpoint = qvec[qvec.size()-1];
        //delq = abs((abs(real(qc1)) + eps) - real(startingpoint))/points1;
        //mom_vector_maker_linear_1(qvec,startingpoint,delq,points1,1);
        mom_vector_maker_along_pminus_xmtoxp(qvec,s,a,m,eps,eps_for_m2k,r,points1);
        
        
        
        /*  circle around the p+                */
        half_circle_drawer_momrep_1(qvec,s,a,m,eps,eps_for_m2k,r,points1/2.0);
        /*  This contour now goes along the pcut up             */
        startingpoint = qvec[qvec.size()-1];
        delq = abs(imag(startingpoint))/points1;
        mom_vector_maker_along_pminus_xptoxm(qvec,s,a,m,eps,eps_for_m2k,r,points1);
        
        /*  This contour now goes right in the real axis                       */
        startingpoint = qvec[qvec.size()-1];
        double retemp = real(startingpoint);
        double imtemp = imag(startingpoint);
        startingpoint = retemp;
        delq = abs(real(max) - real(startingpoint))/points1;
        mom_vector_maker_linear_1(qvec,startingpoint,delq,points1/2.0,1);
        //cout<<"qvec after fifth part:"<<qvec.size()-tmp1<<endl;
        
        
        /* print the size of the created sigma_vector */
        cout<<"qvec created with size = "<<qvec.size()<<endl;

    }
    else if(real(rad)>0.0 && real(rad)<0.05)
    {
        qvec.push_back(min);
        mom_vector_maker_along_pminus_xptoxm(qvec,s,a,m,eps,eps_for_m2k,r,2.0*points1-1.0);
        //cout<<"here"<<endl;
        /*  This contour now goes right in the real axis                       */
        comp startingpoint = qvec[qvec.size()-1];
        double retemp = real(startingpoint);
        double imtemp = imag(startingpoint);
        startingpoint = retemp;
        delq = abs(real(max) - real(startingpoint))/(2.0*points1);
        mom_vector_maker_linear_1(qvec,startingpoint,delq,2.0*points1,1);

        
        /* print the size of the created sigma_vector */
        cout<<"qvec created with size = "<<qvec.size()<<endl;

    }
    else
    {
        double length = real(max - min);
        comp delsig = length/(4.0*points1) ;

        int totpoints = 4.0*points1;
        for(int i=0;i<totpoints;++i)
        {
            comp qmom = min + ii*eps + ((comp)i)*delsig;

            qvec.push_back(qmom);
        }
        cout<<"qvec created with size = "<<qvec.size()<<endl;
    }

}


void check_which_pcut(  comp s,
                        comp q,
                        double m,
                        double eps,
                        int &flag   )
{
    //int flag;

    double xinitial = -1.0; //these are the 'x' values of the pcut
    double xfinal = 1.0;
    double points = 250.0;
    double delx = abs(xinitial - xfinal)/points;

    double tempxval1, tempxval2; //these are x axis values for the two cuts
    double xval1, xval2;
    int startcount = 0;
    for(double x=xinitial;x<=xfinal;x=x+delx)
    {
        comp pcutp = pcut_plus_comp(x,s,q,m,eps);
        comp pcutm = pcut_minus_comp(x,s,q,m,eps);

        if(startcount==0)
        {
            xval1 = (double) real(pcutp);
            xval2 = (double) real(pcutm);
            startcount = 1;
        }

        tempxval1 = (double) real(pcutp);
        tempxval2 = (double) real(pcutm);

        if(tempxval1<xval1) xval1 = tempxval1;
        if(tempxval2<xval2) xval2 = tempxval2;
    }

    if(xval1<xval2) flag = 0;
    else flag = 1;

}

void mom_vector_maker_3(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    comp delq = {0.0,0.0};

    int flag; 
    check_which_pcut(s,q,m,eps,flag);
    if(flag==0)
    {
        double xinitial = -1.0;
        double xfinal = 1.0;
        double xpoints = 250;
        double delx = abs(xinitial-xfinal)/xpoints;
        double xval1,yval1;
        int somecount=0;

        for(double x=xinitial;x<=xfinal;x=x+delx)
        {
            comp pcutm = pcut_minus(x,s,q,m,eps);
            if(somecount==0)
            {
                xval1=(double) real(pcutm);
                yval1=(double) imag(pcutm);
                somecount=1;
            }
            double tempxval1 = (double) real(pcutm);
            double tempyval1 = (double) imag(pcutm);
            if(tempxval1<xval1) xval1 = tempxval1;
            if(tempyval1<yval1) yval1 = tempyval1;

        }

        if(xval1<=1.0e-8 && yval1<=1.0e-3)
        {
            comp startingpoint = min;
            comp delq = abs(startingpoint - r)/(points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint) - yval1 + r)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,1);

        }
        else if(xval1>1.0e-8 && yval1<=1.0e-3)
        {
            comp startingpoint = min;
            delq = abs(imag(startingpoint) - yval1 + r)/(points + points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points + points/10.0,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,1);
        }
        else if(yval1>1.0e-3)
        {
            comp startingpoint = min;
            double totpoints = 4.0*points + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);
        }
    }
    else 
    {
        double xinitial = -1.0;
        double xfinal = 1.0;
        double xpoints = 250;
        double delx = abs(xinitial-xfinal)/xpoints;
        double xval1,yval1;
        int somecount=0;

        for(double x=xinitial;x<=xfinal;x=x+delx)
        {
            comp pcutp = pcut_plus(x,s,q,m,eps); 
            if(somecount==0)
            {
                xval1=(double) real(pcutp);
                yval1=(double) imag(pcutp);
                somecount=1;
            }
            double tempxval1 = (double) real(pcutp);
            double tempyval1 = (double) imag(pcutp);
            if(tempxval1<xval1) xval1 = tempxval1;
            if(tempyval1<yval1) yval1 = tempyval1;

        }

        if(xval1<=1.0e-8 && yval1<=1.0e-3)
        {
            comp startingpoint = min;
            comp delq = abs(startingpoint - r)/(points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint) - yval1 + r)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,1);

        }
        else if(xval1>1.0e-8 && yval1<=1.0e-3)
        {
            comp startingpoint = min;
            delq = abs(imag(startingpoint) - yval1 + r)/(points + points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points + points/10.0,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max)/points;
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,points,1);
        }
        else if(yval1>1.0e-3)
        {
            comp startingpoint = min;
            double totpoints = 4.0*points + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);
        }
    }
}


void closest_pcuts_bottom(  comp s,
                            comp q,
                            double m,
                            double eps,
                            comp &pcutp_close,
                            comp &pcutm_close,
                            double &x_pcutp,
                            double &x_pcutm )
{
    comp ii = {0.0,1.0};
    
    double initialx = 0.0;
    double finalx = 1.0;
    double size = 25000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }

    double tempdiff = 1.0e+16;
    double tempxp = 0.0;
    double tempxm = 0.0;

    for(int i=0;i<pcutpvec.size();++i)
    {
        for(int j=0;j<pcutmvec.size();++j)
        {
            double diff = abs(pcutpvec[i] - pcutmvec[j]);

            //cout<<"xp = "<<xvec[i]<<'\t'<<"xm = "<<xvec[j]<<endl;

            if(diff<=tempdiff)
            {
                //cout<<"diff minimum then previous"<<endl;
                //cout<<"diff = "<<diff<<endl;
                //if(imag(pcutpvec[i])<0.0 && imag(pcutmvec[j])<0.0)
                {
                    tempdiff = diff;
                    tempxp = xvec[i];
                    tempxm = xvec[j];
                }
            }
            
        }
    }

    pcutp_close = pcut_plus(tempxp,s,q,m,eps);
    pcutm_close = pcut_minus(tempxm,s,q,m,eps);
    x_pcutp = tempxp;
    x_pcutm = tempxm;
    cout<<"diff found in pcuts = "<<abs(pcutp_close-pcutm_close)<<endl;
    
}

void slanted_line_maker(    vector<comp> &qvec,
                            comp startingpoint,
                            comp endingpoint,
                            double offset,
                            double points       )
{
    comp ii = {0.0,1.0};
    double x1 = real(startingpoint);
    double y1 = imag(startingpoint);
    double x2 = real(endingpoint);
    double y2 = imag(endingpoint);
    double m = (y2 - y1)/(x2 - x1);
    m = - abs(m);
    //cout<<"m:"<<m<<endl;
    double length = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
    double r = length + 4.0*offset;
    double theta = atan(m);

    double delr = r/points;

    for(int i=1;i<points+1;++i)
    {
        double rval = i*delr;

        comp qval = startingpoint + rval*(cos(theta) + ii*sin(theta));

        qvec.push_back(qval);
    }

    
}

void mom_vector_maker_4(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    comp delq = {0.0,0.0};

    comp pcutp_close;
    comp pcutm_close;
    double x_pcutp;
    double x_pcutm;
    closest_pcuts_bottom(s,q,m,eps,pcutp_close,pcutm_close,x_pcutp,x_pcutm);


    int flag; 
    check_which_pcut(s,q,m,eps,flag);
    if(flag==0)
    {
        double xinitial = -1.0;
        double xfinal = 1.0;
        double xpoints = 250;
        double delx = abs(xinitial-xfinal)/xpoints;
        double xval1,yval1;
        int somecount=0;



        for(double x=xinitial;x<=xfinal;x=x+delx)
        {
            comp pcutm = pcut_minus(x,s,q,m,eps);
            if(somecount==0)
            {
                xval1=(double) real(pcutm);
                yval1=(double) imag(pcutm);
                somecount=1;
            }
            double tempxval1 = (double) real(pcutm);
            double tempyval1 = (double) imag(pcutm);
            if(tempxval1<xval1) xval1 = tempxval1;
            if(tempyval1<yval1) yval1 = tempyval1;

        }



        if(xval1<=1.0e-8 && yval1<=1.0e-5)
        {
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/10.0);
            //startingpoint = startingpoint - ii*r;
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close + pcutm_close)/2.0;
            slanted_line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,r,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max - startingpoint)/(3.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,3.0*points/10.0,1);

        }
        else if(xval1>1.0e-8 && yval1<=1.0e-5)
        {
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/10.0);
            //startingpoint = startingpoint - ii*r;
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close + pcutm_close)/2.0;
            slanted_line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,r,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max - startingpoint)/(3.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,3.0*points/10.0,1);
        }
        else if(yval1>=0.0)
        {
            comp startingpoint = min;
            double totpoints = points;// + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);
        }
    }
    else 
    {
        double xinitial = -1.0;
        double xfinal = 1.0;
        double xpoints = 250;
        double delx = abs(xinitial-xfinal)/xpoints;
        double xval1,yval1;
        int somecount=0;

        for(double x=xinitial;x<=xfinal;x=x+delx)
        {
            comp pcutp = pcut_plus(x,s,q,m,eps); 
            if(somecount==0)
            {
                xval1=(double) real(pcutp);
                yval1=(double) imag(pcutp);
                somecount=1;
            }
            double tempxval1 = (double) real(pcutp);
            double tempyval1 = (double) imag(pcutp);
            if(tempxval1<xval1) xval1 = tempxval1;
            if(tempyval1<yval1) yval1 = tempyval1;

        }

        if(xval1<=1.0e-8 && yval1<=1.0e-5)
        {
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/10.0);
            //startingpoint = startingpoint - ii*r;
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close + pcutm_close)/2.0;
            slanted_line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,r,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max - startingpoint)/(3.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,3.0*points/10.0,1);

        }
        else if(xval1>1.0e-8 && yval1<=1.0e-5)
        {
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/10.0);
            //startingpoint = startingpoint - ii*r;
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close + pcutm_close)/2.0;
            slanted_line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,r,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(real(startingpoint) - rad - r)/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(imag(startingpoint))/(2.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,2.0*points/10.0,2);

            startingpoint = qvec[qvec.size()-1];
            delq = abs(max - startingpoint)/(3.0*points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,3.0*points/10.0,1);
        }
        else if(yval1>=0.0)
        {
            comp startingpoint = min;
            double totpoints = points;// + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);
        }
    }
}



void check_which_pcut_updown(   comp s,
                                comp q,
                                double m,
                                double eps,
                                int &flag   )
{
    //int flag;

    double xinitial = -1.0; //these are the 'x' values of the pcut
    double xfinal = 1.0;
    double points = 250.0;
    double delx = abs(xinitial - xfinal)/points;

    comp pcutp_ini = pcut_plus(xinitial,s,q,m,eps);
    comp pcutm_ini = pcut_minus(xinitial,s,q,m,eps);
    comp pcutp_fin = pcut_plus(xfinal,s,q,m,eps);
    comp pcutm_fin = pcut_minus(xfinal,s,q,m,eps);

    comp select_pcutp, select_pcutm;
    if(imag(pcutp_ini)>imag(pcutp_fin)) select_pcutp = pcutp_ini;
    else select_pcutp = pcutp_fin;

    if(imag(pcutm_ini)>imag(pcutm_fin)) select_pcutm = pcutm_ini;
    else select_pcutm = pcutm_fin;

    if(imag(select_pcutp)>imag(select_pcutm)) flag = 1;
    else flag = 0;

    //if flag = 0 then pcutm is above and pcutp is in the bottom
    //if flag = 1 then pcutp is above and pcutm is in the bottom


}

void line_maker(    vector<comp> &qvec, 
                    comp a,
                    comp b,
                    double points   )
{
    double zinitial = 0.0;
    double zfinal = 1.0;
    double delz = abs(zfinal-zinitial)/points;
    for(int i=1;i<points+1;++i)
    {
        double z = zinitial + i*delz;
        comp x_of_z = (b-a)*z + a;

        qvec.push_back(x_of_z);
    }
}

//this one is used for checking the anomalous behavious
//in Mphib for Im(s)>0 
void mom_vector_maker_6(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double r,
                            int &tag1,
                            int &tag2     )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    comp delq = {0.0,0.0};

    comp pcutp_close;
    comp pcutm_close;
    double x_pcutp;
    double x_pcutm;
    //closest_pcuts_bottom(s,q,m,eps,pcutp_close,pcutm_close,x_pcutp,x_pcutm);


    int flag; 
    if(real(s)<8.925)
    check_which_pcut(s,q,m,eps,flag);
    else
    check_which_pcut_updown(s,q,m,eps,flag);
    //flag = 0;
    if(flag==0)
    {
        double xval = 0;

        comp pcutm = pcut_minus(xval+0.4,s,q,m,eps);
        comp pcutp = pcut_plus(xval-0.45,s,q,m,eps);

        comp startingpoint = min;
        double thetaoffset = acos(-1.0);

        comp endingpoint_for_slanted_maker = pcutp;

        line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
        //slanted_line_maker_2(qvec,startingpoint,endingpoint_for_slanted_maker,0.0,points,0.0);//offset set to zero
        tag1 = qvec.size() - 1;
        
        startingpoint = qvec[qvec.size()-1];
        endingpoint_for_slanted_maker = pcutm;
        line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
        //slanted_line_maker_2(qvec,startingpoint,endingpoint_for_slanted_maker,0.0,points,thetaoffset);//offset set to zero
        tag2 = qvec.size() - 1;

        startingpoint = qvec[qvec.size()-1];
        endingpoint_for_slanted_maker = max;
        line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
        //slanted_line_maker_2(qvec,startingpoint,endingpoint_for_slanted_maker,0.0,points,0.0);//offset set to zero

    }
    else 
    {
        double xval = 0;

        comp pcutm = pcut_minus(xval+0.4,s,q,m,eps);
        comp pcutp = pcut_plus(xval-0.45,s,q,m,eps);

        comp startingpoint = min;
        double thetaoffset = acos(-1.0);

        comp endingpoint_for_slanted_maker = pcutp;
        line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
        //slanted_line_maker_2(qvec,startingpoint,endingpoint_for_slanted_maker,0.0,points,0.0);//offset set to zero
        tag1 = qvec.size() - 1;
        
        startingpoint = qvec[qvec.size()-1];
        endingpoint_for_slanted_maker = pcutm;
        line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
        //slanted_line_maker_2(qvec,startingpoint,endingpoint_for_slanted_maker,0.0,points,thetaoffset);//offset set to zero
        tag2 = qvec.size() - 1;

        startingpoint = qvec[qvec.size()-1];
        endingpoint_for_slanted_maker = max;
        line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
        //slanted_line_maker_2(qvec,startingpoint,endingpoint_for_slanted_maker,0.0,points,0.0);//offset set to zero

    }
}



void pcut_fixer(    vector<comp> &pcutp,
                    vector<comp> &pcutm  )
{
    int size = pcutp.size();
    for(int i=0;i<size;++i)
    {
        if(i==0)
        {
            comp checkval_pcutp = pcutp[i];
            comp checkprevval_pcutp = pcutp[i+1];
            comp checkval_pcutm = pcutm[i];

            double checkdiff_with_pcutp = abs(checkval_pcutp - checkprevval_pcutp);
            double checkdiff_with_pcutm = abs(checkval_pcutm - checkprevval_pcutp);

            if(checkdiff_with_pcutm<checkdiff_with_pcutp)
            {
                comp temp;
                temp = pcutp[i];
                pcutp[i] = pcutm[i];
                pcutm[i] = temp;
            }
            else continue;
        }
        comp checkval_pcutp = pcutp[i];
        comp checkprevval_pcutp = pcutp[i-1];
        comp checkval_pcutm = pcutm[i];

        double checkdiff_with_pcutp = abs(checkval_pcutp - checkprevval_pcutp);
        double checkdiff_with_pcutm = abs(checkval_pcutm - checkprevval_pcutp);

        if(checkdiff_with_pcutm<checkdiff_with_pcutp)
        {
            comp temp;
            temp = pcutp[i];
            pcutp[i] = pcutm[i];
            pcutm[i] = temp;
        }
        else continue;
    }
}


//this function takes the pcutpvec, pcutmvec
//and finds the closest gap between the two 
//cuts and returns the closest pcutp1 and pcutm1,
//along with the second pair of these values 
void closest_pcuts(             vector<double> xvec,
                                vector<comp> pcutpvec,
                                vector<comp> pcutmvec,
                                comp s,
                                comp q,
                                double m,
                                double eps,
                                int &first_i,            //i is for pcutp index
                                int &first_j,            //j is for pcutm index
                                int &second_i,
                                int &second_j      )
{
    comp ii = {0.0,1.0};
    
    

    double tempdiff = 1.0e+16;
    double tempxp = 0.0;
    double tempxm = 0.0;
    
    for(int i=0;i<pcutpvec.size();++i)
    {
        for(int j=0;j<pcutmvec.size();++j)
        {
            double diff = abs(pcutpvec[i] - pcutmvec[j]);

            //cout<<"xp = "<<xvec[i]<<'\t'<<"xm = "<<xvec[j]<<endl;

            if(diff<=tempdiff)
            {
                //cout<<"diff minimum then previous"<<endl;
                //cout<<"diff = "<<diff<<endl;
                //if(imag(pcutpvec[i])<0.0 && imag(pcutmvec[j])<0.0)
                {
                    tempdiff = diff;
                    tempxp = xvec[i];
                    tempxm = xvec[j];
                    first_i = i;
                    first_j = j;
                }
            }
            
        }
    }

    comp pcutp_close1, pcutm_close1;
    pcutp_close1 = pcutpvec[first_i];//pcut_plus(tempxp,s,q,m,eps);
    pcutm_close1 = pcutmvec[first_j];//pcut_minus(tempxm,s,q,m,eps);
    
    //cout<<"first_i = "<<first_i<<endl;
    //cout<<"first_j = "<<first_j<<endl;

    //cout<<"pcutp_close1 = "<<pcutp_close1<<endl;
    //cout<<"pcutm_close1 = "<<pcutm_close1<<endl;
    
    tempdiff = 1.0e+16;

    for(int i=0;i<pcutpvec.size();++i)
    {
        for(int j=0;j<pcutmvec.size();++j)
        {
            double diff = abs(pcutpvec[i] - pcutmvec[j]);

            //cout<<"xp = "<<xvec[i]<<'\t'<<"xm = "<<xvec[j]<<endl;
            if(i==first_i && j==first_j)
            {
                
            }
            else
            {
                if(diff<=tempdiff)
                {
                    //cout<<"diff minimum then previous"<<endl;
                    //cout<<"diff = "<<diff<<endl;
                    //if(imag(pcutpvec[i])<0.0 && imag(pcutmvec[j])<0.0)
                
                    tempdiff = diff;
                    tempxp = xvec[i];
                    tempxm = xvec[j];
                    second_i = i;
                    second_j = j;
                
                }
            }
            
        }
    }

    comp pcutp_close2, pcutm_close2;
    pcutp_close2 = pcutpvec[second_i];//pcut_plus(tempxp,s,q,m,eps);
    pcutm_close2 = pcutmvec[second_j];//pcut_minus(tempxm,s,q,m,eps);

    //cout<<"second_i = "<<second_i<<endl;
    //cout<<"second_j = "<<second_j<<endl;

    //cout<<"pcutp_close2 = "<<pcutp_close2<<endl;
    //cout<<"pcutm_close2 = "<<pcutm_close2<<endl;
    
    cout<<"diff found in pcuts1 = "<<abs(pcutp_close1-pcutm_close1)<<endl;
    cout<<"diff found in pcuts2 = "<<abs(pcutp_close2-pcutm_close2)<<endl;
    
}



//this function tells which pcutvec is to the right 
//and which one is to the left. flag=0 means pcutpvec is 
//to the right and flag=1 means pcutmvec is to the right
void check_pcut_leftright(  vector<comp> pcutpvec,
                            vector<comp> pcutmvec,
                            int &flag   )
{
    //int flag;

    

    double tempxval1, tempxval2; //
    double xval1, xval2;
    int startcount = 0;
    //for(double x=xinitial;x<=xfinal;x=x+delx)
    for(int i=0;i<pcutpvec.size();++i)
    {
        comp pcutp = pcutpvec[i];//pcut_plus(x,s,q,m,eps);
        comp pcutm = pcutmvec[i];//pcut_minus(x,s,q,m,eps);

        if(startcount==0)
        {
            xval1 = (double) real(pcutp);
            xval2 = (double) real(pcutm);
            startcount = 1;
        }

        tempxval1 = (double) real(pcutp);
        tempxval2 = (double) real(pcutm);

        if(tempxval1>xval1) xval1 = tempxval1;
        if(tempxval2>xval2) xval2 = tempxval2;
    }

    if(xval1<xval2) flag = 1;
    else flag = 0;

}



//this function gives out the lowest y-values of 
//pcutpvec and pcutmvec

void lowest_yval_pcuts( vector<comp> pcutpvec,
                        vector<comp> pcutmvec,
                        double &yvalp,
                        double &yvalm,
                        int &lowest_p_index,
                        int &lowest_m_index       )
{
    double tempyvalp, tempyvalm;
    for(int i=0;i<pcutpvec.size();++i)
    {
        if(i==0)
        {
            yvalp = imag(pcutpvec[i]);
            yvalm = imag(pcutmvec[i]);
            lowest_p_index = i;
            lowest_m_index = i;
            //continue;
        }

        tempyvalp = imag(pcutpvec[i]);
        tempyvalm = imag(pcutmvec[i]);

        if(tempyvalp<yvalp)
        {
            yvalp = tempyvalp;
            lowest_p_index = i;
        }
        if(tempyvalm<yvalm)
        {
            yvalm = tempyvalm;
            lowest_m_index = i;
        }
    }
    //cout<<"yvalp = "<<yva
}

void highest_yval_pcuts( vector<comp> &pcutpvec,
                        vector<comp> &pcutmvec,
                        double &yvalp,
                        double &yvalm,
                        int &lowest_p_index,
                        int &lowest_m_index       )
{
    double tempyvalp, tempyvalm;
    for(int i=0;i<pcutpvec.size();++i)
    {
        if(i==0)
        {
            yvalp = imag(pcutpvec[i]);
            yvalm = imag(pcutmvec[i]);
            lowest_p_index = i;
            lowest_m_index = i;
            //continue;
        }

        tempyvalp = imag(pcutpvec[i]);
        tempyvalm = imag(pcutmvec[i]);

        if(tempyvalp>yvalp)
        {
            yvalp = tempyvalp;
            lowest_p_index = i;
        }
        if(tempyvalm>yvalm)
        {
            yvalm = tempyvalm;
            lowest_m_index = i;
        }
    }
    //cout<<"yvalp = "<<yva
}



void highest_xval_pcuts(    vector<comp> &pcutpvec,
                            vector<comp> &pcutmvec,
                            double &xvalp,
                            double &xvalm,
                            int &highest_p_index,
                            int &highest_m_index       )
{
    double tempxvalp, tempxvalm;
    for(int i=0;i<pcutpvec.size();++i)
    {
        if(i==0)
        {
            xvalp = real(pcutpvec[i]);
            xvalm = real(pcutmvec[i]);
            highest_p_index = i;
            highest_m_index = i;
            //continue;
        }

        tempxvalp = real(pcutpvec[i]);
        tempxvalm = real(pcutmvec[i]);

        if(tempxvalp>xvalp)
        {
            xvalp = tempxvalp;
            highest_p_index = i;
        }
        if(tempxvalm>xvalm)
        {
            xvalm = tempxvalm;
            highest_m_index = i;
        }
    }
    //cout<<"yvalp = "<<yva
}



int sign_checker(   double a    )
{
    if(a>=0.0) return 1;
    else return -1;
}

void pcutvec_realaxis_crossing_checker( vector<comp> pcutpvec,
                                        vector<comp> pcutmvec,
                                        int &flag   )
{
    double some_index1;
    double some_index2;
    double tempyvalp, tempyvalm, yvalp, yvalm;
    int somecountp = 0;
    int somecountm = 0;
    for(int i=0;i<pcutpvec.size();++i)
    {
        comp pcutp = pcutpvec[i];
        comp pcutm = pcutmvec[i];

        if(real(pcutp)>0.0)
        {
            if(somecountp==0)
            {
                yvalp = imag(pcutp);
                somecountp = 1;
                
            }

            tempyvalp = imag(pcutp);
            if(tempyvalp<yvalp) yvalp = tempyvalp;

            
        }

        if(real(pcutm)>0.0)
        {
            if(somecountm==0)
            {
                yvalp = imag(pcutm);
                somecountm = 1;
            }

            tempyvalm = imag(pcutm);
            if(tempyvalm<yvalm) yvalm = tempyvalm;

        }
    }

    if(yvalp*yvalm>0.0) flag = 0; //flag=0 means deform the contour
                                  //flag=1 means take a straight line through the real axis;
    else flag = 1;


}

void pcutvec_realaxis_crossing_checker_1(   vector<comp> pcutpvec,
                                            vector<comp> pcutmvec,
                                            int &flag   )
{
    double some_index1;
    double some_index2;
    double tempyvalp, tempyvalm, yvalp, yvalm;
    int somecountp = 0;
    int somecountm = 0;
    int prevpsign;
    int prevmsign;
    int signchanged_for_pcutp = 0;
    int signchanged_for_pcutm = 0;
    for(int i=0;i<pcutpvec.size();++i)
    {
        comp pcutp = pcutpvec[i];
        comp pcutm = pcutmvec[i];

        if(real(pcutp)>0.0)
        {
            int signp = sign_checker(imag(pcutp));

            if(somecountp==0)
            {
                prevpsign = signp;
                somecountp = 1;
                //continue;
                
            }
            else 
            {
                if(prevpsign==signp)
                {
                    signchanged_for_pcutp=0; 
                } 
                else
                {
                    signchanged_for_pcutp=1;
                    break;
                }
                prevpsign = signp;
            }
            //tempyvalp = imag(pcutp);
            //if(tempyvalp<yvalp) yvalp = tempyvalp;
            
            
        }

        if(real(pcutm)>0.0)
        {
            int signm = sign_checker(imag(pcutm));

            if(somecountm==0)
            {
                prevmsign = signm;
                somecountm = 1;
                //continue;
                
            }
            else 
            {

            //signchanged_for_pcut = 0 means sign didn't change
            //signchanged_for_pcut = 1 means sign did change 
                if(prevmsign==signm)
                {
                    signchanged_for_pcutm=0; 
                } 
                else
                {
                    signchanged_for_pcutm=1;
                    break;
                }
                prevmsign = signm;
            }

            
            //tempyvalp = imag(pcutp);
            //if(tempyvalp<yvalp) yvalp = tempyvalp;

        }
    }

    if(signchanged_for_pcutp==0 && signchanged_for_pcutm==0) flag = 1;  
    else if(signchanged_for_pcutp==1 && signchanged_for_pcutm==0) flag = 0;
    else if(signchanged_for_pcutp==0 && signchanged_for_pcutm==1) flag = 0;
    else if(signchanged_for_pcutp==1 && signchanged_for_pcutm==1) flag = 0;
    //flag=0 means deform the contour
    //flag=1 means take a straight line through the real axis;

}



//This is the same vector maker as maker_4 but 
//only difference is the which pcut is to the 
//right has been changed, along with a pcut fixer 
//has been used and slanted line makers have been 
//changed with just line_maker which is a better 
//line maker maker
void mom_vector_maker_41(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    double qvec_r = 0.01;

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=0 means deform the contour first backward and then all the other steps
    //problem_index_flag=1 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);
    
    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/4.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/4.0,1);

            startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/4.0);

            startingpoint = qvec[qvec.size()-1];
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[lowest_p_index] - ii*0.01;
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[lowest_m_index] - ii*0.01;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/4.0);

            startingpoint = qvec[qvec.size()-1];
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/4.0);
        }
        else
        {
            cout<<"deforming the contour downward first"<<endl;
            //comp startingpoint = min;
            //comp delq = abs(startingpoint - rad/2.0)/(points);
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,1);

            //startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            //comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
            }
            else
            {
                endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/2.0);

            startingpoint = qvec[qvec.size()-1];
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/2.0);
        }

    }
    else
    {
        cout<<"taking the straight line"<<endl;
        //cout<<"prob here with axis_flag="<<axis_flag<<endl;
        comp startingpoint = min;
        double totpoints = 4.0*points;// + points/10.0;
        delq = abs(real(min) - real(max))/(totpoints);
        mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);
    }



    
}


//this is the same as maker_41, only change that 
//has been included is that if the midpoint of the 
//difference of the cuts is lower than the lowest y
//value of the cut to the right, then the contour 
//goes straight to the other endpoint 'max'. Also
//if difference between the x axis and the lowest 
//yvalue of the right cut is less than 1.0e-2, we will
//still deform the contour downward.
void mom_vector_maker_42(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    double qvec_r = 0.01;

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=0 means deform the contour first backward and then all the other steps
    //problem_index_flag=1 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);
    
    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/4.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/4.0,1);

            startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/4.0);

            startingpoint = qvec[qvec.size()-1];
            double lowest_index_yval = 0.0;
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[lowest_p_index] - ii*0.01;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[lowest_m_index] - ii*0.01;
            }

            if(lowest_index_yval<=imag(startingpoint))
            {
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/4.0);

                startingpoint = qvec[qvec.size()-1];
                endingpoint_for_slanted_maker = max;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/4.0);
            }
            else
            {
                startingpoint = qvec[qvec.size()-1];
                endingpoint_for_slanted_maker = max;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/2.0);
            }
        }
        else
        {
            cout<<"deforming the contour downward first"<<endl;
            //comp startingpoint = min;
            //comp delq = abs(startingpoint - rad/2.0)/(points);
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,1);

            //startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            //comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
            }
            else
            {
                endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/2.0);

            startingpoint = qvec[qvec.size()-1];
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/2.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            lowest_index_yval = imag(pcutpvec[lowest_p_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
        }
        else
        {
            lowest_index_yval = imag(pcutmvec[lowest_m_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
        }


        if(abs(lowest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour downward first"<<endl;
            comp startingpoint = min;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/2.0);

            endingpoint_for_slanted_maker = max;
            startingpoint = qvec[qvec.size()-1];

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/2.0);


        }
        else 
        {
            cout<<"taking the straight line"<<endl;
            //cout<<"prob here with axis_flag="<<axis_flag<<endl;
            comp startingpoint = min;
            double totpoints = points;// + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);

        }
    }



    
}


//this is the same as maker_42, with added
//qvec that goes to right, to qmax 
//again 
void mom_vector_maker_43(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);
            
            cout<<"turningpoint 0 = "<<qvec[0]<<endl;
            startingpoint = qvec[qvec.size()-1];
            cout<<"turningpoint 1 = "<<startingpoint<<endl;
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            cout<<"turningpoint 2 = "<<startingpoint<<endl;
            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[lowest_p_index] - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[lowest_m_index] - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            if(lowest_index_yval<=imag(startingpoint))
            {
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 3 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0); 

                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 4 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = max;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
                cout<<"end point 5 = "<<qvec[qvec.size() - 1]<<endl;
            }
            else
            {
                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 3 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
                
                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 4 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = max;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
                cout<<"end point 5 = "<<qvec[qvec.size() - 1]<<endl;
            }
        }
        else
        {
            cout<<"deforming the contour downward first"<<endl;
            //comp startingpoint = min;
            //comp delq = abs(startingpoint - rad/2.0)/(points);
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,1);

            //startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            //comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/5.0);


            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            lowest_index_yval = imag(pcutpvec[lowest_p_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
            highest_index_xval = real(pcutpvec[highest_p_index]);
        }
        else
        {
            lowest_index_yval = imag(pcutmvec[lowest_m_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
            highest_index_xval = real(pcutpvec[highest_m_index]);
        }


        if(abs(lowest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour downward first"<<endl;
            comp startingpoint = min;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);


        }
        else 
        {
            cout<<"taking the straight line"<<endl;
            //cout<<"prob here with axis_flag="<<axis_flag<<endl;
            comp startingpoint = min;
            double totpoints = points;// + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);

        }
    }



    
}


//this is the same as maker_42, with added
//qvec that goes to right, to qmax 
//again, gives the weight vector back,
//piecewise linear line maker used 
void mom_vector_maker_43_with_weights(      vector<comp> &qvec, 
                                            vector<comp> &weights,
                                            comp s,
                                            comp min,
                                            comp max,
                                            double a,
                                            double m,
                                            double eps,
                                            double eps_for_m2k,
                                            double points,
                                            double qvec_r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    
    cout<<"axis_flag = "<<axis_flag<<endl;
    if(axis_flag==0)
    cout<<"cut to the right crosses real line"<<endl;
    else
    cout<<"cut to the right doesnot cross real line"<<endl;

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);
            line_maker_with_weights(qvec,weights,startingpoint,-rad/2.0,points/5.0);
            
            cout<<"turningpoint 0 = "<<qvec[0]<<endl;
            startingpoint = qvec[qvec.size()-1];
            cout<<"turningpoint 1 = "<<startingpoint<<endl;
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            cout<<"turningpoint 2 = "<<startingpoint<<endl;
            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[lowest_p_index] - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[lowest_m_index] - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            if(lowest_index_yval<=imag(startingpoint))
            {
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 3 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0); 

                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 4 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = max;
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);
                cout<<"end point 5 = "<<qvec[qvec.size() - 1]<<endl;
            }
            else
            {
                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 3 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
                
                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 4 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = max;
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
                cout<<"end point 5 = "<<qvec[qvec.size() - 1]<<endl;
            }
        }
        else
        {
            cout<<"deforming the contour downward first"<<endl;
            //comp startingpoint = min;
            //comp delq = abs(startingpoint - rad/2.0)/(points);
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,1);

            //startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            //comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.01*abs(imag(pcutpvec[lowest_p_index]));
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.01*abs(imag(pcutpvec[lowest_p_index]));
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,2.0*points/5.0);


            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            lowest_index_yval = imag(pcutpvec[lowest_p_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.0001;
            highest_index_xval = real(pcutpvec[highest_p_index]);
        }
        else
        {
            lowest_index_yval = imag(pcutmvec[lowest_m_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.0001;
            highest_index_xval = real(pcutpvec[highest_m_index]);
        }


        if(abs(lowest_index_yval)<1.0e-3)
        {
            cout<<"still deforming the contour downward first"<<endl;
            comp startingpoint = min;

            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,2.0*points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;

            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);


        }
        else 
        {
            cout<<"taking the straight line"<<endl;
            //cout<<"prob here with axis_flag="<<axis_flag<<endl;
            comp startingpoint = min;
            comp endingpoint = max;
            double totpoints = points;// + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            //mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint,points);
            

        }
    }



    
}



void mom_vector_maker_43_2( double somepoint,   
                            vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - somepoint*rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);
            
            cout<<"turningpoint 0 = "<<qvec[0]<<endl;
            startingpoint = qvec[qvec.size()-1];
            cout<<"turningpoint 1 = "<<startingpoint<<endl;
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            cout<<"turningpoint 2 = "<<startingpoint<<endl;
            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[lowest_p_index] - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[lowest_m_index] - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            if(lowest_index_yval<=imag(startingpoint))
            {
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 3 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0); 

                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 4 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = max;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
                cout<<"end point 5 = "<<qvec[qvec.size() - 1]<<endl;
            }
            else
            {
                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 3 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
                
                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 4 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = max;
                line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
                cout<<"end point 5 = "<<qvec[qvec.size() - 1]<<endl;
            }
        }
        else
        {
            cout<<"deforming the contour downward first"<<endl;
            //comp startingpoint = min;
            //comp delq = abs(startingpoint - rad/2.0)/(points);
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,1);

            //startingpoint = qvec[qvec.size()-1];
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            //comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/5.0);


            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            lowest_index_yval = imag(pcutpvec[lowest_p_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
            highest_index_xval = real(pcutpvec[highest_p_index]);
        }
        else
        {
            lowest_index_yval = imag(pcutmvec[lowest_m_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
            highest_index_xval = real(pcutpvec[highest_m_index]);
        }


        if(abs(lowest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour downward first"<<endl;
            comp startingpoint = min;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);


        }
        else 
        {
            cout<<"taking the straight line"<<endl;
            //cout<<"prob here with axis_flag="<<axis_flag<<endl;
            comp startingpoint = min;
            double totpoints = points;// + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);

        }
    }



    
}


//This one has tags for taking the contour 
//under the second sheet of OPE
void mom_vector_maker_44(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[7500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[7500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[7500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[7500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"taking the straight line"<<endl;
            //cout<<"prob here with axis_flag="<<axis_flag<<endl;
            comp startingpoint = min;
            double totpoints = points;// + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);

        }
    }



    
}


//This one has tags for taking the contour 
//under the second sheet of OPE, 
//with 2 tags that tells where ope should be taken
//to the unphysical sheet
void mom_vector_maker_45(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r,
                            int &tag1,
                            int &tag2    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"taking the straight line"<<endl;
            //cout<<"prob here with axis_flag="<<axis_flag<<endl;
            comp startingpoint = min;
            double totpoints = points;// + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);
            tag1 = 0;
            tag2 = 0;

        }
    }



    
}


//This one has tags for taking the contour 
//under the second sheet of OPE, 
//with 2 tags that tells where ope should be taken
//to the unphysical sheet, the cut still goes through
//the second sheet of ope even when the ope cuts 
//are not crossing the real q axis.
void mom_vector_maker_46(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r,
                            int &tag1,
                            int &tag2    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[6000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[6000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        

        }
    }



    
}

void tag_fixer( Eigen::VectorXcd &Gvec,
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

    for(int i=tag_initial;i<tag_final;++i)
    {
        comp p = qvec[i];
        comp Gvecval = Gvec[i];
        comp Gvecvalnext = Gvec[i+1]; 
        comp Gs1stsheet = GS_pk(s,p,k,m,eps);
        comp Gs2ndsheet = GS_pk_withtags(s,p,k,m,eps,i,i,1,qvec.size());

        comp pnext = qvec[i+1];
        comp Gs1stsheetnext = GS_pk(s,pnext,k,m,eps);
        comp Gs2ndsheetnext = GS_pk_withtags(s,pnext,k,m,eps,i,i,1,qvec.size());
        double real_G1st = real(Gs1stsheet);
        double real_G1stnext = real(Gs1stsheetnext);
        double real_G2nd = real(Gs2ndsheet);
        double real_G2ndnext = real(Gs2ndsheetnext);

        double firstsheet_secondsheetnext = abs(real_G1st - real_G2ndnext);
        double firstsheet_firstsheetnext = abs(real_G1st - real_G1stnext);

        if(firstsheet_firstsheetnext<firstsheet_secondsheetnext)
        {
            
        }
        else
        {
            Gvec[i] = Gs1stsheet;
            Gvec[i+1] = Gs2ndsheetnext; 
            if(somecount==0)
            {
                temp_tag = i+1;
                somecount = 1;
                break;
            }   
        }
    }

    for(int i=temp_tag;i<tag_final;++i)
    {
        comp p = qvec[i];
        comp Gs2ndsheet = GS_pk_withtags(s,p,k,m,eps,i,i,1,qvec.size());
        Gvec[i] = Gs2ndsheet;
    }

    problem_tag = temp_tag;
}

//this gvec_fixer does not work
void Gvec_fixer(    Eigen::VectorXcd &Gvec,
                    Eigen::VectorXcd &fixed_Gvec,
                    vector<comp> qvec,
                    comp s,
                    comp k,
                    double m,
                    double eps    )
{
    Eigen::VectorXcd temp_Gvec(qvec.size());
    temp_Gvec[0] = Gvec[0];
    /*for(int i=0;i<qvec.size()-1;++i)
    {
        comp p = qvec[i];
        comp next_p = qvec[i+1];
        comp this_gvec = Gvec[i];
        comp GS_firstsheet = GS_pk(s,next_p,k,m,eps);
        comp GS_secondsheet = GS_pk_secondsheet(s,next_p,k,m,eps);
        double realGS_firstsheet = real(GS_firstsheet);
        double realGS_secondsheet = real(GS_secondsheet);
        double this_real_gvec = real(this_gvec);
        if(abs(this_real_gvec-realGS_firstsheet)<abs(this_real_gvec-realGS_secondsheet))
        {
            temp_Gvec[i+1] = GS_firstsheet;
        }
        else
        {
            temp_Gvec[i+1] = GS_secondsheet;
        }

        fixed_Gvec = temp_Gvec;


    }
    */

    double past_difference = 0.0;
    vector<int> close_indices;
    vector<double> gap;
    for(int i=0;i<qvec.size()-1;++i)
    {
        comp p = qvec[i];
        comp future_p = qvec[i+1];
        comp GS_firstsheet = GS_pk(s,p,k,m,eps);
        comp GS_secondsheet = GS_pk_secondsheet(s,p,k,m,eps);
        comp future_GS_firstsheet = GS_pk(s,future_p,k,m,eps);
        comp future_GS_secondsheet = GS_pk_secondsheet(s,future_p,k,m,eps);
        double realGS_firstsheet = real(GS_firstsheet);
        double realGS_secondsheet = real(GS_secondsheet);
        double future_realGS_firstsheet = real(future_GS_firstsheet);
        double future_realGS_secondsheet = real(future_GS_secondsheet);
        
        double present_difference = abs(realGS_firstsheet - realGS_secondsheet);
        double future_difference = abs(future_realGS_firstsheet - future_realGS_secondsheet);

        
        if(i!=0)
        {
            if(past_difference>present_difference && present_difference<future_difference)
            {
                close_indices.push_back(i);
                gap.push_back(present_difference);
            }
        }
        past_difference = present_difference;
    }

    double temp_min_gap = gap[0];
    double temp_min = 0.01;
    vector<int> selected_indices;
    vector<double> selected_gap;
    for(int i=0;i<gap.size();++i)
    {
        if(gap[i]<temp_min)
        {
            selected_gap.push_back(gap[i]);
            selected_indices.push_back(close_indices[i]);
        }
    }
    int highest_index;
    int lowest_index;
    int tmp_ind_max = 0;
    int tmp_ind_min = INT_MAX;

    for(int i=0;i<selected_indices.size();++i)
    {
        if(selected_indices[i]>tmp_ind_max)
        {
            tmp_ind_max = selected_indices[i];
        }
        if(selected_indices[i]<tmp_ind_min)
        {
            tmp_ind_min = selected_indices[i];
        }

    }
    cout<<"close indices found at = "<<endl;
    for(int i=0;i<selected_indices.size();++i) cout<<selected_indices[i]<<'\t';
    
    cout<<endl;
    cout<<"min ind = "<<tmp_ind_min<<'\t'<<" max ind = "<<tmp_ind_max<<endl;
    for(int i=0;i<qvec.size();++i)
    {
        comp p = qvec[i];
        comp GS_firstsheet = GS_pk(s,p,k,m,eps);
        comp GS_secondsheet = GS_pk_secondsheet(s,p,k,m,eps);

        if(i>=tmp_ind_min && i<=tmp_ind_max)
            fixed_Gvec[i] = GS_secondsheet;
        else    
            fixed_Gvec[i] = GS_firstsheet;
        
    }
}

//this gvec_fixer_1 works
void Gvec_fixer_1(    Eigen::VectorXcd &Gvec,
                    Eigen::VectorXcd &fixed_Gvec,
                    vector<comp> qvec,
                    comp s,
                    comp k,
                    double m,
                    double eps,
                    int &tag1,
                    int &tag2    )
{

    
    Eigen::VectorXcd temp_Gvec(qvec.size());
    temp_Gvec[0] = Gvec[0];

    double past_difference = 0.0;
    vector<int> close_indices;
    vector<double> gap;
    for(int i=0;i<qvec.size()-1;++i)
    {
        comp p = qvec[i];
        comp future_p = qvec[i+1];
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
        
    }

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

    double temp_max2 = 0.0;
    int temp_ind2 = 0.0;

    for(int i=0;i<gap.size();++i)
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
    }

    if(temp_ind1>temp_ind2)
    {
        int sometemp = temp_ind2;
        temp_ind2 = temp_ind1;
        temp_ind1 = sometemp;
    }
    
    
    
    cout<<"min ind = "<<temp_ind1<<'\t'<<" max ind = "<<temp_ind2<<endl;
    cout<<"gap 1 = "<<temp_max1<<'\t'<<" gap 2 = "<<temp_max2<<endl;
    for(int i=0;i<qvec.size();++i)
    {
        comp p = qvec[i];
        comp GS_firstsheet = GS_pk(s,p,k,m,eps);
        comp GS_secondsheet = GS_pk_secondsheet(s,p,k,m,eps);

        if(i>=temp_ind1 && i<=temp_ind2)
            fixed_Gvec[i] = GS_secondsheet;
        else    
            fixed_Gvec[i] = GS_firstsheet;
        
    }

    tag1 = temp_ind1;
    tag2 = temp_ind2;
}


//This one has tags for taking the contour 
//under the second sheet of OPE, 
//with 2 tags that tells where ope should be taken
//to the unphysical sheet, the cut still goes through
//the second sheet of ope even when the ope cuts 
//are not crossing the real q axis. Also this 
//one goes through the pcuts only once 
void mom_vector_maker_47(    vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r,
                            int &tag1,
                            int &tag2    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    //int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
        

        }
    }



    
}


void mom_vector_maker_47_with_weights(    vector<comp> &qvec,
                                          vector<comp> &weights,
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r,
                            int &tag1,
                            int &tag2    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    //int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);
            line_maker_with_weights(qvec,weights,startingpoint,-delq,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);


        }
    }



    
}


//Same as contour_47 but we can take 
//the first line to left as much as we want
void mom_vector_maker_47_1( double somepoint,
                            vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r,
                            int &tag1,
                            int &tag2    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    //int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - somepoint*rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = somepoint*pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = somepoint*pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = somepoint*pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = somepoint*pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = somepoint*pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = somepoint*pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
        

        }
    }



    
}

//Same as contour_47 but we can take 
//the line that goes through the second OPE 
//sheet to bend to left as much as we want
void mom_vector_maker_47_2( double somepoint,
                            vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r,
                            int &tag1,
                            int &tag2    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    //int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/10.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[9700];//pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[9700];//pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            comp end_temp = endingpoint_for_slanted_maker;
            startingpoint = qvec[qvec.size()-1];
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = -somepoint*1.5*pcutpvec[1000];
            }
            else 
            {
                endingpoint_for_slanted_maker = -somepoint*1.5*pcutmvec[1000];
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,4.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            endingpoint_for_slanted_maker = end_temp;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,4.0*points/10.0);

            //startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker = end_temp;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);



            /////////////////////////////////////////
            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
        

        }
    }



    
}

void mom_vector_maker_47_2_1( double somepoint,
                            vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r,
                            int &tag1,
                            int &tag2    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    //int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                //endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;
                endingpoint_for_slanted_maker = pcutmvec[8000];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                //endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;
                endingpoint_for_slanted_maker = pcutpvec[8000];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1 + 1 ;

            comp temp_real = real(startingpoint);

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[9700];;//pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[9700];//pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            comp temp_imag = imag(endingpoint_for_slanted_maker) + somepoint;
            comp end_temp = endingpoint_for_slanted_maker;
            endingpoint_for_slanted_maker = temp_real + ii*temp_imag;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            endingpoint_for_slanted_maker = end_temp;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            //startingpoint = qvec[qvec.size()-1];
            //if(flag_pcut==0)
            //{
            //    endingpoint_for_slanted_maker = -somepoint*1.5*pcutpvec[1000];
            //}
            //else 
            //{
            //    endingpoint_for_slanted_maker = -somepoint*1.5*pcutmvec[1000];
            //}
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            //startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker = end_temp;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            //startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker = end_temp;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);



            /////////////////////////////////////////
            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
        

        }
    }



    
}

void mom_vector_maker_47_2_2( double somepoint,
                            vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r,
                            int &tag1,
                            int &tag2    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    //int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                //endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;
                endingpoint_for_slanted_maker = pcutmvec[8000];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                //endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;
                endingpoint_for_slanted_maker = pcutpvec[8000];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size();

            double temp_real = real(startingpoint) - somepoint;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[9700];;//pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[9700];//pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            double temp_imag = imag(endingpoint_for_slanted_maker);
            comp end_temp = endingpoint_for_slanted_maker;
            endingpoint_for_slanted_maker = temp_real + ii*temp_imag;
        
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            endingpoint_for_slanted_maker = end_temp;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            //startingpoint = qvec[qvec.size()-1];
            //if(flag_pcut==0)
            //{
            //    endingpoint_for_slanted_maker = -somepoint*1.5*pcutpvec[1000];
            //}
            //else 
            //{
            //    endingpoint_for_slanted_maker = -somepoint*1.5*pcutmvec[1000];
            //}
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            //startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker = end_temp;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            //startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker = end_temp;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);



            /////////////////////////////////////////
            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
        

        }
    }



    
}

void mom_vector_maker_47_2_3( double somepoint,
                            vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r,
                            int &tag1,
                            int &tag2    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    //int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {   
            comp endingpoint_for_slanted_maker;
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            if(flag_pcut==0)
            {
                //lowest_index_yval = imag(pcutmvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutmvec[pcutmvec.size()/2];;//pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                //highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                //lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutpvec[pcutpvec.size()/2];//pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                //highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            } 
            //comp delq = abs(startingpoint - rad)/(points/5.0);
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size();
            delq = abs(startingpoint - rad/4.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            
            double temp_real = real(startingpoint);// - somepoint;
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            //(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);
            
            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[9590];;//pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[9590];//pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            double temp_imag = imag(endingpoint_for_slanted_maker);
            comp end_temp = endingpoint_for_slanted_maker;
            endingpoint_for_slanted_maker = temp_real + ii*temp_imag;
        
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            endingpoint_for_slanted_maker = end_temp;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            /////////////////////////////////////////
            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
        

        }
    }



    
}



//Same as contour_47 but we can take 
//the OFFSET_R to be as long as we want
void mom_vector_maker_47_3( double somepoint,
                            vector<comp> &qvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double eps_for_m2k,
                            double points,
                            double qvec_r,
                            int &tag1,
                            int &tag2    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    //int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/10.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);

            //comp end_temp = endingpoint_for_slanted_maker;
            //startingpoint = qvec[qvec.size()-1];
            //if(flag_pcut==0)
            //{
            //    endingpoint_for_slanted_maker = -somepoint*1.5*pcutpvec[500];
            //}
            //else 
            //{
            //    endingpoint_for_slanted_maker = -somepoint*1.5*pcutmvec[500];
            //}
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            //startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker = end_temp;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            //startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker = end_temp;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);



            /////////////////////////////////////////
            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + somepoint*qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,4.0*points/10.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,3.0*points/10.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/10.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,2.0*points/10.0);
        

        }
    }



    
}


void mom_vector_maker_47_2_A1(  int somepoint,  
                                vector<comp> &qvec,   
                                comp s,
                                comp min,
                                comp max,
                                double a,
                                double m,
                                double eps,
                                double eps_for_m2k,
                                double points,
                                double qvec_r,
                                int &tag1,
                                int &tag2    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;
    //int tag1=0,tag2=0;  // these are added here because 
                        // this OPE function does not need the tags 
                        // but the vector making function still needed them.

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double highestyvalp,highestyvalm;
    int highest_py_index, highest_my_index;
    highest_yval_pcuts(pcutpvec,pcutmvec,highestyvalp,highestyvalm,highest_py_index,highest_my_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            cout<<"deforming the contour backward first"<<endl;
            comp startingpoint = min;
            comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
            mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            //delq = abs(imag(startingpoint) - yval1 + r)/points;
            //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
            comp endingpoint_for_slanted_maker;//(pcutp_close_selected + pcutm_close_selected)/2.0;
            //line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points);

            //startingpoint = qvec[qvec.size()-1];

            //cout<<"last point after step = "<<startingpoint<<endl;
            double lowest_index_yval = 0.0;
            
            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                endingpoint_for_slanted_maker = pcutpvec[selected_p_index - somepoint];//pcutmvec[lowest_m_index];// - ii*0.01;

                highest_index_xval = real(pcutpvec[highest_p_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            else
            {
                lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                endingpoint_for_slanted_maker = pcutmvec[selected_m_index - somepoint];//pcutpvec[lowest_p_index];// - ii*0.01;

                highest_index_xval = real(pcutmvec[highest_m_index]);
                //cout<<"highest x val = "<<highest_index_xval<<endl;
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            tag2 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
        }
        else
        {
            cout<<"deforming the contour to left cut 7500 index first"<<endl;
            

            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            tag2 = qvec.size() - 1;
            
            startingpoint = qvec[qvec.size()-1];    
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        double highest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            highest_index_yval = imag(pcutpvec[highest_py_index]);
            
        }
        else
        {
            highest_index_yval = imag(pcutmvec[highest_py_index]);
        }


        if(abs(highest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[8500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[8500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        
        }
        else 
        {
            cout<<"still still deforming the contour to left cut 7500 index first"<<endl;
            
            comp startingpoint = min;//qvec[qvec.size()-1];
            comp endingpoint_for_slanted_maker;
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5000]/2.0;//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5000]/2.0;//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            
            
            startingpoint = qvec[qvec.size()-1];
            
            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutmvec[5500];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutpvec[5500];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                //highest_index_xval = real(pcutpvec[highest_m_index]);
            }
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);


            startingpoint = qvec[qvec.size()-1];
            tag1 = qvec.size() - 1;

            endingpoint_for_slanted_maker = ii*imag(startingpoint) - 4.5*abs(real(qc1));
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];

            if(flag_pcut==0)
            {
                endingpoint_for_slanted_maker = pcutpvec[7000];//ii*imag(pcutpvec[lowest_p_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_p_index]);
            }
            else
            {
                endingpoint_for_slanted_maker = pcutmvec[7000];//ii*imag(pcutmvec[lowest_m_index]) - ii*0.01;
                highest_index_xval = real(pcutpvec[highest_m_index]);
            }

            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
            
            startingpoint = qvec[qvec.size()-1];    
            tag2 = qvec.size() - 1;
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);

            startingpoint = qvec[qvec.size()-1];
            //cout<<"last point after step = "<<startingpoint<<endl;
            endingpoint_for_slanted_maker = max;
            line_maker(qvec,startingpoint,endingpoint_for_slanted_maker,points/5.0);
        

        }
    }



    
}


//these are initial sigvec makers, these chose the //
//contour by hand mostly using 5 different for     //
//loops. maker and maker_1 are of these kind.      //
void sigma_vector_maker(    vector<comp> &sigvec,
                            comp s,
                            comp sigmin,
                            comp sigmax,
                            double a,
                            double m,
                            double points,
                            double eps   )
{
    comp ii = {0.0,1.0};
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);

    if(real(rad)>0.0)
    {
        //double eps = 1.0e-5;
        double length = + real(sigmac1 - 2.0*eps) 
                        - real(sigmin)
                        + real(sigmac2 + 2.0*eps)
                        - real(sigmax)
                        + 4.0*real(rad + 2.0*eps);
        comp delsig = length/points;

        for(double tempsig=real(sigmin);tempsig<real(sigmac1 - 2.0*eps);tempsig=tempsig+real(delsig))
        {
            comp sigma = tempsig + ii*eps;

            sigvec.push_back(sigma);
        }
        for(double tempsig=0.0;tempsig<real(rad+2.0*eps);tempsig=tempsig+real(delsig))
        {
            comp sigma = sigmac1 - 2.0*eps + ii*tempsig + ii*eps;

            sigvec.push_back(sigma);
        }
        for(double tempsig=0.0;tempsig<2.0*real(rad+2.0*eps);tempsig=tempsig+real(delsig))
        {
            comp sigma = sigmac1 - 2.0*eps + ii*(rad+2.0*eps) + tempsig + ii*eps;

            sigvec.push_back(sigma);
        }
        for(double tempsig=0.0;tempsig<real(rad+2.0*eps);tempsig=tempsig+real(delsig))
        {
            comp sigma = sigmac1 - 2.0*eps + ii*(rad+2.0*eps) - ii*tempsig + 2.0*(rad+2.0*eps) + ii*eps;

            sigvec.push_back(sigma);
        }
        for(double tempsig=real(sigmac2 + 2.0*eps);tempsig>=real(sigmax);tempsig=tempsig-real(delsig))
        {
            comp sigma = tempsig + ii*eps;

            sigvec.push_back(sigma);
        }
    }
    else
    {
        //double eps = 1.0e-5;
        double length = real(sigmax - sigmin);
        comp delsig = length/points;

        for(double tempsig=real(sigmin);tempsig<=real(sigmax);tempsig=tempsig+real(delsig))
        {
            comp sigma = tempsig + ii*eps;

            sigvec.push_back(sigma);
        }
    }

    cout<<"sigvec created with size = "<<sigvec.size()<<endl;
}


void sigma_vector_maker_1(    vector<comp> &sigvec,
                            comp s,
                            comp sigmin,
                            comp sigmax,
                            double a,
                            double m,
                            double points,
                            double eps   )
{
    comp ii = {0.0,1.0};
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);

    cout<<"sigmin:"<<sigmin<<endl;
    cout<<"sigmax:"<<sigmax<<endl;
    cout<<"rad:"<<rad<<endl;
    cout<<"sigmac1:"<<sigmac1<<endl;
    cout<<"sigmac2:"<<sigmac2<<endl;

    if(real(rad)>0.0)
    {
        //double eps = 1.0e-5;
        double length = + real(sigmac1 - 20.0*eps) 
                        - real(sigmin)
                        + real(sigmac2 + 20.0*eps)
                        - real(sigmax)
                        + 4.0*real(rad + 20.0*eps);
        comp delsig = (comp)length/points;

        for(double tempsig=real(sigmin);tempsig<real(sigmac1 - 20.0*eps);tempsig=tempsig+real(delsig))
        {
            comp sigma = tempsig + ii*eps;

            sigvec.push_back(sigma);
        }
        for(double tempsig=0.0 + real(delsig);tempsig<real(rad+20.0*eps);tempsig=tempsig+real(delsig))
        {
            comp sigma = sigmac1 - 20.0*eps + ii*tempsig + ii*eps;

            sigvec.push_back(sigma);
        }
        for(double tempsig=0.0 + real(delsig);tempsig<2.0*real(rad+20.0*eps);tempsig=tempsig+real(delsig))
        {
            comp sigma = sigmac1 - 20.0*eps + ii*(rad+20.0*eps) + tempsig + ii*eps;

            sigvec.push_back(sigma);
        }
        for(double tempsig=0.0 + real(delsig);tempsig<real(rad+20.0*eps);tempsig=tempsig+real(delsig))
        {
            comp sigma = sigmac1 - 20.0*eps + ii*(rad+20.0*eps) - ii*tempsig + 2.0*(rad+20.0*eps) + ii*eps;

            sigvec.push_back(sigma);
        }
        for(double tempsig=real(sigmac2 + 20.0*eps);tempsig>=real(sigmax);tempsig=tempsig-real(delsig))
        {
            comp sigma = tempsig + ii*eps;

            sigvec.push_back(sigma);
        }
    }
    else
    {
        //double eps = 1.0e-5;
        double length = real(sigmax - sigmin);
        comp delsig = length/points;

        for(double tempsig=real(sigmin);tempsig<=real(sigmax);tempsig=tempsig+real(delsig))
        {
            comp sigma = tempsig + ii*eps;

            sigvec.push_back(sigma);
        }
    }

    cout<<"sigvec created with size = "<<sigvec.size()<<endl;
}


//This one makes the sigvec above threshold     //
//It makes a uniform sigvec between sigmin and  //
//sigmax                                        //
void sigma_vector_maker_abovethreshold( vector<comp> &sigvec,
                                        comp s,
                                        comp sigmamin,
                                        comp sigmamax,
                                        double N            )
{
    double length = real(sigmamax - sigmamin);
    comp delsig = (comp)length/N;

    for(double tempsig=real(sigmamin);tempsig<=real(sigmamax);tempsig=tempsig+real(delsig))
    {
        comp sigma = tempsig;
        sigvec.push_back(sigma);
    }

    cout<<"sigvec created with size = "<<sigvec.size()<<endl;
}

//this one makes specific lines in real or imag //
//axis based on the 'ud' number. It takes the   //
//max and min values, divide it by points to get//
//delsig, which it then adds from min towards   //
//max. This might be a little buggy. So will be //
//building a new one.                           //

void sigma_vector_maker_linear( vector<comp> &sigvec,
                                comp min,
                                comp max,
                                double points,
                                double eps,
                                int ud   )
{
    double delsig = abs(max - min)/points;
    comp ii = {0.0,1.0};
    if(ud==1)
    {
        for(int i=1;i<=(int)points;++i)
        {
            comp sigma = min + (comp)i*delsig + ii*eps ;
            //if(sigma<=max)
            //{
                sigvec.push_back(sigma);
            //}
            //else break;
        } 
    }
    else if(ud==2)
    {
        for(int i=1;i<=(int)points;++i)
        {
            comp sigma = min + (comp)i*ii*delsig + ii*eps;
            //if(sigma<=max)
            //{
                sigvec.push_back(sigma);
            //}
            //else break;
        }

    }
    else if(ud==3)
    {
        for(int i=1;i<=(int)points;++i)
        {
            comp sigma = min - (comp)i*ii*delsig + ii*eps;
            //if(sigma<=max)
            //{
                sigvec.push_back(sigma);
            //}
            //else break;
        }

    }
    else if(ud==4)
    {
        for(int i=1;i<=(int)points;++i)
        {
            //cout<<"im here"<<endl;
            comp sigma = min - (comp)i*delsig + ii*eps;
            //if(sigma<=max)
            //{
                sigvec.push_back(sigma);
            //}
            //else break;
        }

    }
}

//This momentum vector makers are defined according
//to sebastians definitions found in the paper
void mom_vector_maker_seba_imsneg(  vector<comp> &qvec,
                                    vector<comp> &weights, 
                                    comp s,
                                    comp kmin,
                                    comp kmax,
                                    double a,
                                    double m,
                                    double eps,
                                    double eps_for_m2k,
                                    double points   )
{
    //initial setup//
    comp ii = {0.0,1.0};
    comp q = pmom(s,sigmab(a,m),m);
    comp qc2 = q_c2(s,a,m,eps);
    comp pplus = pcut_plus(1.0,s,q,m,eps);
    comp pminus = pcut_plus(-1.0,s,q,m,eps);
    comp qc1 = q_c1(real(s),a,m,eps);
    double rad = real(qc1);
    //cout<<"radius = "<<rad<<endl;
    
    if(rad>0.0)
    {
        //cout<<"contour chosen"<<endl;
        //firstnode - secondnode
        comp firstnode = kmin;
        comp secondnode = -2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2);
        line_maker_with_weights(qvec,weights,firstnode,secondnode,2.0*points/10.0);

        //secondnode - thirdnode
        comp thirdnode = 0.5*real(pplus - pminus) - abs(qc2)*ii;
        line_maker_with_weights(qvec,weights,secondnode,thirdnode,2.0*points/10.0);

        //thirdnode - fourthnode 
        comp fourthnode = (2.0/3.0)*(1.0 - 2.0*ii)*abs(qc2);
        line_maker_with_weights(qvec,weights,thirdnode,fourthnode,2.0*points/10.0);

        //fourthnode - fifthnode
        comp fifthnode = (3.0/2.0)*abs(qc2);
        line_maker_with_weights(qvec,weights,fourthnode,fifthnode,2.0*points/10.0);

        //fifthnode - sixthnode 
        comp sixthnode = kmax;
        line_maker_with_weights(qvec,weights,fifthnode,sixthnode,2.0*points/10.0);
    }
    else
    {
        //cout<<"straight line chosen"<<endl;
        line_maker_with_weights(qvec,weights,kmin,kmax,points);     
    }

}

//This momentum vector makers are defined according
//to sebastians definitions found in the paper
void mom_vector_maker_seba_imspos(  vector<comp> &qvec,
                                    vector<comp> &weights, 
                                    comp s,
                                    comp kmin,
                                    comp kmax,
                                    double a,
                                    double m,
                                    double eps,
                                    double eps_for_m2k,
                                    double points,
                                    int &tag1,
                                    int &tag2,
                                    int &switch_for_gvec_fixer   )
{
    //initial setup//
    comp ii = {0.0,1.0};
    comp q = pmom(s,sigmab(a,m),m);
    comp qc2 = q_c2(s,a,m,eps);
    comp pplus = pcut_plus(1.0,s,q,m,eps);
    comp pminus = pcut_plus(-1.0,s,q,m,eps);
    comp qc1 = q_c1(real(s),a,m,eps);
    comp sigb = sigmab(a,m);
    double rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    
    if(rad>0.0)
    {
        cout<<"contour chosen"<<endl;
        //firstnode - secondnode
        comp firstnode = kmin;
        comp secondnode = -2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2);
        line_maker_with_weights(qvec,weights,firstnode,secondnode,4.0*points/20.0);

        //secondnode - thirdnode
        comp sigmac2 = sigc2(s,sigb,m);
        comp p = pmom(s,sigmac2,m);
        comp k = q;
        comp alphapk = ((sqrt(s)-omega_comp(p,m)-omega_comp(k,m))*(sqrt(s)-omega_comp(p,m)-omega_comp(k,m)) - p*p - k*k - m*m)/(2.0*p*k);
        //comp zpk = 2.0*s*(sigb + ii*eps) - (s+sigb - m*m)*(s+m*m-sigmac2);
        //comp zpk = 2.0*s*(sigmac2 + ii*eps) - (s+sigmac2 - m*m)*(s+m*m-sigb);
        
        double x0 = -1.0*real(alphapk);
        //double x0 = real(zpk);
        cout<<"sigma c2 = "<<sigmac2<<endl;
        
        cout<<"x0 = "<<x0<<endl;
        comp thirdnode = -pcut_minus(-x0,s,q,m,eps);
        line_maker_with_weights(qvec,weights,secondnode,thirdnode,4.0*points/20.0);

        tag1 = qvec.size() - 1;
        //thirdnode - fourthnode 
        comp fourthnode = pcut_minus(x0,s,q,m,eps);
        line_maker_with_weights(qvec,weights,thirdnode,fourthnode,1.0*points/20.0);

        tag2 = qvec.size() - 1;
        //fourthnode - fifthnode
        comp fifthnode = (2.0/3.0)*(1.0-2.0*ii)*abs(qc2);
        line_maker_with_weights(qvec,weights,fourthnode,fifthnode,3.0*points/20.0);

        //fifthnode - sixthnode 
        comp sixthnode = (3.0/2.0)*abs(qc2);
        line_maker_with_weights(qvec,weights,fifthnode,sixthnode,4.0*points/20.0);

        //sixthnode - seventhnode
        comp seventhnode = kmax;
        line_maker_with_weights(qvec,weights,sixthnode,seventhnode,4.0*points/20.0);

        switch_for_gvec_fixer = 0;

    }
    else
    {
        cout<<"straight line chosen"<<endl;
        line_maker_with_weights(qvec,weights,kmin,kmax,points);   
        switch_for_gvec_fixer = 1;  
    }

}

//instead of going through the lowest point between
//the two cuts, it goes through the cut from little bit
//above

void mom_vector_maker_seba_imspos_1(  vector<comp> &qvec,
                                    vector<comp> &weights, 
                                    comp s,
                                    comp kmin,
                                    comp kmax,
                                    double a,
                                    double m,
                                    double eps,
                                    double eps_for_m2k,
                                    double points,
                                    int &tag1,
                                    int &tag2,
                                    int &switch_for_gvec_fixer   )
{
    //initial setup//
    comp ii = {0.0,1.0};
    comp q = pmom(s,sigmab(a,m),m);
    comp qc2 = q_c2(s,a,m,eps);
    comp pplus = pcut_plus(1.0,s,q,m,eps);
    comp pminus = pcut_plus(-1.0,s,q,m,eps);
    comp qc1 = q_c1(real(s),a,m,eps);
    comp sigb = sigmab(a,m);
    double rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    
    if(rad>0.0)
    {
        switch_for_gvec_fixer = 0;
        cout<<"contour chosen"<<endl;
        //firstnode - secondnode
        comp firstnode = kmin;
        comp secondnode = -2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2);
        line_maker_with_weights(qvec,weights,firstnode,secondnode,2.0*points/10.0);

        //secondnode - thirdnode
        comp sigmac2 = sigc2(s,sigb,m);
        comp p = pmom(s,sigmac2,m);
        comp k = q;
        comp alphapk = ((sqrt(s)-omega_comp(p,m)-omega_comp(k,m))*(sqrt(s)-omega_comp(p,m)-omega_comp(k,m)) - p*p - k*k - m*m)/(2.0*p*k);
        //comp zpk = 2.0*s*(sigb + ii*eps) - (s+sigb - m*m)*(s+m*m-sigmac2);
        //comp zpk = 2.0*s*(sigmac2 + ii*eps) - (s+sigmac2 - m*m)*(s+m*m-sigb);
        
        double x0 = -1.0*real(alphapk);
        //double x0 = real(zpk);
        cout<<"sigma c2 = "<<sigmac2<<endl;
        
        cout<<"x0 = "<<x0<<endl;
        comp thirdnode = -pcut_minus((-x0-0.97)/2.0,s,q,m,eps);
        line_maker_with_weights(qvec,weights,secondnode,thirdnode,2.0*points/10.0);

        tag1 = qvec.size() - 1;
        //thirdnode - fourthnode 
        comp fourthnode = pcut_minus(2.0*x0/3.0,s,q,m,eps);
        line_maker_with_weights(qvec,weights,thirdnode,fourthnode,2.0*points/10.0);

        tag2 = qvec.size() - 1;
        //fourthnode - fifthnode
        comp fifthnode = real(qc1) + ii*imag(fourthnode);//(2.0/3.0)*(1.0-2.0*ii)*abs(qc2);
        line_maker_with_weights(qvec,weights,fourthnode,fifthnode,2.0*points/10.0);

        //fifthnode - sixthnode 
        comp sixthnode = (3.0/2.0)*abs(qc2);
        line_maker_with_weights(qvec,weights,fifthnode,sixthnode,1.0*points/10.0);

        //sixthnode - seventhnode
        comp seventhnode = kmax;
        line_maker_with_weights(qvec,weights,sixthnode,seventhnode,1.0*points/10.0);

    }
    else
    {
        switch_for_gvec_fixer = 1;
        cout<<"straight line chosen"<<endl;
        line_maker_with_weights(qvec,weights,kmin,kmax,points);     
    }

}

//this one is build finally by looking at a=16 kernel
//performed glockle test and it succeeded
void mom_vector_maker_seba_imspos_2(  vector<comp> &qvec,
                                    vector<comp> &weights, 
                                    comp s,
                                    comp kmin,
                                    comp kmax,
                                    double a,
                                    double m,
                                    double eps,
                                    double eps_for_m2k,
                                    double points,
                                    int &tag1,
                                    int &tag2,
                                    int &switch_for_gvec_fixer   )
{
    //initial setup//
    comp ii = {0.0,1.0};
    comp q = pmom(s,sigmab(a,m),m);
    comp qc2 = q_c2(s,a,m,eps);
    comp pplus = pcut_plus(1.0,s,q,m,eps);
    comp pminus = pcut_plus(-1.0,s,q,m,eps);
    comp qc1 = q_c1(real(s),a,m,eps);
    comp sigb = sigmab(a,m);
    double rad = abs(real(qc1));
    cout<<"radius = "<<rad<<endl;
    
    if(rad>0.0)
    {
        switch_for_gvec_fixer = 0;
        cout<<"contour chosen"<<endl;
        //firstnode - secondnode
        comp firstnode = kmin;
        comp secondnode = -rad/8.0;//-2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2);
        //1
        line_maker_with_weights(qvec,weights,firstnode,secondnode,1.0*points/4.0);

        //secondnode - thirdnode
        comp sigmac2 = sigc2(s,sigb,m);
        comp p = (q_c1(s,a,m,eps) + qc2)/2.0;//pmom(s,sigmac2,m);
        comp k = q;
        comp alphapk = ((sqrt(s)-omega_comp(p,m)-omega_comp(k,m))*(sqrt(s)-omega_comp(p,m)-omega_comp(k,m)) - p*p - k*k - m*m)/(2.0*p*k);
        //comp zpk = 2.0*s*(sigb + ii*eps) - (s+sigb - m*m)*(s+m*m-sigmac2);
        //comp zpk = 2.0*s*(sigmac2 + ii*eps) - (s+sigmac2 - m*m)*(s+m*m-sigb);
        
        double x0 = -1.0*real(alphapk);
        //double x0 = real(zpk);
        cout<<"sigma c2 = "<<sigmac2<<endl;
        
        cout<<"x0 = "<<x0<<endl;
        //comp thirdnode = -pcut_minus((-x0-0.77)/2.0,s,q,m,eps);
        
        
        //tag1 = qvec.size() - 1;
        //thirdnode - fourthnode 
        comp f1 = pcut_minus((x0+1.0)/2.0,s,q,m,eps);
        
        comp f2 = 0.5*(real(pplus - pminus)) - ii*abs(qc2);

        comp temp_f3_1 = pcut_minus(-1,s,q,m,eps);
        comp temp_f3_2 = pcut_minus(1,s,q,m,eps);
        comp temp_f3_3 = pcut_plus(-1,s,q,m,eps);
        comp temp_f3_4 = pcut_plus(1,s,q,m,eps);
        comp f3 = 0.5*(real(pplus - pminus)) - ii*5.0/4.0*abs(qc2);;
        /*if(imag(temp_f3_1)<imag(temp_f3_2) && imag(temp_f3_1)<imag(temp_f3_3) && imag(temp_f3_1)<imag(temp_f3_4))
        {
            f3 = pcut_minus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_2)<imag(temp_f3_1) && imag(temp_f3_2)<imag(temp_f3_3) && imag(temp_f3_2)<imag(temp_f3_4))
        {
            f3 = pcut_minus(1,s,q,m,eps);
        }
        else if(imag(temp_f3_3)<imag(temp_f3_1) && imag(temp_f3_3)<imag(temp_f3_2) && imag(temp_f3_3)<imag(temp_f3_4))
        {
            f3 = pcut_plus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_4)<imag(temp_f3_1) && imag(temp_f3_4)<imag(temp_f3_2) && imag(temp_f3_4)<imag(temp_f3_3))
        {
            f3 = pcut_plus(1,s,q,m,eps);
        }*/

        comp thirdnode = real(secondnode) + ii*imag((f2 + f3)/2.0);//-pcut_minus(-x0,s,q,m,eps);
        //2
        line_maker_with_weights(qvec,weights,secondnode,thirdnode,1.0*points/4.0);
        tag1 = qvec.size() - 1;
        
        //comp f3 = pcut_minus(1,s,q,m,eps);
        //comp fourthnode = real(thirdnode) + ii*imag((f2 + f3)/2.0);//real(thirdnode) + ii*imag(f1)*0.75;//pcut_minus(2.0*x0/3.0,s,q,m,eps);
        //3
        //line_maker_with_weights(qvec,weights,thirdnode,fourthnode,2.0*points/10.0);

        //tag2 = qvec.size() - 1;
        //fourthnode - fifthnode
        //comp fifthnode = (f2 + f3)/2.0;//real(f1) + ii*imag(f1)*0.75;//real(qc1) + ii*imag(fourthnode);//(2.0/3.0)*(1.0-2.0*ii)*abs(qc2);
        //4
        //line_maker_with_weights(qvec,weights,fourthnode,fifthnode,2.0*points/10.0);

        //fifthnode - sixthnode 
        comp sixthnode = real(q_c1(s,a,m,eps)) + (ii*imag((f2+f3)/2.0));//(3.0/2.0)*abs(qc2);
        //5
        //line_maker_with_weights(qvec,weights,fifthnode,sixthnode,1.0*points/4.0);
        line_maker_with_weights(qvec,weights,thirdnode,sixthnode,1.0*points/4.0);
        tag2 = qvec.size() - 1;
        //sixthnode - seventhnode
        comp seventhnode = kmax;
        //6
        line_maker_with_weights(qvec,weights,sixthnode,seventhnode,1.0*points/4.0);

    }
    else
    {
        switch_for_gvec_fixer = 1;
        cout<<"straight line chosen"<<endl;
        line_maker_with_weights(qvec,weights,kmin,kmax,points);     
    }

}


//this one builds upon contour43, uses seba contour imsneg 
//for the first part of contour making in a=16, build for 
//mphib 3d plots
void mom_vector_maker_43_with_weights_with_seba_imsneg(      vector<comp> &qvec, 
                                            vector<comp> &weights,
                                            comp s,
                                            comp min,
                                            comp max,
                                            double a,
                                            double m,
                                            double eps,
                                            double eps_for_m2k,
                                            double points,
                                            double qvec_r    )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    comp qc2 = q_c2(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};
    //double qvec_r = 0.01;

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    
    cout<<"axis_flag = "<<axis_flag<<endl;
    if(axis_flag==0)
    cout<<"cut to the right crosses real line"<<endl;
    else
    cout<<"cut to the right doesnot cross real line"<<endl;

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;


    if(axis_flag==0)
    {
        //if(real(rad)>=0.0505)
        if(problem_index_flag==1)
        {
            if(abs(imag(s))<0.10)
            {
                cout<<"deforming the contour backward first 1 "<<endl;
                comp q = pmom(s,sigmab(a,m),m);
                comp qc2 = q_c2(s,a,m,eps);
                comp pplus = pcut_plus(-1.0,s,q,m,eps);
                comp pminus = pcut_plus(+1.0,s,q,m,eps);
                cout<<"pplus = "<<pplus<<'\t'<<"pminus = "<<pminus<<endl;
                cout<<"pplus-pminus = "<<0.5*real(pplus - pminus)<<endl;
                comp qc1 = q_c1(s,a,m,eps);
                double rad = real(qc1);

                comp firstnode = min;
                comp secondnode = -2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2*qc2/qc1);
                line_maker_with_weights(qvec,weights,firstnode,secondnode,2.0*points/10.0);

                //secondnode - thirdnode
                comp thirdnode = 0.5*real(pplus - pminus) - abs(qc2)*ii;
                line_maker_with_weights(qvec,weights,secondnode,thirdnode,2.0*points/10.0);

                //thirdnode - fourthnode 
                comp fourthnode = (2.0/3.0)*(1.0 - 1.68*ii)*abs(qc2);//(2.0/3.0)*(1.0 - 2.0*ii)*abs(qc2);
                line_maker_with_weights(qvec,weights,thirdnode,fourthnode,2.0*points/10.0);

                //fourthnode - fifthnode
                comp fifthnode = (3.0/2.0)*abs(qc2);
                line_maker_with_weights(qvec,weights,fourthnode,fifthnode,2.0*points/10.0);

                //fifthnode - sixthnode 
                comp sixthnode = max;
                line_maker_with_weights(qvec,weights,fifthnode,sixthnode,2.0*points/10.0);
                //////////////////////////////////////////////////////
            }
            else 
            {
                cout<<"deforming contour background first 2"<<endl;
                comp startingpoint = min;
                comp delq = abs(startingpoint - rad/2.0)/(points/5.0);
                //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points/5.0,1);
                line_maker_with_weights(qvec,weights,startingpoint,-rad/2.0,points/5.0);
            
                cout<<"turningpoint 0 = "<<qvec[0]<<endl;
                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 1 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                //delq = abs(imag(startingpoint) - yval1 + r)/points;
                //mom_vector_maker_linear_1(qvec,startingpoint,-delq,points,2);
                comp endingpoint_for_slanted_maker = (pcutp_close_selected + pcutm_close_selected)/2.0;
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

                startingpoint = qvec[qvec.size()-1];
                cout<<"turningpoint 2 = "<<startingpoint<<endl;
                //cout<<"last point after step = "<<startingpoint<<endl;
                double lowest_index_yval = 0.0;
            
                if(flag_pcut==0)
                {
                    lowest_index_yval = imag(pcutpvec[lowest_p_index]);
                    endingpoint_for_slanted_maker = pcutpvec[lowest_p_index] - ii*0.1*abs(imag(pcutpvec[lowest_p_index]));

                    highest_index_xval = real(pcutpvec[highest_p_index]);
                    //cout<<"highest x val = "<<highest_index_xval<<endl;
                }
                else
                {
                    lowest_index_yval = imag(pcutmvec[lowest_m_index]);
                    endingpoint_for_slanted_maker = pcutmvec[lowest_m_index] - ii*0.1*abs(imag(pcutmvec[lowest_m_index]));

                    highest_index_xval = real(pcutmvec[highest_m_index]);
                    //cout<<"highest x val = "<<highest_index_xval<<endl;
                }

                if(lowest_index_yval<=imag(startingpoint))
                {
                    line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);

                    startingpoint = qvec[qvec.size()-1];
                    cout<<"turningpoint 3 = "<<startingpoint<<endl;
                    //cout<<"last point after step = "<<startingpoint<<endl;
                    endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                    line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0); 

                    startingpoint = qvec[qvec.size()-1];
                    cout<<"turningpoint 4 = "<<startingpoint<<endl;
                    //cout<<"last point after step = "<<startingpoint<<endl;
                    endingpoint_for_slanted_maker = max;
                    line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points/5.0);
                    cout<<"end point 5 = "<<qvec[qvec.size() - 1]<<endl;
                }
                else
                {
                    startingpoint = qvec[qvec.size()-1];
                    cout<<"turningpoint 3 = "<<startingpoint<<endl;
                    //cout<<"last point after step = "<<startingpoint<<endl;
                    endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                    line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
                
                    startingpoint = qvec[qvec.size()-1];
                    cout<<"turningpoint 4 = "<<startingpoint<<endl;
                    //cout<<"last point after step = "<<startingpoint<<endl;
                    endingpoint_for_slanted_maker = max;
                    line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
                    cout<<"end point 5 = "<<qvec[qvec.size() - 1]<<endl;
                }
            }
            
        }
        else
        {
            if(abs(imag(s))>0.0230)
            {    
                cout<<"deforming the contour downward first"<<endl;
            

                comp startingpoint = min;//qvec[qvec.size()-1];
                comp endingpoint_for_slanted_maker;
                if(flag_pcut==0)
                {
                    endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.35*abs(imag(qc2));//abs(imag(pcutpvec[lowest_p_index]));
                    highest_index_xval = real(pcutpvec[highest_p_index]);
                }
                else
                {
                    endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.35*abs(imag(qc2));//*abs(imag(pcutmvec[lowest_m_index]));
                    highest_index_xval = real(pcutpvec[highest_m_index]);
                }
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,2.0*points/5.0);


                startingpoint = qvec[qvec.size()-1];
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);

                startingpoint = qvec[qvec.size()-1];
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = max;
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
            }
            else
            {
                cout<<"deforming the contour backward first 1 "<<endl;
                comp q = pmom(s,sigmab(a,m),m);
                comp qc2 = q_c2(s,a,m,eps);
                comp pplus = pcut_plus(-1.0,s,q,m,eps);
                comp pminus = pcut_plus(+1.0,s,q,m,eps);
                cout<<"pplus = "<<pplus<<'\t'<<"pminus = "<<pminus<<endl;
                cout<<"pplus-pminus = "<<0.5*real(pplus - pminus)<<endl;
                comp qc1 = q_c1(s,a,m,eps);
                double rad = real(qc1);

                comp firstnode = min;
                comp secondnode = -2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2*qc2/qc1);
                line_maker_with_weights(qvec,weights,firstnode,secondnode,2.0*points/10.0);

                //secondnode - thirdnode
                comp thirdnode = 0.5*real(pplus - pminus) - abs(qc2)*ii;
                line_maker_with_weights(qvec,weights,secondnode,thirdnode,2.0*points/10.0);

                //thirdnode - fourthnode 
                comp fourthnode = (2.0/3.0)*(1.0 - 1.68*ii)*abs(qc2);//(2.0/3.0)*(1.0 - 2.0*ii)*abs(qc2);
                line_maker_with_weights(qvec,weights,thirdnode,fourthnode,2.0*points/10.0);

                //fourthnode - fifthnode
                comp fifthnode = (3.0/2.0)*abs(qc2);
                line_maker_with_weights(qvec,weights,fourthnode,fifthnode,2.0*points/10.0);

                //fifthnode - sixthnode 
                comp sixthnode = max;
                line_maker_with_weights(qvec,weights,fifthnode,sixthnode,2.0*points/10.0);
                //////////////////////////////////////////////////////
            }
        }

    }
    else
    {
        double lowest_index_yval = 0.0;
        comp endingpoint_for_slanted_maker;
        if(flag_pcut==0)
        {
            lowest_index_yval = imag(pcutpvec[lowest_p_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutpvec[lowest_p_index]) - ii*0.35*abs(imag(qc2));//*abs(imag(pcutpvec[lowest_p_index]));
            highest_index_xval = real(pcutpvec[highest_p_index]);
        }
        else
        {
            lowest_index_yval = imag(pcutmvec[lowest_m_index]);
            endingpoint_for_slanted_maker = ii*imag(pcutmvec[lowest_m_index]) - ii*0.35*abs(imag(qc2));//*abs(imag(pcutmvec[lowest_m_index]));
            highest_index_xval = real(pcutpvec[highest_m_index]);
        }


        if(abs(lowest_index_yval)<1.0e-2)
        {
            cout<<"still deforming the contour downward first"<<endl;
            comp startingpoint = min;

            if(imag(endingpoint_for_slanted_maker)<0.0)
            {
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,2.0*points/5.0);

                startingpoint = qvec[qvec.size()-1];
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = ii*imag(startingpoint) + abs(highest_index_xval) + qvec_r;
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);

            
                startingpoint = qvec[qvec.size()-1];
                //cout<<"last point after step = "<<startingpoint<<endl;
                endingpoint_for_slanted_maker = max;

                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,1.5*points/5.0);
            }
            else 
            {
                endingpoint_for_slanted_maker = max;
                line_maker_with_weights(qvec,weights,startingpoint,endingpoint_for_slanted_maker,points);
            }


        }
        else 
        {
            cout<<"taking the straight line"<<endl;
            //cout<<"prob here with axis_flag="<<axis_flag<<endl;
            comp startingpoint = min;
            comp endingpoint = max;
            double totpoints = points;// + points/10.0;
            delq = abs(real(min) - real(max))/(totpoints);
            //mom_vector_maker_linear_1(qvec,startingpoint,+delq,totpoints,1);
            line_maker_with_weights(qvec,weights,startingpoint,endingpoint,points);
            

        }
    }



    
}

void mom_vector_maker_seba_imspos_2_with_contour47(  vector<comp> &qvec,
                                    vector<comp> &weights, 
                                    comp s,
                                    comp kmin,
                                    comp kmax,
                                    double a,
                                    double m,
                                    double eps,
                                    double eps_for_m2k,
                                    double points,
                                    int &tag1,
                                    int &tag2,
                                    int &switch_for_gvec_fixer   )
{
    //initial setup//
    comp ii = {0.0,1.0};
    comp q = pmom(s,sigmab(a,m),m);
    comp qc2 = q_c2(s,a,m,eps);
    comp pplus = pcut_plus(1.0,s,q,m,eps);
    comp pminus = pcut_plus(-1.0,s,q,m,eps);
    comp qc1 = q_c1(real(s),a,m,eps);
    comp sigb = sigmab(a,m);
    double rad = abs(real(qc1));
    cout<<"radius = "<<rad<<endl;

    //this section is only for looking into if the 
    //pcuts have crossed the real axis or not

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    if(axis_flag==0)
    //if(rad>0.0)
    {
        switch_for_gvec_fixer = 0;
        cout<<"contour chosen"<<endl;
        //firstnode - secondnode
        comp firstnode = kmin;
        comp secondnode = -rad/3.0;//-2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2);
        //1
        line_maker_with_weights(qvec,weights,firstnode,secondnode,1.0*points/4.0);

        //secondnode - thirdnode
        comp sigmac2 = sigc2(s,sigb,m);
        comp p = (q_c1(s,a,m,eps) + qc2)/2.0;//pmom(s,sigmac2,m);
        comp k = q;
        comp alphapk = ((sqrt(s)-omega_comp(p,m)-omega_comp(k,m))*(sqrt(s)-omega_comp(p,m)-omega_comp(k,m)) - p*p - k*k - m*m)/(2.0*p*k);
        //comp zpk = 2.0*s*(sigb + ii*eps) - (s+sigb - m*m)*(s+m*m-sigmac2);
        //comp zpk = 2.0*s*(sigmac2 + ii*eps) - (s+sigmac2 - m*m)*(s+m*m-sigb);
        
        double x0 = -1.0*real(alphapk);
        //double x0 = real(zpk);
        cout<<"sigma c2 = "<<sigmac2<<endl;
        
        cout<<"x0 = "<<x0<<endl;
        //comp thirdnode = -pcut_minus((-x0-0.77)/2.0,s,q,m,eps);
        
        
        //tag1 = qvec.size() - 1;
        //thirdnode - fourthnode 
        comp f1 = pcut_minus((x0+1.0)/2.0,s,q,m,eps);
        
        comp f2 = 0.5*(real(pplus - pminus)) - ii*abs(qc2);

        comp temp_f3_1 = pcut_minus(-1,s,q,m,eps);
        comp temp_f3_2 = pcut_minus(1,s,q,m,eps);
        comp temp_f3_3 = pcut_plus(-1,s,q,m,eps);
        comp temp_f3_4 = pcut_plus(1,s,q,m,eps);
        comp f3 = 0.5*(real(pplus - pminus)) - ii*5.0/4.0*abs(qc2);;
        /*if(imag(temp_f3_1)<imag(temp_f3_2) && imag(temp_f3_1)<imag(temp_f3_3) && imag(temp_f3_1)<imag(temp_f3_4))
        {
            f3 = pcut_minus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_2)<imag(temp_f3_1) && imag(temp_f3_2)<imag(temp_f3_3) && imag(temp_f3_2)<imag(temp_f3_4))
        {
            f3 = pcut_minus(1,s,q,m,eps);
        }
        else if(imag(temp_f3_3)<imag(temp_f3_1) && imag(temp_f3_3)<imag(temp_f3_2) && imag(temp_f3_3)<imag(temp_f3_4))
        {
            f3 = pcut_plus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_4)<imag(temp_f3_1) && imag(temp_f3_4)<imag(temp_f3_2) && imag(temp_f3_4)<imag(temp_f3_3))
        {
            f3 = pcut_plus(1,s,q,m,eps);
        }*/

        comp thirdnode = real(secondnode) + ii*imag((f2 + f3)/2.0);//-pcut_minus(-x0,s,q,m,eps);
        //2
        line_maker_with_weights(qvec,weights,secondnode,thirdnode,1.0*points/4.0);
        tag1 = qvec.size() - 1;
        
        //comp f3 = pcut_minus(1,s,q,m,eps);
        //comp fourthnode = real(thirdnode) + ii*imag((f2 + f3)/2.0);//real(thirdnode) + ii*imag(f1)*0.75;//pcut_minus(2.0*x0/3.0,s,q,m,eps);
        //3
        //line_maker_with_weights(qvec,weights,thirdnode,fourthnode,2.0*points/10.0);

        //tag2 = qvec.size() - 1;
        //fourthnode - fifthnode
        //comp fifthnode = (f2 + f3)/2.0;//real(f1) + ii*imag(f1)*0.75;//real(qc1) + ii*imag(fourthnode);//(2.0/3.0)*(1.0-2.0*ii)*abs(qc2);
        //4
        //line_maker_with_weights(qvec,weights,fourthnode,fifthnode,2.0*points/10.0);

        //fifthnode - sixthnode 
        comp sixthnode = real(q_c1(s,a,m,eps)) + (ii*imag((f2+f3)/2.0));//(3.0/2.0)*abs(qc2);
        //5
        //line_maker_with_weights(qvec,weights,fifthnode,sixthnode,1.0*points/4.0);
        line_maker_with_weights(qvec,weights,thirdnode,sixthnode,1.0*points/4.0);
        tag2 = qvec.size() - 1;
        //sixthnode - seventhnode
        comp seventhnode = kmax;
        //6
        line_maker_with_weights(qvec,weights,sixthnode,seventhnode,1.0*points/4.0);

    }
    else
    {
        switch_for_gvec_fixer = 1;
        cout<<"straight line chosen"<<endl;
        line_maker_with_weights(qvec,weights,kmin,kmax,points);     
    }

}


void rightmost_left_pcut(   vector<comp> &pcutpvec,
                            vector<comp> &pcutmvec,
                            int flag_pcut,
                            comp &value    )
{
    comp ii = {0.0,1.0};
    double real_part_pcut = 0.0;
    double imag_part_pcut = 0.0;
    double temp_real_pcut = 0.0;
    double temp_imag_pcut = 0.0;
    double actual_real_pcut = 0.0;
    double actual_imag_pcut = 0.0;
    //if flag_pcut==0, then pcutm is to the left
    //if flag_pcut==1, then pcutp is to the left
    int somecounter = 0;
    if(flag_pcut==0)
    {
        for(int i=0;i<pcutmvec.size();++i)
        {
            real_part_pcut = real(pcutmvec[i]);
            imag_part_pcut = imag(pcutmvec[i]);

            if(imag_part_pcut<=0.0)
            {
                if(somecounter==0)
                {
                    actual_real_pcut = real(pcutmvec[i]);
                    actual_imag_pcut = imag(pcutmvec[i]);
                    somecounter=1;
                    continue;
                }

                if(real_part_pcut>actual_real_pcut)
                {
                    actual_real_pcut = real_part_pcut;
                    actual_imag_pcut = imag_part_pcut;
                }
                
            }
        }
    }
    else
    {
        for(int i=0;i<pcutpvec.size();++i)
        {
            real_part_pcut = real(pcutpvec[i]);
            imag_part_pcut = imag(pcutpvec[i]);

            if(imag_part_pcut<=0.0)
            {
                if(somecounter==0)
                {
                    actual_real_pcut = real(pcutpvec[i]);
                    actual_imag_pcut = imag(pcutpvec[i]);
                    somecounter=1;
                    continue;
                }

                if(real_part_pcut>actual_real_pcut)
                {
                    actual_real_pcut = real_part_pcut;
                    actual_imag_pcut = imag_part_pcut;
                }

            }
        }
    }

    value = actual_real_pcut + ii*actual_imag_pcut ;
}


void leftmost_left_pcut(   vector<comp> &pcutpvec,
                            vector<comp> &pcutmvec,
                            int flag_pcut,
                            comp &value    )
{
    comp ii = {0.0,1.0};
    double real_part_pcut = 0.0;
    double imag_part_pcut = 0.0;
    double temp_real_pcut = 0.0;
    double temp_imag_pcut = 0.0;
    double actual_real_pcut = 0.0;
    double actual_imag_pcut = 0.0;
    //if flag_pcut==0, then pcutm is to the left
    //if flag_pcut==1, then pcutp is to the left
    int somecounter = 0;
    if(flag_pcut==0)
    {
        for(int i=0;i<pcutmvec.size();++i)
        {
            real_part_pcut = real(pcutmvec[i]);
            imag_part_pcut = imag(pcutmvec[i]);


            if(somecounter==0)
            {
                actual_real_pcut = real_part_pcut;
                actual_imag_pcut = imag_part_pcut;
                somecounter=1;
                continue;
            }

            if(real_part_pcut<actual_real_pcut)
            {
                actual_real_pcut = real_part_pcut;
                actual_imag_pcut = imag_part_pcut;
            }
            
            
        }
    }
    else
    {
        for(int i=0;i<pcutpvec.size();++i)
        {
            real_part_pcut = real(pcutpvec[i]);
            imag_part_pcut = imag(pcutpvec[i]);


            if(somecounter==0)
            {
                actual_real_pcut = real_part_pcut;
                actual_imag_pcut = imag_part_pcut;
                somecounter=1;
                continue;
            }

            if(real_part_pcut<actual_real_pcut)
            {
                actual_real_pcut = real_part_pcut;
                actual_imag_pcut = imag_part_pcut;
            }
        }
    }

    value = actual_real_pcut + ii*actual_imag_pcut ;
}

void lowest_left_pcut(   vector<comp> &pcutpvec,
                            vector<comp> &pcutmvec,
                            int flag_pcut,
                            comp &value    )
{
    comp ii = {0.0,1.0};
    double real_part_pcut = 0.0;
    double imag_part_pcut = 0.0;
    double temp_real_pcut = 0.0;
    double temp_imag_pcut = 0.0;
    double actual_real_pcut = 0.0;
    double actual_imag_pcut = 0.0;
    //if flag_pcut==0, then pcutm is to the left
    //if flag_pcut==1, then pcutp is to the left
    int somecounter = 0;
    if(flag_pcut==0)
    {
        for(int i=0;i<pcutmvec.size();++i)
        {
            real_part_pcut = real(pcutmvec[i]);
            imag_part_pcut = imag(pcutmvec[i]);


            if(somecounter==0)
            {
                actual_real_pcut = real_part_pcut;
                actual_imag_pcut = imag_part_pcut;
                somecounter=1;
                continue;
            }

            if(imag_part_pcut<actual_imag_pcut)
            {
                actual_real_pcut = real_part_pcut;
                actual_imag_pcut = imag_part_pcut;
            }
            
            
        }
    }
    else
    {
        for(int i=0;i<pcutpvec.size();++i)
        {
            real_part_pcut = real(pcutpvec[i]);
            imag_part_pcut = imag(pcutpvec[i]);


            if(somecounter==0)
            {
                actual_real_pcut = real_part_pcut;
                actual_imag_pcut = imag_part_pcut;
                somecounter=1;
                continue;
            }

            if(imag_part_pcut<actual_imag_pcut)
            {
                actual_real_pcut = real_part_pcut;
                actual_imag_pcut = imag_part_pcut;
            }
        }
    }

    value = actual_real_pcut + ii*actual_imag_pcut ;
}


//the new vector maker to solve for mphib 3d for ims>0
void mom_vector_maker_seba_imspos_3(  vector<comp> &qvec,
                                    vector<comp> &weights, 
                                    comp s,
                                    comp kmin,
                                    comp kmax,
                                    double a,
                                    double m,
                                    double eps,
                                    double eps_for_m2k,
                                    double points,
                                    int &tag1,
                                    int &tag2,
                                    int &switch_for_gvec_fixer   )
{
    comp ii = {0.0,1.0};
    comp sigb = sigmab(a,m);
    comp q = pmom(s,sigb-ii*eps_for_m2k,m);
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp qpls = q_plus(s,a,m,eps);
    comp qmns = q_minus(s,a,m,eps);
    comp diffqplsqmns = abs(imag(-qmns) - imag(qpls))/2.0;
    comp qc1 = q_c1(s,a,m,eps);
    comp qc2 = q_c2(s,a,m,eps);
    rad = real(qc1);
    cout<<"radius = "<<rad<<endl;
    comp delq = {0.0,0.0};

    comp pcutp_close1;
    comp pcutm_close1;
    comp pcutp_close2;
    comp pcutm_close2;
    double x_pcutp;
    double x_pcutm;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int pcutp_i_index1;
    int pcutm_j_index1;
    int pcutp_i_index2;
    int pcutm_j_index2;
    
    
    closest_pcuts(xvec,pcutpvec,pcutmvec, s, q, m, eps, pcutp_i_index1, pcutm_j_index1, pcutp_i_index2, pcutm_j_index2);
    

    cout<<"pcutp_i_index1 = "<<pcutp_i_index1<<endl;
    cout<<"pcutm_j_index1 = "<<pcutm_j_index1<<endl;
    cout<<"pcutp_i_index2 = "<<pcutp_i_index2<<endl;
    cout<<"pcutm_j_index2 = "<<pcutm_j_index2<<endl;

    pcutp_close1 = pcutpvec[pcutp_i_index1];
    pcutm_close1 = pcutmvec[pcutm_j_index1];
    pcutp_close2 = pcutpvec[pcutp_i_index2];
    pcutm_close2 = pcutmvec[pcutm_j_index2];

    int select_pcut; //this selects which of the pcut_close to be chosen, 
                     //0 means the first set
                     //1 means the second set
    if(imag(pcutp_close1)<imag(pcutp_close2)) select_pcut = 0;
    else select_pcut = 1;
    
    comp pcutp_close_selected;
    comp pcutm_close_selected;
    int selected_p_index;
    int selected_m_index;

    if(select_pcut==0)
    {
        pcutp_close_selected = pcutp_close1;
        pcutm_close_selected = pcutm_close1;
        selected_p_index = pcutp_i_index1;
        selected_m_index = pcutm_j_index1;
    }
    else 
    {
        pcutp_close_selected = pcutp_close2;
        pcutm_close_selected = pcutm_close2;
        selected_p_index = pcutp_i_index2;
        selected_m_index = pcutm_j_index2;
    }

    cout<<"selected = "<<pcutp_close_selected<<endl;
    cout<<"selected = "<<pcutm_close_selected<<endl;
    cout<<"selected p index = "<<selected_p_index<<endl;
    cout<<"selected m index = "<<selected_m_index<<endl;

    vector<int> problem_index;

    for(int i=0;i<5;++i)
    {
        problem_index.push_back(i);
        problem_index.push_back(pcutpvec.size()-i-1);
    }

    int problem_index_flag_p;
    //problem_index_flag_p=0 means this index is the end point or near to end point of pcut
    //problem_index_flag_p=1 means this index is not at the end point
    int problem_index_flag_m;

    int problem_index_flag;
    //problem_index_flag=1 means deform the contour first backward and then all the other steps
    //problem_index_flag=0 means deform the contour first downward then all the other steps

    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_p_index==problem_index[i])
        {
            problem_index_flag_p = 0;
            break;
        } 
        else problem_index_flag_p = 1;
    }
    for(int i=0;i<problem_index.size();++i)
    {
        if(selected_m_index==problem_index[i])
        {
            problem_index_flag_m = 0;
            break;
        } 
        else problem_index_flag_m = 1;
    }

    if(problem_index_flag_p==0 && problem_index_flag_m==0) problem_index_flag = 0;
    else problem_index_flag = 1;

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);
    
    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    double yvalp,yvalm;
    int lowest_p_index, lowest_m_index;
    lowest_yval_pcuts(pcutpvec,pcutmvec,yvalp,yvalm,lowest_p_index,lowest_m_index);

    double xvalp,xvalm;
    int highest_p_index, highest_m_index;
    highest_xval_pcuts(pcutpvec,pcutmvec,xvalp,xvalm,highest_p_index,highest_m_index);

    double highest_index_xval = 0.0;

    if(flag_pcut==0)
    cout<<"pcutp is to the right"<<endl;
    else 
    cout<<"pcutm is to the right"<<endl;
    
    cout<<"axis_flag = "<<axis_flag<<endl;
    if(axis_flag==0)
    cout<<"cut to the right crosses real line"<<endl;
    else
    cout<<"cut to the right doesnot cross real line"<<endl;

    cout<<"pcutp highest xval index = "<<highest_p_index<<endl;
    cout<<"pcutp highest xval = "<<pcutpvec[highest_p_index]<<endl;
    cout<<"pcutm highest xval index = "<<highest_m_index<<endl;
    cout<<"pcutm highest xval = "<<pcutmvec[highest_m_index]<<endl;

    if(axis_flag==0)
    {
        cout<<"deforming the contour backward first 1 "<<endl;
        comp q = pmom(s,sigmab(a,m),m);
        comp qc2 = q_c2(s,a,m,eps);
        comp pplus = pcut_plus(-1.0,s,q,m,eps);
        comp pminus = pcut_plus(+1.0,s,q,m,eps);
        cout<<"pplus = "<<pplus<<'\t'<<"pminus = "<<pminus<<endl;
        cout<<"pplus-pminus = "<<0.5*real(pplus - pminus)<<endl;
        comp qc1 = q_c1(s,a,m,eps);
        double rad = real(qc1);

        comp firstnode = kmin;
        comp secondnode;
        if(flag_pcut==0)
        {
            secondnode = pcutmvec[pcutmvec.size()/2];
        }
        else
        {
            secondnode = pcutpvec[pcutpvec.size()/2];
        }

        line_maker_with_weights(qvec,weights,firstnode,secondnode/2.0,points);

        comp thirdnode;
        if(flag_pcut==0)
        {
            thirdnode = pcutmvec[9*pcutmvec.size()/10];
        }
        else
        {
            thirdnode = pcutpvec[9*pcutpvec.size()/10];
        }

        line_maker_with_weights(qvec,weights,qvec[qvec.size()-1],thirdnode,points);

        tag1 = qvec.size() - 1;

        comp tempforthnode1;
        if(flag_pcut==0)
        {
            tempforthnode1 = pcutmvec[lowest_m_index];
        }
        else
        {
            tempforthnode1 = pcutpvec[lowest_p_index];
        }

        comp tempforthnode2;
        if(flag_pcut==0)
        {
            tempforthnode2 = pcutpvec[lowest_p_index];
        }
        else
        {
            tempforthnode2 = pcutmvec[lowest_m_index];
        }

        //comp forthnode = real(thirdnode) + ii*0.45*imag((tempforthnode1 + tempforthnode2)/2.0);
        comp forthnode = real(thirdnode) + ii*imag((tempforthnode1));// + tempforthnode2)/2.0);

        line_maker_with_weights(qvec,weights,thirdnode,forthnode,points);

        comp fifthnode;

        comp pth = phibthreshold(a,m);
        double multiplier = (1.0 - real(s)/real(pth));

        if(flag_pcut==0)
        {
            fifthnode = pcutpvec[multiplier*9.5*pcutpvec.size()/10];
        }
        else
        {
            fifthnode = pcutmvec[multiplier*9.5*pcutmvec.size()/10];
        }

        line_maker_with_weights(qvec,weights,forthnode,fifthnode,points);

        tag2 = qvec.size() - 1;

        comp tempsixthnode;

        if(flag_pcut==0)
        {
            tempsixthnode = pcutpvec[highest_p_index];
        }
        else
        {
            tempsixthnode = pcutmvec[highest_m_index];
        }

        comp sixthnode = real(tempsixthnode) + ii*imag(fifthnode);
        line_maker_with_weights(qvec,weights,fifthnode,sixthnode,points);

        comp seventhnode = kmax;

        line_maker_with_weights(qvec,weights,sixthnode,seventhnode,points);


         




    }
}

void mom_vector_maker_seba_imspos_4(  vector<comp> &qvec,
                                    vector<comp> &weights, 
                                    comp s,
                                    comp kmin,
                                    comp kmax,
                                    double a,
                                    double m,
                                    double eps,
                                    double eps_for_m2k,
                                    double points,
                                    int &tag1,
                                    int &tag2,
                                    int &switch_for_gvec_fixer   )
{
    //initial setup//
    comp ii = {0.0,1.0};
    comp q = pmom(s,sigmab(a,m),m);
    comp qc2 = q_c2(s,a,m,eps);
    comp pplus = pcut_plus(1.0,s,q,m,eps);
    comp pminus = pcut_plus(-1.0,s,q,m,eps);
    comp qc1 = q_c1(real(s),a,m,eps);
    comp sigb = sigmab(a,m);
    double rad = abs(real(qc1));
    cout<<"radius = "<<rad<<endl;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);

    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    if(axis_flag==0)
    {
        switch_for_gvec_fixer = 0;
        cout<<"contour chosen"<<endl;
        
        //firstnode - secondnode
        comp firstnode = kmin;

        comp radend;
        comp radstart;
        rightmost_left_pcut(pcutpvec,pcutmvec,flag_pcut,radend);
        leftmost_left_pcut(pcutpvec,pcutmvec,flag_pcut,radstart);

        double newrad = (real(radend) + real(radstart))/2.0;

        comp secondnode = newrad;//-rad/3.0;//-2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2);
        //1
        line_maker_with_weights(qvec,weights,firstnode,secondnode,1.0*points/4.0);

        //secondnode - thirdnode
        comp sigmac2 = sigc2(s,sigb,m);
        comp p = (q_c1(s,a,m,eps) + qc2)/2.0;//pmom(s,sigmac2,m);
        comp k = q;
        comp alphapk = ((sqrt(s)-omega_comp(p,m)-omega_comp(k,m))*(sqrt(s)-omega_comp(p,m)-omega_comp(k,m)) - p*p - k*k - m*m)/(2.0*p*k);
        //comp zpk = 2.0*s*(sigb + ii*eps) - (s+sigb - m*m)*(s+m*m-sigmac2);
        //comp zpk = 2.0*s*(sigmac2 + ii*eps) - (s+sigmac2 - m*m)*(s+m*m-sigb);
        
        double x0 = -1.0*real(alphapk);
        //double x0 = real(zpk);
        cout<<"sigma c2 = "<<sigmac2<<endl;
        
        cout<<"x0 = "<<x0<<endl;
        //comp thirdnode = -pcut_minus((-x0-0.77)/2.0,s,q,m,eps);
        
        
        //tag1 = qvec.size() - 1;
        //thirdnode - fourthnode 
        comp f1 = pcut_minus((x0+1.0)/2.0,s,q,m,eps);
        
        comp f2 = 0.5*(real(pplus - pminus)) - ii*abs(qc2);

        comp temp_f3_1 = pcut_minus(-1,s,q,m,eps);
        comp temp_f3_2 = pcut_minus(1,s,q,m,eps);
        comp temp_f3_3 = pcut_plus(-1,s,q,m,eps);
        comp temp_f3_4 = pcut_plus(1,s,q,m,eps);
        comp f3 = 0.5*(real(pplus - pminus)) - ii*5.0/4.0*abs(qc2);;
        /*if(imag(temp_f3_1)<imag(temp_f3_2) && imag(temp_f3_1)<imag(temp_f3_3) && imag(temp_f3_1)<imag(temp_f3_4))
        {
            f3 = pcut_minus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_2)<imag(temp_f3_1) && imag(temp_f3_2)<imag(temp_f3_3) && imag(temp_f3_2)<imag(temp_f3_4))
        {
            f3 = pcut_minus(1,s,q,m,eps);
        }
        else if(imag(temp_f3_3)<imag(temp_f3_1) && imag(temp_f3_3)<imag(temp_f3_2) && imag(temp_f3_3)<imag(temp_f3_4))
        {
            f3 = pcut_plus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_4)<imag(temp_f3_1) && imag(temp_f3_4)<imag(temp_f3_2) && imag(temp_f3_4)<imag(temp_f3_3))
        {
            f3 = pcut_plus(1,s,q,m,eps);
        }*/

        comp thirdnode = real(secondnode) + ii*0.85*imag((f2 + f3)/2.0);//-pcut_minus(-x0,s,q,m,eps);
        //2
        line_maker_with_weights(qvec,weights,secondnode,thirdnode,1.0*points/4.0);
        tag1 = qvec.size() - 1;
        
        //comp f3 = pcut_minus(1,s,q,m,eps);
        //comp fourthnode = real(thirdnode) + ii*imag((f2 + f3)/2.0);//real(thirdnode) + ii*imag(f1)*0.75;//pcut_minus(2.0*x0/3.0,s,q,m,eps);
        //3
        //line_maker_with_weights(qvec,weights,thirdnode,fourthnode,2.0*points/10.0);

        //tag2 = qvec.size() - 1;
        //fourthnode - fifthnode
        //comp fifthnode = (f2 + f3)/2.0;//real(f1) + ii*imag(f1)*0.75;//real(qc1) + ii*imag(fourthnode);//(2.0/3.0)*(1.0-2.0*ii)*abs(qc2);
        //4
        //line_maker_with_weights(qvec,weights,fourthnode,fifthnode,2.0*points/10.0);

        //fifthnode - sixthnode 
        comp sixthnode = real(q_c1(s,a,m,eps)) + (ii*imag((f2+f3)/2.0));//(3.0/2.0)*abs(qc2);
        //5
        //line_maker_with_weights(qvec,weights,fifthnode,sixthnode,1.0*points/4.0);
        line_maker_with_weights(qvec,weights,thirdnode,sixthnode,1.0*points/4.0);
        tag2 = qvec.size() - 1;
        //sixthnode - seventhnode
        comp seventhnode = kmax;
        //6
        line_maker_with_weights(qvec,weights,sixthnode,seventhnode,1.0*points/4.0);

    }
    else
    {
        switch_for_gvec_fixer = 1;
        cout<<"straight line chosen"<<endl;
        line_maker_with_weights(qvec,weights,kmin,kmax,points);     
    }

}

void mom_vector_maker_seba_imspos_5(  vector<comp> &qvec,
                                    vector<comp> &weights, 
                                    comp s,
                                    comp kmin,
                                    comp kmax,
                                    double a,
                                    double m,
                                    double eps,
                                    double eps_for_m2k,
                                    double points,
                                    int &tag1,
                                    int &tag2,
                                    int &switch_for_gvec_fixer   )
{
    //initial setup//
    comp ii = {0.0,1.0};
    comp q = pmom(s,sigmab(a,m),m);
    comp qc2 = q_c2(s,a,m,eps);
    comp pplus = pcut_plus(1.0,s,q,m,eps);
    comp pminus = pcut_plus(-1.0,s,q,m,eps);
    comp qc1 = q_c1(real(s),a,m,eps);
    comp sigb = sigmab(a,m);
    double rad = abs(real(qc1));
    //cout<<"radius = "<<rad<<endl;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);

    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    if(axis_flag==0)
    {
        switch_for_gvec_fixer = 0;
        //cout<<"contour chosen"<<endl;
        
        //firstnode - secondnode
        comp firstnode = kmin;

        comp radend;
        comp radstart;
        rightmost_left_pcut(pcutpvec,pcutmvec,flag_pcut,radend);
        leftmost_left_pcut(pcutpvec,pcutmvec,flag_pcut,radstart);

        double newrad = (real(radend) + real(radstart))/2.0;

        comp secondnode = newrad;//-rad/3.0;//-2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2);
        //1
        line_maker_with_weights(qvec,weights,firstnode,secondnode,1.0*points/4.0);

        //secondnode - thirdnode
        comp sigmac2 = sigc2(s,sigb,m);
        comp p = (q_c1(s,a,m,eps) + qc2)/2.0;//pmom(s,sigmac2,m);
        comp k = q;
        comp alphapk = ((sqrt(s)-omega_comp(p,m)-omega_comp(k,m))*(sqrt(s)-omega_comp(p,m)-omega_comp(k,m)) - p*p - k*k - m*m)/(2.0*p*k);
        //comp zpk = 2.0*s*(sigb + ii*eps) - (s+sigb - m*m)*(s+m*m-sigmac2);
        //comp zpk = 2.0*s*(sigmac2 + ii*eps) - (s+sigmac2 - m*m)*(s+m*m-sigb);
        
        double x0 = -1.0*real(alphapk);
        //double x0 = real(zpk);
        //cout<<"sigma c2 = "<<sigmac2<<endl;
        
        //cout<<"x0 = "<<x0<<endl;
        //comp thirdnode = -pcut_minus((-x0-0.77)/2.0,s,q,m,eps);
        
        
        //tag1 = qvec.size() - 1;
        //thirdnode - fourthnode 
        comp f1 = pcut_minus((x0+1.0)/2.0,s,q,m,eps);
        
        comp f2 = 0.5*(real(pplus - pminus)) - ii*abs(qc2);

        comp temp_f3_1 = pcut_minus(-1,s,q,m,eps);
        comp temp_f3_2 = pcut_minus(1,s,q,m,eps);
        comp temp_f3_3 = pcut_plus(-1,s,q,m,eps);
        comp temp_f3_4 = pcut_plus(1,s,q,m,eps);
        comp f3 = 0.5*(real(pplus - pminus)) - ii*5.0/4.0*abs(qc2);;
        /*if(imag(temp_f3_1)<imag(temp_f3_2) && imag(temp_f3_1)<imag(temp_f3_3) && imag(temp_f3_1)<imag(temp_f3_4))
        {
            f3 = pcut_minus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_2)<imag(temp_f3_1) && imag(temp_f3_2)<imag(temp_f3_3) && imag(temp_f3_2)<imag(temp_f3_4))
        {
            f3 = pcut_minus(1,s,q,m,eps);
        }
        else if(imag(temp_f3_3)<imag(temp_f3_1) && imag(temp_f3_3)<imag(temp_f3_2) && imag(temp_f3_3)<imag(temp_f3_4))
        {
            f3 = pcut_plus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_4)<imag(temp_f3_1) && imag(temp_f3_4)<imag(temp_f3_2) && imag(temp_f3_4)<imag(temp_f3_3))
        {
            f3 = pcut_plus(1,s,q,m,eps);
        }*/

        comp lowest_left;
        lowest_left_pcut(pcutpvec,pcutmvec,flag_pcut,lowest_left);

        //comp thirdnode = real(secondnode) + ii*0.85*imag((f2 + f3)/2.0);//-pcut_minus(-x0,s,q,m,eps);
        comp thirdnode = real(secondnode) + ii*imag(lowest_left) - ii*0.01;//-pcut_minus(-x0,s,q,m,eps);
        
        
        //2
        line_maker_with_weights(qvec,weights,secondnode,thirdnode,1.0*points/4.0);
        tag1 = qvec.size() - 1;
        
        //comp f3 = pcut_minus(1,s,q,m,eps);
        //comp fourthnode = real(thirdnode) + ii*imag((f2 + f3)/2.0);//real(thirdnode) + ii*imag(f1)*0.75;//pcut_minus(2.0*x0/3.0,s,q,m,eps);
        //3
        //line_maker_with_weights(qvec,weights,thirdnode,fourthnode,2.0*points/10.0);

        //tag2 = qvec.size() - 1;
        //fourthnode - fifthnode
        //comp fifthnode = (f2 + f3)/2.0;//real(f1) + ii*imag(f1)*0.75;//real(qc1) + ii*imag(fourthnode);//(2.0/3.0)*(1.0-2.0*ii)*abs(qc2);
        //4
        //line_maker_with_weights(qvec,weights,fourthnode,fifthnode,2.0*points/10.0);

        //fifthnode - sixthnode 
        //comp sixthnode = real(q_c1(s,a,m,eps)) + (ii*imag((f2+f3)/2.0));//(3.0/2.0)*abs(qc2);
        comp sixthnode = real(q_c1(s,a,m,eps)) + (ii*imag(thirdnode));//(3.0/2.0)*abs(qc2);
        
        
        //5
        //line_maker_with_weights(qvec,weights,fifthnode,sixthnode,1.0*points/4.0);
        line_maker_with_weights(qvec,weights,thirdnode,sixthnode,1.0*points/4.0);
        tag2 = qvec.size() - 1;
        //sixthnode - seventhnode
        comp seventhnode = kmax;
        //6
        line_maker_with_weights(qvec,weights,sixthnode,seventhnode,1.0*points/4.0);

    }
    else
    {
        switch_for_gvec_fixer = 1;
        //cout<<"straight line chosen"<<endl;
        line_maker_with_weights(qvec,weights,kmin,kmax,points);     
    }

}

//this is extension of the imspos_5, this works for large a's,
//and should work for low a's aswell. We changed the portion 
//where the second line comes down instead to an arbitrary 0.01 imag unit 
//to the midpoint of the lowest branchpoint and the lower point between
//the two cuts

void mom_vector_maker_seba_imspos_5_ext1(  vector<comp> &qvec,
                                    vector<comp> &weights, 
                                    comp s,
                                    comp kmin,
                                    comp kmax,
                                    double a,
                                    double m,
                                    double eps,
                                    double eps_for_m2k,
                                    double points,
                                    int &tag1,
                                    int &tag2,
                                    int &switch_for_gvec_fixer   )
{
    //initial setup//
    comp ii = {0.0,1.0};
    comp q = pmom(s,sigmab(a,m),m);
    comp qc2 = q_c2(s,a,m,eps);
    comp pplus = pcut_plus(1.0,s,q,m,eps);
    comp pminus = pcut_plus(-1.0,s,q,m,eps);
    comp qc1 = q_c1(real(s),a,m,eps);
    comp sigb = sigmab(a,m);
    double rad = abs(real(qc1));
    //cout<<"radius = "<<rad<<endl;

    double initialx = -1.0;
    double finalx = 1.0;
    double size = 10000.0;
    double delx = abs(initialx-finalx)/size;

    vector<double> xvec;
    vector<comp> pcutpvec;
    vector<comp> pcutmvec;

    for(int i=0;i<(int) size+1;++i)
    {
        double x = initialx + i*delx;

        xvec.push_back(x);
        comp pcutp = pcut_plus(x,s,q,m,eps);
        comp pcutm = pcut_minus(x,s,q,m,eps);

        pcutpvec.push_back(pcutp);
        pcutmvec.push_back(pcutm);
    }
    pcut_fixer(pcutpvec,pcutmvec);

    int flag_pcut; 
    
    check_pcut_leftright(pcutpvec,pcutmvec,flag_pcut);

    //flag_pcut=0 means pcutp to the right
    //flag_pcut=1 means pcutm to the right

    int axis_flag;
    //axis_flag=0 means the cut crosses the real axis
    //axis_flag=1 means the cut does not cross the real axis

    pcutvec_realaxis_crossing_checker_1(pcutpvec,pcutmvec,axis_flag);
    
    if(axis_flag==0)
    {
        switch_for_gvec_fixer = 0;
        //cout<<"contour chosen"<<endl;
        
        //firstnode - secondnode
        comp firstnode = kmin;

        comp radend;
        comp radstart;
        rightmost_left_pcut(pcutpvec,pcutmvec,flag_pcut,radend);
        leftmost_left_pcut(pcutpvec,pcutmvec,flag_pcut,radstart);

        double newrad = (real(radend) + real(radstart))/2.0;

        comp secondnode = newrad;//-rad/3.0;//-2.0/(3.0*sqrt(2))*(1.0 + ii)*abs(qc2);
        //1
        line_maker_with_weights(qvec,weights,firstnode,secondnode,1.0*points/4.0);

        //secondnode - thirdnode
        comp sigmac2 = sigc2(s,sigb,m);
        comp p = (q_c1(s,a,m,eps) + qc2)/2.0;//pmom(s,sigmac2,m);
        comp k = q;
        comp alphapk = ((sqrt(s)-omega_comp(p,m)-omega_comp(k,m))*(sqrt(s)-omega_comp(p,m)-omega_comp(k,m)) - p*p - k*k - m*m)/(2.0*p*k);
        //comp zpk = 2.0*s*(sigb + ii*eps) - (s+sigb - m*m)*(s+m*m-sigmac2);
        //comp zpk = 2.0*s*(sigmac2 + ii*eps) - (s+sigmac2 - m*m)*(s+m*m-sigb);
        
        double x0 = -1.0*real(alphapk);
        //double x0 = real(zpk);
        //cout<<"sigma c2 = "<<sigmac2<<endl;
        
        //cout<<"x0 = "<<x0<<endl;
        //comp thirdnode = -pcut_minus((-x0-0.77)/2.0,s,q,m,eps);
        
        
        //tag1 = qvec.size() - 1;
        //thirdnode - fourthnode 
        comp f1 = pcut_minus((x0+1.0)/2.0,s,q,m,eps);
        
        comp f2 = 0.5*(real(pplus - pminus)) - ii*abs(qc2);

        comp temp_f3_1 = pcut_minus(-1,s,q,m,eps);
        comp temp_f3_2 = pcut_minus(1,s,q,m,eps);
        comp low_f3_1;
        if(imag(temp_f3_1)<imag(temp_f3_2))
        {
            low_f3_1 = temp_f3_1;
        }
        else 
        {
            low_f3_1 = temp_f3_2;
        }
        comp temp_f3_3 = pcut_plus(-1,s,q,m,eps);
        comp temp_f3_4 = pcut_plus(1,s,q,m,eps);
        comp low_f3_2;
        if(imag(temp_f3_3)<imag(temp_f3_4))
        {
            low_f3_2 = temp_f3_3;
        }
        else 
        {
            low_f3_2 = temp_f3_4;
        }
        comp low_f3;
        if(imag(low_f3_1)<imag(low_f3_2))
        {
            low_f3 = low_f3_1;
        }
        else 
        {
            low_f3 = low_f3_2;
        }

        //cout<<"branchcut endpoints 1 = "<<temp_f3_1<<endl;
        //cout<<"branchcut endpoints 2 = "<<temp_f3_2<<endl;
        //cout<<"branchcut endpoints 3 = "<<temp_f3_3<<endl;
        //cout<<"branchcut endpoints 4 = "<<temp_f3_4<<endl;
        //cout<<"chosen = "<<low_f3<<endl;

        comp f3 = 0.5*(real(pplus - pminus)) - ii*5.0/4.0*abs(qc2);;
        /*if(imag(temp_f3_1)<imag(temp_f3_2) && imag(temp_f3_1)<imag(temp_f3_3) && imag(temp_f3_1)<imag(temp_f3_4))
        {
            f3 = pcut_minus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_2)<imag(temp_f3_1) && imag(temp_f3_2)<imag(temp_f3_3) && imag(temp_f3_2)<imag(temp_f3_4))
        {
            f3 = pcut_minus(1,s,q,m,eps);
        }
        else if(imag(temp_f3_3)<imag(temp_f3_1) && imag(temp_f3_3)<imag(temp_f3_2) && imag(temp_f3_3)<imag(temp_f3_4))
        {
            f3 = pcut_plus(-1,s,q,m,eps);
        }
        else if(imag(temp_f3_4)<imag(temp_f3_1) && imag(temp_f3_4)<imag(temp_f3_2) && imag(temp_f3_4)<imag(temp_f3_3))
        {
            f3 = pcut_plus(1,s,q,m,eps);
        }*/

        comp lowest_left;
        lowest_left_pcut(pcutpvec,pcutmvec,flag_pcut,lowest_left);

        comp low_qc1;
        if(imag(qc1)<imag(-qc1))
        {
            low_qc1 = qc1;
        }
        else 
        {
            low_qc1 = -qc1;
        }
        
        comp low_qc2;
        if(imag(qc2)<imag(-qc2))
        {
            low_qc2 = qc2;
        }
        else 
        {
            low_qc2 = -qc2;
        }

        //cout<<"qc 1 = "<<qc1<<endl;
        //cout<<"low qc1 = "<<low_qc1<<endl;

        //cout<<"qc 2 = "<<qc2<<endl;
        //cout<<"low qc2 = "<<low_qc2<<endl;
        //comp thirdnode = real(secondnode) + ii*0.85*imag((f2 + f3)/2.0);//-pcut_minus(-x0,s,q,m,eps);
        comp thirdnode = real(secondnode) + ii*imag((low_f3 + low_qc2)/2.0);//real(secondnode) + ii*imag(lowest_left) - ii*0.01;//-pcut_minus(-x0,s,q,m,eps);
        
        
        //2
        line_maker_with_weights(qvec,weights,secondnode,thirdnode,1.0*points/4.0);
        tag1 = qvec.size() - 1;
        
        //comp f3 = pcut_minus(1,s,q,m,eps);
        //comp fourthnode = real(thirdnode) + ii*imag((f2 + f3)/2.0);//real(thirdnode) + ii*imag(f1)*0.75;//pcut_minus(2.0*x0/3.0,s,q,m,eps);
        //3
        //line_maker_with_weights(qvec,weights,thirdnode,fourthnode,2.0*points/10.0);

        //tag2 = qvec.size() - 1;
        //fourthnode - fifthnode
        //comp fifthnode = (f2 + f3)/2.0;//real(f1) + ii*imag(f1)*0.75;//real(qc1) + ii*imag(fourthnode);//(2.0/3.0)*(1.0-2.0*ii)*abs(qc2);
        //4
        //line_maker_with_weights(qvec,weights,fourthnode,fifthnode,2.0*points/10.0);

        //fifthnode - sixthnode 
        //comp sixthnode = real(q_c1(s,a,m,eps)) + (ii*imag((f2+f3)/2.0));//(3.0/2.0)*abs(qc2);
        comp sixthnode = real(q_c1(s,a,m,eps)) + (ii*imag(thirdnode));//(3.0/2.0)*abs(qc2);
        
        
        //5
        //line_maker_with_weights(qvec,weights,fifthnode,sixthnode,1.0*points/4.0);
        line_maker_with_weights(qvec,weights,thirdnode,sixthnode,1.0*points/4.0);
        tag2 = qvec.size() - 1;
        //sixthnode - seventhnode
        comp seventhnode = kmax;
        //6
        line_maker_with_weights(qvec,weights,sixthnode,seventhnode,1.0*points/4.0);

    }
    else
    {
        switch_for_gvec_fixer = 1;
        //cout<<"straight line chosen"<<endl;
        line_maker_with_weights(qvec,weights,kmin,kmax,points);     
    }

}




//this gvec_fixer_2 works, this rescales the tags for 
//any N based on the temp_qvec size of 5000 it generates 
//inside and looks for the largest gap for N of 5000 
//and then rescales the tags based on this
void Gvec_fixer_2(    Eigen::VectorXcd &Gvec,
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
    int points1 = 5000;
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

    cout<<"ran here"<<endl;

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

    double temp_max2 = 0.0;
    int temp_ind2 = 0.0;

    for(int i=0;i<gap.size();++i)
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
    }

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

//this gvec_fixer is based on contour_imspos_2_contour47
//where we have one discontinuity to left 
//and one similarity to right
void Gvec_fixer_4(    Eigen::VectorXcd &Gvec,
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
    //mom_vector_maker_seba_imspos_2(temp_qvec,temp_weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,temp_tag1,temp_tag2,switch_for_gvec_fixer);
    mom_vector_maker_seba_imspos_2_with_contour47(temp_qvec,temp_weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,tag1,tag2,switch_for_gvec_fixer);
                                          
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

//this gvec_fixer is based on contour_imspos_5
//where we have one discontinuity to left 
//and one similarity to right
void Gvec_fixer_5(    Eigen::VectorXcd &Gvec,
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
    //mom_vector_maker_seba_imspos_2(temp_qvec,temp_weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,temp_tag1,temp_tag2,switch_for_gvec_fixer);
    mom_vector_maker_seba_imspos_5(temp_qvec,temp_weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,temp_tag1,temp_tag2,switch_for_gvec_fixer);
                                    
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

    //cout<<"tags before rescaling = "<<temp_ind1<<'\t'<<temp_ind2<<endl;

    comp temp_gvec_tag1_val = temp_Gvec[temp_ind1];
    comp temp_gvec_tag2_val = temp_Gvec[temp_ind2];
    
    int rescaled_tag1 = ceil(temp_ind1*qvec.size()/temp_qvec.size());// + 1 ;
    int rescaled_tag2 = ceil(temp_ind2*qvec.size()/temp_qvec.size());// - 1;
    //int rescaled_tag1 = floor(temp_ind1*qvec.size()/temp_qvec.size());
    //int rescaled_tag2 = floor(temp_ind2*qvec.size()/temp_qvec.size());
    
    
    //cout<<"min ind = "<<rescaled_tag1<<'\t'<<" max ind = "<<rescaled_tag2<<endl;
    //cout<<"gap 1 = "<<temp_max1<<'\t'<<" gap 2 = "<<temp_max2<<endl;
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


//this gvec_fixer is based on contour_imspos_5
//where we have one discontinuity to left 
//and one similarity to right
//this should work without tag fixer 
void Gvec_fixer_6(    Eigen::VectorXcd &Gvec,
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
    int points1 = qvec.size();
    int temp_tag1,temp_tag2;
    //mom_vector_maker_seba_imspos(temp_qvec,temp_weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,temp_tag1,temp_tag2,switch_for_gvec_fixer);
    //mom_vector_maker_47_with_weights(temp_qvec,temp_weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,0.00101,temp_tag1,temp_tag2);
    //mom_vector_maker_seba_imspos_2(temp_qvec,temp_weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,temp_tag1,temp_tag2,switch_for_gvec_fixer);
    mom_vector_maker_seba_imspos_5(temp_qvec,temp_weights,s,0.0,kmax,a,m,eps,eps_for_m2k,(double)points1,temp_tag1,temp_tag2,switch_for_gvec_fixer);
                                    
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

    cout<<"temp ind1 found = "<<temp_ind1<<endl;

    

    double temp_max2 = 0;
    int temp_ind2 = 0;

    for(int i=0;i<gap.size();++i)
    {
        if(gap[i]>temp_max2 && i!=temp_ind1)
        {
            temp_max2 = gap[i];
            temp_ind2 = i;
        }
    }

    cout<<"temp ind2 found = "<<temp_ind2<<endl;

    temp_ind1 = temp_ind1 + 1;


    if(temp_ind1>temp_ind2)
    {
        int sometemp = temp_ind2;
        temp_ind2 = temp_ind1;
        temp_ind1 = sometemp;
    }

    cout<<"tags before rescaling = "<<temp_ind1<<'\t'<<temp_ind2<<endl;

    comp temp_gvec_tag1_val = temp_Gvec[temp_ind1];
    comp temp_gvec_tag2_val = temp_Gvec[temp_ind2];
    
    int rescaled_tag1 = ceil(temp_ind1*qvec.size()/temp_qvec.size());// + 1 ;
    int rescaled_tag2 = ceil(temp_ind2*qvec.size()/temp_qvec.size());// - 1;
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

void tag_fixer_for_gvec_fixer5( Eigen::VectorXcd &Gvec,
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
        
        double first_difference = abs(real(Gvecval) - real(Gs1stsheetnext));
        double second_difference = abs(real(Gvecval) - real(Gs2ndsheetnext));
        
        if(first_difference<=second_difference)
        {
            Gvec[i+1] = GS_pk(s,pnext,k,m,eps);
        }
        else
        {
            Gvec[i+1] = GS_pk_secondsheet(s,pnext,k,m,eps);
            problem_tag = i+1;

        }
        
        
    }

    
}



void contour_for_resonance( vector<comp> &qvec,
                            vector<comp> &weights,
                            comp kmin, 
                            comp kmax,
                            double shift,
                            double qvecpoints   )
{
    
    comp ii = {0.0,1.0};
    comp midpoint = (kmax + kmin)/2.0;
    midpoint = midpoint + ii*shift;
    
    
    line_maker_with_weights(qvec,weights,kmin,midpoint,qvecpoints/2.0);
    line_maker_with_weights(qvec,weights,midpoint,kmax,qvecpoints/2.0);
}


void contour_for_resonance_1( vector<comp> &qvec,
                            vector<comp> &weights,
                            comp q,
                            comp kmin, 
                            comp kmax,
                            double shift,
                            double qvecpoints   )
{
    
    comp ii = {0.0,1.0};
    comp midpoint = (kmax + kmin)/2.0;
    midpoint = midpoint + ii*shift;
    
    comp firstpoint = real(q) - abs(shift);
    //cout<<kmin<<'\t'<<firstpoint<<endl;
    comp secondpoint = real(firstpoint) + ii*imag(q) + ii*shift;
    comp thirdpoint = real(secondpoint) + 2.0*abs(shift) + ii*imag(secondpoint);
    comp forthpoint = real(thirdpoint);


    if(imag(q)<=0.0)
    {
        line_maker_with_weights(qvec,weights,kmin,firstpoint,qvecpoints/5.0);
        line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints/5.0);
        line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints/5.0);
        line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints/5.0);
        line_maker_with_weights(qvec,weights,forthpoint,kmax,qvecpoints/5.0);
    }
    else
    {
        line_maker_with_weights(qvec,weights,kmin,kmax,qvecpoints);
    }
}

void contour_for_resonance_2( vector<comp> &qvec,
                            vector<comp> &weights,
                            comp kmin, 
                            double reshift,
                            double imshift,
                            comp kmax,
                            double qvecpoints   )
{
    
    comp ii = {0.0,1.0};

    comp shift = reshift + ii*imshift;
    
    
    line_maker_with_weights(qvec,weights,kmin,shift,qvecpoints/2.0);
    line_maker_with_weights(qvec,weights,shift,kmax,qvecpoints/2.0);
}


//this should be used with the rotated definition of m2k
void contour_for_resonance_3( vector<comp> &qvec,
                            vector<comp> &weights,
                            comp s,
                            double m,
                            comp kmin, 
                            comp kmax,
                            double eps,
                            double shift,
                            double qvecpoints   )
{
    
    comp ii = {0.0,1.0};


    comp m2k_cut_1 = M2kbranchcut_right_momrep_plus_eps(s,m,eps);
    comp m2k_cut_2 = M2kbranchcut_right_momrep_minus_eps(s,m,eps);
    comp selected_m2k_cut;
    if(real(m2k_cut_1)>real(m2k_cut_2))
    {
        selected_m2k_cut = m2k_cut_1;
    }
    else
    {
        selected_m2k_cut = m2k_cut_2;
    }

    
    
    if(imag(selected_m2k_cut)<0.0)
    {
        comp firstpoint = kmin;
        comp secondpoint = ii*imag(selected_m2k_cut) + ii*shift;

        line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints/4.0);

        comp thirdpoint = ii*imag(secondpoint) + real(selected_m2k_cut) + real(selected_m2k_cut)/2.0;

        line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints/2.0);

        comp forthpoint = kmax;

        line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints/4.0);
    }
    else 
    {
        comp firstpoint = kmin;
        comp secondpoint = ii*shift; 
        
        line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints/4.0);
        
        comp thirdpoint = ii*imag(secondpoint) + real(selected_m2k_cut) + real(selected_m2k_cut)/2.0;

        line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints/2.0);
        
        comp forthpoint = kmax;
        
        line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints/4.0);

    }

}

//this contour is for re(s)<9 and im(s)<0 
void contour_for_resonance_4( vector<comp> &qvec,
                            vector<comp> &weights,
                            comp s,
                            double m,
                            comp kmin, 
                            comp kmax,
                            double eps,
                            double qvecpoints,
                            int &tag_for_m2k   )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);

    vector<comp> m2kcut;
    for(int i=0;i<500;++i)
    {
        double delsigk = abs(1.0)/500;
        double delsig = i*delsigk;
        comp sigk = 4.0 + (comp)delsig;

        comp m2kcut_in_p = pmom(s,sigk,m);
        m2kcut.push_back(m2kcut_in_p);
    }

    comp firstpoint = kmin;
    comp secondpoint = -0.15 -0.35*ii;
    line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints/4.0);

    comp thirdpoint = m2kcut[2*m2kcut.size()/3];
    line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints/4.0);

    tag_for_m2k = qvec.size() - 1;

    comp forthpoint = kmax; 
    line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints/2.0);


}

//this qvec creator is for 1st resonances whose poles 
//lie below real(s)<6 

void contour_for_resonance_5( vector<comp> &qvec,
                            vector<comp> &weights,
                            comp s,
                            double m,
                            comp kmin, 
                            comp kmax,
                            double eps,
                            double qvecpoints,
                            int &tag_for_m2k   )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);

    vector<comp> m2kcut;
    for(int i=0;i<500;++i)
    {
        double delsigk = abs(1.0)/500;
        double delsig = i*delsigk;
        comp sigk = 4.0 + (comp)delsig;

        comp m2kcut_in_p = pmom(s,sigk,m);
        if(imag(m2kcut_in_p)>0.0)
        {
            comp temp = -real(m2kcut_in_p) - ii*imag(m2kcut_in_p);
            m2kcut_in_p = temp;
        }
        m2kcut.push_back(m2kcut_in_p);
    }

    comp firstpoint = kmin;
    comp secondpoint = -0.05 -0.05*ii;
    line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints/4.0);

    comp thirdpoint = secondpoint;//real(secondpoint) + ii*imag(m2kcut[m2kcut.size()/2]);
    //line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints/1.0);
    
    comp forthpoint = m2kcut[m2kcut.size()/2];
    line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints/4.0);
    
    tag_for_m2k = qvec.size() - 1;

    comp fifthpoint = ii*imag(qvec[qvec.size()-1]) + real(kmax);
    line_maker_with_weights(qvec,weights,forthpoint,fifthpoint,qvecpoints/4.0);

    comp sixthpoint = kmax; 
    line_maker_with_weights(qvec,weights,fifthpoint,sixthpoint,qvecpoints/4.0);


}

//this contour is for second sheet of dS the top part is s plane
//this will work with simag>0.0 part only
//does not fully work for all s
void contour_for_resonance_6( vector<comp> &qvec,
                            vector<comp> &weights,
                            double a,
                            comp s,
                            double m,
                            comp kmin,
                            comp kmax,
                            double eps,
                            double qvecpoints,
                            int &tag_for_m2k   )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);

    vector<comp> m2kcut;
    for(int i=0;i<500;++i)
    {
        double delsigk = abs(1.0)/500;
        double delsig = i*delsigk;
        comp sigk = 4.0 + (comp)delsig;

        comp m2kcut_in_p = -pmom(s,sigk,m);
        m2kcut.push_back(m2kcut_in_p);
    }

    comp sigb = sigmab(a,m);
    comp qpole = -pmom(s,sigb,m);
    comp m2kbc = -M2kbranchcut_right_momrep_plus_eps(s,m,eps);
    
    double shiftreal = 0.05;
    //if(real(s)>9.0)
    int set = 0;

    if(set==0)
    {
        comp firstpoint = kmin;
        comp secondpoint = real(m2kbc) + 0.005*ii;//-0.005 + 0.005*ii;
        //line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints*2.0/10.0);

        //comp thirdpoint = real(qpole) - shiftreal + ii*imag(qpole);
        
        comp thirdpoint = (qpole + m2kbc)/2.0;
        
        
        //line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints*2.0/10.0);
        line_maker_with_weights(qvec,weights,firstpoint,thirdpoint,qvecpoints*2.0/10.0);

        //comp forthpoint = real(thirdpoint) + ii*imag(m2kcut[2*m2kcut.size()/3])/2.0;
        comp forthpoint = real(thirdpoint) + ii*imag(m2kcut[2*m2kcut.size()/3]);
        line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints*2.0/10.0);

        comp fifthpoint = m2kcut[2*m2kcut.size()/3];
        line_maker_with_weights(qvec,weights,forthpoint,fifthpoint,qvecpoints*2.0/10.0);

        tag_for_m2k = qvec.size() - 1;

        comp sixthpoint = real(kmax);
        line_maker_with_weights(qvec,weights,fifthpoint,sixthpoint,qvecpoints*2.0/10.0);

        comp seventhpoint = kmax;
        line_maker_with_weights(qvec,weights,sixthpoint,seventhpoint,qvecpoints*2.0/10.0);


    }
    else 
    {
        
        comp firstpoint = kmin;
        comp secondpoint = real(m2kbc);// + 0.005*ii;//-0.005 + 0.005*ii;
        line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints*2.0/10.0);

        //comp thirdpoint = real(qpole) - shiftreal + ii*imag(qpole);
        
        comp thirdpoint = (qpole + m2kbc)/2.0;
        
        
        line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints*1.0/10.0);
        //line_maker_with_weights(qvec,weights,firstpoint,thirdpoint,qvecpoints*2.0/10.0);

        //comp forthpoint = real(thirdpoint) + ii*imag(m2kcut[2*m2kcut.size()/3])/2.0;
        comp forthpoint = real(thirdpoint) + ii*imag(m2kcut[2*m2kcut.size()/3]);
        line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints*1.0/10.0);

        comp fifthpoint = m2kcut[2*m2kcut.size()/3];
        line_maker_with_weights(qvec,weights,forthpoint,fifthpoint,qvecpoints*2.0/10.0);

        tag_for_m2k = qvec.size() - 1;

        comp sixthpoint = real(kmax);
        line_maker_with_weights(qvec,weights,fifthpoint,sixthpoint,qvecpoints*2.0/10.0);

        comp seventhpoint = kmax;
        line_maker_with_weights(qvec,weights,sixthpoint,seventhpoint,qvecpoints*2.0/10.0);


    }
    
    
    

}


//this contour is made for the third sheet of dS
//the bottom simag panel 
void contour_for_resonance_7( vector<comp> &qvec,
                            vector<comp> &weights,
                            double a,
                            comp s,
                            double m,
                            comp kmin,
                            comp kmax,
                            double eps,
                            double qvecpoints,
                            int &tag_for_m2k_1,
                            int &tag_for_m2k_2   )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);
    double shiftreal = 0.05;
    double shiftcontour_from_m2kbc = 0.015;

    vector<comp> m2kcut;
    for(int i=0;i<500;++i)
    {
        double delsigk = abs(1.0)/500;
        double delsig = i*delsigk;
        comp sigk = 4.0 + (comp)delsig;

        comp m2kcut_in_p;
        if(imag(s)<0.0)
        m2kcut_in_p = pmom(s,sigk,m);
        else 
        m2kcut_in_p = -pmom(s,sigk,m);

        m2kcut.push_back(m2kcut_in_p);
    }

    comp sigb = sigmab(a,m);
    comp qpole;
    comp m2kbc;

    if(imag(s)<0.0)
    {
        qpole = pmom(s,sigb,m);
        m2kbc = M2kbranchcut_right_momrep_plus_eps(s,m,eps);
    }
    else
    {
        qpole = -pmom(s,sigb,m);
        m2kbc = -M2kbranchcut_right_momrep_plus_eps(s,m,eps);
    }
    //if(real(s)>9.0)
    {
        comp firstpoint = kmin;
        comp secondpoint = -0.05 - 0.05*ii;
        line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints*1.0/10.0);

        
        comp thirdpoint = m2kcut[1*m2kcut.size()/10];//(qpole + m2kbc)/2.0;
        
        line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints*2.0/10.0);
        tag_for_m2k_1 = qvec.size() - 1;

        comp forthpoint = real((qpole + m2kbc)/2.0) + ii*imag(thirdpoint);
        
        line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints*1.0/10.0);

        comp fifthpoint = (qpole + m2kbc)/2.0;
        
        line_maker_with_weights(qvec,weights,forthpoint,fifthpoint,qvecpoints*1.0/10.0);

        comp sixthpoint = real(m2kbc) + ii*shiftcontour_from_m2kbc;
        //line_maker_with_weights(qvec,weights,fifthpoint,sixthpoint,qvecpoints*1.0/10.0);

        comp seventhpoint = (-fifthpoint + m2kbc)/2.0;//(kmin + fifthpoint)/2.0;
        
        //line_maker_with_weights(qvec,weights,sixthpoint,seventhpoint,qvecpoints*1.0/10.0);

        comp eighthpoint = thirdpoint;
        //line_maker_with_weights(qvec,weights,seventhpoint,eighthpoint,qvecpoints*1.0/10.0);
        
        tag_for_m2k_2 = qvec.size() - 1;

        comp ninthpoint = real(kmax) + ii*imag(seventhpoint);
        //line_maker_with_weights(qvec,weights,eighthpoint,ninthpoint,qvecpoints*1.0/10.0);
        
        comp tenthpoint = kmax;
        //line_maker_with_weights(qvec,weights,ninthpoint,tenthpoint,qvecpoints*1.0/10.0);
        

    }

}

//this contour is for sheet -1 to look for mirror poles 
void contour_for_resonance_for_sheet_minus1_imsneg( vector<comp> &qvec,
                            vector<comp> &weights,
                            double a,
                            comp s,
                            double m,
                            comp kmin,
                            comp kmax,
                            double eps,
                            double qvecpoints,
                            int &tag_for_m2k   )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);

    vector<comp> m2kcut;
    for(int i=0;i<500;++i)
    {
        double delsigk = abs(1.0)/500;
        double delsig = i*delsigk;
        comp sigk = 4.0 + (comp)delsig;

        comp m2kcut_in_p = -pmom(real(s)-ii*imag(s),sigk,m);
        m2kcut.push_back(m2kcut_in_p);
    }

    comp sigb = sigmab(a,m);
    comp qpole = -pmom(s,sigb,m);
    comp m2kbc = -M2kbranchcut_right_momrep_plus_eps(s,m,eps);
    
    double shiftreal = 0.05;
    //if(real(s)>9.0)
    int set = 0;

    if(set==0)
    {
        comp firstpoint = kmin;
        comp secondpoint = real(m2kbc) + 0.005*ii;//-0.005 + 0.005*ii;
        //line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints*2.0/10.0);

        //comp thirdpoint = real(qpole) - shiftreal + ii*imag(qpole);
        
        comp thirdpoint = (qpole + m2kbc)/2.0;
        
        
        //line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints*2.0/10.0);
        line_maker_with_weights(qvec,weights,firstpoint,thirdpoint,qvecpoints*2.0/10.0);

        //comp forthpoint = real(thirdpoint) + ii*imag(m2kcut[2*m2kcut.size()/3])/2.0;
        comp forthpoint = real(thirdpoint) + ii*imag(m2kcut[2*m2kcut.size()/3]);
        line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints*2.0/10.0);

        comp fifthpoint = m2kcut[2*m2kcut.size()/3];
        line_maker_with_weights(qvec,weights,forthpoint,fifthpoint,qvecpoints*2.0/10.0);

        tag_for_m2k = qvec.size() - 1;

        comp sixthpoint = real(kmax);
        line_maker_with_weights(qvec,weights,fifthpoint,sixthpoint,qvecpoints*2.0/10.0);

        comp seventhpoint = kmax;
        line_maker_with_weights(qvec,weights,sixthpoint,seventhpoint,qvecpoints*2.0/10.0);


    }
    else 
    {
        
        comp firstpoint = kmin;
        comp secondpoint = real(m2kbc);// + 0.005*ii;//-0.005 + 0.005*ii;
        line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints*2.0/10.0);

        //comp thirdpoint = real(qpole) - shiftreal + ii*imag(qpole);
        
        comp thirdpoint = (qpole + m2kbc)/2.0;
        
        
        line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints*1.0/10.0);
        //line_maker_with_weights(qvec,weights,firstpoint,thirdpoint,qvecpoints*2.0/10.0);

        //comp forthpoint = real(thirdpoint) + ii*imag(m2kcut[2*m2kcut.size()/3])/2.0;
        comp forthpoint = real(thirdpoint) + ii*imag(m2kcut[2*m2kcut.size()/3]);
        line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints*1.0/10.0);

        comp fifthpoint = m2kcut[2*m2kcut.size()/3];
        line_maker_with_weights(qvec,weights,forthpoint,fifthpoint,qvecpoints*2.0/10.0);

        tag_for_m2k = qvec.size() - 1;

        comp sixthpoint = real(kmax);
        line_maker_with_weights(qvec,weights,fifthpoint,sixthpoint,qvecpoints*2.0/10.0);

        comp seventhpoint = kmax;
        line_maker_with_weights(qvec,weights,sixthpoint,seventhpoint,qvecpoints*2.0/10.0);


    }
    
    
    

}


//this contour is for sheet -1 of the three body amplitude
//this is to look for mirror poles 
void contour_for_resonance_for_sheet_minus1_imspos( vector<comp> &qvec,
                            vector<comp> &weights,
                            comp s,
                            double m,
                            comp kmin, 
                            comp kmax,
                            double eps,
                            double qvecpoints,
                            int &tag_for_m2k   )
{
    comp ii = {0.0,1.0};
    double pi = acos(-1.0);

    vector<comp> m2kcut;
    for(int i=0;i<500;++i)
    {
        double delsigk = abs(1.0)/500;
        double delsig = i*delsigk;
        comp sigk = 4.0 + (comp)delsig;

        comp m2kcut_in_p = pmom(real(s) - ii*imag(s),sigk,m);
        m2kcut.push_back(m2kcut_in_p);
    }

    comp firstpoint = kmin;
    comp secondpoint = -0.15 -0.35*ii;
    line_maker_with_weights(qvec,weights,firstpoint,secondpoint,qvecpoints/4.0);

    comp thirdpoint = m2kcut[2*m2kcut.size()/3];
    line_maker_with_weights(qvec,weights,secondpoint,thirdpoint,qvecpoints/4.0);

    tag_for_m2k = qvec.size() - 1;

    comp forthpoint = kmax; 
    line_maker_with_weights(qvec,weights,thirdpoint,forthpoint,qvecpoints/2.0);


}


//this one uses the sigma_vector_maker_linear   //
//and makes the contour. This uses fixed gaps   //
//of 0.01 from the branching points on the real //
//and imaginary axis.                           // 

void sigma_vector_maker_2(  vector<comp> &sigvec,
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double point1,
                            double point2,
                            double eps   )
{
    comp ii = {0.0,1.0};
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);

    if(real(rad)>0.0)
    {
    

    cout<<"sigmin:"<<min<<endl;
    cout<<"sigmax:"<<max<<endl;
    cout<<"rad:"<<rad<<endl;
    cout<<"sigmac1:"<<sigmac1<<endl;
    cout<<"sigmac2:"<<sigmac2<<endl;

    sigma_vector_maker_linear(sigvec,min,sigmac1-0.01,point1,eps,1); // the last 1 is for making vec on real axis
    sigma_vector_maker_linear(sigvec,sigmac1-0.01,sigmac1+rad+0.001,point2,eps,2);
    sigma_vector_maker_linear(sigvec,sigmac1-0.01+ii*(rad+0.001),sigmac1+2.0*(rad+0.001)+ii*(rad+0.001),2.0*point2,eps,1);
    sigma_vector_maker_linear(sigvec,sigmac1+2.0*(rad+0.001)+ii*(rad+0.001),sigmac2+0.001 ,point2,eps,3);
    sigma_vector_maker_linear(sigvec,sigmac2+0.001,max,point2,eps,4);

    cout<<"sigvec created with size = "<<sigvec.size()<<endl;
    }
    else
    {
        double length = real(max - min);
        comp delsig = length/(point1 + 5.0*point2) ;

        for(double tempsig=real(min);tempsig<=real(max);tempsig=tempsig+real(delsig))
        {
            comp sigma = tempsig + ii*eps;

            sigvec.push_back(sigma);
        }
    }
    

}

//this one uses the sigma_vector_maker_linear   //
//and makes the contour. This uses fixed gaps   //
//of 0.0 from the branching points on the real  //
//and imaginary axis.                           // 


void sigma_vector_maker_3(  vector<comp> &sigvec,
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double point1,
                            double point2,
                            double eps   )
{
    comp ii = {0.0,1.0};
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);

    if(real(rad)>0.0)
    {
    

    cout<<"sigmin:"<<min<<endl;
    cout<<"sigmax:"<<max<<endl;
    cout<<"rad:"<<rad<<endl;
    cout<<"sigmac1:"<<sigmac1<<endl;
    cout<<"sigmac2:"<<sigmac2<<endl;

    sigma_vector_maker_linear(sigvec,min,sigmac1,point1,eps,1); // the last 1 is for making vec on real axis
    sigma_vector_maker_linear(sigvec,sigmac1,sigmac1+rad,point2,eps,2);
    sigma_vector_maker_linear(sigvec,sigmac1+ii*rad,sigmac1+2.0*(rad)+ii*(rad),2.0*point2,eps,1);
    sigma_vector_maker_linear(sigvec,sigmac1+2.0*(rad)+ii*(rad),sigmac2 ,point2,eps,3);
    sigma_vector_maker_linear(sigvec,sigmac2,max,point2,eps,4);

    cout<<"sigvec created with size = "<<sigvec.size()<<endl;
    }
    else
    {
        double length = real(max - min);
        comp delsig = length/(point1 + 5.0*point2) ;

        for(double tempsig=real(min);tempsig<=real(max);tempsig=tempsig+real(delsig))
        {
            comp sigma = tempsig + ii*eps;

            sigvec.push_back(sigma);
        }
    }
    

}

//this one uses the sigma_vector_maker_linear   //
//and makes the contour. This uses fixed gaps   //
//of 100*eps from the branching points on the   //
//real  and imaginary axis.                     // 

void sigma_vector_maker_4(  vector<comp> &sigvec,
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double point1,
                            double point2,
                            double eps   )
{
    comp ii = {0.0,1.0};
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);

    if(real(rad)>0.0)
    {
    

    cout<<"sigmin:"<<min<<endl;
    cout<<"sigmax:"<<max<<endl;
    cout<<"rad:"<<rad<<endl;
    cout<<"sigmac1:"<<sigmac1<<endl;
    cout<<"sigmac2:"<<sigmac2<<endl;

    sigma_vector_maker_linear(sigvec,min,sigmac1-100*eps,point1,eps,1); // the last 1 is for making vec on real axis
    sigma_vector_maker_linear(sigvec,sigmac1-100*eps,sigmac1+rad,point2,eps,2);
    sigma_vector_maker_linear(sigvec,sigmac1-100*eps+ii*(rad+100*eps),sigmac1+2.0*(rad+100*eps)+ii*(rad+100*eps),2.0*point2,eps,1);
    sigma_vector_maker_linear(sigvec,sigmac1+2.0*(rad+100*eps)+ii*(rad+100*eps),sigmac2 + 2.0*(100*eps) ,point2,eps,3);
    sigma_vector_maker_linear(sigvec,sigmac2+100*eps,max,point2,eps,4);

    cout<<"sigvec created with size = "<<sigvec.size()<<endl;
    }
    else
    {
        double length = real(max - min);
        comp delsig = length/(point1 + 5.0*point2) ;

        for(double tempsig=real(min);tempsig<=real(max);tempsig=tempsig+real(delsig))
        {
            comp sigma = tempsig + ii*eps;

            sigvec.push_back(sigma);
        }
    }
    

}


//We will make a linear sigma vector. It will have a    //
//starting point, we will provide with the delsig and   //
//number of points. The delsig will have a plus or minus//
//sign with it which tells it to go up or down and left //
//right to form the contour in the complex plane. It    //
//will also have a int,'1' sets it to form the contour  //
//along the real axis and '2' sets it to form the       //
//contour in the imaginary axis.                        //

void sigma_vector_maker_linear_1(   vector<comp> &sigvec,
                                    comp startingpoint,
                                    comp delsig,
                                    double points,
                                    int a               )
{
    comp ii = {0.0,1.0};
    if(a==1)
    {
        for(int i=1;i<=(int)points;++i)
        {
            comp tempsig = startingpoint + ((comp)i)*delsig;
            sigvec.push_back(tempsig);
        }
    }
    else if(a==2)
    {
        for(int i=1;i<=(int)points;++i)
        {
            comp tempsig = startingpoint + ii*((comp)i)*delsig;
            sigvec.push_back(tempsig);
        }
    }

}

//this one makes the contour for each s, we give it an  //
//epsratio, this determines how large the box is around //
//circular cut. The created sigvec do exclude the       //
//minimum value                                         //
void sigma_vector_maker_5(  vector<comp> &sigvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double points1,
                            double points2,
                            double epsratio     )
{
    comp ii = {0.0,1.0};
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp delsig = {0.0,0.0};

    if(real(rad)>0.0)
    {
        /*  this is the initial straight line along the real axis */
        comp startingpoint = min + ii*eps;
        delsig = abs((sigmac1-epsratio*eps)-min)/points1;
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,points1,1);
        
        /*  this is the first contour up in the imaginary axis   */
        startingpoint = sigvec[sigvec.size()-1];
        delsig = abs(rad+((comp)epsratio*eps))/points2;
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,points2,2);

        /*  this is again parallel to the real axis, to the right, this 
            goes for 2*radius and for 2*points2                         */
        startingpoint = sigvec[sigvec.size()-1];
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,2.0*points2,1);

        /*  the contour now comes down, so we flip the sign of the delsig */
        startingpoint = sigvec[sigvec.size()-1];
        sigma_vector_maker_linear_1(sigvec,startingpoint,-delsig,points2,2);

        /*  this is the end part, the contour goes to the left to max */
        startingpoint = sigvec[sigvec.size()-1];
        delsig = abs(real(startingpoint) - max)/points2;
        sigma_vector_maker_linear_1(sigvec,startingpoint,-delsig,points2,1);

        /* print the size of the created sigma_vector */
        cout<<"sigvec created with size = "<<sigvec.size()<<endl;

    }
    else
    {
        double length = real(max - min);
        comp delsig = length/(points1 + 5.0*points2) ;

        int totpoints = points1 + 5*points2;
        for(int i=1;i<=totpoints;++i)
        {
            comp sigma = min + ii*eps + ((comp)i)*delsig;

            sigvec.push_back(sigma);
        }
    }

}

void sigma_vector_maker_6(  vector<comp> &sigvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double points1,
                            double points2,
                            double epsratio     )
{
    comp ii = {0.0,1.0};
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp delsig = {0.0,0.0};

    if(real(rad)>0.0)
    {
        /*  this is the initial straight line along the real axis */
        comp startingpoint = min + ii*eps;
        delsig = abs((sigmac1-epsratio*eps)-min)/points1;
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,points1,1);
        
        /*  this is the first contour up in the imaginary axis   */
        startingpoint = sigvec[sigvec.size()-1];
        delsig = abs(rad+((comp)epsratio*eps))/points2;
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,points2,2);

        /*  this is again parallel to the real axis, to the right, this 
            goes for 2*radius and for 2*points2                         */
        startingpoint = sigvec[sigvec.size()-1];
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,2.0*points2,1);

        /*  the contour now comes down, so we flip the sign of the delsig */
        startingpoint = sigvec[sigvec.size()-1];
        sigma_vector_maker_linear_1(sigvec,startingpoint,-delsig,points2,2);

        /*  this is the end part, the contour goes to the left to max */
        startingpoint = sigvec[sigvec.size()-1];
        delsig = abs(real(startingpoint) - max)/(2.0*points2);
        sigma_vector_maker_linear_1(sigvec,startingpoint,-delsig,2.0*points2,1);

        /* print the size of the created sigma_vector */
        cout<<"sigvec created with size = "<<sigvec.size()<<endl;

    }
    else
    {
        double length = real(max - min);
        comp delsig = length/(points1 + 6.0*points2) ;

        int totpoints = points1 + 6*points2;
        for(int i=1;i<=totpoints;++i)
        {
            comp sigma = min + ii*eps + ((comp)i)*delsig;

            sigvec.push_back(sigma);
        }
    }

}


void sigma_vector_maker_7(  vector<comp> &sigvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double points1,
                            double points2,
                            double epsratio     )
{
    comp ii = {0.0,1.0};
    comp rad = radius(s,sigmab(a,m),m);
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp delsig = {0.0,0.0};

    if(real(rad)>0.0)
    {
        /*  this is the initial straight line along the real axis */
        comp startingpoint = min + ii*eps;
        delsig = abs((sigmac1-epsratio*eps)-min)/points1;
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,points1,1);
        
        /*  this is the first contour up in the imaginary axis   */
        startingpoint = sigvec[sigvec.size()-1];
        delsig = abs(rad+((comp)epsratio*eps))/points2;
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,points2,2);

        /*  this is again parallel to the real axis, to the right, this 
            goes for 2*radius and for 2*points2                         */
        startingpoint = sigvec[sigvec.size()-1];
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,2.0*points2,1);

        /*  the contour now comes down, so we flip the sign of the delsig */
        startingpoint = sigvec[sigvec.size()-1];
        sigma_vector_maker_linear_1(sigvec,startingpoint,-delsig,points2,2);

        /*  this is the end part, the contour goes to the left to max */
        startingpoint = sigvec[sigvec.size()-1];
        delsig = abs(real(startingpoint) - max)/(6.0*points2);
        sigma_vector_maker_linear_1(sigvec,startingpoint,-delsig,6.0*points2,1);

        /* print the size of the created sigma_vector */
        cout<<"sigvec created with size = "<<sigvec.size()<<endl;

    }
    else
    {
        double length = real(max - min);
        comp delsig = length/(points1 + 10.0*points2) ;

        int totpoints = points1 + 10*points2;
        for(int i=1;i<=totpoints;++i)
        {
            comp sigma = min + ii*eps + ((comp)i)*delsig;

            sigvec.push_back(sigma);
        }
    }

}


//this is the sigma_vec_maker that uses complex energy  //
//with scomp = Re s - I eps, we also choose that sigma  //
//goes between Im(sigma_+) and Im((sqrt(s) - m)^2)      //
void sigma_vector_maker_8(  vector<comp> &sigvec,   
                            comp s,
                            comp min,
                            comp max,
                            double a,
                            double m,
                            double eps,
                            double points1,
                            double points2,
                            double epsratio     )
{
    comp ii = {0.0,1.0};
    comp rad = abs(radius(s,sigmab(a,m),m));
    comp sigmac1 = sigc1(s,sigmab(a,m),m);
    comp sigmac2 = sigc2(s,sigmab(a,m),m);
    comp sigb = sigmab(a,m);
    comp delsig = {0.0,0.0};

    if(real(rad)>0.0)
    {
        /*  this is the initial straight line along the real axis */
        comp startingpoint = min + ii*eps;
        delsig = abs((sigmac1-epsratio*eps)-min)/points1;
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,points1,1);
        
        /*  this is the first contour up in the imaginary axis   */
        startingpoint = sigvec[sigvec.size()-1];
        delsig = abs(rad+((comp)epsratio*eps))/points2;
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,points2,2);

        /*  this is again parallel to the real axis, to the right, this 
            goes for 2*radius and for 2*points2                         */
        startingpoint = sigvec[sigvec.size()-1];
        sigma_vector_maker_linear_1(sigvec,startingpoint,delsig,2.0*points2,1);

        /*  the contour now comes down, so we flip the sign of the delsig */
        startingpoint = sigvec[sigvec.size()-1];
        delsig = abs(rad + abs(imag(sigpplus_witheps(s,sigb,m,0.0)))+ ((comp)epsratio*eps + abs(imag(max) - imag(sigpplus_witheps(s,sigb,m,0.0)) )/2.0 ) )/points2;
        sigma_vector_maker_linear_1(sigvec,startingpoint,-delsig,points2,2);

        /*  this is the end part, the contour goes to the left to max */
        startingpoint = sigvec[sigvec.size()-1];
        delsig = abs(real(startingpoint) - real(max))/(6.0*points2);
        sigma_vector_maker_linear_1(sigvec,startingpoint,-delsig,6.0*points2,1);

        /* print the size of the created sigma_vector */
        cout<<"sigvec created with size = "<<sigvec.size()<<endl;

    }
    else
    {
        double length = real(max - min);
        comp delsig = length/(points1 + 10.0*points2) ;

        int totpoints = points1 + 10*points2;
        for(int i=1;i<=totpoints;++i)
        {
            comp sigma = min + ii*eps + ((comp)i)*delsig;

            sigvec.push_back(sigma);
        }
    }

}


#endif