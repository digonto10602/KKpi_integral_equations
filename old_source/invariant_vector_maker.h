#ifndef INVARIANT_VECTOR_MAKER_H
#define INVARIANT_VECTOR_MAKER_H

#include<bits/stdc++.h>


using namespace std;

typedef complex<double> comp;

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