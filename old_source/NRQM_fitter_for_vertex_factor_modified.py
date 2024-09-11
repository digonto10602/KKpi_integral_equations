#this is the working code cell. The function takes the input file name as a string
#fits the data and makes the plot

import numpy as np 
import math
import matplotlib as mpl 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable


def eps_linear_fitting_function(eps, A, B):
    result = A + eps*B 
    return result

def N_linear_fitting_function(N, A, B):
    result = A + B/N 
    return result

def s_linear_fitting_function(s, m, c):
    result = m*s + c 
    return result 

def NRQM_fitting_function(k, Asq):
    spole = 8.9999991916618836
    Eb = math.sqrt(spole)
    kappasq = 3.0 - math.sqrt(spole)
    kappa = math.sqrt(kappasq)
    modc = 96.351
    s0 = 1.00624
    pi = math.pi
    m = 1.0
    firstterm = modc*Asq

    secondterm = 256.0*pow(pi,5.0/2.0)/pow(3,1.0/4.0)
    thirdterm = m*m*kappa*kappa/(k*k*(kappa*kappa + 3.0*k*k/4.0))
    forthterm = pow(np.sin(s0*np.arcsinh(np.sqrt(3.0)*k/(2.0*kappa))),2.0)/pow(np.sinh(pi*s0/2.0),2.0)

    result = firstterm*secondterm*thirdterm*forthterm 
    return result 

def NRQM_fitting_function_modified(k, Asq, B):
    spole = 8.99999922238197#8.999792091437811#8.895133263228455
    Eb = math.sqrt(spole)
    kappasq = 3.0 - math.sqrt(spole)
    kappa = math.sqrt(kappasq)
    modc = 96.351
    s0 = 1.00624
    pi = math.pi
    m = 1.0
    firstterm = modc*Asq

    secondterm = 256.0*pow(pi,5.0/2.0)/pow(3,1.0/4.0)
    thirdterm = m*m*kappa*kappa/(k*k*(kappa*kappa + 3.0*k*k/4.0))
    forthterm = pow(np.sin(B*s0*np.arcsinh(np.sqrt(3.0)*k/(2.0*kappa))),2.0)/pow(np.sinh(pi*s0/2.0),2.0)

    result = firstterm*secondterm*thirdterm*forthterm 
    return result 

def fitting_code_eps(epsfilename, reguessA, reguessB, imguessA, imguessB):
    #fout<<setprecision(16)<<s<<'\t'<<eps<<'\t'<<size<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
    reguess = [reguessA,reguessB]
    imguess = [imguessA,imguessB]
    (s, eps, N, reRes, imRes) = np.genfromtxt(epsfilename,unpack=True)
    reParameters, reCovariance = curve_fit(eps_linear_fitting_function,eps,reRes, p0=reguess) 
    imParameters, imCovariance = curve_fit(eps_linear_fitting_function,eps,imRes, p0=imguess) 
    reFitA = reParameters[0]
    reFitB = reParameters[1]
    reSE = np.sqrt(np.diag(reCovariance))
    SE_reFitA = reSE[0]
    SE_reFitB = reSE[1]
    imFitA = imParameters[0]
    imFitB = imParameters[1]
    imSE = np.sqrt(np.diag(imCovariance))
    SE_imFitA = imSE[0]
    SE_imFitB = imSE[1]
    return reFitA, reFitB, SE_reFitA, SE_reFitB, imFitA, imFitB, SE_imFitA, SE_imFitB

def fitting_code_N(epsfilename, reguessA, reguessB, imguessA, imguessB):
    #fout<<setprecision(16)<<s<<'\t'<<eps<<'\t'<<size<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
    reguess = [reguessA,reguessB]
    imguess = [imguessA,imguessB]
    (s, eps, N, reRes, imRes) = np.genfromtxt(epsfilename,unpack=True)
    reParameters, reCovariance = curve_fit(N_linear_fitting_function,N,reRes, p0=reguess) 
    imParameters, imCovariance = curve_fit(N_linear_fitting_function,N,imRes, p0=imguess) 
    reFitA = reParameters[0]
    reFitB = reParameters[1]
    reSE = np.sqrt(np.diag(reCovariance))
    SE_reFitA = reSE[0]
    SE_reFitB = reSE[1]
    imFitA = imParameters[0]
    imFitB = imParameters[1]
    imSE = np.sqrt(np.diag(imCovariance))
    SE_imFitA = imSE[0]
    SE_imFitB = imSE[1]
    return reFitA, reFitB, SE_reFitA, SE_reFitB, imFitA, imFitB, SE_imFitA, SE_imFitB

def fitting_code_vertex(vertexfilename, reguessA, reguessB, imguessA, imguessB):
    #fout<<setprecision(16)<<s<<'\t'<<eps<<'\t'<<size<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
    reguess = [reguessA,reguessB]
    imguess = [imguessA,imguessB]
    (s, sigk, reqpole, imqpole, reRes, imRes, reInvRes, imInvRes, ressbRes, imssbRes) = np.genfromtxt(vertexfilename,unpack=True)
    reParameters, reCovariance = curve_fit(s_linear_fitting_function,s,reInvRes, p0=reguess) 
    imParameters, imCovariance = curve_fit(s_linear_fitting_function,s,imInvRes, p0=imguess) 
    reFitA = reParameters[0]
    reFitB = reParameters[1]
    reSE = np.sqrt(np.diag(reCovariance))
    SE_reFitA = reSE[0]
    SE_reFitB = reSE[1]
    imFitA = imParameters[0]
    imFitB = imParameters[1]
    imSE = np.sqrt(np.diag(imCovariance))
    SE_imFitA = imSE[0]
    SE_imFitB = imSE[1]
    return reFitA, reFitB, SE_reFitA, SE_reFitB, imFitA, imFitB, SE_imFitA, SE_imFitB

def fitting_code_NRQMvertex(vertexfilename, reguessA):
    #fout<<setprecision(16)<<s<<'\t'<<eps<<'\t'<<size<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
    reguess = [reguessA]
    
    (sigk, k, vertexres, SEvertex, Bval, SEBval, ksqvertexres, sbfix) = np.genfromtxt(vertexfilename,unpack=True)
    reParameters, reCovariance = curve_fit(NRQM_fitting_function,k,vertexres, p0=reguess) 
    #imParameters, imCovariance = curve_fit(s_linear_fitting_function,s,imInvRes, p0=imguess) 
    reFitA = reParameters[0]
    #reFitB = reParameters[1]
    reSE = np.sqrt(np.diag(reCovariance))
    SE_reFitA = reSE[0]
    #SE_reFitB = reSE[1]
    #imFitA = imParameters[0]
    #imFitB = imParameters[1]
    #imSE = np.sqrt(np.diag(imCovariance))
    #SE_imFitA = imSE[0]
    #SE_imFitB = imSE[1]
    return reFitA, SE_reFitA


def fitting_code_NRQMvertex_modified(vertexfilename, reguessA,reguessB):
    #fout<<setprecision(16)<<s<<'\t'<<eps<<'\t'<<size<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
    reguess = [reguessA,reguessB]
    
    (k, m2kressq, vertexressq, res) = np.genfromtxt(vertexfilename,unpack=True)



    reParameters, reCovariance = curve_fit(NRQM_fitting_function_modified,k,res, p0=reguess) 
    #imParameters, imCovariance = curve_fit(s_linear_fitting_function,s,imInvRes, p0=imguess) 
    reFitA = reParameters[0]
    reFitB = reParameters[1]
    reSE = np.sqrt(np.diag(reCovariance))
    SE_reFitA = reSE[0]
    SE_reFitB = reSE[1]
    #imFitA = imParameters[0]
    #imFitB = imParameters[1]
    #imSE = np.sqrt(np.diag(imCovariance))
    #SE_imFitA = imSE[0]
    #SE_imFitB = imSE[1]
    return reFitA, SE_reFitA, reFitB, SE_reFitB



def eps_plot_function(epsfilename, a, reguessA, reguessB, imguessA, imguessB):
    (s, eps, N, reRes, imRes) = np.genfromtxt(epsfilename,unpack=True)
    print("eps file loaded: ",epsfilename)
    fig, ax = plt.subplots(2,figsize = (12,10))
    
    sval = '%0.3f'%s[0]
    Nval = N[0]

    title_str = "$am = $" + str(a) + ", $s/m^2 = $" + str(sval) + ", $N = $" + str(Nval)
    fig.suptitle(title_str,fontsize=20)

    reFitA, reFitB, SE_reFitA, SE_reFitB, imFitA, imFitB, SE_imFitA, SE_imFitB = fitting_code_eps(epsfilename, reguessA, reguessB, imguessA, imguessB)
    print("reFitA = ",reFitA," reFitB = ", reFitB, " SE_reFitA = ",SE_reFitA, " SE_reFitB = ", SE_reFitB)
    print("imFitA = ",imFitA," imFitB = ", imFitB, " SE_imFitA = ",SE_imFitA, " SE_imFitB = ", SE_imFitB)

    ax[0].set_xlabel("$ \\epsilon$",fontsize=20)
    ax[0].set_ylabel("Re $\\mathcal{M}_{\\phi b}$",fontsize=20)
    ax[1].set_xlabel("$ \\epsilon$",fontsize=20)
    ax[1].set_ylabel("Im $\\mathcal{M}_{\\phi b}$",fontsize=20)

    ax[0].tick_params(axis='both', which='major', labelsize=20)
    ax[0].tick_params(axis='both', which='minor', labelsize=20)
    
    ax[1].tick_params(axis='both', which='major', labelsize=20)
    ax[1].tick_params(axis='both', which='minor', labelsize=20)
    
    
    ax[0].set_xscale('log')
    ax[1].set_xscale('log')
    ax[0].plot(eps,reRes,linestyle="solid",linewidth=2,color="darkred", label="data")
    ax[0].plot(eps,eps_linear_fitting_function(eps,reFitA,reFitB),linestyle="dashed",linewidth=2,color="red", label="fit")
    ax[0].plot(eps,eps_linear_fitting_function(eps,reFitA,0),linestyle="dashed",linewidth=2,color="grey", label="$\\epsilon \\rightarrow 0$")
    
    
    ax[1].plot(eps,imRes,linestyle="solid",linewidth=2,color="darkred", label="data")
    ax[1].plot(eps,eps_linear_fitting_function(eps,imFitA,imFitB),linestyle="dashed",linewidth=2,color="red", label="fit")
    ax[1].plot(eps,eps_linear_fitting_function(eps,imFitA,0),linestyle="dashed",linewidth=2,color="grey", label="$\\epsilon \\rightarrow 0$")
    
    ax[0].legend(loc='best')
    ax[1].legend(loc='best')
    outputfile_str = epsfilename + ".pdf"

    print("output file: ",outputfile_str)
    
    fig.tight_layout()
    
    plt.savefig(outputfile_str)
    plt.close()

def N_plot_function(Nfilename, a, reguessA, reguessB, imguessA, imguessB):
    (s, eps, N, reRes, imRes) = np.genfromtxt(Nfilename,unpack=True)
    print("N file loaded: ",Nfilename)
    fig, ax = plt.subplots(2,figsize = (12,10))
    
    sval = '%0.3f'%s[0]
    Nval = N[0]
    epsval = eps[0]

    title_str = "$am = $" + str(a) + ", $s/m^2 = $" + str(sval) + ", $\\epsilon = $" + str(epsval)
    fig.suptitle(title_str,fontsize=20)

    reFitA, reFitB, SE_reFitA, SE_reFitB, imFitA, imFitB, SE_imFitA, SE_imFitB = fitting_code_N(Nfilename, reguessA, reguessB, imguessA, imguessB)
    print("reFitA = ",reFitA," reFitB = ", reFitB, " SE_reFitA = ",SE_reFitA, " SE_reFitB = ", SE_reFitB)
    print("imFitA = ",imFitA," imFitB = ", imFitB, " SE_imFitA = ",SE_imFitA, " SE_imFitB = ", SE_imFitB)

    ax[0].set_xlabel("$N$",fontsize=20)
    ax[0].set_ylabel("Re $\\mathcal{M}_{\\phi b}$",fontsize=20)
    ax[1].set_xlabel("$N$",fontsize=20)
    ax[1].set_ylabel("Im $\\mathcal{M}_{\\phi b}$",fontsize=20)

    ax[0].tick_params(axis='both', which='major', labelsize=20)
    ax[0].tick_params(axis='both', which='minor', labelsize=20)
    ax[1].tick_params(axis='both', which='major', labelsize=20)
    ax[1].tick_params(axis='both', which='minor', labelsize=20)
    

    ax[0].plot(N,reRes,linestyle="solid",linewidth=2,color="darkred", label="data")
    ax[0].plot(N,N_linear_fitting_function(N,reFitA,reFitB),linestyle="dashed",linewidth=2,color="red", label="fit")
    ax[0].plot(N,N_linear_fitting_function(N,reFitA,0),linestyle="dashed",linewidth=2,color="grey", label="$\\epsilon \\rightarrow 0$")
    
    ax[1].plot(N,imRes,linestyle="solid",linewidth=2,color="darkred", label="data")
    ax[1].plot(N,N_linear_fitting_function(N,imFitA,imFitB),linestyle="dashed",linewidth=2,color="red", label="fit")
    ax[1].plot(N,N_linear_fitting_function(N,imFitA,0),linestyle="dashed",linewidth=2,color="grey", label="$\\epsilon \\rightarrow 0$")
    
    ax[0].legend(loc='best')
    ax[1].legend(loc='best')
    outputfile_str = Nfilename + ".pdf"

    print("output file: ",outputfile_str)
    
    fig.tight_layout()
    plt.savefig(outputfile_str)
    plt.close()



def multiN_plot_function(Nfilename1,Nfilename2,Nfilename3,Nfilename4, a, reguessA, reguessB, imguessA, imguessB):
    (s1, eps1, N1, reRes1, imRes1) = np.genfromtxt(Nfilename1,unpack=True)
    (s2, eps2, N2, reRes2, imRes2) = np.genfromtxt(Nfilename2,unpack=True)
    (s3, eps3, N3, reRes3, imRes3) = np.genfromtxt(Nfilename3,unpack=True)
    (s4, eps4, N4, reRes4, imRes4) = np.genfromtxt(Nfilename4,unpack=True)
    print("N file loaded: ",Nfilename1)
    print("N file loaded: ",Nfilename2)
    print("N file loaded: ",Nfilename3)
    print("N file loaded: ",Nfilename4)
    fig, ax = plt.subplots(2,figsize = (12,10))
    
    sval = '%0.3f'%s1[0]
    Nval = N1[0]
    epsval1 = eps1[0]
    epsval2 = eps2[0]
    epsval3 = eps3[0]
    epsval4 = eps4[0]

    title_str = "$am = $" + str(a) + ", $s/m^2 = $" + str(sval) 
    fig.suptitle(title_str,fontsize=20)

    reFitA1, reFitB1, SE_reFitA1, SE_reFitB1, imFitA1, imFitB1, SE_imFitA1, SE_imFitB1 = fitting_code_N(Nfilename1, reguessA, reguessB, imguessA, imguessB)
    print("reFitA = ",reFitA1," reFitB = ", reFitB1, " SE_reFitA = ",SE_reFitA1, " SE_reFitB = ", SE_reFitB1)
    print("imFitA = ",imFitA1," imFitB = ", imFitB1, " SE_imFitA = ",SE_imFitA1, " SE_imFitB = ", SE_imFitB1)

    reFitA2, reFitB2, SE_reFitA2, SE_reFitB2, imFitA2, imFitB2, SE_imFitA2, SE_imFitB2 = fitting_code_N(Nfilename2, reguessA, reguessB, imguessA, imguessB)
    print("reFitA = ",reFitA2," reFitB = ", reFitB2, " SE_reFitA = ",SE_reFitA2, " SE_reFitB = ", SE_reFitB2)
    print("imFitA = ",imFitA2," imFitB = ", imFitB2, " SE_imFitA = ",SE_imFitA2, " SE_imFitB = ", SE_imFitB2)

    reFitA3, reFitB3, SE_reFitA3, SE_reFitB3, imFitA3, imFitB3, SE_imFitA3, SE_imFitB3 = fitting_code_N(Nfilename3, reguessA, reguessB, imguessA, imguessB)
    print("reFitA = ",reFitA3," reFitB = ", reFitB3, " SE_reFitA = ",SE_reFitA3, " SE_reFitB = ", SE_reFitB3)
    print("imFitA = ",imFitA3," imFitB = ", imFitB3, " SE_imFitA = ",SE_imFitA3, " SE_imFitB = ", SE_imFitB3)

    reFitA4, reFitB4, SE_reFitA4, SE_reFitB4, imFitA4, imFitB4, SE_imFitA4, SE_imFitB4 = fitting_code_N(Nfilename4, reguessA, reguessB, imguessA, imguessB)
    print("reFitA = ",reFitA4," reFitB = ", reFitB4, " SE_reFitA = ",SE_reFitA4, " SE_reFitB = ", SE_reFitB4)
    print("imFitA = ",imFitA4," imFitB = ", imFitB4, " SE_imFitA = ",SE_imFitA4, " SE_imFitB = ", SE_imFitB4)

    ax[0].set_xlabel("$N$",fontsize=20)
    ax[0].set_ylabel("Re $\\mathcal{M}_{\\phi b}$",fontsize=20)
    ax[1].set_xlabel("$N$",fontsize=20)
    ax[1].set_ylabel("Im $\\mathcal{M}_{\\phi b}$",fontsize=20)

    ax[0].tick_params(axis='both', which='major', labelsize=20)
    ax[0].tick_params(axis='both', which='minor', labelsize=20)
    ax[1].tick_params(axis='both', which='major', labelsize=20)
    ax[1].tick_params(axis='both', which='minor', labelsize=20)
    

    #file 1 
    ax[0].scatter(N1,reRes1, label="data, $\\epsilon = $" + str(epsval1))
    ax[0].plot(N1,N_linear_fitting_function(N1,reFitA1,reFitB1),linestyle="solid",linewidth=2, label="fit")
    ax[0].plot(N1,N_linear_fitting_function(N1,reFitA1,0),linestyle="dashed",linewidth=2, label="$N \\rightarrow \\infty$")
    
    ax[1].scatter(N1,imRes1, label="data, $\\epsilon = $"+ str(epsval1))
    ax[1].plot(N1,N_linear_fitting_function(N1,imFitA1,imFitB1),linestyle="solid",linewidth=2, label="fit")
    ax[1].plot(N1,N_linear_fitting_function(N1,imFitA1,0),linestyle="dashed",linewidth=2, label="$N \\rightarrow \\infty$")
    
    #file 2 
    ax[0].scatter(N2,reRes2, label="data, $\\epsilon = $"+ str(epsval2))
    ax[0].plot(N2,N_linear_fitting_function(N2,reFitA2,reFitB2),linestyle="solid",linewidth=2, label="fit")
    ax[0].plot(N1,N_linear_fitting_function(N2,reFitA2,0),linestyle="dashed",linewidth=2, label="$N \\rightarrow \\infty$")
    
    ax[1].scatter(N2,imRes2, label="data, $\\epsilon = $"+ str(epsval2))
    ax[1].plot(N2,N_linear_fitting_function(N2,imFitA2,imFitB2),linestyle="solid",linewidth=2, label="fit")
    ax[1].plot(N2,N_linear_fitting_function(N2,imFitA2,0),linestyle="dashed",linewidth=2, label="$N \\rightarrow \\infty$")
    
    #file 3 
    ax[0].scatter(N3,reRes3, label="data, $\\epsilon = $"+ str(epsval3))
    ax[0].plot(N3,N_linear_fitting_function(N3,reFitA3,reFitB3),linestyle="solid",linewidth=2, label="fit")
    ax[0].plot(N3,N_linear_fitting_function(N3,reFitA3,0),linestyle="dashed",linewidth=2, label="$N \\rightarrow \\infty$")
    
    ax[1].scatter(N3,imRes3, label="data, $\\epsilon = $"+ str(epsval3))
    ax[1].plot(N3,N_linear_fitting_function(N3,imFitA3,imFitB3),linestyle="solid",linewidth=2, label="fit")
    ax[1].plot(N3,N_linear_fitting_function(N3,imFitA3,0),linestyle="dashed",linewidth=2, label="$N \\rightarrow \\infty$")
    
    #file 4
    #ax[0].scatter(N4,reRes4, label="data, $\\epsilon = $"+ str(epsval4))
    #ax[0].plot(N4,N_linear_fitting_function(N4,reFitA4,reFitB4),linestyle="solid",linewidth=2, label="fit")
    #ax[0].plot(N1,N_linear_fitting_function(N4,reFitA4,0),linestyle="dashed",linewidth=2, label="$N \\rightarrow \\infty$")
    
    #ax[1].scatter(N4,imRes4, label="data, $\\epsilon = $"+ str(epsval4))
    #ax[1].plot(N4,N_linear_fitting_function(N4,imFitA4,imFitB4),linestyle="solid",linewidth=2, label="fit")
    #ax[1].plot(N4,N_linear_fitting_function(N4,imFitA4,0),linestyle="dashed",linewidth=2, label="$N \\rightarrow \\infty$")
    
    ax[0].legend(loc='best')
    ax[1].legend(loc='best')
    #ax[0].set_yscale('log')
    #ax[1].set_yscale('log')
    outputfile_str = Nfilename1 + "+" + Nfilename2 + "+" + Nfilename3  + ".pdf"

    print("output file: ",outputfile_str)
    
    fig.tight_layout()
    plt.savefig(outputfile_str)
    plt.close()



def vertexfactor_making_function(vertexfilename, reguessA, reguessB, imguessA, imguessB ):
    (s, reqpole, imqpole, sigk, reRes, imRes, reInvRes, imInvRes, ressbRes, imssbRes) = np.genfromtxt(vertexfilename,unpack=True)
    print("vertex file loaded: ",vertexfilename)
    

    reFitA, reFitB, SE_reFitA, SE_reFitB, imFitA, imFitB, SE_imFitA, SE_imFitB = fitting_code_vertex(vertexfilename, reguessA, reguessB, imguessA, imguessB)
    print("reFitA = ",reFitA," reFitB = ", reFitB, " SE_reFitA = ",SE_reFitA, " SE_reFitB = ", SE_reFitB)
    print("imFitA = ",imFitA," imFitB = ", imFitB, " SE_imFitA = ",SE_imFitA, " SE_imFitB = ", SE_imFitB)

    vertexres = -1.0/reFitA 
    sb = -reFitB/reFitA
    return vertexres, sb, reqpole[0], reFitA, reFitB  
    

#epsfilename = "Mphib_momrep_a_16_N_5000_contour4.dat"
def sigk_vertexfilemaking_func():
    a = 16
    reguessA = 0
    reguessB = 0
    imguessA = 0
    imguessB = 0

    f = open("vertexfactor_vs_sigk_file.dat", "w")
    #f.write("Now the file has more content!")


    for i in range(0,201):
        spole = 8.782848835105487
        vertexfilename = "vertexfactor_momrep_a_16_N_500_eps_0_sigkcount_" + str(i) + ".dat"
        (s, sigk, reqpole1, imqpole1, reRes, imRes, reInvRes, imInvRes, ressbRes, imssbRes) = np.genfromtxt(vertexfilename,unpack=True)
        vertexres, sb, reqpole, m, c = vertexfactor_making_function(vertexfilename,reguessA,reguessB,imguessA,imguessB)
        print(sigk[0],vertexres,sb)
        sigkval = '%5.16f'%sigk[0]
        vertexval = '%5.16f'%vertexres
        sbval = '%5.16f'%sb
        spoleval = '%5.16f'%spole
        qpoleval = '%5.16f'%reqpole 
        mval = '%5.16f'%m
        cval = '%5.16f'%c 
        f.write(sigkval  + "\t" + vertexval + "\t" + sbval + "\t" + spoleval + "\t" + qpoleval + "\t" + mval + "\t" + cval +  "\n" )
    f.close()

reguessA = 1.0
reguessB = 1.0
a = 10000
N = 7000
BSnumber = 3
filename = "vertexfactor_vs_k_BS"+str(BSnumber)+"_for_a_"+str(a)+"_N_" + str(N) + ".dat"
#(sigk, k, vertexres, SEvertex, Bval, SEBval, ksqvertexres, sbfix) = np.genfromtxt(filename,unpack=True)
(k, m2kressq, vertexressq, res) = np.genfromtxt(filename,unpack=True)

Asqfit, Asqerr, B, BERR = fitting_code_NRQMvertex_modified(filename,reguessA,reguessB)

print("Asq = ",Asqfit," AsqERR = ",Asqerr)
print("B = ",B," Berr = ",BERR)

outputfilename = "modified_NRQM_with_" + filename 
fittedfile = "fitted_NRQM_with_" + filename 
f = open(outputfilename, "w")
#f1 = open(fittedfile,"w")

for i in range(0,len(k)):
    nrqmvertex = NRQM_fitting_function_modified(k[i],Asqfit,B)
    
    f.write( str(k[i]) + "\t" + str(nrqmvertex)  +  "\n" )
    
f.close()


