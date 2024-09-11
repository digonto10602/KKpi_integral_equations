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
    (s, sigk, reRes, imRes, reInvRes, imInvRes, ressbRes, imssbRes) = np.genfromtxt(vertexfilename,unpack=True)
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
    (s, sigk, reRes, imRes, reInvRes, imInvRes, ressbRes, imssbRes) = np.genfromtxt(vertexfilename,unpack=True)
    print("vertex file loaded: ",vertexfilename)
    

    reFitA, reFitB, SE_reFitA, SE_reFitB, imFitA, imFitB, SE_imFitA, SE_imFitB = fitting_code_vertex(vertexfilename, reguessA, reguessB, imguessA, imguessB)
    print("reFitA = ",reFitA," reFitB = ", reFitB, " SE_reFitA = ",SE_reFitA, " SE_reFitB = ", SE_reFitB)
    print("imFitA = ",imFitA," imFitB = ", imFitB, " SE_imFitA = ",SE_imFitA, " SE_imFitB = ", SE_imFitB)

    vertexres = -1.0/reFitA 
    sb = -reFitB/reFitA
    return vertexres, sb, reFitA, reFitB  
    

#epsfilename = "Mphib_momrep_a_16_N_5000_contour4.dat"


def NRQM_with_vertexfunction_plot_function(a,filename1, filename2, filename3, 
                                             nrqmfile1, nrqmfile2, nrqmfile3,
                                             fitfile1,  fitfile2,  fitfile3  ):
    #kval + "\t" + vertexval + "\t" + nrqmvertexval + "\t" + sigkval + "\t" + sbfixval +  "\n"
    plt.rcParams.update({'font.size': 16})
    plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)
    (k1, vertexres1, nrqmvertexres1, sigk1, sbfix1) = np.genfromtxt(filename1,unpack=True)
    (k2, vertexres2, nrqmvertexres2, sigk2, sbfix2) = np.genfromtxt(filename2,unpack=True)
    (k3, vertexres3, nrqmvertexres3, sigk3, sbfix3) = np.genfromtxt(filename3,unpack=True)
    (k4, vertexres4, nrqmvertexres4, sigk4, sbfix4) = np.genfromtxt(nrqmfile1,unpack=True)
    (k5, vertexres5, nrqmvertexres5, sigk5, sbfix5) = np.genfromtxt(nrqmfile2,unpack=True)
    (k6, vertexres6, nrqmvertexres6, sigk6, sbfix6) = np.genfromtxt(nrqmfile3,unpack=True)
    (fitk7, fitvertexres7, Asq7, AsqERR7, B7, BERR7 ) = np.genfromtxt(fitfile1,unpack=True)
    (fitk8, fitvertexres8, Asq8, AsqERR8, B8, BERR8 ) = np.genfromtxt(fitfile2,unpack=True)
    (fitk9, fitvertexres9, Asq9, AsqERR9, B9, BERR9 ) = np.genfromtxt(fitfile3,unpack=True)
    
    print("file1 loaded: ",filename1)
    print("file2 loaded: ",filename2)
    print("file3 loaded: ",filename3)
    fig, ax = plt.subplots(figsize = (12,10))
    
    #sval = '%0.3f'%s[0]
    #Nval = N[0]
    #epsval = eps[0]
    A_BS1sq = Asq7[0]
    A_BS1sq_ERR = AsqERR7[0]
    A_BS2sq = Asq8[0]
    A_BS2sq_ERR = AsqERR8[0]
    A_BS3sq = Asq9[0]
    A_BS3sq_ERR = AsqERR9[0]
    
    B_fit_BS1 = B7[0]
    B_fit_BS1_ERR = BERR7[0]
    B_fit_BS2 = B8[0]
    B_fit_BS2_ERR = BERR8[0]
    B_fit_BS3 = B9[0]
    B_fit_BS3_ERR = BERR9[0]

    A_BS1sqval = '%1.3f'%A_BS1sq
    A_BS1sq_ERRval = '%1.3f'%A_BS1sq_ERR
    A_BS2sqval = '%1.3f'%A_BS2sq
    A_BS2sq_ERRval = '%1.3f'%A_BS2sq_ERR
    A_BS3sqval = '%1.3f'%A_BS3sq
    A_BS3sq_ERRval = '%1.3f'%A_BS3sq_ERR

    Bfit_BS1_val = '%1.3f'%B_fit_BS1
    Bfit_BS1_ERR_val = '%1.3f'%B_fit_BS1_ERR
    Bfit_BS2_val = '%1.3f'%B_fit_BS2
    Bfit_BS2_ERR_val = '%1.3f'%B_fit_BS2_ERR
    Bfit_BS3_val = '%1.3f'%B_fit_BS3
    Bfit_BS3_ERR_val = '%1.3f'%B_fit_BS3_ERR

    title_str = "$am = $" + str(a) 
    fig.suptitle(title_str,fontsize=20)


    ax.set_xlabel("$k/m$",fontsize=20)
    ax.set_ylabel("$|\\Gamma^{(u)}(k)|^2$",fontsize=20)
    ax.set_yscale('log')
    
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.tick_params(axis='both', which='minor', labelsize=20)
    
    nrqmlabel1 = "NRQM BS1,$|A|^2=$" +  A_BS1sqval + "$\\pm$" + A_BS1sq_ERRval + ", $B=$" + Bfit_BS1_val +"$\pm $" + Bfit_BS1_ERR_val
    nrqmlabel2 = "NRQM BS2,$|A|^2=$" +  A_BS2sqval + "$\\pm$" + A_BS2sq_ERRval + ", $B=$" + Bfit_BS2_val +"$\pm $" + Bfit_BS2_ERR_val
    nrqmlabel3 = "NRQM BS3,$|A|^2=$" +  A_BS3sqval + "$\\pm$" + A_BS3sq_ERRval + ", $B=$" + Bfit_BS3_val +"$\pm $" + Bfit_BS3_ERR_val

    ax.scatter(k1,vertexres1, s=120, facecolor="none", edgecolor="blue", linewidth=2, label="Rel. BS1")
    ax.plot(k4,nrqmvertexres4,linestyle="solid",linewidth=2,color="darkblue", label=nrqmlabel1)
    ax.scatter(fitk7,fitvertexres7, s=120, marker="s",facecolor="white", edgecolor="lightskyblue", linewidth=2)
    ax.scatter(k2,vertexres2, s=120, facecolor="none", edgecolor="red", linewidth=2, label="Rel. BS2")
    ax.plot(k5,nrqmvertexres5,linestyle="solid",linewidth=2,color="darkred", label=nrqmlabel2)
    ax.scatter(fitk8,fitvertexres8, s=120, marker="s",facecolor="white", edgecolor="darkorange", linewidth=2)
    ax.scatter(k3,vertexres3, s=120, facecolor="none", edgecolor="green", linewidth=2, label="Rel. BS3")
    ax.plot(k6,nrqmvertexres6,linestyle="solid",linewidth=2,color="darkgreen", label=nrqmlabel3)
    ax.scatter(fitk9,fitvertexres9, s=120, marker="s",facecolor="white", edgecolor="turquoise", linewidth=2)
    
    
    ax.legend(loc='best')


    outputfile_str = "modified_vertexfactor_comparison_vs_k_with_extended_NRQM_result.pdf"

    print("output file: ",outputfile_str)
    
    plt.legend()
    fig.tight_layout()
    plt.savefig(outputfile_str)
    plt.close()

def NRQM_with_vertexfunction_by_kappasq_plot_function(a,filename1,filename2,filename3, nrqmfile1, nrqmfile2, nrqmfile3):
    #kval + "\t" + vertexval + "\t" + nrqmvertexval + "\t" + sigkval + "\t" + sbfixval +  "\n"
    plt.rcParams.update({'font.size': 16})
    plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)
    (k1, vertexres1, nrqmvertexres1, sigk1, sbfix1) = np.genfromtxt(filename1,unpack=True)
    (k2, vertexres2, nrqmvertexres2, sigk2, sbfix2) = np.genfromtxt(filename2,unpack=True)
    (k3, vertexres3, nrqmvertexres3, sigk3, sbfix3) = np.genfromtxt(filename3,unpack=True)
    (k4, vertexres4, nrqmvertexres4, sigk4, sbfix4) = np.genfromtxt(nrqmfile1,unpack=True)
    (k5, vertexres5, nrqmvertexres5, sigk5, sbfix5) = np.genfromtxt(nrqmfile2,unpack=True)
    (k6, vertexres6, nrqmvertexres6, sigk6, sbfix6) = np.genfromtxt(nrqmfile3,unpack=True)
    print("file1 loaded: ",filename1)
    print("file2 loaded: ",filename2)
    print("file3 loaded: ",filename3)
    fig, ax = plt.subplots(figsize = (12,10))
    
    kappasq1 = 3.0 - math.sqrt(sbfix1[0])
    kappasq2 = 3.0 - math.sqrt(sbfix2[0])
    kappasq3 = 3.0 - math.sqrt(sbfix3[0])

    
    #sval = '%0.3f'%s[0]
    #Nval = N[0]
    #epsval = eps[0]

    title_str = "$am = $" + str(a) 
    fig.suptitle(title_str,fontsize=20)


    ax.set_xlabel("$k/m$",fontsize=20)
    ax.set_ylabel("$|\\Gamma(k)|^2/\\kappa^2$",fontsize=20)
    ax.set_yscale('log')
    
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.tick_params(axis='both', which='minor', labelsize=20)

    A_BS1sq = 0.9883543584593427
    A_BS1sq_ERR = 0.0005988691319459065
    A_BS2sq = 1.011825249531947
    A_BS2sq_ERR = 0.00045955262275461205
    A_BS3sq = 0.7867626730028752
    A_BS3sq_ERR = 0.006908180622985164

    A_BS1sqval = '%1.3f'%A_BS1sq
    A_BS1sq_ERRval = '%1.3f'%A_BS1sq_ERR
    A_BS2sqval = '%1.3f'%A_BS2sq
    A_BS2sq_ERRval = '%1.3f'%A_BS2sq_ERR
    A_BS3sqval = '%1.3f'%A_BS3sq
    A_BS3sq_ERRval = '%1.3f'%A_BS3sq_ERR

    nrqmlabel1 = "NRQM BS1,$|A|^2=$" +  A_BS1sqval + "$\\pm$" + A_BS1sq_ERRval
    nrqmlabel2 = "NRQM BS2,$|A|^2=$" +  A_BS2sqval + "$\\pm$" + A_BS2sq_ERRval
    nrqmlabel3 = "NRQM BS3,$|A|^2=$" +  A_BS3sqval + "$\\pm$" + A_BS3sq_ERRval
    

    ax.scatter(k1,vertexres1/kappasq1, s=120, facecolor="none", edgecolor="blue", linewidth=2, label="Rel. BS1")
    ax.plot(k4,nrqmvertexres4/kappasq1,linestyle="solid",linewidth=2,color="darkblue", label=nrqmlabel1)
    ax.scatter(k2,vertexres2/kappasq2, s=120, facecolor="none", edgecolor="red", linewidth=2, label="Rel. BS2")
    ax.plot(k5,nrqmvertexres5/kappasq2,linestyle="solid",linewidth=2,color="darkred", label=nrqmlabel2)
    ax.scatter(k3,vertexres3/kappasq3, s=120, facecolor="none", edgecolor="green", linewidth=2, label="Rel. BS3")
    ax.plot(k6,nrqmvertexres6/kappasq3,linestyle="solid",linewidth=2,color="darkgreen", label=nrqmlabel3)
    
    
    ax.legend(loc='best')


    outputfile_str = "kappasq_vertexfactor_comparison_vs_k_with_Extended_NRQM_result.pdf"

    print("output file: ",outputfile_str)
    
    plt.legend()
    fig.tight_layout()
    plt.savefig(outputfile_str)
    plt.close()

def NRQM_with_vertexfunction_by_kappasq_plot_function_1(a,filename1, filename2, filename3, 
                                             nrqmfile1, nrqmfile2, nrqmfile3,
                                             fitfile1,  fitfile2,  fitfile3):
    #kval + "\t" + vertexval + "\t" + nrqmvertexval + "\t" + sigkval + "\t" + sbfixval +  "\n"
    plt.rcParams.update({'font.size': 16})
    plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)
    (k1, vertexres1, nrqmvertexres1, sigk1, sbfix1) = np.genfromtxt(filename1,unpack=True)
    (k2, vertexres2, nrqmvertexres2, sigk2, sbfix2) = np.genfromtxt(filename2,unpack=True)
    (k3, vertexres3, nrqmvertexres3, sigk3, sbfix3) = np.genfromtxt(filename3,unpack=True)
    (k4, vertexres4, nrqmvertexres4, sigk4, sbfix4) = np.genfromtxt(nrqmfile1,unpack=True)
    (k5, vertexres5, nrqmvertexres5, sigk5, sbfix5) = np.genfromtxt(nrqmfile2,unpack=True)
    (k6, vertexres6, nrqmvertexres6, sigk6, sbfix6) = np.genfromtxt(nrqmfile3,unpack=True)
    (fitk7, fitvertexres7, Asq7, AsqERR7, B7, BERR7 ) = np.genfromtxt(fitfile1,unpack=True)
    (fitk8, fitvertexres8, Asq8, AsqERR8, B8, BERR8 ) = np.genfromtxt(fitfile2,unpack=True)
    (fitk9, fitvertexres9, Asq9, AsqERR9, B9, BERR9 ) = np.genfromtxt(fitfile3,unpack=True)
    print("file1 loaded: ",filename1)
    print("file2 loaded: ",filename2)
    print("file3 loaded: ",filename3)
    fig, ax = plt.subplots(figsize = (12,10))
    
    kappasq1 = 3.0 - math.sqrt(sbfix1[0])
    kappasq2 = 3.0 - math.sqrt(sbfix2[0])
    kappasq3 = 3.0 - math.sqrt(sbfix3[0])

    
    #sval = '%0.3f'%s[0]
    #Nval = N[0]
    #epsval = eps[0]

    title_str = "$am = $" + str(a) 
    fig.suptitle(title_str,fontsize=20)


    ax.set_xlabel("$k/m$",fontsize=20)
    ax.set_ylabel("$|\\Gamma^{(u)}(k)|^2/\\frac{\\kappa^2}{k^2(\\kappa^2 + 3k^2/4)}$",fontsize=20)
    ax.set_yscale('log')
    ax.set_ylim(1000,100000)
    
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.tick_params(axis='both', which='minor', labelsize=20)

    A_BS1sq = Asq7[0]
    A_BS1sq_ERR = AsqERR7[0]
    A_BS2sq = Asq8[0]
    A_BS2sq_ERR = AsqERR8[0]
    A_BS3sq = Asq9[0]
    A_BS3sq_ERR = AsqERR9[0]
    
    B_fit_BS1 = B7[0]
    B_fit_BS1_ERR = BERR7[0]
    B_fit_BS2 = B8[0]
    B_fit_BS2_ERR = BERR8[0]
    B_fit_BS3 = B9[0]
    B_fit_BS3_ERR = BERR9[0]

    A_BS1sqval = '%1.3f'%A_BS1sq
    A_BS1sq_ERRval = '%1.3f'%A_BS1sq_ERR
    A_BS2sqval = '%1.3f'%A_BS2sq
    A_BS2sq_ERRval = '%1.3f'%A_BS2sq_ERR
    A_BS3sqval = '%1.3f'%A_BS3sq
    A_BS3sq_ERRval = '%1.3f'%A_BS3sq_ERR

    Bfit_BS1_val = '%1.3f'%B_fit_BS1
    Bfit_BS1_ERR_val = '%1.3f'%B_fit_BS1_ERR
    Bfit_BS2_val = '%1.3f'%B_fit_BS2
    Bfit_BS2_ERR_val = '%1.3f'%B_fit_BS2_ERR
    Bfit_BS3_val = '%1.3f'%B_fit_BS3
    Bfit_BS3_ERR_val = '%1.3f'%B_fit_BS3_ERR

    nrqmlabel1 = "NRQM BS1,$|A|^2=$" +  A_BS1sqval + "$\\pm$" + A_BS1sq_ERRval + ", $B=$" + Bfit_BS1_val +"$\pm $" + Bfit_BS1_ERR_val
    nrqmlabel2 = "NRQM BS2,$|A|^2=$" +  A_BS2sqval + "$\\pm$" + A_BS2sq_ERRval + ", $B=$" + Bfit_BS2_val +"$\pm $" + Bfit_BS2_ERR_val
    nrqmlabel3 = "NRQM BS3,$|A|^2=$" +  A_BS3sqval + "$\\pm$" + A_BS3sq_ERRval + ", $B=$" + Bfit_BS3_val +"$\pm $" + Bfit_BS3_ERR_val


    ax.scatter(k1,vertexres1/(kappasq1/(np.square(k1)*(kappasq1 + 3.0*np.square(k1)/4.0))), s=120, facecolor="none", edgecolor="blue", linewidth=2, label="Rel. BS1",alpha=0.2)
    ax.plot(k4,nrqmvertexres4/(kappasq1/(np.square(k4)*(kappasq1 + 3.0*np.square(k4)/4.0))),linestyle="solid",linewidth=2,color="darkblue", label=nrqmlabel1)
    ax.scatter(k2,vertexres2/(kappasq2/(np.square(k2)*(kappasq2 + 3.0*np.square(k2)/4.0))), s=120, facecolor="none", edgecolor="red", linewidth=2, label="Rel. BS2",alpha=0.2)
    ax.plot(k5,nrqmvertexres5/(kappasq2/(np.square(k5)*(kappasq2 + 3.0*np.square(k5)/4.0))),linestyle="solid",linewidth=2,color="darkred", label=nrqmlabel2)
    ax.scatter(k3,vertexres3/(kappasq3/(np.square(k3)*(kappasq3 + 3.0*np.square(k3)/4.0))), s=120, facecolor="none", edgecolor="green", linewidth=2, label="Rel. BS3",alpha=0.2)
    ax.plot(k6,nrqmvertexres6/(kappasq3/(np.square(k6)*(kappasq3 + 3.0*np.square(k6)/4.0))),linestyle="solid",linewidth=2,color="darkgreen", label=nrqmlabel3)
    
    
    ax.legend(loc='best')


    outputfile_str = "modified_kappasq_vertexfactor_comparison_vs_k_with_Extended_NRQM_result.pdf"

    print("output file: ",outputfile_str)
    
    plt.legend()
    fig.tight_layout()
    plt.savefig(outputfile_str)
    plt.close()


a = 100000
filename1 = "modified_NRQM_with_vertexfactor_vs_sigk_BS1_file_a="+str(a)+"_extrapolated.dat"
filename2 = "modified_NRQM_with_vertexfactor_vs_sigk_BS2_file_a="+str(a)+"_extrapolated.dat"
filename3 = "modified_NRQM_with_vertexfactor_vs_sigk_BS3_file_a="+str(a)+"_extrapolated.dat"
nrqmfile1 = "modified_OnlyNRQM_for_vertexfactor_vs_sigk_BS1_file_a="+str(a)+"_extrapolated.dat"
nrqmfile2 = "modified_OnlyNRQM_for_vertexfactor_vs_sigk_BS2_file_a="+str(a)+"_extrapolated.dat"
nrqmfile3 = "modified_OnlyNRQM_for_vertexfactor_vs_sigk_BS3_file_a="+str(a)+"_extrapolated.dat"
fitfile1 = "fitted_NRQM_with_vertexfactor_vs_sigk_BS1_file_a="+str(a)+"_extrapolated.dat"
fitfile2 = "fitted_NRQM_with_vertexfactor_vs_sigk_BS2_file_a="+str(a)+"_extrapolated.dat"
fitfile3 = "fitted_NRQM_with_vertexfactor_vs_sigk_BS3_file_a="+str(a)+"_extrapolated.dat"


NRQM_with_vertexfunction_plot_function(a,filename1,filename2,filename3,nrqmfile1,nrqmfile2,nrqmfile3,fitfile1,fitfile2,fitfile3)
NRQM_with_vertexfunction_by_kappasq_plot_function_1(a,filename1,filename2,filename3,nrqmfile1,nrqmfile2,nrqmfile3,fitfile1,fitfile2,fitfile3)
#NRQM_with_vertexfunction_by_kappasq_plot_function(a,filename1,filename2,filename3,nrqmfile1,nrqmfile2,nrqmfile3)