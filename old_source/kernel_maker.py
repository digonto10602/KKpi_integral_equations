#this is the working code cell. The function takes the input file name as a string
#makes the contour plot and save it in a pdf file

import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os.path
from scipy.optimize import curve_fit
import scipy.interpolate


def N_linear_fitting_function(N, A, B):
    result = A + B/N 
    return result

def N_quadratic_fitting_function(N, A, B, C):
    result = A + B/N + C/(N*N)
    return result 

def fitting_code_N_linear(filename, reguessA, reguessB, imguessA, imguessB):
    #fout<<setprecision(16)<<s<<'\t'<<eps<<'\t'<<size<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
    reguess = [reguessA,reguessB]
    imguess = [imguessA,imguessB]
    (res, ims, N, reRes, imRes) = np.genfromtxt(filename,unpack=True)
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

def fitting_code_N_linear_seba(filename, reguessA, reguessB, imguessA, imguessB):
    #fout<<setprecision(16)<<s<<'\t'<<eps<<'\t'<<size<<'\t'<<real(result)<<'\t'<<imag(result)<<endl;
    reguess = [reguessA,reguessB]
    imguess = [imguessA,imguessB]
    (N, reRes, imRes) = np.genfromtxt(filename,unpack=True)
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

def fitting_code_N_quadratic(filename, reguessA, reguessB, reguessC, imguessA, imguessB, imguessC):
    reguess = [reguessA, reguessB, reguessC]
    imguess = [imguessA, imguessB, imguessC]
    (res, ims, N, reRes, imRes) = np.genfromtxt(filename, unpack=True)
    reParameters, reCovariance = curve_fit(N_quadratic_fitting_function,N,reRes, p0=reguess)
    imParameters, imCovariance = curve_fit(N_quadratic_fitting_function,N,imRes, p0=imguess)
    reFitA = reParameters[0]
    reFitB = reParameters[1]
    reFitC = reParameters[2]
    reSE = np.sqrt(np.diag(reCovariance))
    SE_reFitA = reSE[0]
    SE_reFitB = reSE[1]
    SE_reFitC = reSE[2]
    imFitA = imParameters[0]
    imFitB = imParameters[1]
    imFitC = imParameters[2]
    imSE = np.sqrt(np.diag(imCovariance))
    SE_imFitA = imSE[0]
    SE_imFitB = imSE[1]
    SE_imFitC = imSE[2]
    return reFitA, reFitB, reFitC, SE_reFitA, SE_reFitB, SE_reFitC, imFitA, imFitB, imFitC, SE_imFitA, SE_imFitB, SE_imFitC 


def plotFit( filename1, filename2 ):
    plt.rcParams.update({'font.size': 12})
    plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)
    (res1, ims1, N1, reRes1, imRes1) = np.genfromtxt(filename1, unpack=True)
    ##(res2, ims2, N2, reRes2, imRes2) = np.genfromtxt(filename2, unpack=True)
    (N2, reRes2, imRes2) = np.genfromtxt(filename2, unpack=True)
    title_str = "$am = $" + str(16) + ", $ s/m^2 = "+str(res1[0]) + "+" +str(ims1[0]) + "\ i$" 

    fig, ax = plt.subplots(2,1,figsize=(12,10))

    fig.suptitle(title_str)

    ax[0].set_xlabel("$N$", fontsize=20)
    ax[0].set_ylabel("Re $(\mathcal{M}_{\phi b})$", fontsize=20)

    ax[1].set_xlabel("$N$", fontsize=20)
    ax[1].set_ylabel("Im $(\mathcal{M}_{\phi b})$", fontsize=20)

    ax[0].tick_params(axis='both', which='major', labelsize=20)
    ax[0].tick_params(axis='both', which='minor', labelsize=20)
    ax[1].tick_params(axis='both', which='major', labelsize=20)
    ax[1].tick_params(axis='both', which='minor', labelsize=20)
    
    (reFitA1, reFitB1, SE_reFitA1, SE_reFitB1, imFitA1, imFitB1, SE_imFitA1, SE_imFitB1) = fitting_code_N_linear(filename1,1.0,1.0,1.0,1.0)
    (reFitA2, reFitB2, SE_reFitA2, SE_reFitB2, imFitA2, imFitB2, SE_imFitA2, SE_imFitB2) = fitting_code_N_linear_seba(filename2,1.0,1.0,1.0,1.0)

    ax[0].scatter(N1,reRes1,s=100, marker="o", color='darkred', alpha=0.55, label="digonto's Data" )
    ax[0].plot(N1,N_linear_fitting_function(N1,reFitA1,reFitB1),linestyle='solid', linewidth=2, color='black' )
    ax[0].axhline(y=reFitA1,linestyle='dashed', linewidth=2, color='red', label="digonto, $N \\rightarrow \\infty$" )
    ax[0].scatter(N2,reRes2,s=100, marker="o", color='purple', alpha=0.55, label="sebastian's Data" )
    ax[0].plot(N2,N_linear_fitting_function(N2,reFitA2,reFitB2),linestyle='solid', linewidth=2, color='black' )
    ax[0].axhline(y=reFitA2,linestyle='dashed', linewidth=2, color='orange', label="sebastian, $N \\rightarrow \\infty$" )

    
    ax[1].scatter(N1,imRes1,s=100, marker="o", color='teal', alpha=0.55, label="digonto's Data" )
    ax[1].plot(N1,N_linear_fitting_function(N1,imFitA1,imFitB1), linestyle='solid', linewidth=2, color='black' )
    ax[1].axhline(y=imFitA1,linestyle='dashed', linewidth=2, color='blue', label="digonto, $N \\rightarrow \\infty$" )
    ax[1].scatter(N2,imRes2,s=100, marker="o", color='darkgreen', alpha=0.55, label="sebastian's Data" )
    ax[1].plot(N2,N_linear_fitting_function(N2,imFitA2,imFitB2), linestyle='solid', linewidth=2, color='black' )
    ax[1].axhline(y=imFitA2,linestyle='dashed', linewidth=2, color='green', label="sebastian, $N \\rightarrow \\infty$" )

    output_file = filename1 + ".pdf" #"Mphib_contour_comparison_linear_fit.pdf"

    fig.tight_layout()
    fig.legend()
    print("output file = ",output_file)
    plt.savefig(output_file)
    plt.close()

def plotFit3( filename1, filename2, filename3 ):
    plt.rcParams.update({'font.size': 12})
    plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)
    (res1, ims1, N1, reRes1, imRes1) = np.genfromtxt(filename1, unpack=True)
    (res2, ims2, N2, reRes2, imRes2) = np.genfromtxt(filename2, unpack=True)
    (res3, ims3, N3, reRes3, imRes3) = np.genfromtxt(filename3, unpack=True)

    title_str = "$am = $" + str(16) + ", $ s/m^2 = 8.6 - 0.05\ i$" 

    fig, ax = plt.subplots(2,1,figsize=(12,10))

    fig.suptitle(title_str)

    ax[0].set_xlabel("$N$", fontsize=20)
    ax[0].set_ylabel("Re $(\mathcal{M}_{\phi b})$", fontsize=20)

    ax[1].set_xlabel("$N$", fontsize=20)
    ax[1].set_ylabel("Im $(\mathcal{M}_{\phi b})$", fontsize=20)

    ax[0].tick_params(axis='both', which='major', labelsize=20)
    ax[0].tick_params(axis='both', which='minor', labelsize=20)
    ax[1].tick_params(axis='both', which='major', labelsize=20)
    ax[1].tick_params(axis='both', which='minor', labelsize=20)
    #ax[2].tick_params(axis='both', which='major', labelsize=20)
    #ax[2].tick_params(axis='both', which='minor', labelsize=20)
    
    (reFitA1, reFitB1, SE_reFitA1, SE_reFitB1, imFitA1, imFitB1, SE_imFitA1, SE_imFitB1) = fitting_code_N_linear(filename1,1.0,1.0,1.0,1.0)
    (reFitA2, reFitB2, SE_reFitA2, SE_reFitB2, imFitA2, imFitB2, SE_imFitA2, SE_imFitB2) = fitting_code_N_linear(filename2,1.0,1.0,1.0,1.0)
    (reFitA3, reFitB3, SE_reFitA3, SE_reFitB3, imFitA3, imFitB3, SE_imFitA3, SE_imFitB3) = fitting_code_N_linear(filename3,1.0,1.0,1.0,1.0)

    ax[0].scatter(N1,reRes1,s=100, marker="o", color='darkred', alpha=0.55, label="digonto's Data" )
    ax[0].plot(N1,N_linear_fitting_function(N1,reFitA1,reFitB1),linestyle='solid', linewidth=2, color='black' )
    ax[0].axhline(y=reFitA1,linestyle='dashed', linewidth=2, color='red', label="digonto, $N \\rightarrow \\infty$" )
    ax[0].scatter(N2,reRes2,s=100, marker="o", color='purple', alpha=0.55, label="sebastian's Data" )
    ax[0].plot(N2,N_linear_fitting_function(N2,reFitA2,reFitB2),linestyle='solid', linewidth=2, color='black' )
    ax[0].axhline(y=reFitA2,linestyle='dashed', linewidth=2, color='orange', label="sebastian, $N \\rightarrow \\infty$" )
    ax[0].scatter(N3,reRes3,s=100, marker="s", color='saddlebrown', alpha=0.55, label="straight cont. Data" )
    ax[0].plot(N3,N_linear_fitting_function(N3,reFitA3,reFitB3),linestyle='solid', linewidth=2, color='black' )
    ax[0].axhline(y=reFitA3,linestyle='dashed', linewidth=2, color='grey', label="straight, $N \\rightarrow \\infty$" )

    
    ax[1].scatter(N1,imRes1,s=100, marker="o", color='teal', alpha=0.55, label="digonto's Data" )
    ax[1].plot(N1,N_linear_fitting_function(N1,imFitA1,imFitB1), linestyle='solid', linewidth=2, color='black' )
    ax[1].axhline(y=imFitA1,linestyle='dashed', linewidth=2, color='blue', label="digonto, $N \\rightarrow \\infty$" )
    ax[1].scatter(N2,imRes2,s=100, marker="o", color='darkgreen', alpha=0.55, label="sebastian's Data" )
    ax[1].plot(N2,N_linear_fitting_function(N2,imFitA2,imFitB2), linestyle='solid', linewidth=2, color='black' )
    ax[1].axhline(y=imFitA2,linestyle='dashed', linewidth=2, color='green', label="sebastian, $N \\rightarrow \\infty$" )
    ax[1].scatter(N3,imRes3,s=100, marker="s", color='paleturquoise', alpha=0.55, label="straight cont. Data" )
    ax[1].plot(N3,N_linear_fitting_function(N3,imFitA3,imFitB3), linestyle='solid', linewidth=2, color='black' )
    ax[1].axhline(y=imFitA3,linestyle='dashed', linewidth=2, color='grey', label="straight, $N \\rightarrow \\infty$" )

    output_file = "Mphib_contour_comparison_linear_fit_3files.pdf"

    fig.tight_layout()
    fig.legend()
    print("output file = ",output_file)
    plt.savefig(output_file)
    plt.close()


def print_contour_files(opefile,qvecfile,pcutfile,thresholdfile,opecuttracerfile,scount,N):
    #a = 1.42
    sval=0
    #print("s in beginning: ",sval)

    

    x, y, z1, z2 = np.genfromtxt(opefile,unpack=True)
    
    print("OPE file loaded: ",opefile)
    
    xval, pcutpR, pcutpI, pcutmR, pcutmI,req,imq,absqminuspcutp,absqminuspcutm = np.genfromtxt(pcutfile,unpack=True)
    qvecind, qvecx, qvecy = np.genfromtxt(qvecfile,unpack=True)

    

    print("pcut file loaded: ",pcutfile)
    print("qvec file loaded: ",qvecfile)

    (reqplusplus1, imqplusplus1,
     reqplusminus1, imqplusminus1,
     reqminusplus1, imqminusplus1,
     reqminusminus1, imqminusminus1,
     req1, imq1, reMbcleftplus, imMbcleftplus, 
     reMbcleftminus, imMbcleftminus,
     reMbcrightplus, imMbcrightplus, 
     reMbcrightminus, imMbcrightminus,
     res3,ims3,eps,a ) = np.genfromtxt(thresholdfile,unpack=True)

    print("threshold file loaded: ",thresholdfile)

    ( xval1, pcutpR1, pcutpI1, pcutmR1, pcutmI1, 
     req2, imq2, absqminusqcutp1, absqminusqcutm1 ) = np.genfromtxt(opecuttracerfile,unpack=True)
    
    print("cuttracer file loaded: ",opecuttracerfile)

    sval = res3 #+ scount*dels3
    sval2 = '%0.3f'%sval 
    #sigkval = '%0.3f'%sigk
    print("Re s = ",sval," Im s = ",str(ims3))
    
    #reqval0 = '%0.3f'%reactualq
    #imqval0 = '%0.3f'%imactualq

    rekval0 = '%0.3f'%req[0]
    imkval0 = '%0.6f'%imq[0]

    
    xi = np.linspace(x.min(), x.max(), N)
    yi = np.linspace(y.min(), y.max(), N)
    Xi,Yi = np.meshgrid(xi,yi)
    zi = scipy.interpolate.griddata((x, y), z2, (Xi, Yi), method='linear')

    fig, ax = plt.subplots(figsize = (12,8))

    ims3val = '%0.3f'%ims3


    if(ims3==0):
        title_str = "am = " + str(a) + " $s/m^2 = $" + str(sval2) + ", $\epsilon =$" + str(eps) #+ ", $\sigma_k =$" + str(sigkval)
    elif(ims3>0):
        title_str = "am = " + str(a) + " $s/m^2 = $" + str(sval2) + " + " + str(ims3val) + "i" + ", $\epsilon =$" + str(eps) #+ ", $\sigma_k =$" + str(sigkval)
    else:
        title_str = "am = " + str(a) + " $s/m^2 = $" + str(sval2) + str(ims3val) + "i" + ", $\epsilon =$" + str(eps) #+ ", $\sigma_k =$" + str(sigkval)
    fig1_str = "Im $K(p,q), q = \pm$" + str(rekval0) + " $\pm i$ " + str(imkval0) #+ ", $q = \pm$" + str(reqval0) + " $\pm i$ " + str(imqval0)

    fig2_str = "$p_{cut}$"
    #qvecxval = '%0.3f'%qvecx[pcount]
    #qvecyval = '%0.3f'%qvecy[pcount]
    #if(qvecy[pcount]==0):
    #    fig2_str = "$p_{cut}$,  $k = $" + str(qvecxval) 
    #elif(qvecy[pcount]>0):
    #    fig2_str = "$p_{cut}$,  $k = $" + str(qvecxval) + " $ + $ " + str(qvecyval) + " i"
    #else:
    #    fig2_str = "$p_{cut}$,  $k = $" + str(qvecxval) + " " + str(qvecyval) + " i"
    
    #title_str = "$p_{cut}$" + ", am = " + str(a) + " $s/m^2 = $" + str(sval) +  ", $\epsilon =$" + str(eps) 
    
    #qmom = qvecx[somenum] + j*qvecy[somenum]
    #qvecxval = '%0.3f'%qvecx[somenum]
    ##qvecyval = '%0.3f'%qvecy[somenum]
    #if qvecy[somenum]>=0.0:
    #    title_str = "Im $K(q,k)$" + ", am = " + str(a) + " $s/m^2 = $" + str(res) + ",$p = $" + str(qvecxval)+ "+" +str(qvecyval)+ "i" + ", $\epsilon =$" + str(eps) 
    #else:
    #    title_str = "Im $K(q,k)$" + ", am = " + str(a) + " $s/m^2 = $" + str(res) + ",$p = $" + str(qvecxval) + str(qvecyval)+ "i" + ", $\epsilon =$" + str(eps) 
    

    

    fig.suptitle(title_str,fontsize=20)
    ax.set_title(fig1_str,fontsize=20)
    #ax[1].set_title(fig2_str,fontsize=20)

    ax2xlimvalmin = req[0] - 0.01
    ax2xlimvalmax = req[0] + 0.01
    ax2ylimvalmin = imq[0] - 0.01
    ax2ylimvalmax = imq[0] + 0.01
    ax.set_xlim([x.min(),x.max()])
    ax.set_ylim([y.min(),y.max()])
    #ax[1].set_xlim([ax2xlimvalmin,ax2xlimvalmax])
    #ax[1].set_ylim([ax2ylimvalmin,ax2ylimvalmax])
    #ax[2].set_xlim([-1,1])
    #ax[2].set_ylim([-0.01,0.01])


    ax.set_xlabel("Re $p$",fontsize=20)
    ax.set_ylabel("Im $p$",fontsize=20)
    #ax[1].set_xlabel("Re $p$",fontsize=20)
    #ax[1].set_ylabel("Im $p$",fontsize=20)
    #ax[2].set_xlabel("$x$",fontsize=20)
    #ax[2].set_ylabel("$|q - p_{cut}(x)|$",fontsize=20)
    #ax.set_zlim(-500,500)
    #tickz = [-500,-400,-300,-200,-100,0,100,200,300,400,500]

    
    
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.tick_params(axis='both', which='minor', labelsize=20)
    #ax[1].tick_params(axis='both', which='major', labelsize=20)
    #ax[1].tick_params(axis='both', which='minor', labelsize=20)
    #ax[2].tick_params(axis='both', which='major', labelsize=20)
    #ax[2].tick_params(axis='both', which='minor', labelsize=20)

    
    #Here we plot the kernel(p,q) with qvec and thresholds
    print("z min = ",zi.min())
    print("z max = ",zi.max())  

    ax.contour(Xi, Yi, zi,levels=np.linspace(zi.min(),zi.max(),num=50), linewidths=0.5, colors='k',linestyles='solid')
    h0 = ax.contourf(Xi, Yi, zi,levels=np.linspace(zi.min(),zi.max(),num=20))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right","5%", pad="3%")
    #cax.yaxis.limit_range_for_scale(-500,500)
    cbar = fig.colorbar(h0, cax=cax, orientation="vertical")
    cax.xaxis.set_ticks_position("default")
    cax.xaxis.set_label_position("bottom")
    
    cax.tick_params(labelsize=18)


    ax.plot(qvecx,qvecy,linestyle="solid",linewidth=2,color="yellow")
    
    #ax.plot(-req,-imq, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="darkviolet")
    #ax.plot(reactualq,imactualq, marker="o", fillstyle="none", markersize=10, linewidth=1.5, markeredgecolor="white", markerfacecolor="darkviolet")
    #ax.plot(-reactualq,-imactualq, marker="o", fillstyle="none", markersize=10, linewidth=1.5, markeredgecolor="white", markerfacecolor="darkviolet")
    
    
    #ax.plot(reqminus,imqminus, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="cyan",label="$q_-$")
    #ax.plot(-reqplus,-imqplus, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="yellow")
    #ax.plot(-reqminus,-imqminus, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="cyan")
    ax.plot(pcutpR1,pcutpI1,linestyle="solid",linewidth=2,color="red")
    ax.plot(pcutmR1,pcutmI1,linestyle="solid",linewidth=2,color="blue")
    ax.plot(reqplusplus1,imqplusplus1, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="red",label="$q_+$")
    ax.plot(reqplusminus1,imqplusminus1, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="red",label="$q_+$")
    ax.plot(reqminusplus1,imqminusplus1, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="blue",label="$q_+$")
    ax.plot(reqminusminus1,imqminusminus1, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="blue",label="$q_+$")
    
    ax.plot(reMbcrightplus,imMbcrightplus, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="orange")
    ax.plot(reMbcrightminus,imMbcrightminus, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="orange")
    ax.plot(req,imq, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="darkviolet", label="$q$")
    ax.plot(-req,-imq, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="darkviolet", label="$q$")
    

    
    


    

    #qplusval = complex(reqplus, imqplus)
    #qminusval = complex(reqminus,imqminus) 
    #qval = complex(req,imq)
    #textstr = '\n'.join((
    #r'$q_+={0.real:.3f} + {0.imag:.3f}i$'.format(qplusval),
    #r'$q_-={0.real:.3f} + {0.imag:.3f}i$'.format(qminusval),
    #r'$-q_+={0.real:.3f} + {0.imag:.3f}i$'.format(-qplusval),
    #r'$-q_-={0.real:.3f} + {0.imag:.3f}i$'.format(-qminusval),
    #r'$q={0.real:.3f} + {0.imag:.3f}i$'.format(qval)))

    

    #ax.plot(bc_left_plus, 0, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="orange")
    #ax.plot(bc_left_minus, 0, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="orange")
    #ax.plot(0,imMbcrightplus, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="orange",label="$\mathcal{M}_2$ branch cut")
    #ax.plot(0,imMbcrightminus, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="orange")
    #ax.plot(req,imq, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="red", label="$q$")
    #ax.plot(-req,-imq, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="red")
    #ax.plot(reqplus,imqplus, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="yellow",label="$q_+$")
    #ax.plot(reqminus,imqminus, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="cyan",label="$q_-$")
    #ax.plot(-reqplus,-imqplus, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="yellow")
    #ax.plot(-reqminus,-imqminus, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="cyan")
    
    #ax[0].plot(pcutpR,pcutpI,color="red")
    #ax[0].plot(pcutmR,pcutmI,color="blue")
    
    #ax[1].plot(pcutpR,pcutpI, marker="o", markersize=2, markeredgecolor="gray", markerfacecolor="red",alpha=0.2)
    #ax[1].plot(pcutmR,pcutmI, marker="o", markersize=2, markeredgecolor="gray", markerfacecolor="blue",alpha=0.2)
   #ax[1].scatter(pcutmR,pcutmI, marker="o",s=2,edgecolors=None, color="blue",alpha=0.1)
   #ax[1].plot(req,imq, marker="o", markersize=3, markeredgecolor="black", markerfacecolor="darkviolet")
    #ax[1].plot(-req,-imq, marker="o", markersize=15, markeredgecolor="black", markerfacecolor="darkviolet")
    #ax[1].plot(qvecx,qvecy,linestyle="solid",marker="o",linewidth=0.5,color="green")
    #ax[1].plot(qvecx[pcount],qvecy[pcount],marker="o", markersize=10, markeredgecolor="orange", markerfacecolor="orange")

    line1x,line1y = [-1,1],[0,0]

    #ax[2].plot(xval,absqminuspcutp, marker="o", markersize=5, markeredgecolor="darkred", markerfacecolor="red",alpha=0.3)
    #ax[2].plot(xval,absqminuspcutm, marker="o", markersize=5, markeredgecolor="darkblue", markerfacecolor="blue",alpha=0.3)

    #ax[2].plot(line1x,line1y, linestyle='solid' )
    #ax.plot(pcutplusvecx,pcutplusvecy, marker="o", markersize=0.5, markeredgecolor="white", markerfacecolor="white")
    #ax.plot(pcutminusvecx,pcutminusvecy, marker="o", markersize=0.5, markeredgecolor="gold", markerfacecolor="white")
    
    #ax.plot(reqc1,imqc1, marker="o", markersize=20, markeredgecolor="orange", markerfacecolor="yellow")
    #ax.plot(reqc2,imqc2, marker="o", markersize=20, markeredgecolor="green", markerfacecolor="cyan")
    
    #ax.plot(sigvecx,sigvecy, marker="o", markersize=1, markeredgecolor="blue", markerfacecolor="blue")
    

    #ax.plot(np.linspace(-0.2,0.2,100),np.linspace(0,0,100), marker="o", markersize=10, markeredgecolor="black", markerfacecolor="white")
    #ax.plot(qvecx,qvecy, marker="o", markersize=5, markeredgecolor="black", markerfacecolor="white")
    
    #ax.plot(0,1.0, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="red")
    #ax.plot(0,-1.0, marker="o", markersize=10, markeredgecolor="black", markerfacecolor="red")
    #props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
    #    verticalalignment='top', bbox=props)

    outputfile_str = "glockletest_" + pcutfile + ".pdf"

    print("output file: ",outputfile_str)
    
    #fig.tight_layout()
    plt.savefig(outputfile_str)
    plt.close()


scountinitial = 0
scountfinal = 1
pcountinitial = 0
pcountfinal = 204



#for scount in range(0,301,1):
    #eps = 1.0e-5
    #simag = 0.0
#    N = 250
#    sval=0
#    opefile = 'kernel_scount_' + str(scount) + '.dat'
    #qvecfile = 'qvec_scount_' + str(scount) + '.dat'
#    qvecfile = 'contour_smooth.dat"
#    pcutfile = 'pcut_scount_' + str(scount) + '.dat'
#    thresholdfile = 'threshold_scount_' + str(scount) + '.dat'
#    opecuttracerfile = 'opecuttracer_scount_' + str(scount) + '.dat'
#    if(os.path.exists(opefile)):
#        print_contour_files(opefile,qvecfile,pcutfile,thresholdfile,opecuttracerfile,scount,N)
#    else:
#        continue 
    #print("for scount ",scount,", files have been created = ",(pcount+1)," out of ",(pcountfinal+1))
    #scount = 0
    #pcount = 11
    #dels = 0.0256555
    #sinitial = 8.72

N = 250
sval=0
opefile = 'kernel_scount_' + str(0) + '.dat'
qvecfile = 'qvec_scount_' + str(0) + '.dat'
#qvecfile = 'contour_smooth_1.dat'
pcutfile = 'pcut_scount_' + str(0) + '.dat'
thresholdfile = 'threshold_scount_' + str(0) + '.dat'
opecuttracerfile = 'opecuttracer_scount_' + str(0) + '.dat'
if(os.path.exists(opefile)):
    print_contour_files(opefile,qvecfile,pcutfile,thresholdfile,opecuttracerfile,0,N)

#filename1 = "Mphib_s=8.600000+0.100000i.dat"
#filename2 = "test_8.6+0.1i.dat"
#filename3 = "straightMphib_s=8.6-0.05i.dat"
#plotFit(filename1, filename2)
#plotFit3(filename1,filename2,filename3)