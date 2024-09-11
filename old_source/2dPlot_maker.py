import numpy as np 
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.optimize import curve_fit
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable 
import os.path


#x = np.linspace(1,10,10)
#y = np.linspace(1,10,10)

#for i in range(0,len(x),1):
#    print(x[i],y[i])


def test_raul_plot(filename):
    plt.rcParams.update({'font.size': 20})
    plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)
    (res,ims,N,ReM,ImM) = np.genfromtxt(filename,unpack=True)
    #(sebares,sebaReM,sebaImM) = np.genfromtxt(sebadata,unpack=True)

    x = np.empty(101)
    y = np.empty(101)
    z1 = np.empty(101)
    z2 = np.empty(101) 

    somecount = 101
    
    fig, ax = plt.subplots(2,1,figsize = (12,10))

    title_str = "am = " + str(16) + ", Im $(s/m^2) = $ " + str(ims[0]) + ", N = " + str(N[0]) 
    fig.suptitle(title_str,fontsize=20)

    ax[0].set_xlabel("$N$")
    ax[0].set_ylabel("Re $\mathcal{M}_{\phi b}$")

    ax[1].set_xlabel("$N$")
    ax[1].set_ylabel("Im $\mathcal{M}_{\phi b}$")

    ax[0].tick_params(axis='both', which='major', labelsize=20)
    ax[0].tick_params(axis='both', which='minor', labelsize=20)

    ax[1].tick_params(axis='both', which='major', labelsize=20)
    ax[1].tick_params(axis='both', which='minor', labelsize=20)
    
    ax[0].scatter(N,ReM,s=100,color='darkred',alpha=0.55,label="data")
    ax[1].scatter(N,ImM,s=100,color='teal',alpha=0.55,label="data")
    #ax.plot(sebares,sebaReM,color='darkred',linewidth=3,label="real, sebastian")
    #ax.plot(sebares,sebaImM,color='teal',linewidth=3,label="imag, sebastian")
    plt.tight_layout()
    plt.legend()

    output_str = filename + ".pdf"
    #plt.show()
    plt.savefig(output_str)
    plt.close()
    print("output file = ",output_str)

 

def test_sebadata_M2d_printer(filename, sebadata):
    plt.rcParams.update({'font.size': 20})
    plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)
    (res,ims,ReM,ImM,Red,Imd,ReM2,ImM2,ReM2denom,ImM2denom,a,N) = np.genfromtxt(filename,unpack=True)
    (sebares,sebaReM,sebaImM) = np.genfromtxt(sebadata,unpack=True)

    x = np.empty(101)
    y = np.empty(101)
    z1 = np.empty(101)
    z2 = np.empty(101) 

    somecount = 101
    
    fig, ax = plt.subplots(figsize = (12,5))

    title_str = "am = " + str(a[0]) + ", Im $(s/m^2) = $ " + str(ims[0]) + ", N = " + str(N[0]) 
    fig.suptitle(title_str,fontsize=20)

    ax.set_xlabel("Re ($s/m^2$)")
    ax.set_ylabel("$d_S$")

    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.tick_params(axis='both', which='minor', labelsize=20)
    
    ax.scatter(res,ReM,color='darkred',alpha=0.25,label="real, digonto")
    ax.scatter(res,ImM,color='teal',alpha=0.25,label="imag, digonto")
    ax.plot(sebares,sebaReM,color='darkred',linewidth=3,label="real, sebastian")
    ax.plot(sebares,sebaImM,color='teal',linewidth=3,label="imag, sebastian")
    plt.tight_layout()
    plt.legend()

    output_str = "sebastian_test_N=500_simag=-0.2.pdf"
    #plt.show()
    plt.savefig(output_str)
    plt.close()
    print("output file = ",output_str)

 

def M2d_printer(filename):
    plt.rcParams.update({'font.size': 20})
    plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)
    (res,ims,ReM,ImM,Red,Imd,ReM2,ImM2,ReM2denom,ImM2denom,a,N) = np.genfromtxt(filename,unpack=True)

    x = np.empty(101)
    y = np.empty(101)
    z1 = np.empty(101)
    z2 = np.empty(101) 

    somecount = 101
    
    fig, ax = plt.subplots(figsize = (12,5))

    title_str = "am = " + str(a[0]) + ", Im $(s/m^2) = $ " + str(ims[0]) + ", N = " + str(N[0]) 
    fig.suptitle(title_str,fontsize=20)

    ax.set_xlabel("Re ($s/m^2$)")
    ax.set_ylabel("$d_S$")

    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.tick_params(axis='both', which='minor', labelsize=20)
    
    ax.plot(res,Red,color='darkred',linewidth=3,label="real")
    ax.plot(res,Imd,color='teal',linewidth=3,label="imag")
    plt.tight_layout()
    plt.legend()

    output_str = filename +  ".pdf"
    #plt.show()
    plt.savefig(output_str)
    plt.close()
    print("output file = ",output_str)

        



def M2d_from_M3d_printer(filename):
    plt.rcParams.update({'font.size': 20})
    plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)
    (res,ims,ReM,ImM,Red,Imd) = np.genfromtxt(filename,unpack=True)

    x = np.empty(101)
    y = np.empty(101)
    z1 = np.empty(101)
    z2 = np.empty(101) 

    somecount = 101
    for i in range(0,101,1):
        for j in range(0,101,1):
            #print(i,j,res[somecount*i + j],ims[somecount*i + j],Red[somecount*i + j],Imd[somecount*i + j])
            x[j] = res[somecount*i + j]
            y[j] = ims[somecount*i + j]

            z1[j] = Red[somecount*i + j]

            z2[j] = Imd[somecount*i + j]

        fig, ax = plt.subplots(figsize = (12,5))

        title_str = "am = " + str(16) + ", Im $(s/m^2) = $ " + str(y[0]) 
        fig.suptitle(title_str,fontsize=20)

        ax.set_xlabel("Re ($s/m^2$)")
        ax.set_ylabel("$d_S$")

        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.tick_params(axis='both', which='minor', labelsize=20)
    
        ax.plot(x,z1,color='darkred',linewidth=3,label="real")
        ax.plot(x,z2,color='teal',linewidth=3,label="imag")
        plt.tight_layout()
        plt.legend()

        output_str = "plot_" + str(i) +  ".pdf"
        #plt.show()
        plt.savefig(output_str)
        plt.close()
        print("output file = ",output_str)


def Re_M3d_printer(filename):
    plt.rcParams.update({'font.size': 20})
    plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)
    (res,ims,ReM,ImM,Red,Imd,sigk) = np.genfromtxt(filename,unpack=True)

    
    x = np.linspace(res.min(),res.max(),100)
    y = np.linspace(ims.min(),ims.max(),100)
    print(len(x),len(y))

    print(res.min(),res.max())
    print(ims.min(),ims.max())


    X1,Y1 = np.meshgrid(x,y)
    X2,Y2 = np.meshgrid(x,y)
    Z1 = scipy.interpolate.griddata((res,ims),Red,(X1,Y1), method='linear')
    #Z2 = scipy.interpolate.griddata((res,ims),ImM,(X2,Y2), method='linear')


    for i in range(0,len(X1)):
        for j in range(0,len(Y1)):
            if Z1[i,j]<-1000:
                Z1[i,j] = 'nan'
            if Z1[i,j]>1000:
                Z1[i,j] = 'nan'
            print(i,j,Z1[i,j])

    #for i in range(0,len(X2)):
    #    for j in range(0,len(Y2)):
    #        if Z2[i,j]<-1000:
    #            Z2[i,j] = 'nan'
    #        if Z2[i,j]>1000:
    #            Z2[i,j] = 'nan'
    #        print(i,j,Z2[i,j])
    #print(z)

    fig, ax = plt.subplots(figsize = (12,8))

    ax = plt.axes(projection='3d')
    #ax[0] = plt.axes(projection='3d')
    #ax[1] = plt.axes(projection='3d')
    #ax.set_zlim(-10000,10000)


    title_str = "am = " + str(16) + ", $\sigma_p = \sigma_k = $" + str(sigk[0])
    fig.suptitle(title_str,fontsize=20)
    #ax.set_title(fig1_str,fontsize=20)

    ax.set_zlim3d(-500,500)
    ax.set_xlabel("Re$(s)$")
    ax.set_ylabel("Im$(s)$")
    ax.set_zlabel("Re ($d_S$)")
    h1 = ax.plot_surface(X1,Y1,Z1, cmap='viridis',lw=0.5,  alpha=0.5)
    h3 = ax.contour(X1,Y1,Z1, colors='k',levels=250,linestyle='solid',linewidth=2.0)
    #h1 = ax.plot_surface(X1,Y1,Z1, cmap='viridis')

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.tick_params(axis='both', which='minor', labelsize=16)
    ax.tick_params(labelsize=12)

    #ax[0].set_zlim3d(-7000,7000)
    #ax[0].set_xlabel("Re$(s)$")
    #ax[0].set_ylabel("Im$(s)$")
    #ax[0].set_zlabel("$\mathcal{M}_{\phi b}$")
    #h1 = ax[0].plot_surface(X1,Y1,Z1, cmap='viridis')

    #ax.set_zlim(-10000,10000)
    #ax[1].set_zlim3d(-7000,7000)
    #ax[1].set_xlabel("Re$(s)$")
    #ax[1].set_ylabel("Im$(s)$")
    #ax[1].set_zlabel("$\mathcal{M}_{\phi b}$")
    #h2 = ax[1].plot_surface(X2,Y2,Z2, cmap='viridis')


    fig.colorbar(h1, shrink=0.5, aspect=8)

    #fig.colorbar(h1, shrink=0.5, aspect=8)
    #fig.colorbar(h2, shrink=0.5, aspect=8)
    plt.tight_layout()

    output_str = "Re_" + filename + ".pdf"
    #plt.show()
    plt.savefig(output_str)
    plt.close()

def Im_M3d_printer(filename):
    plt.rcParams.update({'font.size': 20})
    plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)
    (res,ims,ReM,ImM,Red,Imd,sigk) = np.genfromtxt(filename,unpack=True)


    x = np.linspace(res.min(),res.max(),100)
    y = np.linspace(ims.min(),ims.max(),100)
    print(len(x),len(y))

    print(res.min(),res.max())
    print(ims.min(),ims.max())


    X1,Y1 = np.meshgrid(x,y)
    X2,Y2 = np.meshgrid(x,y)
    #Z1 = scipy.interpolate.griddata((res,ims),Red,(X1,Y1), method='linear')
    Z2 = scipy.interpolate.griddata((res,ims),Imd,(X2,Y2), method='linear')


    #for i in range(0,len(X1)):
    #    for j in range(0,len(Y1)):
    #        if Z1[i,j]<-1000:
    #            Z1[i,j] = 'nan'
    #        if Z1[i,j]>1000:
    #            Z1[i,j] = 'nan'
    #        print(i,j,Z1[i,j])

    for i in range(0,len(X2)):
        for j in range(0,len(Y2)):
            if Z2[i,j]<-1000:
                Z2[i,j] = 'nan'
            if Z2[i,j]>1000:
                Z2[i,j] = 'nan'
            print(i,j,Z2[i,j])
    #print(z)

    fig, ax = plt.subplots(figsize = (12,8))

    ax = plt.axes(projection='3d')
    #ax[0] = plt.axes(projection='3d')
    #ax[1] = plt.axes(projection='3d')
    #ax.set_zlim(-10000,10000)


    title_str = "am = " + str(16) + ", $\sigma_p = \sigma_k = $" + str(sigk[0])
    fig.suptitle(title_str,fontsize=20)
    #ax.set_title(fig1_str,fontsize=20)

    ax.set_zlim3d(-500,500)
    ax.set_xlabel("Re$(s)$")
    ax.set_ylabel("Im$(s)$")
    ax.set_zlabel("Im ($d_S$)")
    h1 = ax.plot_surface(X1,Y1,Z2, cmap='viridis',lw=0.5,  alpha=0.5)
    h3 = ax.contour(X1,Y1,Z2, colors='k',levels=250,linestyle='solid',linewidth=2.0)
    #h1 = ax.plot_surface(X1,Y1,Z1, cmap='viridis')

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.tick_params(axis='both', which='minor', labelsize=16)
    ax.tick_params(labelsize=12)

    #ax[0].set_zlim3d(-7000,7000)
    #ax[0].set_xlabel("Re$(s)$")
    #ax[0].set_ylabel("Im$(s)$")
    #ax[0].set_zlabel("$\mathcal{M}_{\phi b}$")
    #h1 = ax[0].plot_surface(X1,Y1,Z1, cmap='viridis')

    #ax.set_zlim(-10000,10000)
    #ax[1].set_zlim3d(-7000,7000)
    #ax[1].set_xlabel("Re$(s)$")
    #ax[1].set_ylabel("Im$(s)$")
    #ax[1].set_zlabel("$\mathcal{M}_{\phi b}$")
    #h2 = ax[1].plot_surface(X2,Y2,Z2, cmap='viridis')


    fig.colorbar(h1, shrink=0.5, aspect=8)

    #fig.colorbar(h1, shrink=0.5, aspect=8)
    #fig.colorbar(h2, shrink=0.5, aspect=8)
    plt.tight_layout()

    output_str = "Im_" + filename + ".pdf"
    #plt.show()
    plt.savefig(output_str)
    plt.close()

#print(x)
#print(X)
filename0 = "fixeds3imag_Mphib_a_16_N_500_contour43_1.dat"
sebadata = "test_boom.dat"
#filename1 = "Msigk2msq_momrep_a_16_N_1000_3d_withoutomp_qss.dat"
filename2 = "fixeds3imag_Mphib_a_16_N_5000_8500.dat"
filename3 = "Msigk3.85msq_momrep_a_16_N_1000_3d_withoutomp_qss.dat"


filename1 = "Mphib_s=8.6-0.1i.dat"
filename2 = "Mphib_s=8.6-0.05i.dat"
filename3 = "Mphib_s=8.6-0.15i.dat"

#test_sebadata_M2d_printer(filename0, sebadata)
#Re_M3d_printer(filename2)
#M2d_printer(filename0)
#M2d_from_M3d_printer(filename2)

test_raul_plot(filename1)
test_raul_plot(filename2)
test_raul_plot(filename3)