import numpy as np 


#filename = "vertexfactor_vs_k_for_a_16_N_1000.dat"
filename = "vertexfactor_vs_sigk_BS3_file_a=10000_N=3000.dat"

(sigkval,vertexval,vertexsqval,sbval,spoleval,qpoleval,mval,cval,Nval) = np.genfromtxt(filename, unpack=True)
#(qpoleval,m2kressq,vertexval,vertexsqval) = np.genfromtxt(filename, unpack=True)

sumval = 0
for i in range(0,len(vertexval),1):
    if(i==0):
        sumval = sumval + vertexval[i]*(qpoleval[i])
    else:
        sumval = sumval + vertexval[i]*(qpoleval[i]-qpoleval[i-1])

#filename1 = "normalized_" + filename
#f = open(filename1, "w")

vertexval = vertexval/sumval;
vertexsqval = vertexsqval/sumval;

for i in range(0,len(vertexval),1):
    #print(sigkval[i]  + "\t" + vertexval[i] + "\t" + vertexsqval[i] + "\t" + sbval[i] + "\t" + spoleval[i] + "\t" + qpoleval[i] + "\t" + mval[i] + "\t" + cval[i] + "\t" + Nval[i] +  "\n" )
    #print(sigkval[i],vertexval[i],vertexsqval[i],sbval[i],spoleval[i],qpoleval[i],mval[i],cval[i],Nval[i])
    #print(qpoleval[i],m2kressq[i],vertexval[i],vertexsqval[i])
    print(qpoleval[i],vertexval[i],vertexsqval[i])
    
#print("sum = ",sumval)
#f.close()