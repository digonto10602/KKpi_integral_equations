import numpy as np 
import math 

a = 10
b = 11
delc = (math.log(a) - math.log(b))/10
newa = a 
for i in range(0,10,1):
    
    newa = newa + math.exp(-delc)
    delc = np.abs((math.log(newa) - math.log(b)))/(10)

    print(newa,delc,math.exp(-delc))
print(delc)