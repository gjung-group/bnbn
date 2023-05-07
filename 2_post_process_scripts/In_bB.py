#!/usr/bin/env python
# coding: utf-8

# In[7]:


import numpy as np
import math
import matplotlib.pyplot as plt

def function_Phi(x,y,A,B,C,ForG):
    alpha=-4.5
    beta=-2.598076211353316
    gamma=4.5
    delta=beta

    if np.abs(B-C) < 0.0000001:
        print('we have a singularity\n')
        D = (A-B)/(10^(-16)) 
    else:
        D = (A-B)/(B-C)
    phi = math.atan((1.0/(delta/beta*D-1.0)*((delta*alpha-beta*gamma)/(beta*delta)))-(gamma)/(delta))

    if np.abs(B-C) < 0.0000001:
        c1 = (10^(-16))/(2.0*(gamma*np.cos(phi)+delta*np.sin(phi)))

    else:
        c1= (B-C)/(2.0*(gamma*np.cos(phi)+delta*np.sin(phi)))

    
    Ax = 0.0;
    Ay = 1/np.sqrt(3);
    G1 = 4.0*np.pi/(np.sqrt(3.0));
    c0 = A - 2.0 * c1 * np.cos(phi - G1 * Ay) - 4.0 * c1 * np.cos(G1 * Ay / 2.0 + phi) * np.cos(np.sqrt(3.0) * G1 * Ax / 2.0)
    aBN = 2.479
    G1 = 4.0*np.pi/(np.sqrt(3.0)*aBN)
    if ForG == 1:
        f1 = 2.0*c1*np.cos(phi-G1*y) + 4.0*c1*np.cos(G1*y/2.0 + phi)*np.cos(np.sqrt(3.0)*G1*x/2.0)
        out = c0 + f1
    else:
        #out = c0*np.cos(np.sqrt(3)*G1*x/2)*np.cos(G1*y/2-phi)-2*c1*np.cos(G1*y+phi)-1j*2*np.sqrt(3)*c1*np.sin(np.sqrt(3)*G1*x/2)*np.sin(G1*y/2-phi)
        out = 2*c1*np.cos(np.sqrt(3)*G1*x/2)*np.cos(G1*y/2-phi)-2*c1*np.cos(G1*y+phi)-1j*2*np.sqrt(3)*c1*np.sin(np.sqrt(3)*G1*x/2)*np.sin(G1*y/2-phi)

    return out

xxx=np.loadtxt('dx_B.dat')
yyy=np.loadtxt('dy_B.dat')
dxv=np.loadtxt('bB',usecols=(0,))
dyv=np.loadtxt('bB',usecols=(1,))
# dxv=np.loadtxt('tN',usecols=(0,))
# dyv=np.loadtxt('tN',usecols=(1,))

#top layer

# B = -2.7001; #%%AA stacking -ab
# C = -2.7161; #% AB stacking-ab
# A = -2.6971; #%BA stacking-ab

# B = 1.7666; #%AA stacking -aa
# C = 1.7169; #% AB stacking-aa
# A = 1.6636; #%BA stacking-aa

#B = -2.1843; #AA stacking -bb
#C = -2.2444; # AB stacking-bb
#A = -2.3393; #BA stacking-bb

#bottom layer

# B = -2.7001; #%%AA stacking -ab
# C = -2.6954; #% AB stacking-ab
# A = -2.7190; #%BA stacking-ab

B = 1.7666; #%AA stacking -aa
C = 1.6664; #% AB stacking-aa
A = 1.7128; #%BA stacking-aa

#B = -2.1843; #AA stacking -bb
#C = -2.3294; # AB stacking-bb
#A = -2.2591; #BA stacking-bb


OnsiteC1 = np.zeros((np.size(xxx),np.size(yyy)))

indexdx = np.size(dxv)

indexdy = np.size(dyv)

OnsiteC1 = function_Phi(yyy,xxx,A,B,C,1)
np.savetxt('onsite_bB',OnsiteC1)
# print(OnsiteC1)
#plt.scatter(dxv,dyv,c=np.real(OnsiteC1),s=55,cmap='RdBu',edgecolor='none')
#plt.colorbar()                   


# In[3]:


#x1 = 0.0
#y1 = 1/np.sqrt(3)
#x2 = 0.0
#y2 = 0.0
#G1 = 4.0*np.pi/(np.sqrt(3.0))
#m1 = np.cos(np.sqrt(3.0)*G1*x1/2.0)
#m2 = np.cos(np.sqrt(3.0)*G1*x2/2.0)
#a1 = np.cos(G1*y1) - np.cos(G1*y2)
#a3 = 2.0*np.cos(G1*y1/2.0) * m1 - 2.0*np.cos(G1*y2/2.0)* m2
#alpha = a1 + a3
#a2 = np.sin(G1*y1) - np.sin(G1*y2)
#a4 = 2.0*np.sin(G1*y1/2.0) * m1 - 2.0*np.sin(G1*y2/2.0)* m2
#beta = a2 - a4
#delta=beta
#gamma=-alpha
#     phi = math.atan((1.0/(delta/beta*D-1.0)*((delta*alpha-beta*gamma)/(beta*delta)))-(gamma)/(delta));

#print(alpha,beta,gamma)


# In[15]:



    


# In[ ]:




