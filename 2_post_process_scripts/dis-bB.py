#!/usr/bin/env python
# coding: utf-8

# In[8]:


import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import griddata
# b_x = np.loadtxt('tB',usecols=(0,))
# b_y = np.loadtxt('tB',usecols=(1,))
# b_z = np.loadtxt('tB',usecols=(2,))

# t_x = np.loadtxt('bB.dat',usecols=(0,))
# t_y = np.loadtxt('bB.dat',usecols=(1,))
# t_z = np.loadtxt('bB.dat',usecols=(2,))

b_x = np.loadtxt('tB',usecols=(0,))
b_y = np.loadtxt('tB',usecols=(1,))
b_z = np.loadtxt('tB',usecols=(2,))

t_x = np.loadtxt('bB',usecols=(0,))
t_y = np.loadtxt('bB',usecols=(1,))
t_z = np.loadtxt('bB',usecols=(2,))

del_x = []
del_y = []
for elx, ely in zip(b_x,b_y):
    del_r_min = 1.43
    for elx2, ely2 in zip(t_x, t_y):
        del_r = np.sqrt((elx-elx2)**2 + (ely-ely2)**2)
        if del_r < del_r_min:
            del_x_min = elx-elx2
            del_y_min = ely-ely2
            del_r_min = del_r
    del_x.append(del_x_min)
    del_y.append(del_y_min)

del_r = np.sqrt(np.asarray(del_x)**2+np.asarray(del_y)**2)

np.savetxt('dx_B.dat',del_x)
np.savetxt('dy_B.dat',del_y)
#np.savetxt('dis_B.dat',del_r)
#plt.scatter(b_x,b_y,c=del_r,s=55,cmap='RdBu',edgecolor='none')
#plt.colorbar()


# In[15]:


# t_y = np.loadtxt('bB.dat',usecols=(1,))
#d_r = np.loadtxt('dis.dat',usecols=(0,))

#dd_r=d_r.flatten()
#np.savetxt('dd.dat',dd_r)
#np.shape(dd_r)


# In[ ]:




