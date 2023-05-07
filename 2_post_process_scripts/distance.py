import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import griddata
t_x = np.loadtxt('top.xyz',usecols=(0,))
b_x = np.loadtxt('bot.xyz',usecols=(0,))
t_y = np.loadtxt('top.xyz',usecols=(1,))
b_y = np.loadtxt('bot.xyz',usecols=(1,))
t_z = np.loadtxt('top.xyz',usecols=(2,))
b_z = np.loadtxt('bot.xyz',usecols=(2,))
order = np.argsort(b_z)
f_interp = scipy.interpolate.griddata((t_x, t_y), t_z, (b_x, b_y),method="linear")
diff = f_interp-b_z
diff = diff[~np.isnan(diff)]
dis=np.average(diff)/2
dis=dis.flatten()
np.savetxt('dis.txt',dis,fmt="%s")
#f=open("dis.txt","w")
#f.writelines(dis)
#f.close()
