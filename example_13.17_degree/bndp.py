
import Quest as Quest
import numpy as np
import matplotlib.pyplot as plt 

a = 10.919041854
T = a*np.array([[1,0],[-0.5,np.sqrt(3)/2]])
point = Quest.rspace(T)
tB_file = np.genfromtxt("tB")
for x,y,z in tB_file:
    point.add_orbit(x=x,y=y,z=z,label="tB")
tN_file = np.genfromtxt("tN")
for x,y,z in tN_file:
    point.add_orbit(x=x,y=y,z=z,label="tN")
bB_file = np.genfromtxt("bB")
for x,y,z in bB_file:
    point.add_orbit(x=x,y=y,z=z,label="bB")
bN_file = np.genfromtxt("bN")
for x,y,z in bN_file:
    point.add_orbit(x=x,y=y,z=z,label="bN")
    
gamma_0_aa=-2.700 #first neighbor interaction
gamma_1_aa=0.8310
gamma_0_aap=-0.2102

gamma_0_bb=-2.700##irst neighbor interaction
gamma_1_bb=0.6602
gamma_0_bbp=-0.2102


gamma_0_ab=-2.7 #first neighbor interaction
gamma_1_ab=0.3989
gamma_0_abp=-0.2041

a=1.43
a_1=3.261
a_0=2.4795

coeff_aa=(np.log(gamma_0_aap/gamma_0_aa))/(a-a_0)
coeff_bb=(np.log(gamma_0_bbp/gamma_0_bb))/(a-a_0)
coeff_ab=(np.log(gamma_0_abp/gamma_0_ab))/(a-a_0)

q_sigma_aa=coeff_aa*a_1
q_pi_aa=coeff_aa*a
q_sigma_bb=coeff_bb*a_1
q_pi_bb=coeff_bb*a
q_sigma_ab=coeff_ab*a_1
q_pi_ab=coeff_ab*a

def t_aa(X,Y,Z):
    return ((a_1/np.sqrt(X**2+Y**2+Z**2))**2)*(gamma_1_aa*np.exp(q_sigma_aa*(1-(np.sqrt(X**2+Y**2+Z**2)/3.261))))+(1-(a_1/np.sqrt(X**2+Y**2+Z**2))**2)*(gamma_0_aa*np.exp(q_pi_aa*(1-(np.sqrt(X**2+Y**2+Z**2)/1.43))))


def t_bb(X,Y,Z):
    return ((a_1/np.sqrt(X**2+Y**2+Z**2))**2)*(gamma_1_bb*np.exp(q_sigma_bb*(1-(np.sqrt(X**2+Y**2+Z**2)/3.261))))+(1-(a_1/np.sqrt(X**2+Y**2+Z**2))**2)*(gamma_0_bb*np.exp(q_pi_bb*(1-(np.sqrt(X**2+Y**2+Z**2)/1.43))))

def t_ab(X,Y,Z):
    return ((a_1/np.sqrt(X**2+Y**2+Z**2))**2)*(gamma_1_ab*np.exp(q_sigma_ab*(1-(np.sqrt(X**2+Y**2+Z**2)/3.261))))+(1-(a_1/np.sqrt(X**2+Y**2+Z**2))**2)*(gamma_0_ab*np.exp(q_pi_ab*(1-(np.sqrt(X**2+Y**2+Z**2)/1.43))))

H_tot = Quest.H_tot(point,True,True,cutoff=9)

H0 = Quest.H0()
H0.set_hoppings('tB','tB',np.array([0.007210,0.022764,-0.048173]),cutoff=3)
H0.set_hoppings('tN','tN',np.array([0.193699,0.019136,-0.037515]),cutoff=3)
H0.set_hoppings('tB','tN',np.array([-2.720361,-0.211012,0.079806]),cutoff=3)
H0.set_hoppings('bB','bB',np.array([0.007210,0.022764,-0.048173]),cutoff=3)
H0.set_hoppings('bN','bN',np.array([0.193699,0.019136,-0.037515]),cutoff=3)
H0.set_hoppings('bB','bN',np.array([-2.720361,-0.211012,0.079806]),cutoff=3)


onsite = [1.753584]*19+[-2.293102]*19+[1.753484]*19+[-2.293102]*19
onsite = onsite
H0.set_onsite(onsite=np.array(onsite))
H0.set_hoppings('tB','bB',t_aa,isFn=True,cutoff=6)
H0.set_hoppings('tB','bN',t_ab,isFn=True,cutoff=6)
H0.set_hoppings('tN','bB',t_ab,isFn=True,cutoff=6)
H0.set_hoppings('tN','bN',t_bb,isFn=True,cutoff=6)
kpoint = Quest.kspace(T)
HSpoint = [0.5,0],  [0,0] , [1/3,1/3] ,[0.5,0] 
N = 60
N_grids = int(1.46*N),int((2.3-1.46)*N),int((4-2.3)*N)
kpoint.gen_bandpath(np.array(HSpoint),np.array(N_grids))

H_tot += H0

Solver = Quest.Solver()
Solver.Solve(H_tot,kpoint)
Quest.Save_object("KPOINT",kpoint)
Quest.Save_object("EIGS",Solver)
