
import QUEST as Quest
import numpy as np

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
def F2G2_inter(x,y,z):
    r = np.sqrt(x**2 + y**2 + z**2)
    a1 = 2.4795/np.sqrt(3)
    a2 = 2.4795*2/np.sqrt(3)
    a3 = 2.4795*np.sqrt(7/3)
    hopping = np.zeros_like(r)
    cond1 = np.abs(r-a1)<0.2
    hopping[cond1] = -2.6971*np.exp(-2.45*(r[cond1]-a1)/a1)
    cond2 = np.abs(r-a2)<0.2
    hopping[cond2] = -0.2207*np.exp(-2.45*(r[cond2]-a2)/a2)
    cond3 = np.abs(r-a3)<0.2
    hopping[cond3] = 0.0779*np.exp(-2.45*(r[cond3]-a3)/a3)
    return hopping

H0.set_hoppings('tB','tB',np.array([0.00530,0.0223,-0.0483,-0.0029,-0.0033,0.0002]),cutoff=6)
H0.set_hoppings('tN','tN',np.array([0.1923,0.01925,-0.0373,-0.0027,0.000,-0.0009]),cutoff=6)
H0.set_hoppings('tB','tN',F2G2_inter,True,cutoff=6)
H0.set_hoppings('tB','tB',np.array([0.00530,0.0223,-0.0483,-0.0029,-0.0033,0.0002]),cutoff=6)
H0.set_hoppings('tN','tN',np.array([0.1923,0.01925,-0.0373,-0.0027,0.000,-0.0009]),cutoff=6)
H0.set_hoppings('bB','bN',F2G2_inter,True,cutoff=6)


H0.set_hoppings('tB','bB',t_aa,isFn=True,cutoff=6)
H0.set_hoppings('tB','bN',t_ab,isFn=True,cutoff=6)
H0.set_hoppings('tN','bB',t_ab,isFn=True,cutoff=6)
H0.set_hoppings('tN','bN',t_bb,isFn=True,cutoff=6)
onsite = np.loadtxt('onsite')
onsite = onsite
H0.set_onsite(onsite=np.array(onsite))

kpoint = Quest.kspace(T)
HSpoint = [1/3,1/3], [0,0] , [1/2,0] ,[2/3,2/3] 
N = 60
N_grids = int(1.46*N),int((2.3-1.46)*N),int((4-2.3)*N)
kpoint.gen_bandpath(np.array(HSpoint),np.array(N_grids))

H_tot += H0
Solver = Quest.Solver()
Solver.Solve(H_tot,kpoint)
Quest.Save_object("KPOINT_on",kpoint)
Quest.Save_object("EIGS_on",Solver)
eigs=Quest.Load_object("./EIGS_on")
kpoint=Quest.Load_object("./KPOINT_on")
np.savetxt('kpt_on',np.hstack(kpoint.point))
np.savetxt('eigvals_on',eigs.eigvals)

