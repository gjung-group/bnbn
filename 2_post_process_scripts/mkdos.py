#!/bin/bash
for i in {-10..10}
do
cp BNBN.cart __tBG_dkl.py ele_$i
cd ele_$i
sed -n '4p' BNBN.cart > a0
a1=`awk '{print $2}' a0`
aa=`echo "scale=8;$a1" | bc -l`
cat > dos.py << EOF
import __tBG_dkl as Quest
import numpy as np
import matplotlib.pyplot as plt 

a = $aa
T = a*np.array([[1,0],[-0.5,np.sqrt(3)/2]])
point = Quest.Position(T)
tB_file = np.genfromtxt("tB")
for x,y,z in tB_file:
    point.add_atom(X=x,Y=y,Z=z,label="tB")
tN_file = np.genfromtxt("tN")
for x,y,z in tN_file:
    point.add_atom(X=x,Y=y,Z=z,label="tN")
bB_file = np.genfromtxt("bB")
for x,y,z in bB_file:
    point.add_atom(X=x,Y=y,Z=z,label="bB")
bN_file = np.genfromtxt("bN")
for x,y,z in bN_file:
    point.add_atom(X=x,Y=y,Z=z,label="bN")
    
gamma_0_aa=-2.700 #first neighbor interaction
gamma_1_aa=0.6924
gamma_0_aap=-0.2102

gamma_0_bb=-2.700##irst neighbor interaction
gamma_1_bb=0.3171
gamma_0_bbp=-0.2102


gamma_0_ab=-2.7 #first neighbor interaction
gamma_1_ab=0.5289
gamma_0_abp=-0.2041

a=1.43
a_1=3.4
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


H0 = Quest.Hamilton(point,True,True)
H0.set_hoppings('tB','tB',np.array([0.007210,0.022764,-0.048173]),arg=3)
H0.set_hoppings('tN','tN',np.array([0.193699,0.019136,-0.037515]),arg=3)
H0.set_hoppings('tB','tN',np.array([-2.720361,-0.211012,0.079806]),arg=3)
H0.set_hoppings('bB','bB',np.array([0.007210,0.022764,-0.048173]),arg=3)
H0.set_hoppings('bN','bN',np.array([0.193699,0.019136,-0.037515]),arg=3)
H0.set_hoppings('bB','bN',np.array([-2.720361,-0.211012,0.079806]),arg=3)

H0.set_onsite(label='tB',onsite=1.753484)
H0.set_onsite(label='tN',onsite=-2.293102)
H0.set_onsite(label='bB',onsite=1.753484)
H0.set_onsite(label='bN',onsite=-2.293102)
H0.set_hoppings('tB','bB',t_aa,isFn=True,arg=6)
H0.set_hoppings('tB','bN',t_ab,isFn=True,arg=6)
H0.set_hoppings('tN','bB',t_ab,isFn=True,arg=6)
H0.set_hoppings('tN','bN',t_bb,isFn=True,arg=6)
kpoint = Quest.Kspace(T)
kpoint.gen_Kpoints([7,7])
kpoint.Weight /= np.sum(kpoint.Weight)
H0.finalize()
Solver = Quest.Solver()
Solver.Solve(H0,kpoint,True)
Quest.Save_object("KPOINT",kpoint)
Quest.Save_object("EIGS",Solver)

eigs = Quest.Load_object("./EIGS")
bound = [-10,10]
energy_resolution = 0.001
epsilon = 0.001
N_energy = int(round((bound[1]-bound[0])/energy_resolution+1))
energies = np.linspace(bound[0],bound[1], num=N_energy)
Quest.Save_object("./E_LDOS",energies)
gauss = lambda Ei : np.exp(-((energies-Ei)**2)/(2*epsilon**2))/(epsilon*np.sqrt(2*np.pi))
LDOS = np.zeros((eigs.eigvals[0].size,N_energy))
for eigvals, eigvecs in eigs.eigs:
    for i_e in range(eigvals.size):
        eigval = eigvals[i_e]
        eigvec = eigvecs[:,i_e]
        gauss_eigs = gauss(eigval)
        LDOS += np.outer(np.real_if_close(eigvec.conj()*eigvec),gauss_eigs)
        Quest.Save_object("./LDOS",LDOS)
EOF
cd ../
done

