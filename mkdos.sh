#!/bin/bash
for i in {-10..10}
do
cd ele_$i
cat > mkdos.py << EOF
import numpy as np
import __tBG_dkl as tBG

eigs = tBG.Load_object("./EIGS")

bound = [-8,8]
energy_resolution = 0.003
epsilon = 0.01

N_energy = int(round((bound[1]-bound[0])/energy_resolution+1))
energies = np.linspace(bound[0],bound[1], num=N_energy)
tBG.Save_object("./E_LDOS",energies)
#exit()
gauss = lambda Ei : np.exp(-((energies-Ei)**2)/(2*epsilon**2))/(epsilon*np.sqrt(2*np.pi))
LDOS = np.zeros((eigs.eigvals[0].size,N_energy))

for eigvals, eigvecs in eigs.eigs:
    for i_e in range(eigvals.size):
        eigval = eigvals[i_e]
        eigvec = eigvecs[:,i_e]
        gauss_eigs = gauss(eigval)
        LDOS += np.outer(np.real_if_close(eigvec.conj()*eigvec),gauss_eigs)

tBG.Save_object("./LDOS",LDOS)
EOF
cd ../
done

