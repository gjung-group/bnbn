import numpy as np
import matplotlib.pyplot as plt
plt.figure(figsize=(5,5))
W_p=np.loadtxt('eigvals_on')
px_lim=W_p.shape[0]
vbm_p=np.max(np.array(W_p[W_p<0]))
cbm_p=np.min(np.array(W_p[W_p>0]))
gap_p=cbm_p-vbm_p
print("band gap is:",gap_p)

plt.plot(W_p-vbm_p,linestyle="-", linewidth=2,color="grey")
plt.xlim(0,px_lim)
plt.ylim(-4,8)
plt.xticks([])
plt.xlabel("  K          $\Gamma$          M        K", fontsize=20)
plt.ylabel("Energy(eV)", fontsize=20)
plt.show()

