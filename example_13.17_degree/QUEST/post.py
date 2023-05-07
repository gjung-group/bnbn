__all__= ["BerryCurv","OrbitalMagnetization", "ChernMarker", "Conductivity_Optical", "Indexing_sparseband"]

#################################################################################
#   QUEST.Post module
#   
#   LIST :
#       BerryCurvature
#       ChernMarker of Resta's formula
#       Optical conductivity
#       Sparse band reindexing
#
#################################################################################
## BerryCurvature

from .solver import Solver
from .base import *
from scipy.special import expit
import numpy as np


def BerryCurv(Hamilton, Kpoints, eigvals=None, eigvecs=None):
    """
    input:
        Hamilton object
        Kpoints object
        (optional) "eigs" corresponding the Hamilton and Kpoints

    return:
        Berry(N_kpoint,N_band)

    ref:
        X. Wang et al., PRB 74, 195118 (2006)
    """
    N_orbit = Hamilton.N_orbit
    delX = Hamilton.delX
    delY = Hamilton.delY
    if(Hamilton.isSpin ==True):
        delX = sp.bmat([[delX,delX],
                        [delX,delX]])
        delY = sp.bmat([[delY,delY],
                        [delY,delY]])
    if(eigvecs is None):
        solver = Solver()
    def eig_call(H, index_kpt):
        if(eigvecs is None):
            eigval, eigvec = solver._solve_matrix(H.toarray())
        else:
            eigval = eigvals[index_kpt]
            eigvec = eigvecs[index_kpt]
        return eigval, eigvec
    Berry = np.zeros((Kpoints.N_kpt,N_orbit))
    for ik,kx,ky,w in Kpoints:
        H = Hamilton(kx,ky)
        dHdkx = -1j*H.multiply(delX)
        dHdky = -1j*H.multiply(delY)
        eigval, eigvec = eig_call(H,ik)
        Vxnm = eigvec.conj().T@(dHdkx.A@eigvec)
        Vynm = eigvec.conj().T@(dHdky.A@eigvec)
        np.fill_diagonal(Vxnm, 0,wrap=True)
        np.fill_diagonal(Vynm, 0,wrap=True)
        Berry[i] = 2*np.imag(np.sum((Vxnm*(Vynm.T))/(((eigval*np.ones((N_orbit,N_orbit))).T-eigval)**2+eps), axis=1))
    return Berry # Not yet checked normailized 

def OrbitalMagnetization(Hamilton, Kpoints, mu=0, temp=1e-5, eigvals=None, eigvecs=None, density=None):
    """
    input:
        Hamilton object
        Kpoints object
        mu : chemical potential
        (optional) "eigs" corresponding the Hamilton and Kpoints

    return:
        Magnetization (float)
    """
    isSolving = False
    if(eigvecs is None):
        isSolving = True
        solver = Solver()
    M = 0.0
    for ik,kx,ky,w in Kpoints:
        if(isSolving):
            eigval, eigvec = solver._solve_matrix(Hamilton(kx,ky,density=density).A)
        else:
            eigval, eigvec = eigvals[ik], eigvecs[ik]
        fermi = expit((mu-eigval)/(K_B*temp))
        dHdkx,dHdky = Hamilton.dHdk(kx,ky)
        Vxnm = eigvec.conj().T@(dHdkx.A@eigvec)
        Vynm = eigvec.conj().T@(dHdky.A@eigvec)
        np.fill_diagonal(Vxnm, 0,wrap=True)
        np.fill_diagonal(Vynm, 0,wrap=True)
        # sum m -> product fermi -> sum n 
        M += w*np.sum(np.imag(np.sum((Vxnm*(Vynm.T))*(np.add.outer(eigval,eigval)-2*mu)/(np.subtract.outer(eigval,eigval)**2+eps), axis=1))*fermi)
    return M # Not yet checked normailized 


## Local Chern marker
def ChernMarker(eigvals, eigvecs, position, Kpoints, mu):
    """
    input :
        "eigs"
        rspace object 
        kspace object
        Target Mu
    
    return:
        RChern(index_site)

    ref:
        Marsal, Quentin, Dániel Varjas, and Adolfo G. Grushin. "Topological Weaire-Thorpe models of amorphous matter." Proceedings of the National Academy of Sciences 117.48 (2020): 30260-30265.
    """
    N_atom = position.size
    N_site = eigvecs[0].shape[1]
    RChern = np.zeros(N_site)
    position = np.append(position,position)
    weight_sum = Kpoints.weight_sum()
    for ik,kx,ky,w in Kpoints:
        eigval = eigvals[ik]
        eigvec = eigvecs[ik]
        P = np.diag(eigval < mu)
        Q = np.eye(N_site) - P
        x = eigvec.T.conj() @ np.diag(position['x']) @ eigvec
        y = eigvec.T.conj() @ np.diag(position['y']) @ eigvec
        PxP, PyP = P @ x @ P, P @ y @ P
        C = -2 * pi * 1j * (PxP @ PyP - PyP @ PxP)
        C = eigvec @ C @ eigvec.T.conj()

        RChern += w*np.diag(C).real/weight_sum
    areas = _voronoi_volumes(np.array(position[['x','y']].tolist()))
    return RChern[0:N_atom]/areas[0:N_atom], RChern[N_atom:]/areas[N_atom:]


def _voronoi_volumes(position):
    from scipy.spatial import Voronoi, ConvexHull
    v = Voronoi(position)
    vol = np.zeros(v.npoints)
    for i, reg_num in enumerate(v.point_region):
        indices = v.regions[reg_num]
        if -1 in indices: # some regions can be opened
            vol[i] = -1
        else:
            vol[i] = ConvexHull(v.vertices[indices]).volume
    vol[vol==-1] = np.average(vol[vol!=-1])
    return vol


def Conductivity_Optical_older(Hamilton, Kpoints, mu, temp, omega , eigvals=None, eigvecs=None):
    """
    input:
        Hamilton object
        Kpoints object
        Mu, Temp : Chemical potential and temperature of fermi-dirac function
        omega : list of target frequency-energies [eV] 
        (optional) "eigs" corresponding the Hamilton and Kpoints
    return:
        Conductivitiy tensor([xx,xy,yx,yy],N_omega)
    """
    assert Hamilton.isSmall == False, "Optical conductivity for small hamiltonian is not implemented"
    N_orbit = Hamilton.N_orbit
    delX = Hamilton.delX
    delY = Hamilton.delY
    gamma = 0.1
    if(Hamilton.isSpin ==True):
        delX = sp.bmat([[delX,delX],
                        [delX,delX]])
        delY = sp.bmat([[delY,delY],
                        [delY,delY]])
    if(eigvecs is None):
        solver = Solver()
    def eig_call(H, index_kpt):
        if(eigvecs is None):
            eigval,eigvec = solver._solve_matrix(H.toarray())
        else:
            eigval = eigvals[index_kpt]
            eigvec = eigvecs[index_kpt]
        return eigval, eigvec
    Conductivity = np.zeros((4,len(omega)))
    for ik,kx,ky,w in Kpoints:
        H = Hamilton(kx,ky)
        dHdkx = -1j*H.multiply(delX)
        dHdky = -1j*H.multiply(delY)
        eigval, eigvec = eig_call(H,ik)
        cond = np.logical_or(eigval<=mu+1.2*max(omega), eigval>=mu-1.2*max(omega))
        if(np.sum(cond)==0):
            continue
        eigval = eigval[cond]
        eigvec = eigvec[:,cond]
        fermi = expit((mu-eigval)/(K_B*temp))
        delE = (eigval - (np.ones((eigval.size,eigval.size))*eigval).T).flatten()
        delFermi = (fermi - (np.ones((eigval.size,eigval.size))*fermi).T).flatten()
        Vxnm = eigvec.conj().T@(dHdkx*eigvec)
        Vynm = eigvec.conj().T@(dHdky*eigvec)
        sigma_xx = delFermi*(Vxnm*Vxnm.T).flatten()
        sigma_xy = delFermi*(Vxnm*Vynm.T).flatten()
        sigma_yx = delFermi*(Vynm*Vxnm.T).flatten()
        sigma_yy = delFermi*(Vynm*Vynm.T).flatten()
        # Lorentzian Convolution 
        for j in range(len(omega)):
            cond_0 = gamma/(2*pi*omega[j]*((omega[j]-delE)**2+(gamma/2)**2))
            Conductivity[0,j] += -w*np.sum(sigma_xx*cond_0)
            Conductivity[1,j] += -w*np.sum(sigma_xy*cond_0)
            Conductivity[2,j] += -w*np.sum(sigma_yx*cond_0)
            Conductivity[3,j] += -w*np.sum(sigma_yy*cond_0)
    return Conductivity / Kpoints.Sum_weight


def Conductivity_Optical(Hamilton, Kpoints, mu, temp, omega , eigvals=None, eigvecs=None):
    gamma = 0.055
    isSolving = False
    if(eigvecs is None):
        isSolving = True
        solver = Solver()
    Conductivity = np.zeros((4,len(omega)),dtype=complex)
    for ik,kx,ky,w in Kpoints:
        if(isSolving):
            eigval, eigvec = solver._solve_matrix(Hamilton(kx,ky).A)
        else:
            eigval, eigvec = eigvals[ik], eigvecs[ik]
        cond = np.logical_or(eigval<=mu+1.2*max(omega), eigval>=mu-1.2*max(omega))
        if(np.sum(cond)==0):
            continue
        eigval = eigval[cond]
        eigvec = eigvec[:,cond]
        fermi = expit((mu-eigval)/(K_B*temp))
        delE = (eigval - (np.ones((eigval.size,eigval.size))*eigval).T).flatten()
        delE[delE==0]=np.inf
        delFermi = (fermi - (np.ones((eigval.size,eigval.size))*fermi).T).flatten()
        dHdkx,dHdky = Hamilton.dHdk(kx,ky)
        Vxnm = eigvec.conj().T@(dHdkx.A@eigvec)
        Vynm = eigvec.conj().T@(dHdky.A@eigvec)
        sigma_xx = delFermi*(Vxnm*Vxnm.T).flatten()
        sigma_xy = delFermi*(Vxnm*Vynm.T).flatten().imag
        sigma_yx = delFermi*(Vynm*Vxnm.T).flatten().imag
        for j in range(len(omega)):
            efactorxx = 1./(delE)/(omega[j]+delE+1.j*gamma)
            efactor = 1./((delE)**2-(omega[j]+1.j*gamma)**2)
            Conductivity[0,j] += w*np.sum(sigma_xx*efactorxx)
            Conductivity[1,j] += w*np.sum(sigma_xy*efactor)
            Conductivity[2,j] += w*np.sum(sigma_yx*efactor)
    return Conductivity / Kpoints.Sum_weight


def Indexing_sparseband(eigvals,plotrange=[-2,2]):
    N_kpt = eigvals.shape[0]
    plt.ylim(plotrange)
    plt.xlim([0,N_kpt-1])
    plt.plot(np.arange(N_kpt),eigvals,"k.",markersize=0.5)
    print(">> Right click : add point \n>> Left click : remove point \n>> Middle click : END")
    clickpoints = np.array(plt.ginput(-1,timeout=0))
    E_ref = interp1d(clickpoints[:,0],clickpoints[:,1],fill_value="extrapolate")(np.arange(N_kpt))

    idx_ref = (1/(eigvals.T-E_ref)).argmax(axis=0)
    idx_sort = (eigvals.T.argsort(axis=0)-idx_ref).T
    min_ofMax = np.min(np.max(idx_sort,axis=1))
    max_ofMin = np.max(np.min(idx_sort,axis=1))
    condition = np.logical_and(idx_sort<=min_ofMax, idx_sort>=max_ofMin)
    print(f"The bands are merged to {min_ofMax-max_ofMin+1} bands.")
    compressed_eigvals = (eigvals[condition]).reshape((N_kpt,min_ofMax-max_ofMin+1))
    plt.close()
    return compressed_eigvals
