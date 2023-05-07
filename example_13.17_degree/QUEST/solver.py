__all__ = ["Solver", "CuSolver","HubbardSolver", "CuCPGF"]

from .debug import *
from .base import *
from .hamilton import *
import numpy as np
import scipy.sparse as sp
# import scipy.sparse.linalg
from scipy.integrate import cumtrapz
from scipy.special import expit
from scipy.optimize import root_scalar
from scipy.fftpack import dct
from scipy.linalg import eigh
from functools import partial
import math
import time

try:    # Cupy check
    import cupy as cp
    import cupyx.scipy.sparse as csp
    import cupyx.scipy.sparse.linalg
except:
    pass

class Solver:
    def __init__(self, isSparse=False,**kwargs_sparse):
        self.isSparse = isSparse
        self.isCupy = False
        self.solver= partial(sp.linalg.eigsh,k=kwargs_sparse['k'],sigma=kwargs_sparse['sigma'],which='LM') if self.isSparse else eigh
        if(self.isSparse):
            self.N_eigval = kwargs_sparse['k']
        self.eigvals = []
        self.eigvecs = []
    def _solve_matrix(self,H_matrix):
        return self.solver(H_matrix)
    def Solve(self,Hamilton:'Callable Hamilton object',Kpoints:'kspace object', save_eigvec = False, density = None, Sxy = None):
        """
        Solve eigval, eigvec for H(KPOINTS)
        """
        N_orbit = Hamilton.N_orbit*2 if Hamilton.isSpin else Hamilton.N_orbit
        if(self.isSparse):
            N_eigval = self.N_eigval
        else:
            N_eigval = N_orbit
        self.eigvals = np.zeros((Kpoints.N_kpt,N_eigval))
        if(save_eigvec):
            self.eigvecs = np.zeros((Kpoints.N_kpt,N_orbit,N_eigval),dtype=np.complex64)
        for ik,kx,ky,w in Kpoints:
            H = Hamilton(kx, ky, density=density, Sxy=Sxy) if self.isSparse else Hamilton(kx, ky, density=density, Sxy=Sxy).toarray()
            eigval, eigvec = self._solve_matrix(H)
            self.eigvals[ik] = eigval
            if(save_eigvec):
                self.eigvecs[ik] = eigvec
    @timerun
    def Spectral_density(self, Hamilton, Kpoints, energies:np.ndarray, epsilon:float=0.01, Operator=None, density = None, Sxy = None, rescale_bound:"Not used for ED solver"=None):
        """
        Solve Spectral density function for identity matrix and Operator you input.
        """
        N_energy = len(energies)
        LDOS = np.zeros((Hamilton.N_orbit*2 if Hamilton.isSpin else Hamilton.N_orbit,N_energy))
        if(Operator!=None):
            SDOS = np.zeros((Hamilton.N_orbit,N_energy),dtype=np.complex64)
        Emax = -np.inf
        Emin = np.inf
        for ik,kx,ky,w in Kpoints:
            H = Hamilton(kx, ky, density=density, Sxy=Sxy) if self.isSparse else Hamilton(kx, ky, density=density, Sxy=Sxy).toarray()
            eigvals, eigvecs = self._solve_matrix(H)    # Calculate self.eigs
            lorentz = epsilon/(np.subtract.outer(energies,eigvals)**2+epsilon**2)/pi
            LDOS += w*np.real_if_close((lorentz@(eigvecs.conj()*eigvecs).T),1e5).T
            if(Operator!=None):
                SDOS += w*np.real_if_close((lorentz@Operator(eigvecs,eigvecs,sum=False).T),1e5).T
            Emax = max(Emax,eigvals.max())
            Emin = min(Emin,eigvals.min())
        print(f"Energy range in Solver = [{Emin},{Emax}]")
        if(Operator!=None):
            return LDOS, energies, SDOS
        return LDOS, energies
    @timerun
    def Unfolded_spectral_density(self, Hamilton, Kpoints, energies:np.ndarray, epsilon:float=0.01, density = None, Sxy = None):
        """
        Solve Unfolded spectral function.
        
        Input
            Hamilton : Supercell (Moire) system Hamiltonian
            Kpoints : k-path in primitive BZ 
            etc

        Output 
            Spectral function matrix, indexed for [label, kpoint, energy]
        """
        N_energy = len(energies)
        Spectral = np.zeros((len(Hamilton.position.label_unique), Kpoints.N_point ,N_energy))
        print(Hamilton.position.label_unique)
        if(Hamilton.isSpin):
            x = np.vstack([Hamilton.position.point['x'], Hamilton.position.point['x']])
            y = np.vstack([Hamilton.position.point['y'], Hamilton.position.point['y']])
            cond = []
            for label in Hamilton.position.label_unique:
                cond.append(np.append([Hamilton.position.label == label, Hamilton.position.label == label]))
        else:
            x = Hamilton.position.point['x']
            y = Hamilton.position.point['y']
            cond = []
            for label in Hamilton.position.label_unique:
                cond.append(Hamilton.position.label == label)
        for ik,kx,ky,w in Kpoints:
            H = Hamilton(kx, ky, density=density, Sxy=Sxy) if self.isSparse else Hamilton(kx, ky, density=density, Sxy=Sxy).toarray()
            eigvals, eigvecs = self._solve_matrix(H)
            eigvecs = np.exp(-1j*(kx*x+ky*y))*eigvecs
            ## divide sublattice
            lorentz = epsilon/(np.subtract.outer(energies,eigvals)**2+epsilon**2)/pi
            for i_label in range(len(Hamilton.position.label_unique)):
                eigvecs_sub = np.sum(eigvecs[cond[i_label]],axis=0)
                Spectral[i_label,ik] = np.real_if_close((lorentz@(eigvecs_sub.conj()*eigvecs_sub)),1e5)
        return Spectral, energies
        
        
class CuSolver(Solver):
    def __init__(self, Device=0, isSparse=False,**kwargs_sparse):
        self.isSparse = isSparse
        self.isCupy = True
        self.solver= partial(csp.linalg.eigsh,k=kwargs_sparse['k'],sigma=kwargs_sparse['sigma'],which='LM') if self.isSparse else cp.linalg.eigh
        self.eigvals = []
        self.eigvecs = []
        if(self.isSparse):
            self.N_eigval = kwargs_sparse['k']
        self.Device = Device
        self.gpu = cp.cuda.Device(self.Device)
        self.gpu.use()
    def _solve_matrix(self,H_matrix):
        """
        Solve eigval, eigvec for a matrix
        """
        with self.gpu:
            eigval, eigvec = self.solver(csp.csr_matrix(H_matrix) if self.isSparse else cp.array(H_matrix))
            return eigval.get(), eigvec.get()

class HubbardSolver:
    def __init__(self, Solver):
        self.Spectral_density = Solver.Spectral_density
        self.temp = 0.0
        self.doping = 0.0
        self.Mz = 0
        self.Mx = None
        self.My = None
        self.Operator = None
        self.isCupy = Solver.isCupy
    def set_initial_M(self,Mx=None,My=None,Mz=0.0):
        self.Mz = Mz
        if((Mx is not None) or (My is not None)):
            print("If Mx or My is not None, the Fock's hamiltonian is considered")
            self.Mx = Mx
            self.My = My
            self.Operator = partial(Spinxy_Operator,sum=False,isCupy=self.isCupy) 
    def set_temp(self,temp):
        assert temp>=0.0
        self.temp = temp
    def set_doping(self,doping:"Number of electrons per Cell"=0.0):
        self.doping = doping
    def Solve(self, Hamilton:'Callable Hamilton object', Kpoints:'kspace object', energies:np.ndarray, init_density=None, ref_density=None ,maxiter=1000, mixing_beta=0.7, tol=1e-5, rescale_bound=[-10,10], epsilon=0.005):
        """
            bound, epsilon, energy_resolution options are all corresponding with Spectral_density function in Solver
        """
        assert Hamilton.isSpin
        self.__Natom = Hamilton.N_orbit
        ## Initialize   (density, energy-bound)
        if(init_density is not None):
            density = init_density
        else:
            results = self.Spectral_density(Hamilton, Kpoints, Operator=self.Operator, energies=energies, epsilon=epsilon, rescale_bound=rescale_bound, density = density-ref_density, Sxy = Sxy)
            density = self.__density_from_result(results)[0]
        if(ref_density is not None):
            assert density.shape == ref_density.shape
        else:
            ref_density = 0.5
        Mixer = MIXER(density, algorithm="modified")
        if((self.Mx is not None) or (self.My is not None)):
            Sxy = self.Mx + i*self.My
            Mixer_Sxy = MIXER(Sxy, algorithm="linear")
        else:
            Sxy = None
        ## Main iteration
        for i in range(maxiter):
            print(f"{i+1}th iteration")
            results = self.Spectral_density(Hamilton, Kpoints, Operator=self.Operator, energies=energies, epsilon=epsilon, rescale_bound=rescale_bound, density = density-ref_density, Sxy = Sxy)
            Densities = self.__density_from_result(results)
            density, distance = Mixer.mix(Densities[0], mixing_beta=mixing_beta)
            if(Sxy is not None):
                Sxy, _ = Mixer_Sxy.mix(Densities[1], mixing_beta=mixing_beta)
            diff = np.abs(density-Densities[0])
            self.E_total += Hamilton._E_const(density,Sxy)
            print(rf"""==> Iteration {i+1} step ==>
    Distance = {distance:.3e}
    Density of 1st atom = {density[0]+density[self.__Natom]:.8f}
    Magnetization of 1st atom = {-density[0]+density[self.__Natom]:.8f}
    Avg_diff = {np.mean(diff):.3e}
    Max_diff = {diff.max():.3e}
    Ef={self.E_fermi:.8f}
    Etot={self.E_total:.8f}
    """)
            if(distance<tol):
                print(f"Convergence Done in {i+1}th-iteration\n(distance={distance})\n")
                return density, Sxy
        print(f"Not converged (distance={distance})\n")
        return density, Sxy
    def __find_fermilevel(self,TDOS,energies):
        N_electron = self.__Natom+self.doping
        cumdensity = cumtrapz(y=TDOS,x=energies,initial=0.0)
        assert abs((cumdensity[-1] - self.__Natom*2)/(self.__Natom*2)) < 1e-5, f"Check the bound, {cumdensity[-1]}"
        index_max = np.max(np.where(cumdensity<N_electron))
        fermi_0 = (energies[index_max] + energies[index_max+1])/2
        f_cost = lambda mu : np.trapz(expit((mu-energies)/(K_B*self.temp))*TDOS,x=energies) - N_electron
        self.E_fermi = root_scalar(f_cost,bracket=[fermi_0-1.0, fermi_0+1.0], method='brentq').root
        self.E_total= np.trapz(energies*expit((self.E_fermi-energies)/(K_B*self.temp))*TDOS,x=energies)
    def __density_from_result(self,results):  # LDOS(index,energy)
        LDOS = results[0]
        energies = results[1]
        if(len(results)==3):
            SDOS = results[2]
        else:
            SDOS = None
        TDOS = np.sum(LDOS,axis=0)  # TDOS(energy)
        self.__find_fermilevel(TDOS,energies)
        FDdist = expit((self.E_fermi-energies)/(K_B*self.temp))
        density = np.trapz(LDOS*FDdist,energies,axis=1)  # density(index)
        Sdensity = np.zeros_like(density)
        if(SDOS is not None):
            Sdensity = np.trapz(SDOS*FDdist,energies,axis=1)
        return density, Sdensity[:self.__Natom]
                # [density_up,density_dn], Sxy


class CuCPGF:
    r"""
        GPU-accelerated Chebyshev-polynomial Green function(CPGF) method solver
    """
    def __init__(self, Device=0, N_multi=100, vector="local", N_moments=300):
        self.isCupy = True
        self.Device = Device
        self.N_moments = N_moments
        self.gpu = cp.cuda.Device(self.Device)
        self.gpu.use()
        self.Vector_generator = Vector_generator(kind=vector,N_multi=N_multi)
    @timerun
    def Spectral_density(self, Hamilton, Kpoints, energies:np.ndarray, epsilon=0.01, rescale_bound:"Energy range of the Hamiltonian"=[-15,15], Operator=None, density = None, Sxy = None):
        rescale_bound = np.sort(rescale_bound)
        self._a = np.abs(rescale_bound[1]-rescale_bound[0])/2.0
        self._b = (rescale_bound[1]+rescale_bound[0])/2.0
        self.epsilon = epsilon
        self.N_energies = len(energies)
        self.N_dim = Hamilton.N_dim
        self.energies_complex_rescaled = (energies+1j*epsilon-self._b)/self._a
        self.Vector_generator.finalize(N_dim=self.N_dim)
        self.N_vector = self.Vector_generator.N_vector
        r_green = -2/np.sqrt(1.0-self.energies_complex_rescaled**2)*1j*np.exp(-1j*np.arange(self.N_moments)[:,None]*np.arccos(self.energies_complex_rescaled))
        r_green[0] /= 2
        LDOS = np.zeros((self.N_vector,self.N_energies))
        isSDOS=False
        if(Operator!=None):
            isSDOS=True
            SDOS = np.zeros((self.N_vector,self.N_energies))
        for ik,kx,ky,w in Kpoints:
            self._H_GPU = csp.csr_matrix(Hamilton(kx, ky, density=density, Sxy=Sxy)-sp.diags(self._b*np.ones(self.N_dim), dtype=np.complex64, format="csr"))
            self._H_GPU = self._H_GPU/self._a
            moments,moments_op = self.__moments(Operator=Operator)
            moments = np.real_if_close(moments)
            LDOS = LDOS+(-moments@r_green.imag)*w/pi/self._a
            if(isSDOS):
                moments_op = np.real_if_close(moments_op)
                SDOS = SDOS+(moments_op@r_green.imag)*w/pi/self._a
        if(isSDOS):
            return LDOS, energies, SDOS
        return LDOS, energies
    @timerun
    def __moments(self,Operator=None):
        if(Operator==None):
            moments = np.zeros((self.N_vector,self.N_moments),dtype=np.complex64)   # moments matrix in CPU
            for i_vec,f_vec, vecs in self.Vector_generator:
                alpha_zero = cp.array(vecs.T,dtype=cp.complex64)
                moment = cp.zeros((f_vec-i_vec,self.N_moments),dtype=cp.complex64) # moment array in GPU
                alpha = alpha_zero.copy()
                alpha_next = self._H_GPU.dot(alpha)
                moment[:,0] = (alpha.conj()*alpha).sum(axis=0)
                moment[:,1] = (alpha.conj()*alpha_next).sum(axis=0)
                for i_mom in cp.arange(1,self.N_moments//2):
                    alpha, alpha_next = alpha_next, 2*self._H_GPU.dot(alpha_next) - alpha
                    moment[:,2*i_mom] = 2* (alpha.conj()*alpha).sum(axis=0) - moment[:,0]
                    moment[:,2*i_mom+1] = 2* (alpha_next.conj()*alpha).sum(axis=0) - moment[:,1]
                if (self.N_moments %2):
                    moment[:,self.N_moments -1] = 2* (alpha_next.conj()*alpha_next).sum(axis=0) - moment[:,0]
                moments[i_vec:f_vec] = moment.get(stream=cp.cuda.get_current_stream())
            return moments,None
        else:
            moments = np.zeros((self.N_vector,self.N_moments),dtype=np.complex64)
            moments_op = np.zeros((self.N_vector,self.N_moments),dtype=np.complex64)
            for i_vec,f_vec, vecs in self.Vector_generator:
                moment = cp.zeros((f_vec-i_vec,self.N_moments),dtype=cp.complex64)
                moment_op = cp.zeros((f_vec-i_vec,self.N_moments),dtype=cp.complex64)
                alpha_zero = cp.array(vecs.T,dtype=cp.complex64)
                alpha = alpha_zero.copy()
                alpha_next = self._H_GPU.dot(alpha)
                moment[:,0] = (alpha.conj()*alpha).sum(axis=0)
                moment[:,1] = (alpha.conj()*alpha_next).sum(axis=0)
                moment_op[:,0] = Operator(alpha,alpha)
                moment_op[:,1] = Operator(alpha,alpha_next)
                for i_mom in cp.arange(1,self.N_moments//2):
                    alpha, alpha_next = alpha_next, 2*self._H_GPU.dot(alpha_next) - alpha
                    moment[:,2*i_mom] = 2* (alpha.conj()*alpha).sum(axis=0) - moment[:,0]
                    moment[:,2*i_mom+1] = 2* (alpha_next.conj()*alpha).sum(axis=0) - moment[:,1]
                    moment_op[:,i_mom+1] = Operator(alpha_zero, alpha_next)
                if (self.N_moments %2):
                    moment[:,self.N_moments -1] = 2* (alpha_next.conj()*alpha_next).sum(axis=0) - moment[:,0]
                for i_mom in cp.arange(self.N_moments//2+1,self.N_moments):
                    alpha, alpha_next = alpha_next, 2*self._H_GPU.dot(alpha_next) - alpha
                    moment_op[:,i_mom] = Operator(alpha_zero, alpha_next)
                moments[i_vec:f_vec] = moment.get(stream=cp.cuda.get_current_stream())
                moments_op[i_vec:f_vec] = moment_op.get(stream=cp.cuda.get_current_stream())
            return moments, moments_op


class _CuKPM_legacy:
    r"""
        Copyright 2020-2021 Dongkyu Lee* and Kwant authors.

        This class is based on Kwant-KPM. It is subject to the license terms  at
        http://kwant-project.org/license.

        *: risinghermes@gmail.com 
    """
    def __init__(self, Device=0, kernel="Jackson", vector="local"):
        self.SAMPLING = 2
        self.kernel=kernel
        self.isCupy = True # for Check this solver using GPU or not
        self.Device = Device
        self.gpu = cp.cuda.Device(self.Device)
        self.gpu.use()
        self.Vector_generator = Vector_generator(kind=vector)
    @timerun    ## debugging
    def Spectral_density(self, Hamilton, Kpoints, Operator=None,bound=[-1,1],epsilon=0.01,energy_resolution=0.01, density = None, Sxy = None):
        energies,self.N_moments,self._a,self._b = self.__energy(bound=bound,epsilon=epsilon,energy_resolution=energy_resolution)
        self.N_dim = Hamilton.N_dim
        self.Vector_generator.finalize(N_dim=self.N_dim)
        self.N_vector = self.Vector_generator.N_vector
        LDOS = np.zeros((self.N_vector,self.SAMPLING*self.N_moments))
        isSDOS=False
        if(Operator!=None):
            isSDOS=True
            SDOS = np.zeros((self.N_vector,self.SAMPLING*self.N_moments))
        for ik,kx,ky,w in Kpoints:
            self._H_GPU = csp.csr_matrix(Hamilton(kx, ky, density=density, Sxy=Sxy)-sp.diags(self._b*np.ones(self.N_dim), dtype=np.complex64, format="csr"))
            self._H_GPU = self._H_GPU/self._a
            moments,moments_op = self.__momentum(Operator=Operator)
            moments = np.real_if_close(moments)
            moments_tilde = self.__conv_kernel(moments)
            rho = self.__dct_fft(moments_tilde)
            LDOS = LDOS+rho*w
            if(isSDOS):
                moments_op = np.real_if_close(moments_op)
                moments_tilde = self.__conv_kernel(moments_op)
                rho = self.__dct_fft(moments_tilde)
                SDOS = SDOS+rho*w
        if(isSDOS):
            return LDOS, energies, SDOS
        return LDOS,energies
    def __energy(self,bound,epsilon,energy_resolution):
        a = np.abs(bound[1]-bound[0])/(2.0-epsilon)
        b = (bound[1]+bound[0])/2.0
        num_moments = math.ceil(1.6*a/energy_resolution)
        raw,step = np.linspace(pi,0,num_moments*self.SAMPLING,endpoint=False,retstep=True)
        Chebyshev_nodes=np.cos(raw+step/2)
        energies = a*Chebyshev_nodes+b
        return energies,num_moments,a,b
    def __momentum(self,Operator=None):
        if(Operator==None):
            moments = np.zeros((self.N_vector,self.N_moments),dtype=np.complex64)
            moment = cp.zeros(self.N_moments,dtype=np.complex64)
            for i_vec, vec in self.Vector_generator:
                alpha_zero = cp.array(vec,dtype=np.complex64)
                alpha = alpha_zero.copy()
                alpha_next = self._H_GPU.dot(alpha)
                moment[0] = cp.vdot(alpha,alpha)
                moment[1] = cp.vdot(alpha,alpha_next)
                for i_mom in cp.arange(1,self.N_moments//2):
                    alpha, alpha_next = alpha_next, 2*self._H_GPU.dot(alpha_next) - alpha
                    moment[2*i_mom] = 2* cp.vdot(alpha,alpha) - moment[0]
                    moment[2*i_mom+1] = 2* cp.vdot(alpha_next,alpha) - moment[1]
                if (self.N_moments %2):
                    moment[self.N_moments -1] = 2* cp.vdot(alpha_next,alpha_next) - moment[0]
                moments[i_vec] = moment.get()
            return moments,None
        else:
            moments = np.zeros((self.N_vector,self.N_moments),dtype=np.complex64)
            moment = cp.zeros(self.N_moments,dtype=np.complex64)
            moments_op = np.zeros((self.N_vector,self.N_moments),dtype=np.complex64)
            moment_op = cp.zeros(self.N_moments,dtype=np.complex64)
            for i_vec, vec in self.Vector_generator:
                alpha_zero = cp.array(vec,dtype=np.complex64)
                alpha = alpha_zero.copy()
                alpha_next = self._H_GPU.dot(alpha)
                moment[0] = cp.vdot(alpha,alpha)
                moment[1] = cp.vdot(alpha,alpha_next)
                moment_op[0] = Operator(alpha,alpha)
                moment_op[1] = Operator(alpha,alpha_next)
                for i_mom in cp.arange(1,self.N_moments//2):
                    alpha, alpha_next = alpha_next, 2*self._H_GPU.dot(alpha_next) - alpha
                    moment[2*i_mom] = 2* cp.vdot(alpha,alpha) - moment[0]
                    moment[2*i_mom+1] = 2* cp.vdot(alpha_next,alpha) - moment[1]
                    moment_op[i_mom+1] = Operator(alpha_zero, alpha_next)
                if (self.N_moments %2):
                    moment[self.N_moments -1] = 2* cp.vdot(alpha_next,alpha_next) - moment[0]
                for i_mom in cp.arange(self.N_moments//2+1,self.N_moments):
                    alpha, alpha_next = alpha_next, 2*self._H_GPU.dot(alpha_next) - alpha
                    moment_op[i_mom] = Operator(alpha_zero, alpha_next)
                moments[i_vec] = moment.get()
                moments_op[i_vec] = moment_op.get()
            return moments, moments_op
    def __conv_kernel(self,moments):
        m= np.arange(self.N_moments)
        if(self.kernel=="Jackson"):
            kernel_array = ((self.N_moments -m +1)*np.cos(pi*m/(self.N_moments+1)) + np.sin(pi*m/(self.N_moments+1))/np.tan(pi/(self.N_moments+1)))/(self.N_moments+1)
        elif(self.kernel=="Lorentz"):
            l=4
            kernel_array = np.sinh(l*(1-m/self.N_moments))/np.sinh(l)
        else:
            raise RuntimeError("The 'kernel' should be 'Jackson' or 'Lorentz'.")
        return moments * kernel_array /self._a
    def __dct_fft(self,moments_tilde):
        N_sampling = self.N_moments*self.SAMPLING
        moments_ext = np.zeros((self.N_vector,N_sampling), dtype=moments_tilde.dtype)
        moments_ext[:,0:self.N_moments]=moments_tilde
        raw,step = np.linspace(pi,0.0,N_sampling,endpoint=False,retstep=True)
        Chebyshev_nodes=np.cos(raw+step/2)
        gammas = np.flip(dct(moments_ext,type=3),axis=1)
        return gammas/(pi*np.sqrt(1-Chebyshev_nodes**2))
