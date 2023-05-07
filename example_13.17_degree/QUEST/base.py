
## Object I/O method
import pickle
def Save_object(filename,obj):
    pickle.dump(obj, file = open(filename, "wb"), protocol=-1)

def Load_object(filename):
    return pickle.load(open(filename, "rb"))

## constants
eps=1e-10

from scipy.constants import value
import scipy.sparse as sp
K_B = value(u'Boltzmann constant in eV/K')
import numpy as np
from numpy import pi
pi = pi


## Hopping equations
def SK_hopping(x,y,z):
    a = 2.46
    r0 = 0.184*a
    a0=a/np.sqrt(3)
    d0=3.35
    Vpppi = -2.7
    Vppsig = 0.48
    d_square = x**2 + y**2 + z**2
    R_square = x**2 + y**2
    t_eq = Vpppi*np.exp(-(np.sqrt(d_square)-a0)/r0)*R_square/d_square \
        + Vppsig*np.exp(-(np.sqrt(d_square)-d0)/r0)*z**2/d_square
    return t_eq


def Spinxy_Operator(alpha,beta,sum=True,isCupy=False):
    # return a constant value of <alpha|O|beta>
    # for O = S+ = [0 1]
    #              [0 0]
    N_half = alpha.size/2
    N_half_test = beta.size/2
    assert N_half==N_half_test,"Two array must be same size"
    if(cupy):
        vdot = cp.vdot
        hstack = cp.hstack
    else:
        vdot = np.vdot
        hstack = np.hstack
    if(sum):
        return vdot(alpha[0:N_half],beta[N_half:])
    else:
        return alpha[0:N_half].conj()*beta[N_half:]
            # [CupCdn,CdnCup]

def Rotation_matrix(theta):
    return np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])


class Vector_generator:
    def __init__(self, kind="local",N_multi=1, isSparse=False):
        self.kind = kind
        self.isSparse = isSparse
        self.isFinalized = False
        self.N_multi=N_multi
        self.N_vector = 0
    def __iter__(self): ## make iterator for "for-loop"
        """
            return index, kx, ky, weight
        """
        assert self.isFinalized, "Finalize first before you use this class for iterator."
        self.__index = 0
        self.__final = self.N_vector
        return self
    def __next__(self): ## make iterator for "for-loop"
        """
            return index, kx, ky, weight
        """
        if self.__index < self.__final:
            if self.__index+self.N_multi<self.__final:
                N_vector = self.N_multi
            else:
                N_vector = self.__final-self.__index
            vectors = self.Vectors[self.__index:self.__index+N_vector]
            if(self.isSparse==False):
                vectors = vectors.A
            self.__index = self.__index+N_vector
            return self.__index-N_vector, self.__index, vectors 
        else:
            raise StopIteration
    def finalize(self, N_dim, locations :"array of indicies of location" = None):
        if(self.kind == "local"):
            if(locations == None):
                self.N_vector = N_dim
                self.Vectors = sp.eye(m=self.N_vector,dtype=np.complex64,format="csr")
            else:
                self.N_vector = len(locations)
                self.Vectors = sp.csr_matrix((np.ones(self.N_vector),(np.arange(self.N_vector),locations)),shape=(self.N_vector,N_dim))
        elif(self.kind == "random"):
            raise NotImplementedError()
        else:
            raise RuntimeError("The 'kind' must be 'local' or 'random'.")
        self.isFinalized = True


class MIXER:
    r"""
        Vector mixing algorithm
        How to use:
            mixer = MIXER(n_init, algoritm ="modified")
            n_next = mixer.mix(n_out)

        Ref: 
            Johnson, Duane D. "Modified Broyden's method for accelerating convergence in self-consistent calculations." Physical Review B 38.18 (1988): 12807.
    """
    def __init__(self, n_init, algorithm="linear"):
        print(f"= The {algorithm} method mixing =")
        if(algorithm=="linear"):
            self.__init_linear__(n_init)
            self.mix = self.__linear_mix
        elif(algorithm=="anderson"):
            self.__init_anderson__(n_init)
            self.mix = self.__anderson_mix
        elif(algorithm=="broyden"):
            assert NotImplementedError()
            self.__init_broyden__(n_init)
            self.mix = self.__broyden_mix
        elif(algorithm=="modified"):
            self.__init_modified__(n_init)
            self.mix = self.__modified_mix
        else:
            assert ValueError()
    def __init_linear__(self, n_init):
        self.n_old = n_init
    def __linear_mix(self, n_out, mixing_beta=0.3):
        distance = self.Distance(self.n_old, n_out)
        self.n_old = (1-mixing_beta)*self.n_old + mixing_beta*n_out 
        return self.n_old, distance
    def __init_anderson__(self, n_init):
        self.n_in_old2 = n_init.copy()
        self.n_in_old = n_init.copy()
        self.n_out_old = n_init.copy() + 1e-5*np.random.randn(*(n_init.shape))
        self.F_old = self.n_out_old - self.n_in_old2
    def __anderson_mix(self, n_out, mixing_beta=0.5):
        distance = self.Distance(self.n_in_old, n_out)
        F = n_out - self.n_in_old
        mixing_alpha = np.dot(F,F-self.F_old)/self.Distance(F,self.F_old)**2
        print(mixing_alpha)
        n_new = (1-mixing_beta)*((1-mixing_alpha)*self.n_in_old + mixing_alpha*self.n_in_old2) + mixing_beta*((1-mixing_alpha)*n_out + mixing_alpha*self.n_out_old) 
        self.n_in_old, self.n_in_old2 = n_new, self.n_in_old
        self.n_out_old = n_out
        self.F_old = F
        return n_new, distance
    def __init_modified__(self, n_init, mixing_beta=0.5):
        # Assume G1 = beta*I
        self.n_old = n_init
        self.F_old = np.zeros_like(n_init)
        self.w = np.array([0.0])
        self.dF = np.empty((0,n_init.size))
        self.dN = np.zeros_like(n_init)
    def __modified_mix(self, n_out, mixing_beta=0.5):
        distance = self.Distance(self.n_old, n_out)
        F = n_out - self.n_old
        delF = (F - self.F_old)/np.linalg.norm(F-self.F_old)
        self.dF = np.vstack([self.dF, delF])
        self.dN[-1] /= np.linalg.norm(F-self.F_old)
        gamma = np.sum(self.w*np.dot(self.dF,F)*np.linalg.inv(1e-4*np.eye(self.dF.shape[0])+np.outer(self.w,self.w)*np.inner(self.dF,self.dF)),axis=1)
        n_new = self.n_old + mixing_beta*F - (self.w*gamma)@(mixing_beta*self.dF+self.dN)
        delN = (n_new - self.n_old)
        self.dN = np.vstack([self.dN, delN])
        if(np.sum(np.isnan(n_new))>0):
            print("Warning : NaN out")
            return self.n_old, 0.0
        self.n_old = n_new
        self.F_old = F
        self.w = np.append(self.w, 1.0)
        return n_new, distance
    def Distance(self,a,b):
        return np.sqrt(np.dot(a-b,a-b))