__all__ = ['H_tot', 'H0', 'H_Haldane', 'H_Rashba', 'H_U', 'H_V']

#################################################################################
#   QUEST.Hamilton module
#   
#   How to use :
#       H0 = QUEST.H0()
#       H0.set_hoppings(~~)
#       H_tot += H0
#       H_haldane = QUEST.H_Haldane()
#       H_haldane.set_haldane(~~)
#       H_tot += H_haldane
#       Then, if you call H_tot(kx,ky,density), it will return scipy.spmatrix of Hamiltonian
#
#################################################################################


from .debug import *
from .base import *
import numpy as np
import scipy.sparse as sp
from copy import deepcopy
from itertools import combinations
import warnings
#warnings.filterwarnings('error')
warnings.filterwarnings(action='ignore')

class H_tot:
    """
        The Bowl of all hamiltonians.

        For Developers,
            This class includes below attributes.
            1. __position_input    : Input orbits' position
            2. position            : Expanded position (if isSmall=True case)
            3. delX, delY, delZ, R : Differences of 'position', and sqrt(dX^2 + dY^2)
            4. NN                  : Order of neighbors, can be used as 
                col, row = labels_to_rowcol(label1, label2)
                orders = NN[(col.row)]
            5. nonzero             : Boolean sparse matrix. If an element position in any attributes has value, the same element in the nonzero matrix has True.  
    """
    def __init__(self, position:"rspace object", period_1=True, period_2=True, cutoff=10.0):
        self.Hamils = []
        self.isSpin = False
        self.isSmall = False
        self.maxorders = []
        self.__position_input = position
        self.__position_input._finalize()
        self.position = 0
        self.period_1 = period_1
        self.period_2 = period_2
        self.cutoff = cutoff
        self.__initialize()
    @property
    def N_orbit(self):
        return self.__position_input.N_orbit
    @property
    def N_dim(self):
        """
             Dimension of a matrix when you call this hamiltonian
        """
        return self.position.N_orbit*2 if self.isSpin else self.position.N_orbit
    def __iadd__(self,other):
        assert callable(other)
        other._finalize(self)   # Finalizations of all Hamiltonian is done in this line
        if(hasattr(other, "isSpin")):
            self.isSpin = self.isSpin or other.isSpin
        self.Hamils.append(other)
        return self
    def __call__(self, kx=0, ky=0, density=None, Sxy=None):
        eikdr = -1j*(kx*self.delX+ky*self.delY).tocsr()
        eikdr.data = np.exp(eikdr.data)
        eikdr = eikdr.todok()
        eikdr += (self.nonzero>(eikdr !=0)).astype(np.complex64)
        H_tot = sp.dok_matrix((self.N_dim,self.N_dim), dtype=np.complex64)
        for H in self.Hamils:
            H_tot += H(eikdr=eikdr, isSpin=self.isSpin, density=density, Sxy=Sxy,kx=kx,ky=ky)
        if(self.isSmall):
            assert (self.N_dim/(self.N_orbit*self.multiple)==1) or (self.N_dim/(self.N_orbit*self.multiple)==2)
            H_small = sp.dok_matrix((self.N_orbit, self.N_orbit), dtype=np.complex64)
            for i in range(self.multiple):
                H_small += H_tot[0:self.N_orbit,i*self.N_orbit:(i+1)*self.N_orbit]
            if(self.isSpin):
                raise NotImplementedError()
                H_small_dn = sp.dok_matrix((self.N_orbit, self.N_orbit), dtype=np.complex64)
                H_small_updn = sp.dok_matrix((self.N_orbit, self.N_orbit), dtype=np.complex64)
                for i in range(self.multiple):
                    H_small_dn += H_tot[(self.multiple)*self.N_orbit:(self.multiple+1)*self.N_orbit,(self.multiple+i)*self.N_orbit:(self.multiple+i+1)*self.N_orbit]
                    H_small_updn += H_tot[0:self.N_orbit,(self.multiple+i)*self.N_orbit:(self.multiple+i+1)*self.N_orbit]
                return sp.bmat([[H_small,H_small_updn],[H_small_updn.conj().T,H_small_dn]])
            return H_small
        else:
            return H_tot
    def _E_const(self,density,Sxy=None):
        """
            The constant part in a total energy of each hamiltonian
        """
        E_const = 0.0
        for H in self.Hamils:
            E_const += H._E_const(isSpin = self.isSpin, density = density, Sxy = Sxy)
        return E_const
    def dHdk(self, kx=0, ky=0, density=None, Sxy=None):
        """
            return dH/dkx, dH/dky
        """
        eikdr = -1j*(kx*self.delX+ky*self.delY).tocsr()
        eikdr.data = np.exp(eikdr.data)
        dHtot_dK = sp.dok_matrix((self.N_dim,self.N_dim), dtype=np.complex64)
        for H in self.Hamils:
            dHtot_dK += H(eikdr, isSpin = self.isSpin, density = density, Sxy = Sxy)
        if(self.isSpin):
            dHtot_dkx = dHtot_dK.multiply(-1j*sp.bmat([[self.delX,self.delX],[self.delX,self.delX]]))
            dHtot_dky = dHtot_dK.multiply(-1j*sp.bmat([[self.delY,self.delY],[self.delY,self.delY]]))
        else:
            dHtot_dkx = dHtot_dK.multiply(-1j*self.delX)
            dHtot_dky = dHtot_dK.multiply(-1j*self.delY)
        if(self.isSmall):
            assert (self.N_dim/(self.N_orbit*self.multiple)==1) or (self.N_dim/(self.N_orbit*self.multiple)==2)
            dHs_dkx = sp.dok_matrix((self.N_orbit, self.N_orbit), dtype=np.complex64)
            dHs_dky = sp.dok_matrix((self.N_orbit, self.N_orbit), dtype=np.complex64)
            for i in range(self.multiple):
                dHs_dkx += dHtot_dkx[0:self.N_orbit,i*self.N_orbit:(i+1)*self.N_orbit]
                dHs_dky += dHtot_dky[0:self.N_orbit,i*self.N_orbit:(i+1)*self.N_orbit]
            if(self.isSpin):
                raise NotImplementedError()
                dHs_dkx_dn = sp.dok_matrix((self.N_orbit, self.N_orbit), dtype=np.complex64)
                dHs_dkx_updn = sp.dok_matrix((self.N_orbit, self.N_orbit), dtype=np.complex64)
                dHs_dky_dn = sp.dok_matrix((self.N_orbit, self.N_orbit), dtype=np.complex64)
                dHs_dky_updn = sp.dok_matrix((self.N_orbit, self.N_orbit), dtype=np.complex64)
                for i in ranage(self.multiple):
                    dHs_dkx_dn += dHtot_dkx[(self.multiple)*self.N_orbit:(self.multiple+1)*self.N_orbit,(self.multiple+i)*self.N_orbit:(self.multiple+i+1)*self.N_orbit]
                    dHs_dkx_updn += dHtot_dkx[0:self.N_orbit,(self.multiple+i)*self.N_orbit:(self.multiple+i+1)*self.N_orbit]
                    dHs_dky_dn += dHtot_dky[(self.multiple)*self.N_orbit:(self.multiple+1)*self.N_orbit,(self.multiple+i)*self.N_orbit:(self.multiple+i+1)*self.N_orbit]
                    dHs_dky_updn += dHtot_dky[0:self.N_orbit,(self.multiple+i)*self.N_orbit:(self.multiple+i+1)*self.N_orbit]
                return sp.bmat([[dHs_dkx,dHs_dkx_updn],[dHs_dkx_updn.conj().T,dHs_dkx_dn]]), sp.bmat([[dHs_dky,dHs_dky_updn],[dHs_dky_updn.conj().T,dHs_dky_dn]]) 
            return dHs_dkx, dHs_dky
        else:
            return dHtot_dkx, dHtot_dky
    def __initialize(self):
        N_T1 = 0
        N_T2 = 0
        if(self.period_1 == True):
            N_T1 = int(np.rint(self.cutoff/np.linalg.norm(self.__position_input.Tvector[0])))
        if(self.period_2 == True):
            N_T2 = int(np.rint(self.cutoff/np.linalg.norm(self.__position_input.Tvector[1])))
        if(N_T1!=0 or N_T2!=0):
            print("Small-cell Hamiltonian")
            self.isSmall = True
            self.multiple = (2*N_T1+1)*(2*N_T2+1)
            self.__exp_position(N_T1,N_T2)
        else:
            self.position = self.__position_input
        self.delX = sp.csr_matrix((self.N_dim,self.N_dim), dtype=float)
        self.delY = sp.csr_matrix((self.N_dim,self.N_dim), dtype=float)
        self.delZ = sp.csr_matrix((self.N_dim,self.N_dim), dtype=float)
        self.nonzero = sp.csr_matrix((self.N_dim,self.N_dim), dtype=float)
        self.R = sp.csr_matrix((self.N_dim,self.N_dim), dtype=float)
        self.NN = sp.csr_matrix((self.N_dim,self.N_dim), dtype=int)  # 1: First-NN 2: Second-NN ...
        self.__cal_delXYZ()
    def __exp_position(self,N_T1,N_T2):
        assert self.isSmall
        self.position = deepcopy(self.__position_input)
        self.position.make_supercell([-N_T1,N_T1],[-N_T2,N_T2])
    def __cal_delXYZ(self):
        Tinv = np.linalg.inv(self.position.Tvector).T
        dx = ((self.position.point['x']*np.ones((self.N_dim,self.N_dim))).T-self.position.point['x']).T
        dy = ((self.position.point['y']*np.ones((self.N_dim,self.N_dim))).T-self.position.point['y']).T
        dz = ((self.position.point['z']*np.ones((self.N_dim,self.N_dim))).T-self.position.point['z']).T
        if(self.period_1):
            U = dx*Tinv[0,0] + dy*Tinv[0,1]
            dx -= np.rint(U-eps)*self.position.Tvector[0,0] 
            dy -= np.rint(U-eps)*self.position.Tvector[0,1]
        if(self.period_2):
            V = dx*Tinv[1,0] + dy*Tinv[1,1]
            dx -= np.rint(V-eps)*self.position.Tvector[1,0] 
            dy -= np.rint(V-eps)*self.position.Tvector[1,1]
        cond_cutoff = (dx**2+dy**2) <= self.cutoff**2
        self.delX = sp.csr_matrix(cond_cutoff*dx)
        self.delY = sp.csr_matrix(cond_cutoff*dy)
        self.delZ = sp.csr_matrix(cond_cutoff*dz)
        self.__cal_RNN()
    def __cal_RNN(self):    # Calculate planer-distance R, nonzero indexes and NN-orders
        self.R += self.delX.power(2) + self.delY.power(2)
        D = self.R +self.delZ.power(2)
        self.nonzero = ((D!=0).astype(float)).todok()
        self.nonzero.setdiag(np.ones(self.position.N_orbit))
        self.R = self.R.sqrt()
        ##  Different label case
        for label_from, label_to in combinations(self.position.label_unique, 2):
            row, col = self.labels_to_rowcol(label_from, label_to)
            part_R = D[(row,col)].tocsr()
            part_R.data /= np.mean(part_R.data[part_R.data <(part_R.data.min()*1.2)])
            _,order_unique = np.unique(np.round(part_R.data,decimals=0),return_inverse=True)
            order_unique += 1
            _NN_label = sp.csr_matrix((order_unique,part_R.indices,part_R.indptr),shape=part_R.shape)
            self.NN[(row,col)] = _NN_label
            self.NN[(col,row)] = _NN_label
            self.maxorders.append(order_unique.max())
        ##  Same label case
        for label_from in self.position.label_unique:
            row, col = self.labels_to_rowcol(label_from, label_from)
            part_R = D[(row,col)].tocsr()
            part_R.data /= np.mean(part_R.data[part_R.data <(part_R.data.min()*1.2)])
            _,order_unique = np.unique(np.round(part_R.data,decimals=0),return_inverse=True)
            order_unique += 1
            _NN_label = sp.csr_matrix((order_unique,part_R.indices,part_R.indptr),shape=part_R.shape)
            self.NN[(row,col)] = _NN_label
            self.maxorders.append(order_unique.max())
    def reset_cutoff(self, cutoff):
        self.cutoff = cutoff
        self.__initialize()
    def labels_to_rowcol(self,label_from,label_to):
        cond_from = self.position.label == label_from
        cond_to = self.position.label == label_to
        col, row = np.meshgrid(np.where(cond_from),np.where(cond_to))
        return row,col



class H0:
    """
        Non-interacting tight-binding Hamiltonian 
    """
    def __init__(self):
        self.isSpin = False
        self.hopping = np.array([],dtype=[('From','<U4'),('To','<U4'),('isFn',bool),('t',object),('cutoff',float)])
        self.onsite = 0
        self.isFinalized = False
        self.H0 = np.empty(0)
    def __call__(self, eikdr, isSpin=None, **kwargs_dummy):
        assert self.isFinalized
        if(isSpin==None):
            isSpin = self.isSpin
        H0_k = self.H0.multiply(eikdr) + sp.diags(self.onsite)
        if(isSpin):
            return sp.bmat([[H0_k,None],[None,H0_k]])
        else:
            return H0_k
    def _E_const(self,**kwargs_dummy):
        E_const = 0.0
        return E_const
    def set_hoppings(self, label_from, label_to, t_hoppings, isFn:bool=False, cutoff=None):
        """
            if "isFn" is True, t_hoppings should be given as "callable function" likes t(delta_x, delta_y, delta_z).
            if "isFn" is False, t_hoppings should be numpy-array for hopping parameters.
            The 't_hoppings' array does not contain on-site energy. 
            Onsite-energy are handled in set_onsite

            % cutoff = 'R_lim' is limit length of hopping for hopping equation.
                  = 'N_order' is maximum hopping order for hopping array.
        """
        condition_same = (self.hopping['From']==label_from)*(self.hopping['To']==label_to)
        if(np.sum(condition_same)!=0):
            print(f"The formal hopping parameter from {label_from} to {label_to} is removed")
            condition_same = np.logical_or(condition_same,(self.hopping['From']==label_to)*(self.hopping['To']==label_from))
            self.hopping = np.compress(np.logical_not(condition_same),self.hopping)
        if(isFn==False):
            if(type(t_hoppings) is not np.ndarray):
                raise RuntimeError("t_hoppings should be numpy.ndarray when isFn==False")
            if(cutoff==None):
                cutoff = t_hoppings.size
        else:
            if(callable(t_hoppings)==False):
                raise RuntimeError("t_hoppings should be callable function when isFn==True")
            assert cutoff !=None, "Cutoff should be float-number when the hopping is a equation"
        self.hopping = np.append(self.hopping,np.array((label_from,label_to,isFn,t_hoppings,cutoff),dtype=self.hopping.dtype))
        if(label_from!=label_to):
            self.hopping = np.append(self.hopping,np.array((label_to,label_from,isFn,t_hoppings,cutoff),dtype=self.hopping.dtype))
    def set_onsite(self, onsite):
        """
            onsite can be a number or ndarray(size of N_orbit)
        """
        self.onsite = onsite
    def _finalize(self, total):
        self.N_orbit = total.position.N_orbit
        self.__finalize_onsite()
        self.__finalize_H0(total)
        self.isFinalized = True
    def __finalize_onsite(self):
        if((type(self.onsite) is float) or (type(self.onsite) is int)):
            self.onsite = np.ones(self.N_orbit) * self.onsite
        elif(self.onsite.size == self.N_orbit):
            pass
        elif(self.N_orbit%self.onsite.size==0): # Small-cell case
            self.onsite = np.tile(self.onsite,reps=int(self.N_orbit/self.onsite.size))
        else:
            raise RuntimeError(f"Check the dimension of the onsite array : {self.onsite.size}")
    def __finalize_H0(self, total):
        # Cutoff test 
        i=0
        ##  Different label case
        for label_from, label_to in combinations(total.position.label_unique, 2):
            hopping = self.hopping[np.logical_and(self.hopping['From']==label_from,self.hopping['To']==label_to)]
            if(hopping.size==0):
                continue
            if(hopping['isFn']):
                assert total.cutoff >= hopping['cutoff'], "Your cutoff is too smaller than hopping range, (H0)"
            else:
                assert total.maxorders[i] >= hopping['cutoff'], "Your cutoff is too smaller than hopping order, (H0)"
            i += 1
        ##  Same label case
        for label_from in total.position.label_unique:
            hopping = self.hopping[np.logical_and(self.hopping['From']==label_from,self.hopping['To']==label_from)]
            if(hopping.size==0):
                continue
            if(hopping['isFn']):
                assert total.cutoff >= hopping['cutoff'], "Your cutoff is too smaller than hopping range, (H0)"
            else:
                assert total.maxorders[i] >= hopping['cutoff'], "Your cutoff is too smaller than hopping order, (H0)"
            i += 1
        # Finalization start
        self.H0 = sp.dok_matrix((self.N_orbit,self.N_orbit),dtype=np.complex64)
        ##  Different label case
        for label_from, label_to in combinations(total.position.label_unique, 2):
            hopping = self.hopping[np.logical_and(self.hopping['From']==label_from,self.hopping['To']==label_to)]
            if(hopping.size!=0):
                hopping = hopping[0]
                row, col = total.labels_to_rowcol(label_from, label_to)
                if(hopping['isFn']):    # Hopping equation case
                    nonzero = total.nonzero[(row,col)]
                    R_cut = total.R[(row,col)] <= hopping['cutoff']
                    cond_mat = nonzero.multiply(R_cut).astype('bool')  #logical_and
                    cond = cond_mat.nonzero()
                    delX = total.delX[(row,col)][cond].A1
                    delY = total.delY[(row,col)][cond].A1
                    delZ = total.delZ[(row,col)][cond].A1
                    _H_label = sp.csr_matrix((hopping['t'](delX,delY,delZ),cond_mat.indices,cond_mat.indptr),shape=cond_mat.shape)
                else:                   # Hopping parameters case
                    NN = total.NN[(row,col)].multiply((total.NN[(row,col)]<=hopping['cutoff']).astype('int'))
                    t = hopping['t'][NN.data-1]
                    _H_label = sp.csr_matrix((t,NN.indices,NN.indptr),shape=NN.shape)
                self.H0[(row,col)] = _H_label
                self.H0[(col,row)] = _H_label.conj()
        ##  Same label case
        for label_from in total.position.label_unique:
            hopping = self.hopping[np.logical_and(self.hopping['From']==label_from,self.hopping['To']==label_from)]
            if(hopping.size!=0):
                hopping = hopping[0]
                row, col = total.labels_to_rowcol(label_from, label_from)
                if(hopping['isFn']):    # Hopping equation case
                    nonzero = total.nonzero[(row,col)]
                    R_cut = total.R[(row,col)] <= hopping['cutoff']
                    cond_mat = nonzero.multiply(R_cut).astype('bool')  #logical_and
                    cond = cond_mat.nonzero()
                    delX = total.delX[(row,col)][cond].A1
                    delY = total.delY[(row,col)][cond].A1
                    delZ = total.delZ[(row,col)][cond].A1
                    _H_label = sp.csr_matrix((hopping['t'](delX,delY,delZ),cond_mat.indices,cond_mat.indptr),shape=cond_mat.shape)
                else:                   # Hopping parameters case
                    NN = total.NN[(row,col)].multiply((total.NN[(row,col)]<=hopping['cutoff']).astype('int'))
                    t = hopping['t'][NN.data-1]
                    _H_label = sp.csr_matrix((t,NN.indices,NN.indptr),shape=NN.shape)
                self.H0[(row,col)] = _H_label


class H_Haldane:
    """
        Haldane(or Kane-mele, intrinsic)-SOC Hamiltonian
            H_uu = 1j * t * nu * c^dagger_i * c_j
            H_dd = -H_uu  if isKM==True else +H_uu
        
        For obtain nu, Two labels about sublattices should be defined 'label' and 'label_sub'
    """
    def __init__(self, label, label_sub, isKM = True):
        self.isSpin = False
        self.isFinalized = False
        self.H_soc = np.empty(0)
        self.t_haldane = 0.0
        self.isKM = isKM
        if(self.isKM):
            self.isSpin = True
        self.label = label
        self.label_sub = label_sub
    def __call__(self, eikdr, isSpin=None, **kwargs_dummy):
        assert self.isFinalized
        if(isSpin==None):
            isSpin = self.isSpin
        H_soc_k = self.H_soc.multiply(eikdr)
        if(isSpin):
            return sp.bmat([[self.t_haldane*H_soc_k,None],[None,-self.t_haldane*H_soc_k if self.isKM else self.t_haldane*H_soc_k]])
        else:
            return H_soc_k
    def _E_const(self,**kwargs_dummy):
        E_const = 0.0
        return E_const
    def set_haldane(self, t_haldane):
        self.t_haldane = t_haldane
    def _finalize(self, total):
        self.N_orbit = total.position.N_orbit
        self.H_soc = sp.csr_matrix((self.N_orbit,self.N_orbit),dtype=np.complex64)
        self.__finalize(total)
        self.isFinalized = True
    def __finalize(self,total):
        cond = total.NN==1
        row, col = total.labels_to_rowcol(self.label, self.label)
        delX = total.delX.multiply(cond)[(row,col)]
        delY = total.delY.multiply(cond)[(row,col)]
        R_inv = total.R.multiply(cond)[(row,col)]
        R_inv.data = 1/R_inv.data
        delX = delX.multiply(R_inv)
        delY = delY.multiply(R_inv)
        # pickup a FNN sub-orbit
        row_sub, col_sub = total.labels_to_rowcol(self.label, self.label_sub)
        delX_sub = total.delX.multiply(cond)[(row_sub,col_sub)]
        delX_one = delX_sub.data[delX_sub.indptr[:-1]]
        delY_sub = total.delY.multiply(cond)[(row_sub,col_sub)]
        delY_one = delY_sub.data[delY_sub.indptr[:-1]]
        #   normalization
        R_one = np.sqrt(delX_one**2 + delY_one**2)
        delX_one /= R_one
        delY_one /= R_one
        # cross product with the FNN sub-orbit
        nu = delX.T.multiply(delY_one).T - delY.T.multiply(delX_one).T
        nu.data = np.where((np.arcsin(nu.data)//(pi/3))%2==1, 1, -1)
        self.H_soc[(row,col)] = 1j*nu

class H_Haldane_strained(H_Haldane):
    """
        Haldane(or Kane-mele, intrinsic)-SOC Hamiltonian
            H_uu = 1j * t * nu * c^dagger_i * c_j
            H_dd = -H_uu  if isKM==True else +H_uu
        
        For obtain nu, Two labels about sublattices should be defined 'label' and 'label_sub'
    """
    def __init__(self, label, label_sub, isKM = True):
        super().__init__(label, label_sub, isKM)
    def __finalize(self,total):
        cond = total.NN==1
        row, col = total.labels_to_rowcol(self.label, self.label)
        delX = total.delX.multiply(cond)[(row,col)]
        delY = total.delY.multiply(cond)[(row,col)]
        R_inv = total.R.multiply(cond)[(row,col)]
        R_inv.data = 1/R_inv.data
        delX = delX.multiply(R_inv)
        delY = delY.multiply(R_inv)
        # pickup a FNN sub-orbit
        row_sub, col_sub = total.labels_to_rowcol(self.label, self.label_sub)
        delX_sub = total.delX.multiply(cond)[(row_sub,col_sub)]
        delX_one = delX_sub.data[delX_sub.indptr[:-1]]
        delY_sub = total.delY.multiply(cond)[(row_sub,col_sub)]
        delY_one = delY_sub.data[delY_sub.indptr[:-1]]
        #   normalization
        R_one = np.sqrt(delX_one**2 + delY_one**2)
        delX_one /= R_one
        delY_one /= R_one
        # cross product with the FNN sub-orbit
        nu = delX.T.multiply(delY_one).T - delY.T.multiply(delX_one).T
        nu.data = np.where((np.arcsin(nu.data)//(pi/3))%2==1, 1, -1)
        # strained effect
        strain = R_one.copy()
        strain.data = np.exp(-3.37*(strain.data-2.46)/2.46)
        self.H_soc[(row,col)] = 1j*nu.multiply(strain)


class H_Rashba:
    """
        Rashba(or PIA)-SOC Hamiltonian
            H_ud = 1j * (s X d_ij)_z * t * c^dagger_i * c_j     where d_ij is a normalized vector pointing from site i to site j
    """
    def __init__(self, label_i, label_j):
        self.isSpin = True
        self.isFinalized = False
        self.H_soc = np.empty(0)
        self.t_rashba = 0.0
        self.label_i = label_i
        self.label_j = label_j
    def __call__(self, eikdr, **kwargs_dummy):
        assert self.isFinalized
        H_soc_k = self.H_soc.multiply(eikdr)
        return sp.bmat([[None,self.t_rashba*H_soc_k],[self.t_rashba*H_soc_k.conj().T,None]])
    def _E_const(self,**kwargs_dummy):
        E_const = 0.0
        return E_const
    def set_rashba(self, t_rashba):
        self.t_rashba = t_rashba
    def _finalize(self, total):
        self.N_orbit = total.position.N_orbit
        self.H_soc = sp.csr_matrix((self.N_orbit,self.N_orbit),dtype=np.complex64)
        self.__finalize(total)
        if(self.label_i!=self.label_j):
            self.label_i, self.label_j = self.label_j, self.label_i
            self.__finalize(total)
        self.isFinalized = True
    def __finalize(self,total):
        cond = total.NN==1
        row, col = total.labels_to_rowcol(self.label_i, self.label_j)
        delX = total.delX.multiply(cond)[(row,col)]
        delY = total.delY.multiply(cond)[(row,col)]
        R_inv = total.R.multiply(cond)[(row,col)]
        R_inv.data = 1/R_inv.data
        delX = delX.multiply(R_inv)
        delY = delY.multiply(R_inv)
        self.H_soc[(row,col)] = delX + 1j*delY
        


class H_Rashba_strained(H_Rashba):
    def __init__(self, label_i, label_j):
        super().__init__(label_i, label_j)
    def __finalize(self,total):
        cond = total.NN==1
        row, col = total.labels_to_rowcol(self.label_i, self.label_j)
        delX = total.delX.multiply(cond)[(row,col)]
        delY = total.delY.multiply(cond)[(row,col)]
        R = total.R.multiply(cond)[(row,col)]
        R_inv = R.copy()
        R_inv.data = 1/R_inv.data
        delX = delX.multiply(R_inv)
        delY = delY.multiply(R_inv)
        # strained effect
        strain = R.copy()
        strain.data = np.exp(-3.37*(strain.data-1.420281662)/1.420281662)
        self.H_soc[(row,col)] = (delX + 1j*delY).multiply(strain)
    


class H_PIA_strained(H_Rashba):
    def __init__(self, label_i, label_j):
        super().__init__(label_i, label_j)
    def __finalize(self,total):
        cond = total.NN==1
        row, col = total.labels_to_rowcol(self.label_i, self.label_j)
        delX = total.delX.multiply(cond)[(row,col)]
        delY = total.delY.multiply(cond)[(row,col)]
        R_inv = total.R.multiply(cond)[(row,col)]
        R_inv.data = 1/R_inv.data
        delX = delX.multiply(R_inv)
        delY = delY.multiply(R_inv)
        # strained effect
        strain = R.copy()
        strain.data = np.exp(-3.37*(strain.data-2.46)/2.46)
        self.H_soc[(row,col)] = (delX + 1j*delY).multiply(strain)
    



class H_U:
    """
        Hamiltonian for Hubbard U correction 
    """
    def __init__(self, isFock=False, label=None):
        self.label=label
        self.isSpin = False
        self.isFock = isFock
        if(self.isFock):
            self.isSpin = True
        self.isFinalized = False
        self.U = 0.0
    def __call__(self, eikdr, isSpin=None, density=0.5, Sxy=0.0,**kwargs_dummy):
        assert self.isFinalized
        if(isSpin==None):
            isSpin = self.isSpin
        if(isSpin):
            if(self.cond_label==1.0):
                cond_label = 1.0
            else:
                cond_label = np.hstack([self.cond_label,self.cond_label])
        else:
            cond_label = self.cond_label
        density = self.__density(isSpin,density*cond_label)
        if(self.isFock):   # Non-collinear case (size of density should be 2*atomic number)
            Sxy = self.__density(False,Sxy*self.cond_label)
            return sp.bmat([[sp.diags(self.U*density[self.N_orbit:]), sp.diags(-self.U*Sxy.conj())],
                            [sp.diags(-self.U*Sxy), sp.diags(self.U*density[:self.N_orbit])]])
        elif(isSpin):   # Spin-degeneracy case (size of density should be 2*atomic number)
            return sp.bmat([[sp.diags(self.U*density[self.N_orbit:]),None],
                            [None,sp.diags(self.U*density[:self.N_orbit])]])
        else:   # No-spin-degeneracy case (size of density should be same with atomic number)
            return sp.diags(self.U*density)
    def _E_const(self,isSpin, density, Sxy=0.0,**kwargs_dummy):
        assert self.isFinalized
        if(isSpin==None):
            isSpin = self.isSpin
        if(self.isFock):
            E_const = -self.U*np.dot(density[:self.N_orbit]*self.cond_label,density[self.N_orbit:]*self.cond_label) + self.U*np.vdot(Sxy*self.cond_label,Sxy*self.cond_label)
        elif(isSpin):
            E_const = -self.U*np.dot(density[:self.N_orbit]*self.cond_label,density[self.N_orbit:]*self.cond_label)
        else:
            print("Spin-degenerate case consider that the spin-up density same with the spin-dn density")
            E_const = -self.U*np.dot(density*self.cond_label,density*self.cond_label)
        return np.real_if_close(E_const)
    def _finalize(self,total):
        self.N_orbit = total.position.N_orbit
        if(self.label==None):
            self.cond_label = 1.0
        else:
            self.cond_label = total.position.label == self.label
        self.isFinalized = True
    def __density(self,isSpin,density):
        if(density is None):
            if(isSpin):
                return np.zeros(self.N_orbit*2)
            else:
                return np.zeros(self.N_orbit)
        elif(type(density) is float):
            if(isSpin):
                return np.ones(self.N_orbit*2)*density
            else:
                return np.ones(self.N_orbit)*density
        elif(density.size == self.N_orbit):
            if(isSpin):
                return np.append(density,density)
            else:
                return density
        elif(density.size == self.N_orbit*2):
            if(isSpin):
                return density
            else:
                raise RuntimeError(f"Check the dimension of density{density.size}")
        elif(self.N_orbit%density.size==0):
            if(isSpin):
                return np.tile(density,reps=int(2*self.N_orbit/density.size))
            else:
                return np.tile(density,reps=int(self.N_orbit/density.size))
        elif(self.N_orbit%int(density.size/2)==0):  # [spinup,spindn] case
            assert isSpin
            return np.append(np.tile(density[:int(density.size/2)],reps=int(self.N_orbit/int(density.size/2))) , np.tile(density[int(density.size/2):],reps=int(self.N_orbit/int(density.size/2))))
        else:
            raise RuntimeError(f"Check the dimension of density{density.size}")
    def set_U(self, U:float=0.0):
        self.U = U


class H_V:
    """
        Hamiltonian for Hubbard V correction (extended Hubbard model)
    """
    def __init__(self, isFock=False):
        self.isFock = isFock
        if(isFock):
            raise NotImplementedError()
        self.isSpin = False
        self.isFinalized = False
        self.V = np.array([],dtype=[('From','<U4'),('To','<U4'),('isFn',bool),('V',object),('cutoff',float)])
        self.H_V = np.empty(0)
    def __call__(self, eikdr, isSpin=None, density=0.5, Sxy=0.0,**kwargs_dummy):
        assert self.isFinalized
        if(isSpin==None):
            isSpin = self.isSpin
        density = self.__density(isSpin,density)
        if(isSpin):
            H_V = sp.bmat([[self.H_V,self.H_V],[self.H_V,self.H_V]])
            H_VN = H_V@density
            assert H_VN.ndim == 1
            return sp.diags(H_VN)
        else:
            H_VN = self.H_V@density
            assert H_VN.ndim == 1
            return sp.diags(2*H_VN)
    def _E_const(self,isSpin, density, Sxy=0.0,**kwargs_dummy):
        if(isSpin==None):
            isSpin = self.isSpin
        if(self.isFock):
            raise NotImplementedError()
            # E_const = 
        elif(isSpin):
            rho = density[:self.N_orbit]+density[self.N_orbit:]
            E_const = -rho@self.H_V@rho/2
        else:
            print("Spin-degeneracy case consider that spin-up density same with spin-dn density")
            rho = 2*density
            E_const = -rho@self.H_V@rho/2
        return np.real_if_close(E_const)
    def __density(self,isSpin,density):
        if(density is None):
            if(isSpin):
                return np.zeros(self.N_orbit*2)
            else:
                return np.zeros(self.N_orbit)
        elif(type(density) is float):
            if(isSpin):
                return np.ones(self.N_orbit*2)*density
            else:
                return np.ones(self.N_orbit)*density
        elif(density.size == self.N_orbit):
            if(isSpin):
                return np.append(density,density)
            else:
                return density
        elif(density.size == self.N_orbit*2):
            if(isSpin):
                return density
            else:
                raise RuntimeError(f"Check the dimension of density{density.size}")
        elif(self.N_orbit%density.size==0):
            if(isSpin):
                return np.tile(density,reps=int(2*self.N_orbit/density.size))
            else:
                return np.tile(density,reps=int(self.N_orbit/density.size))
        elif(self.N_orbit%int(density.size/2)==0):  # [spinup,spindn] case
            assert isSpin
            return np.append(np.tile(density[:int(density.size/2)],reps=int(self.N_orbit/int(density.size/2))) , np.tile(density[int(density.size/2):],reps=int(self.N_orbit/int(density.size/2))))
        else:
            raise RuntimeError(f"Check the dimension of density{density.size}")
    def set_V(self,label_from, label_to, V, isFn:bool=False, cutoff=1.0):
        """
            Usage : Same with set_hoppings in H0
        """
        condition_same = (self.V['From']==label_from)*(self.V['To']==label_to)
        if(np.sum(condition_same)!=0):
            print(f"The formal V parameter from {label_from} to {label_to} is removed")
            condition_same = np.logical_or(condition_same,(self.V['From']==label_to)*(self.V['To']==label_from))
            self.V = np.compress(np.logical_not(condition_same),self.V)
        if(isFn==False):
            if(type(V) is not np.ndarray):
                raise RuntimeError("V should be numpy.ndarray when isFn==False")
            if(cutoff==None):
                cutoff = V.size
        else:
            if(callable(V)==False):
                raise RuntimeError("V should be callable function when isFn==True")
            assert cutoff !=None, "Cutoff should be float-number when the V is a equation"
        self.V = np.append(self.V,np.array((label_from,label_to,isFn,V,cutoff),dtype=self.V.dtype))
        if(label_from!=label_to):
            self.V = np.append(self.V,np.array((label_to,label_from,isFn,V,cutoff),dtype=self.V.dtype))
    def _finalize(self,total):
        self.N_orbit = total.position.N_orbit
        # Cutoff test 
        i=0
        ##  Different label case
        for label_from, label_to in combinations(total.position.label_unique, 2):
            V = self.V[np.logical_and(self.V['From']==label_from,self.V['To']==label_to)]
            if(V.size==0):
                continue
            if(V['isFn']):
                assert total.cutoff >= V['cutoff'], "Your cutoff is too smaller than V range, (HV)"
            else:
                assert total.maxorders[i] >= V['cutoff'], "Your cutoff is too smaller than V order, (HV)"
            i += 1
        ##  Same label case
        for label_from in total.position.label_unique:
            V = self.V[np.logical_and(self.V['From']==label_from,self.V['To']==label_from)]
            if(V.size==0):
                continue
            if(V['isFn']):
                assert total.cutoff >= V['cutoff'], "Your cutoff is too smaller than V range, (HV)"
            else:
                assert total.maxorders[i] >= V['cutoff'], "Your cutoff is too smaller than V order, (HV)"
            i += 1
        # Finalization start
        self.H_V = sp.csr_matrix((self.N_orbit,self.N_orbit),dtype=np.complex64)
        ##  Different label case
        for label_from, label_to in combinations(total.position.label_unique, 2):
            V = self.V[np.logical_and(self.V['From']==label_from,self.V['To']==label_to)]
            if(V.size!=0):
                V = V[0]
                row, col = total.labels_to_rowcol(label_from, label_to)
                if(V['isFn']):    # V equation case
                    nonzero = total.nonzero[(row,col)]
                    R_cut = total.R[(row,col)] <= V['cutoff']
                    cond_mat = nonzero.multiply(R_cut).astype('bool')  #logical_and
                    cond = cond_mat.nonzero()
                    delX = total.delX[(row,col)][cond].A1
                    delY = total.delY[(row,col)][cond].A1
                    delZ = total.delZ[(row,col)][cond].A1
                    _H_label = sp.csr_matrix((V['V'](delX,delY,delZ),cond_mat.indices,cond_mat.indptr),shape=cond_mat.shape)
                else:                   # V parameters case
                    NN = total.NN[(row,col)].multiply((total.NN[(row,col)]<=V['cutoff']).astype('int'))
                    t = V['V'][NN.data-1]
                    _H_label = sp.csr_matrix((t,NN.indices,NN.indptr),shape=NN.shape)
                self.H_V[(row,col)] = _H_label
                self.H_V[(col,row)] = _H_label
        ##  Same label case
        for label_from in total.position.label_unique:
            V = self.V[np.logical_and(self.V['From']==label_from,self.V['To']==label_from)]
            if(V.size!=0):
                V = V[0]
                row, col = total.labels_to_rowcol(label_from, label_from)
                if(V['isFn']):    # V equation case
                    nonzero = total.nonzero[(row,col)]
                    R_cut = total.R[(row,col)] <= V['cutoff']
                    cond_mat = nonzero.multiply(R_cut).astype('bool')  #logical_and
                    cond = cond_mat.nonzero()
                    delX = total.delX[(row,col)][cond].A1
                    delY = total.delY[(row,col)][cond].A1
                    delZ = total.delZ[(row,col)][cond].A1
                    _H_label = sp.csr_matrix((V['V'](delX,delY,delZ),cond_mat.indices,cond_mat.indptr),shape=cond_mat.shape)
                else:                   # V parameters case
                    NN = total.NN[(row,col)].multiply((total.NN[(row,col)]<=V['cutoff']).astype('int'))
                    t = V['V'][NN.data-1]
                    _H_label = sp.csr_matrix((t,NN.indices,NN.indptr),shape=NN.shape)
                self.H_V[(row,col)] = _H_label
        self.isFinalized = True

class H_ex:
    """
        Applying exchange field 
    """
    def __init__(self, lambda_ex=None):
        self.isSpin = True
        self.isFinalized = False
        self.lambda_ex = lambda_ex
    def __call__(self, **kwargs_dummy):
        assert self.isFinalized
        H_ex = sp.diags(np.ones(self.N_orbit)*self.lambda_ex)
        return sp.bmat([[H_ex, None],[None,-H_ex]])
    def _E_const(self,**kwargs_dummy):
        E_const = 0.0
        return E_const
    def _finalize(self,total):
        assert self.lambda_ex is not None
        self.N_orbit = total.position.N_orbit
        self.isFinalized = True
    def set_lambda_ex(self, lambda_ex):
        self.lambda_ex = lambda_ex

