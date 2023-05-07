__all__ = ['rspace', 'graphene', 'kspace']

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from .debug import *
from .base import *


class _SPACE:   # Mechanism for points in space
    def __init__(self):
        self.point = []
        self.isFinalized = False
        self.dtype = [('x',float),('y',float),('z',float)]
    @property
    def N_point(self):
        return len(self.point)
    @property
    def Label_unique(self):
        return np.unique(position.label)
    def _append(self,point):
        if(self.isFinalized):
            self.point = np.vstack((self.point,point))
        else:
            if(isinstance(point, np.ndarray)):
                self.point.append(point.tolist())
            else:
                self.point.append(point)
    def _add_point(self,x,y,z):
        if(self.isFinalized):
            self.point = np.append(self.point,np.asarray((float(x),float(y),float(z)),dtype=self.dtype))
        else:
            self.point.append([float(x),float(y),float(z)])
    def _finalize(self):
        if(not self.isFinalized):
            self.point = np.asarray(self.point)
            self.point.dtype = self.dtype
            self.isFinalized = True
    @staticmethod
    def __isSameSide(point,point_ref,path1,path2):
        L1 = np.cross((path2-path1),(point_ref-path1))
        L2 = np.cross((path2-path1),(point-path1))
        return L1*L2 >= 0
    def _isInside(self,point,point_ref,Path):
        isInside = np.ones(point.shape[0])
        for i in range(Path.shape[0]):
            isInside = isInside*self.__isSameSide(point,point_ref,Path[i-1],Path[i])
        return isInside
    def _isSamepoint_XY(self,point):
        return (abs((point[0]-self.point['x'])**2+(point[1]-self.point['y'])**2)<eps).flatten()


class rspace(_SPACE):   # x, y, z, label
    def __init__(self, Tvector):
        super().__init__()
        self.label = np.array([],dtype='<U4')
        self.Tvector = np.asarray(Tvector)
    def __add__(self,other):
        self.__iadd__(ohter)
        return self
    def __iadd__(self,other):
        super()._append(other.point)
        self.label = np.append(self.label, other.label)
        return self
    def set_orbit_array(self,points,labels):
        assert isinstance(points, np.ndarray) and isinstance(labels, np.ndarray), "Inputs should be ndarray"
        self.point = points
        self.label = np.array(labels,dtype=self.label.dtype)
        super()._finalize()
    def add_orbit(self,x,y,z,label):
        super()._add_point(x,y,z)
        self.label = np.append(self.label,label)
    @property
    def N_orbit(self):
        return super().N_point
    @property
    def label_unique(self):
        return np.unique(self.label)
    def make_supercell(self,range_T1=[0,0],range_T2=[0,0]):
        super()._finalize()
        assert (isinstance(range_T1[0], int)) and (isinstance(range_T1[1], int)) and (range_T1[1]>=0) and (range_T1[0]<=0)
        assert (isinstance(range_T2[0], int)) and (isinstance(range_T2[1], int)) and (range_T2[1]>=0) and (range_T2[0]<=0)
        copy_point = self.point.copy()
        copy_label = self.label.copy()
        for N1 in range(range_T1[0],range_T1[1]+1):
            for N2 in range(range_T2[0],range_T2[1]+1):
                if(N1==0 and N2==0):
                    continue
                append_point = copy_point.copy()
                append_point['x'] += N1*self.Tvector[0,0] + N2*self.Tvector[1,0]
                append_point['y'] += N1*self.Tvector[0,1] + N2*self.Tvector[1,1]
                super()._append(append_point)
                self.label = np.append(self.label, copy_label)
        self.Tvector[0] = self.Tvector[0]*(range_T1[1]-range_T1[0]+1)
        self.Tvector[1] = self.Tvector[1]*(range_T2[1]-range_T2[0]+1)
    def sliding(self,label=None,delX=0,delY=0):
        '''
            If the 'label' is None, all atoms will be slided.
        '''
        super()._finalize()
        cond = self.__condition_label(label)
        self.point['x'][cond] += delX
        self.point['y'][cond] += delY
    def rotating(self,label=None,theta:'unit of degree'=0):
        '''
            If the 'label' is None, all atoms will be rotated.
        '''
        super()._finalize()
        cond = self.__condition_label(label)
        theta_rad = np.deg2rad(theta)
        X_tmp = self.point[cond]['x']*np.cos(theta_rad) - self.point['y']*np.sin(theta_rad)
        T_tmp = self.point[cond]['x']*np.sin(theta_rad) + self.point['y']*np.cos(theta_rad)
        self.point['x'][cond] = X_tmp
        self.point['y'][cond] = Y_tmp
    def rotating_Tvec(self,theta:'unit of degree'=0):
        theta_rad = np.deg2rad(theta)
        X_tmp = self.Tvector[:,0]*np.cos(theta_rad) - self.Tvector[:,1]*np.sin(theta_rad)
        Y_tmp = self.Tvector[:,0]*np.sin(theta_rad) + self.Tvector[:,1]*np.cos(theta_rad)
        self.Tvector[:,0]=X_tmp
        self.Tvector[:,1]=Y_tmp
    def slicing_rectangle(self,label=None,axis='x',bound=[0,10]):
        super()._finalize()
        assert axis == 'x' or axis =='y', "axis should be 'x' or 'y'"
        cond_label = self.__condition_label(label)
        cond_axis = np.logical_and(self.point[axis] > bound[0], self.point[axis]<bount[1])
        self.compress(np.logical_or(cond_axis,np.logical_not(cond_label)))
    def slicing_unitcell(self):
        super()._finalize()
        cond = super()._isInside(self.point[['x','y']], (self.Tvector[0]+self.Tvector[1])/2, self.__cell_path())
        self.compress(cond)
    def compress(self,condition):
        super()._finalize()
        self.point = np.compress(condition,self.point)
        self.label = np.compress(condition,self.label)
    def __condition_label(self,label):
        if(label==None):
            return np.ones_like(self.label,dtype=bool)
        else:
            if(label in self.label):
                return self.label == label
            else:
                raise Warning("The label is not finded in label-list")
                return np.zeros_like(self.label,dtype=bool)
    def __cell_path(self):
        return np.array([[0,0],
                        self.Tvector[0],
                        (self.Tvector[0]+self.Tvector[1]),
                        self.Tvector[1],
                        [0,0]])-eps
    def plot(self):
        super()._finalize()
        for label in np.unique(self.label):
            cond = self.__condition_label(label)
            plt.scatter(self.point[cond]['x'], self.point[cond]['y'], label=label)
        Cellpath = self.__cell_path()
        plt.plot(Cellpath[:,0],Cellpath[:,1],'r:')
        plt.axhline(y=0,linewidth=1, color='k',linestyle='--')
        plt.axvline(x=0,linewidth=1, color='k',linestyle='--')
        plt.legend()
        plt.show()
        


class kspace(_SPACE):   # x, y, z, weight
    """
        2D kspace point generator
    """
    def __init__(self, Tvector):
        super().__init__()
        self.weight = np.array([],dtype=float)
        self.Tvector = Tvector[0:2,0:2]
        self.Gvector = 2*pi*np.linalg.inv(self.Tvector.T)
    def __call__(self, index=None):
        self._finalize()
        if(index == None):
            return self.point, self.weight
        else:
            return self.point[index], self.weight[index]
    def __iter__(self): ## make iterator for "for-loop"
        """
            return index, kx, ky, weight
        """
        self._finalize()
        self.__index =0
        self.__final = self.N_kpt
        return self
    def __next__(self): ## make iterator for "for-loop"
        """
            return index, kx, ky, weight
        """
        if self.__index < self.__final:
            kpoint = self.point[self.__index]
            weight = self.weight[self.__index] 
            self.__index = self.__index+1
            return self.__index-1, float(kpoint['x']), float(kpoint['y']), weight 
        else:
            raise StopIteration
    @property
    def N_kpt(self):
        return super().N_point
    @property
    def Sum_weight(self):
        return np.sum(self.weight)
    def weight_normalization(self):
        self.weight /= self.Sum_weight
    def add_kpoint(self,kx,ky,weight=1.0):
        if(abs(kx)<1e-10 and abs(ky)<1e-10):
            kx += eps
            ky += eps
        super()._add_point(kx,ky,0.0)
        self.weight = np.append(self.weight,weight)
    def gen_MonkhorstPack(self,dim):
        '''
            Monkhorst-Pack mesh generation.
            If you want consider K and M point, the dimension of mesh must be mod3.
        eg)
            dim = [9,9] for 9X9 mesh
        '''
        for i in range(dim[0]):
            for j in range(dim[1]):
                K = i/dim[0]*self.Gvector[0] + j/dim[1]*self.Gvector[1]
                self.add_kpoint(K[0], K[1], 1.0)
    def find_b_GRgrid(self,dim,N_sym=1, g_add=None):
        '''
            Print the b values of symmetry-preserving generalized regular k-grid algorithm.
            INPUT :
                dim : Dimension for grid mesh
                N_sym : Symmetry number of major symmetry(eg. C3 sym -> 3)
                g_add(optional) : List of addtional symmetry matrix (must be 2X2 shape)
                      ( eg. for inversion symmetry, g_add = [np.array([[-1,0],[0,-1]])] )
        '''
        G = self.Gvector[0:2,0:2]
        T = self.Tvector[0:2,0:2]
        g_list = self.__g_list_symmetry(N_sym, g_add)
        for b in range(0,dim[1]):
            isSym =  True
            H = np.array([[dim[0], b],
                        [0, dim[1]]])
            for g in g_list:
                X = (T) @ g @ np.linalg.inv(T)
                assert abs(abs(np.linalg.det(X))-1)<eps
                M = H@X@np.linalg.inv(H)
                isSym = isSym*np.all(np.abs(np.rint(M)-M)<eps)
            if(isSym):
                print(f"b = {b} is symmetry preserving")
                K = np.linalg.inv(H).T @ G
                kpoints = []
                for i in range(-2,3):
                    for j in range(-2,3):
                        if(i==0 and j==0): continue
                        kpoints.append(np.array([i,j])@K)
                distance, counts = np.unique(np.round(np.linalg.norm(np.asarray(kpoints),axis=1),decimals=10), return_counts=True)
                print(f" => Distance = { distance[0:3]}, Number of vectors = {counts[0:3]}")
    def gen_IrreducibleGeneralizedRegular(self,dim,N_sym=1,b=0, g_add=None):
        '''
            Irreducible generalized regular mesh generation.
            INPUT :
                dim : Dimension for grid mesh
                b : selected value from "find_b_GRgrid" function
                N_sym : Symmetry number of major symmetry(eg. C3 sym -> 3)
                g_add(optional) : List of addtional symmetry matrix (must be 2X2 shape)
                      ( eg. for inversion symmetry, g_add = [np.array([[-1,0],[0,-1]])] )
        '''
        G = self.Gvector[0:2,0:2]
        T = self.Tvector[0:2,0:2]
        H = np.array([[dim[0], b],
                    [0, dim[1]]])
        K = np.linalg.inv(H).T @ G
        kpoints = np.array([[i,j] for i in range(dim[0]) for j in range(dim[1])]) @ K
        kpoints -= np.floor((kpoints @ np.linalg.inv(G))+eps)@G-eps
        kpoints_frac = kpoints @ np.linalg.inv(G)
        g_list = self.__g_list_symmetry(N_sym, g_add)
        hashtable = np.zeros(dim[0]*dim[0]*dim[1]*dim[1])
        weight = []
        First = []
        uniquecount = 0
        for i in range(kpoints.shape[0]):
            index = H@kpoints_frac[i]
            index = np.mod(index,dim)
            index = np.rint(index[1]*dim[0]+index[0]).astype(int)
            if(hashtable[index]==0):
                uniquecount += 1
                hashtable[index] = uniquecount
                First.append(i)
                weight.append(1)
                for g in g_list:
                    g_frac = G@g@np.linalg.inv(G)
                    kpt_rot = kpoints_frac[i]@g_frac
                    kpt_rot = np.around((kpt_rot-np.floor(kpt_rot))%(1-1e-7),decimals=7)
                    assert np.all(np.logical_and(kpt_rot<1,kpt_rot>=0))
                    index = H @ kpt_rot
                    index = np.mod(index,dim)
                    assert index[1]<dim[1]
                    index = np.rint(index[1]*dim[0]+index[0]).astype(int)
                    if(hashtable[index]==0):
                        hashtable[index] = uniquecount
                        weight[-1] += 1
        assert np.sum(weight)==kpoints.shape[0],f"{np.sum(weight)}!={kpoints.shape[0]}"
        print(f"{len(kpoints)} kpoints reduced to {len(First)}!")
        for [kx,ky],w in zip(kpoints[First],weight):
            self.add_kpoint(kx, ky, w)
        self.weight_normalization()

        fig, ax = plt.subplots(figsize=(8, 12))
        ax.scatter(kpoints[:,0],kpoints[:,1],c='gray',s=10)
        ax.add_patch(
            patches.Polygon(
                # xy is a numpy array with shape Nx2.
                np.array([[0, 0], G[0], G[0]+G[1], G[1]]), # xy
                closed=True,
                edgecolor = 'black',
                linestyle = 'dashdot', 
                fill = False))
        ax.scatter(kpoints[First,0],kpoints[First,1],c=weight)
        plt.show() 
    def __g_list_symmetry(self,N_sym,g_add):
        g_list = []
        degree_0 = 2*pi/N_sym
        for n_degree in range(N_sym):
            degree = n_degree*degree_0
            g_list.append(np.array([[np.cos(degree),-np.sin(degree)],[np.sin(degree),np.cos(degree)]]))
        if(g_add!=None):
            g_list += g_add
        return g_list
    def gen_Densemesh(self,point,dim,scale):    # Should Change Later!!!!!!!!!!!!!!!!!!!!!
        ''' 
            Generation dense k-mesh near the Kpoint
            The dense mesh replace the points
        '''
        super()._finalize()
        condition_same = super()._isSamepoint_XY(point)
        if( np.sum(condition_same) != 1):
            raise RuntimeError(f"No kpoint at ({point[0]}, {point[1]}) point ")
        else:
            ## dense-point generation
            weight_point = self.weight[condition_same]
            points = np.empty(((2*dim[0]-1)*(2*dim[1]-1),2))
            densevector = self.Gvector*(scale/np.sqrt(np.sum(self.Gvector[0]**2)))
            index = 0
            for n_i in range(-dim[0]+1,dim[0]):
                for n_j in range(-dim[1]+1,dim[1]):
                    K_xy = n_i/dim[0]*densevector[0] + n_j/dim[1]*densevector[1]
                    points[index] = K_xy
                    index = index +1
            ## slice to hexagonal
            slice_point = np.array([(0.5+eps)*densevector[0],
                                (0.5+eps)*densevector[1],
                                (0.5+eps)*(densevector[1]-densevector[0]),
                                -(0.5+eps)*densevector[0],
                                -(0.5+eps)*densevector[1],
                                (0.5+eps)*(densevector[0]-densevector[1])])
            points = np.compress(super()._isInside(point=points,point_ref=np.array([0,0]),Path=slice_point),points,axis=0)
            N_dense = points.shape[0]
            self.point = np.compress(np.logical_not(condition_same),self.point)
            self.weight = np.compress(np.logical_not(condition_same),self.weight)
            for i in range(N_dense):
                self.add_kpoint(points[i][0]+point[0],points[i][1]+point[1],weight_point/N_dense)
    def gen_bandpath(self,HSpoint,N_grids):
        HSpoint = np.array(HSpoint)
        N_HSpoint = HSpoint.shape[0]
        for i_HS in range(N_HSpoint-1):
            for i_grid in range(N_grids[i_HS]):
                point = np.dot(HSpoint[i_HS] + (HSpoint[i_HS+1]-HSpoint[i_HS])*i_grid/N_grids[i_HS],self.Gvector)
                self.add_kpoint(kx=point[0],ky=point[1])
        point = np.dot(HSpoint[-1],self.Gvector)
        self.add_kpoint(kx=point[0],ky=point[1])
    @property
    def BZpath(self):
        raise NotImplementedError()
    def move_inBZ(self):
        raise NotImplementedError()
    def __cell_path(self):
        return np.array([[0,0],
                        self.Gvector[0],
                        (self.Gvector[0]+self.Gvector[1]),
                        self.Gvector[1],
                        [0,0]])-eps
    def plot(self):
        super()._finalize()
        plt.scatter(self.point['x'], self.point['y'])
        Cellpath = self.__cell_path()
        plt.plot(Cellpath[:,0],Cellpath[:,1],'r:')
        plt.axhline(y=0,linewidth=1, color='k',linestyle='--')
        plt.axvline(x=0,linewidth=1, color='k',linestyle='--')
        plt.show()


class graphene(rspace):
    def __init__(self, flag="G", N_x=1, N_y=1, a=2.46, z=0.0):
        assert isinstance(N_x,int) and isinstance(N_y,int)
        assert (N_x > 1) and (N_y > 1)
        Tvector = np.array([[],[]]) # should be changed
        super().__init__(Tvector)
        self.a = a
        self.__gen_graphene(N_x, N_y, z)
        raise NotImplementedError()
    def __gen_graphene(self, N_x, N_y, z):
        tan30 = np.tan(pi/6)
        N_atom = 2*N_x*N_y
        for i in range(N_atom):
            y_order, x_order = divmod(i,2*N_x)
            x = (0.5*x_order)*self.a
            y = (0.5*tan30*((x_order+y_order)%2) - y_order*1.5*tan30)*self.a
            label = label_flag + ('B'  if (x_order+y_order)%2 else 'A')
            super().add_orbit(label=label,x=x,y=y,z=z)
