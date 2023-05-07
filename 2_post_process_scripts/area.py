
import numpy as np
import matplotlib.pyplot as plt
aG = 2.505
acc = aG/np.sqrt(3.0)
def crossProduct(ourpoint, point1, point2):
    x1 = point1[0]
    y1 = point1[1]
    x2 = point2[0]
    y2 = point2[1]
    xA = ourpoint[0]
    yA = ourpoint[1]
    v1 = (x2-x1, y2-y1)   # Vector 1
    v2 = (x2-xA, y2-yA)   # Vector 1
    crossProd = v1[0]*v2[1] - v1[1]*v2[0] 
    return crossProd
def insideAAbis(ourpoint,rAAp,rAApp):
    A = np.array([rAAp,0])
    B = np.array([rAAp/2.0,rAApp])
    C = np.array([-rAAp/2.0,rAApp])
    D = np.array([-rAAp,0])
    E = np.array([-rAAp/2.0,-rAApp])
    F = np.array([rAAp/2.0,-rAApp])
    if np.logical_and(np.less_equal(crossProduct(ourpoint,A,B),0.0), np.less_equal(crossProduct(ourpoint,B,C),0.0)):
        if np.logical_and(np.less_equal(crossProduct(ourpoint,C,D),0.0), np.less_equal(crossProduct(ourpoint,D,E),0.0)):
            if np.logical_and(np.less_equal(crossProduct(ourpoint,E,F),0.0), np.less_equal(crossProduct(ourpoint,F,A),0.0)):
                return np.array([True])
            else:
                return np.array([False])
        else:
            return np.array([False])
    else:
        return np.array([False])
    
def insideABbis(ourpointO,rAAp,rAApp):
    #ourpoint = np.array([ourpoint[0] - aCC,ourpoint[1]])
    d1 = np.sqrt(acc**2 - (acc/2.0)**2)
    d2 = 3.0/4.0*acc
    A = np.array([rAAp,0])
    B = np.array([rAAp/2.0,rAApp])
    C = np.array([-rAAp/2.0,rAApp])
    D = np.array([-rAAp,0])
    E = np.array([-rAAp/2.0,-rAApp])
    F = np.array([rAAp/2.0,-rAApp])
    #ourpointVec = [np.array([ourpointO[0]- acc ,ourpointO[1]]),
    #                 np.array([ourpointO[0]+acc/2.0,ourpointO[1]+aG/2]),
    #                 np.array([ourpointO[0]+acc/2.0,ourpointO[1]-aG/2])]
    ourpointVec = [np.array([ourpointO[0] ,ourpointO[1]- acc]),
                     np.array([ourpointO[0]+aG/2,ourpointO[1]+acc/2.0]),
                     np.array([ourpointO[0]-aG/2,ourpointO[1]+acc/2.0])]
    for ourpoint in ourpointVec:
        if np.logical_and(np.less_equal(crossProduct(ourpoint,A,B),0.0), np.less_equal(crossProduct(ourpoint,B,C),0.0)):
            if np.logical_and(np.less_equal(crossProduct(ourpoint,C,D),0.0), np.less_equal(crossProduct(ourpoint,D,E),0.0)):
                if np.logical_and(np.less_equal(crossProduct(ourpoint,E,F),0.0), np.less_equal(crossProduct(ourpoint,F,A),0.0)):
                    status = np.array([True])
                    break
                else:
                    status = np.array([False])
            else:
                status = np.array([False])
        else:
            status = np.array([False])
            
    return status
    
    
def insideBAbis(ourpointO,rAAp,rAApp):
    #ourpoint = np.array([ourpoint[0] - aCC,ourpoint[1]])
    d1 = np.sqrt(acc**2 - (acc/2.0)**2)
    d2 = 3.0/4.0*acc
    A = np.array([rAAp,0])
    B = np.array([rAAp/2.0,rAApp])
    C = np.array([-rAAp/2.0,rAApp])
    D = np.array([-rAAp,0])
    E = np.array([-rAAp/2.0,-rAApp])
    F = np.array([rAAp/2.0,-rAApp])

    ourpointVec = [np.array([ourpointO[0],ourpointO[1]+ acc ]),
                     np.array([ourpointO[0]+aG/2,ourpointO[1]-acc/2.0]),
                     np.array([ourpointO[0]-aG/2,ourpointO[1]-acc/2.0])]
    for ourpoint in ourpointVec:
        if np.logical_and(np.less_equal(crossProduct(ourpoint,A,B),0.0), np.less_equal(crossProduct(ourpoint,B,C),0.0)):
            if np.logical_and(np.less_equal(crossProduct(ourpoint,C,D),0.0), np.less_equal(crossProduct(ourpoint,D,E),0.0)):
                if np.logical_and(np.less_equal(crossProduct(ourpoint,E,F),0.0), np.less_equal(crossProduct(ourpoint,F,A),0.0)):
                    status = np.array([True])
                    break
                else:
                    status = np.array([False])
            else:
                status = np.array([False])
        else:
            status = np.array([False])
    
    return status
        
        
def insideSP1bis(ourpointO,rAA,w1,rAAp,rAApp):
    #ourpoint = np.array([ourpoint[0] - aCC,ourpoint[1]])
    A0 = np.array([rAA/2,0])
    B0 = np.array([rAAp/2.0,w1/2])
    C0 = np.array([-rAAp/2.0,w1/2])
    D0 = np.array([-rAA/2,0])
    E0 = np.array([-rAAp/2.0,-w1/2])
    F0 = np.array([rAAp/2.0,-w1/2])
    theta = -60*np.pi/180
    R = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    AR1 = np.matmul(R,A0)
    BR1 = np.matmul(R,B0)
    CR1 = np.matmul(R,C0)
    DR1 = np.matmul(R,D0)
    ER1 = np.matmul(R,E0)
    FR1 = np.matmul(R,F0)
    theta = 60*np.pi/180
    R = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    AR2 = np.matmul(R,A0)
    BR2 = np.matmul(R,B0)
    CR2 = np.matmul(R,C0)
    DR2 = np.matmul(R,D0)
    ER2 = np.matmul(R,E0)
    FR2 = np.matmul(R,F0)
    ourpointVec = [np.array([ourpointO[0] ,ourpointO[1]- acc/2.0]),
                     np.array([ourpointO[0]+aG/4,ourpointO[1]+acc/4.0]),
                     np.array([ourpointO[0]-aG/4,ourpointO[1]+acc/4.0])]
    for ourpoint,pointvec in zip(ourpointVec, [[A0,B0,C0,D0,E0,F0],[AR1,BR1,CR1,DR1,ER1,FR1],[AR2,BR2,CR2,DR2,ER2,FR2]]):
        #A,B,C,D,E,F = pointvec
        A,B,C,D,E,F = pointvec
        if np.logical_and(np.less_equal(crossProduct(ourpoint,A,B),0.0), np.less_equal(crossProduct(ourpoint,B,C),0.0)):
            if np.logical_and(np.less_equal(crossProduct(ourpoint,C,D),0.0), np.less_equal(crossProduct(ourpoint,D,E),0.0)):
                if np.logical_and(np.less_equal(crossProduct(ourpoint,E,F),0.0), np.less_equal(crossProduct(ourpoint,F,A),0.0)):
                    status = np.array([True])
                    break
                else:
                    status = np.array([False])
            else:
                status = np.array([False])
        else:
            status = np.array([False])
            
    return status
        
def insideSP2bis(ourpointO,rAA,w1,rAAp,rAApp):
    A0 = np.array([rAA/2,0])
    B0 = np.array([rAAp/2.0,w1/2])
    C0 = np.array([-rAAp/2.0,w1/2])
    D0 = np.array([-rAA/2,0])
    E0 = np.array([-rAAp/2.0,-w1/2])
    F0 = np.array([rAAp/2.0,-w1/2])
    theta = 60*np.pi/180
    R = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    AR1 = np.matmul(R,A0)
    BR1 = np.matmul(R,B0)
    CR1 = np.matmul(R,C0)
    DR1 = np.matmul(R,D0)
    ER1 = np.matmul(R,E0)
    FR1 = np.matmul(R,F0)
    theta = -60*np.pi/180
    R = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    AR2 = np.matmul(R,A0)
    BR2 = np.matmul(R,B0)
    CR2 = np.matmul(R,C0)
    DR2 = np.matmul(R,D0)
    ER2 = np.matmul(R,E0)
    FR2 = np.matmul(R,F0)
    #ourpointVec = [np.array([ourpointO[0]- acc/2.0 ,ourpointO[1]]),
    #                 np.array([ourpointO[0]+acc/4.0,ourpointO[1]+aG/4]),
    #                 np.array([ourpointO[0]+acc/4.0,ourpointO[1]-aG/4])]
    ourpointVec = [np.array([ourpointO[0] ,ourpointO[1]+ acc/2.0]),
                     np.array([ourpointO[0]+aG/4,ourpointO[1]-acc/4.0]),
                     np.array([ourpointO[0]-aG/4,ourpointO[1]-acc/4.0])]
    for ourpoint,pointvec in zip(ourpointVec, [[A0,B0,C0,D0,E0,F0],[AR1,BR1,CR1,DR1,ER1,FR1],[AR2,BR2,CR2,DR2,ER2,FR2]]):
        #A,B,C,D,E,F = pointvec
        A,B,C,D,E,F= pointvec    
        if np.logical_and(np.less_equal(crossProduct(ourpoint,A,B),0.0), np.less_equal(crossProduct(ourpoint,B,C),0.0)):
            if np.logical_and(np.less_equal(crossProduct(ourpoint,C,D),0.0), np.less_equal(crossProduct(ourpoint,D,E),0.0)):
                if np.logical_and(np.less_equal(crossProduct(ourpoint,E,F),0.0), np.less_equal(crossProduct(ourpoint,F,A),0.0)):
                    status = np.array([True])
                    break
                else:
                    status = np.array([False])
            else:
                status = np.array([False])
        else:
            status = np.array([False])
            
    return status
        
def insideSP3bis(ourpointO,rAA,w1,rAAp,rAApp):
    #ourpoint = np.array([ourpoint[0] - aCC,ourpoint[1]])
    A0 = np.array([rAA/2,0])
    B0 = np.array([rAAp/2.0,w1/2])
    C0 = np.array([-rAAp/2.0,w1/2])
    D0 = np.array([-rAA/2,0])
    E0 = np.array([-rAAp/2.0,-w1/2])
    F0 = np.array([rAAp/2.0,-w1/2])
    theta = 60*np.pi/180
    R = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    AR1 = np.matmul(R,A0)
    BR1 = np.matmul(R,B0)
    CR1 = np.matmul(R,C0)
    DR1 = np.matmul(R,D0)
    ER1 = np.matmul(R,E0)
    FR1 = np.matmul(R,F0)
    theta = -60*np.pi/180
    R = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    AR2 = np.matmul(R,A0)
    BR2 = np.matmul(R,B0)
    CR2 = np.matmul(R,C0)
    DR2 = np.matmul(R,D0)
    ER2 = np.matmul(R,E0)
    FR2 = np.matmul(R,F0)
    theta = -180*np.pi/180

    
    ourpointVec = [np.array([ourpointO[0]-aG/2,ourpointO[1]]),
                   np.array([ourpointO[0]+aG/2,ourpointO[1]]),
                   np.array([ourpointO[0]-aG/4,ourpointO[1]-acc*3.0/4.0]),
                    np.array([ourpointO[0]+aG/4,ourpointO[1]+acc*3.0/4.0]),
                    np.array([ourpointO[0]+aG/4,ourpointO[1]-acc*3.0/4.0]),
                     np.array([ourpointO[0]-aG/4,ourpointO[1]+acc*3.0/4.0])]
    for ourpoint,pointvec in zip(ourpointVec, [[A0,B0,C0,D0,E0,F0],[A0,B0,C0,D0,E0,F0],
                                              [AR1,BR1,CR1,DR1,ER1,FR1],[AR1,BR1,CR1,DR1,ER1,FR1],
                                              [AR2,BR2,CR2,DR2,ER2,FR2],[AR2,BR2,CR2,DR2,ER2,FR2]]):
                                              #[AR3,BR3,CR3,DR3,ER3,FR3],[AR4,BR4,CR4,DR4,ER4,FR4],[AR5,BR5,CR5,DR5,ER5,FR5]]):
        A,B,C,D,E,F = pointvec 

        if np.logical_and(np.less_equal(crossProduct(ourpoint,A,B),0.0), np.less_equal(crossProduct(ourpoint,B,C),0.0)):
            if np.logical_and(np.less_equal(crossProduct(ourpoint,C,D),0.0), np.less_equal(crossProduct(ourpoint,D,E),0.0)):
                if np.logical_and(np.less_equal(crossProduct(ourpoint,E,F),0.0), np.less_equal(crossProduct(ourpoint,F,A),0.0)):
                    status = np.array([True])
                    break
                else:
                    status = np.array([False])
            else:
                status = np.array([False])
        else:
            status = np.array([False])
            
    return status
thetaRAAVec = []
ratioVec = []
ratioVecSP = []

# data = np.genfromtxt("generate05.xyz",skip_header=4)
# data = np.genfromtxt("generate02.xyz",skip_header=4)
data = np.genfromtxt("generateWithBN.xyz",skip_header=4)
#data = np.genfromtxt("aap10.xyz",skip_header=4)

# data = np.genfromtxt("generate02p.xyz",skip_header=4)
# data = np.genfromtxt("generate05p.xyz",skip_header=4)

dataI = np.genfromtxt("generate.xyz",skip_header=4)

xdata = data[:,1]
ydata = data[:,2]
zdata = data[:,3]

xdataI = dataI[:,1]
ydataI = dataI[:,2]
zdataI = dataI[:,3]

cutoff = 10000

xdataC = xdata[np.logical_and(xdata<cutoff, ydata<cutoff)]
ydataC = ydata[np.logical_and(xdata<cutoff, ydata<cutoff)]
zdataC = zdata[np.logical_and(xdata<cutoff, ydata<cutoff)]

xdataIC = xdataI[np.logical_and(xdata<cutoff, ydata<cutoff)]
ydataIC = ydataI[np.logical_and(xdata<cutoff, ydata<cutoff)]
zdataIC = zdataI[np.logical_and(xdata<cutoff, ydata<cutoff)]

xdataC_top = xdataC[zdataIC>18]
ydataC_top = ydataC[zdataIC>18]
zdataC_top = zdataC[zdataIC>18]

xdataIC_top = xdataIC[zdataIC>18]
ydataIC_top = ydataIC[zdataIC>18]
zdataIC_top = zdataIC[zdataIC>18]

xdataC_middle = xdataC[np.logical_and(zdataIC<18, zdataIC>14)]
ydataC_middle = ydataC[np.logical_and(zdataIC<18, zdataIC>14)]
zdataC_middle = zdataC[np.logical_and(zdataIC<18, zdataIC>14)]

xdataIC_middle = xdataIC[np.logical_and(zdataIC<18, zdataIC>14)]
ydataIC_middle = ydataIC[np.logical_and(zdataIC<18, zdataIC>14)]
zdataIC_middle = zdataIC[np.logical_and(zdataIC<18, zdataIC>14)]

xdataC_bottom = xdataC[zdataIC<14]
ydataC_bottom = ydataC[zdataIC<14]
zdataC_bottom = zdataC[zdataIC<14]

xdataIC_bottom = xdataIC[zdataIC<14]
ydataIC_bottom = ydataIC[zdataIC<14]
zdataIC_bottom = zdataIC[zdataIC<14]

b_x = xdataC_middle
b_y = ydataC_middle
b_z = zdataC_middle

t_x = xdataC_top
t_y = ydataC_top
t_z = zdataC_top

b_xI = xdataIC_middle
b_yI = ydataIC_middle
b_zI = zdataIC_middle

t_xI = xdataIC_top
t_yI = ydataIC_top
t_zI = zdataIC_top

print(len(b_x),len(t_x))


layerIndex = np.genfromtxt("generate.e")[:,2]
sublattices = np.genfromtxt("generate.e")[:,1]
    
    
    
f = open("v","r")
g = open("v2","w")
    
neighVec = []
    
i = 0
for line in f:
    a = line.split()
    i = i + 1
    if (i%2==0):
        g.write(line)
        neighVec.append(a)
            
    
            
f.close()
g.close()
    
neighVecBN = np.asarray(neighVec)[zdataIC>18]
sublatticesBN = sublattices[zdataIC>18]


del_r = []
del_x = []
del_y = []
    
    #print(xdataBN)
for x1, y1, z1, x1I, y1I, z1I, neigh, SL in zip(t_x,t_y,t_z,t_xI,t_yI,t_zI,neighVecBN, sublattices):
    diffMin = 100
    dxMin = 100
    dyMin = 100
    for el in neigh:
        x2I = xdataIC[int(el)-1]
        y2I = ydataIC[int(el)-1]
        z2I = zdataIC[int(el)-1]
        x2 = xdataC[int(el)-1]
        y2 = ydataC[int(el)-1]
        z2 = zdataC[int(el)-1]
        diff1 = np.sqrt((x2I-x1I)**2 + (y2I-y1I)**2)# + (z2I-z1I)**2)\
        dx = x2-x1
        if (diff1 < diffMin and layerIndex[int(el)-1] == 1 and SL == sublattices[int(el)-1]):
            #print("hi")
            diffMin = diff1
            dxMin = x2-x1
            dyMin = y2-y1

    del_r.append(diffMin)
    del_x.append(dxMin)
    del_y.append(dyMin)

cutoffbis = 10000

with open('displacementsCarrFormat.txt', 'w') as f1:
    for x, y, z, disp, dxVal, dyVal  in zip(t_x, t_y, t_z, del_r, del_x, del_y):#, dzTop):
        if (x<cutoffbis and y<cutoffbis):
            f1.write(str(x) + "   " + str(y) + "   " + str(z) + "  " + str(disp)  + "   " + str(dxVal) + "   " + str(dyVal) +"\n")
            
            
Factor = 6
#FF = 1
rAA = 2*acc*np.sqrt(1/12)
w1 = acc/Factor
rAApp = (acc-w1)/2.0
rAAp = 2/np.sqrt(3)*rAApp
w2 = acc-2*rAApp
w3 = acc-2*rAApp

data = np.genfromtxt("displacementsCarrFormat.txt")#, dtype=(np.str, np.float, np.float, np.float))

xdata = data[:,0]
ydata = data[:,1]
zdataDX = data[:,4]
zdataDY = data[:,5]                    

#cp prepareSL.py mm_$i
num = 1200
x = np.linspace(min(xdata),max(xdata),num)
y = np.linspace(min(ydata),max(ydata),num)

from scipy.interpolate import griddata

zdataDXInterp=griddata((xdata,ydata),zdataDX,(x[None,:], y[:,None]), method='linear')
zdataDYInterp=griddata((xdata,ydata),zdataDY,(x[None,:], y[:,None]), method='linear')

zdataDXInterp1D = zdataDXInterp.flatten()
zdataDYInterp1D = zdataDYInterp.flatten()

zdataDXInterp1Dnan = zdataDXInterp1D[np.logical_not(np.logical_or(np.isnan(zdataDXInterp1D),np.isnan(zdataDYInterp1D)))]
zdataDYInterp1Dnan = zdataDYInterp1D[np.logical_not(np.logical_or(np.isnan(zdataDXInterp1D),np.isnan(zdataDYInterp1D)))]

Surface = [0 if insideAAbis(np.array([dx,dy]),rAAp,rAApp)
             else 3 if insideABbis(np.array([dx,dy]),rAAp,rAApp)
             else 4 if insideBAbis(np.array([dx,dy]),rAAp,rAApp)
             else 1 if insideSP1bis(np.array([dx,dy]),rAA,w1,rAAp,rAApp)
            #else 0 if np.logical_and(np.logical_and(insideSP3(np.array([dx,dy]),w3),insideSP2(np.array([dx,dy]),w2)),insideSP1(np.array([dx,dy]),w1))
             else 1 if insideSP2bis(np.array([dx,dy]),rAA,w1,rAAp,rAApp)
             else 1 if insideSP3bis(np.array([dx,dy]),rAA,w1,rAAp,rAApp)
             else 2
             for dx,dy in zip(zdataDXInterp1Dnan,zdataDYInterp1Dnan)] # Changed role of dx and dy here to match the conventions from sketch with the ones coming out of code

Sur = np.asarray(Surface)
AA=100*np.sum(Sur==0)/np.sum(Sur)
SP=100*np.sum(Sur==1)/np.sum(Sur)
EL=100*np.sum(Sur==2)/np.sum(Sur)
AB=100*np.sum(Sur==3)/np.sum(Sur)
BA=100*np.sum(Sur==4)/np.sum(Sur)
print(AA,AB,BA)
np.savetxt('ref6',(AA,AB,BA))

