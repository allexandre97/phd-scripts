import numpy as np 
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter as sf


with open('properties.xvg', 'r') as f:
    lines = f.readlines()
    t = 4*np.array([float(_.split()[0]) for _ in lines[45:]])
    X = np.array([float(_.split()[1]) for _ in lines[45:]])
    Y = np.array([float(_.split()[2]) for _ in lines[45:]])
    Z = np.array([float(_.split()[3]) for _ in lines[45:]])
    Vxx = np.array([float(_.split()[4]) for _ in lines[45:]])
    Vxy = np.array([float(_.split()[5]) for _ in lines[45:]])
    Vxz = np.array([float(_.split()[6]) for _ in lines[45:]])
    Vyx = np.array([float(_.split()[7]) for _ in lines[45:]])
    Vyy = np.array([float(_.split()[8]) for _ in lines[45:]])
    Vyz = np.array([float(_.split()[9]) for _ in lines[45:]])
    Vzx = np.array([float(_.split()[10]) for _ in lines[45:]])
    Vzy = np.array([float(_.split()[11]) for _ in lines[45:]])
    Vzz = np.array([float(_.split()[12]) for _ in lines[45:]])
    Pxx = np.array([float(_.split()[13]) for _ in lines[45:]])
    Pxy = np.array([float(_.split()[14]) for _ in lines[45:]])
    Pxz = np.array([float(_.split()[15]) for _ in lines[45:]])
    Pyx = np.array([float(_.split()[16]) for _ in lines[45:]])
    Pyy = np.array([float(_.split()[17]) for _ in lines[45:]])
    Pyz = np.array([float(_.split()[18]) for _ in lines[45:]])
    Pzx = np.array([float(_.split()[19]) for _ in lines[45:]])
    Pzy = np.array([float(_.split()[20]) for _ in lines[45:]])
    Pzz = np.array([float(_.split()[21]) for _ in lines[45:]])
    T = np.array([float(_.split()[22]) for _ in lines[45:]])
f.close()


BOX = np.zeros((t.shape[0], 3))
BOX[:,0] += X
BOX[:,1] += Y
BOX[:,2] += Z

PRE = np.zeros((t.shape[0],3,3))
PRE[:,0,0] += Pxx
PRE[:,0,1] += Pxy
PRE[:,0,2] += Pxz
PRE[:,1,0] += Pyx
PRE[:,1,1] += Pyy
PRE[:,1,2] += Pyz
PRE[:,2,0] += Pzx
PRE[:,2,1] += Pzy
PRE[:,2,2] += Pzz

VIR = np.zeros((t.shape[0],3,3))
VIR[:,0,0] += Vxx
VIR[:,0,1] += Vxy
VIR[:,0,2] += Vxz
VIR[:,1,0] += Vyx
VIR[:,1,1] += Vyy
VIR[:,1,2] += Vyz
VIR[:,2,0] += Vzx
VIR[:,2,1] += Vzy
VIR[:,2,2] += Vzz

POW = np.zeros_like(VIR)
POW[:-2] += (VIR[2:] - VIR[:-2])/(2*(t[1]-t[0]))

n = 50

BO = sf(BOX, 80*n, 2, axis = 0)
PR = sf(PRE, 80*n, 2, axis = 0)
VR = sf(VIR, 80*n, 2, axis = 0)
TM = sf(T, 80*n, 2, axis = 0)

VL = BO[:,0]*BO[:,1]*BO[:,2]

KE = np.zeros_like(VR)

for i in range(KE.shape[0]):
    KE[i] += VL[i]*PR[i] + VR[i]

GM = BO[:,2]/2 + (PR[:,2,2] - (PR[:,0,0]+PR[:,1,1])/2)

loc = 'upper right'

fig, ax = plt.subplots(2, 2, figsize = (14,14), sharex = True)

a = ax[0,0].plot(t/1000, VR[:,0,0], label = 'Virial xx', c = 'royalblue')
b = ax[0,0].plot(t/1000, VR[:,1,1], label = 'Virial yy', c = 'seagreen')
c = ax[0,0].plot(t/1000, VR[:,2,2], label = 'Virial zz', c = 'firebrick')
ax[0,0].set_xlabel('Time (ns)', fontsize = 20, fontweight = 'bold')
ax[0,0].set_ylabel('Virial\n(kJ/mol)', fontsize = 20, fontweight = 'bold')
ax[0,0].tick_params(labelsize = 18)
ax2 = ax[0,0].twinx()
d = ax2.plot(t/1000, PR[:,2,2], label = 'P', c = 'k')
ax2.set_ylabel('Pressure (bar)', fontsize = 20, fontweight = 'bold')
ax2.tick_params(labelsize = 18)
abcd = a+b+c+d
labs = [_.get_label() for _ in abcd]
ax[0,0].legend(abcd, labs, fontsize=12, loc = loc)

a = ax[0,1].plot(t/1000, VR[:,0,1], label = 'Virial xy/yx', c = 'royalblue')
b = ax[0,1].plot(t/1000, VR[:,0,2], label = 'Virial xz/zx', c = 'seagreen')
c = ax[0,1].plot(t/1000, VR[:,1,2], label = 'Virial yz/zy', c = 'firebrick')
ax[0,1].set_xlabel('Time (ns)', fontsize = 20, fontweight = 'bold')
ax[0,1].set_ylabel('Virial\n(kJ/mol)', fontsize = 20, fontweight = 'bold')
ax[0,1].tick_params(labelsize = 18)
ax2 = ax[0,1].twinx()
d = ax2.plot(t/1000, PR[:,2,2], label = 'P', c = 'k')
ax2.set_ylabel('Pressure (bar)', fontsize = 20, fontweight = 'bold')
ax2.tick_params(labelsize = 18)
abcd = a+b+c+d
labs = [_.get_label() for _ in abcd]
ax[0,1].legend(abcd, labs, fontsize=12, loc = loc)


a = ax[1,0].plot(t/1000, KE[:,0,0], label = 'Ek xx', c = 'royalblue')
b = ax[1,0].plot(t/1000, KE[:,1,1], label = 'Ek yy', c = 'seagreen')
c = ax[1,0].plot(t/1000, KE[:,2,2], label = 'Ek zz', c = 'firebrick')
ax[1,0].set_xlabel('Time (ns)', fontsize = 20, fontweight = 'bold')
ax[1,0].set_ylabel('Kinetic Energy\n(kJ/mol)', fontsize = 20, fontweight = 'bold')
ax[1,0].tick_params(labelsize = 18)
ax2 = ax[1,0].twinx()
d = ax2.plot(t/1000, PR[:,2,2], label = 'P', c = 'k')
ax2.set_ylabel('Pressure (bar)', fontsize = 20, fontweight = 'bold')
ax2.tick_params(labelsize = 18)
abcd = a+b+c+d
labs = [_.get_label() for _ in abcd]
ax[1,0].legend(abcd, labs, fontsize=12, loc = loc)



a = ax[1,1].plot(t/1000, KE[:,0,1], label = 'Ek xy/yx', c = 'royalblue')
b = ax[1,1].plot(t/1000, KE[:,0,2], label = 'Ek xz/zx', c = 'seagreen')
c = ax[1,1].plot(t/1000, KE[:,1,2], label = 'Ek yz/zy', c = 'firebrick')
ax[1,1].set_xlabel('Time (ns)', fontsize = 20, fontweight = 'bold')
ax[1,1].set_ylabel('Kinetic Energy\n(kJ/mol)', fontsize = 20, fontweight = 'bold')
ax[1,1].tick_params(labelsize = 18)
ax2 = ax[1,1].twinx()
d = ax2.plot(t/1000, PR[:,2,2], label = 'P', c = 'k')
ax2.set_ylabel('Pressure (bar)', fontsize = 20, fontweight = 'bold')
ax2.tick_params(labelsize = 18)
abcd = a+b+c+d
labs = [_.get_label() for _ in abcd]
ax[1,1].legend(abcd, labs, fontsize=12, loc = loc)

plt.tight_layout()
plt.show()
