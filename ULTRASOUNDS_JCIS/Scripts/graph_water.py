import numpy as np
from scipy.signal import savgol_filter as sf
import pickle
import matplotlib.pyplot as plt


'''Load pressure tensor from properties xvg file'''
with open('properties.xvg', 'r') as f:
    lines = f.readlines()
    t = 4*np.array([float(_.split()[0]) for _ in lines[45:]])/1000
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

'''Store pressure tensor in one array'''
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

'''Savitzsky-Golay filter on pressure to substract noise'''
PR = sf(PRE, 80*50+1, 2, axis = 0)

'''Minimum and maximum pressure values (plus a little bit of room around them)'''
pmin = np.min(PR[:,2,2])*(1+0.1)
pmax = np.max(PR[:,2,2])*(1+0.1)



with open('water_contact.pickle', 'rb') as f:    
    NC = pickle.load(f)
f.close()

N = int(np.ceil((500*NC.shape[0])/16666))

t = np.linspace(0, 500, N)

win = 8*50+1

fig, ax = plt.subplots(1, 2, figsize = (12, 5))

ax[0].plot(t, sf(NC[:N,0,0], win, 2), color = 'royalblue', label = 'Bead 1')
ax[0].plot(t, sf(NC[:N,1,0], win, 2), color = 'firebrick', label = 'Bead 2')
ax[0].plot(t, sf(NC[:N,2,0], win, 2), color = 'seagreen', label = 'Bead 3')
ax[0].plot(t, sf(NC[:N,3,0], win, 2), color = 'gold', label = 'Bead 4')

axP = ax[0].twinx()
axP.plot(t, PR, color = 'k')
axP.set_ylabel('P / bar', fontsize = 20, fontweight = 'bold')

ax[0].set_title('SN1 tail', fontsize = 22, fontweight = 'bold')
ax[0].set_xlabel('Time / ns', fontsize = 20, fontweight = 'bold')
ax[0].set_ylabel('# Contacts', fontsize = 20, fontweight = 'bold')
ax[0].tick_params(labelsize = 18)
ax[0].legend(fontsize = 18)

ax[1].plot(t, sf(NC[:N,0,1], win, 2), color = 'royalblue', label = 'Bead 1')
ax[1].plot(t, sf(NC[:N,1,1], win, 2), color = 'firebrick', label = 'Bead 2')
ax[1].plot(t, sf(NC[:N,2,1], win, 2), color = 'seagreen', label = 'Bead 3')
ax[1].plot(t, sf(NC[:N,3,1], win, 2), color = 'gold', label = 'Bead 4')

axP = ax[1].twinx()
axP.plot(t, PR, color = 'k')
axP.set_ylabel('P / bar', fontsize = 20, fontweight = 'bold')

ax[1].set_title('SN2 tail', fontsize = 22, fontweight = 'bold')
ax[1].set_xlabel('Time / ns', fontsize = 20, fontweight = 'bold')
ax[1].set_ylabel('# Contacts', fontsize = 20, fontweight = 'bold')
ax[1].tick_params(labelsize = 18)
ax[1].legend(fontsize = 18)

fig.tight_layout()

plt.savefig('water_contact.png', dpi = 300)
plt.close()
