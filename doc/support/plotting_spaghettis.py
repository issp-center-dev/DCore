from matplotlib import *
import matplotlib.pyplot as plt
import numpy as np

filename = 'spaghettis_to_plot.dat'      # Name of file
emin = -1.0                             # Minimum of energy range to plot
emax = 1.5                              # Maximum of energy range to plot
zmin = 0.0 #z.min()                     # Minimum of colour range
zmax = 5.0 #z.max()                     # Maximum of colour range
kpos = [0, 50, 100, 150, 200]           # Position of high-symmetry k-point 
kname = ['M', 'G', 'X', 'Z', 'M']       # Name of high-symmetry k-point

b = np.loadtxt(filename)
n_k = int(b[:,0][-1] + 1)
mesh_size = len(b[:,0])/n_k
klist = b[:,0][::mesh_size]
elist = b[:,1][0:mesh_size]
x,y = np.meshgrid(klist,elist)
z = b[:,2].reshape(n_k,mesh_size)

fig, ax = plt.subplots()
p = ax.pcolormesh(x,y,z.T, cmap=cm.Blues, vmin=zmin, vmax=zmax)
cb = fig.colorbar(p, ax=ax)
ax.set_xlim(klist.min(),klist.max())
ax.set_ylim(emin,emax)
ax.hlines(0.0,0,n_k-1)
plt.title(filename,fontsize=24)
for i in kpos: ax.vlines(i,emin,emax,alpha=0.5)
plt.xticks(kpos,kname,fontsize=16)
plt.ylabel('Energy (%s)'%'eV', fontsize=18)
plt.yticks(fontsize=16)
plt.grid()
plt.savefig(''.join(kname)+'_'+filename.replace('dat','png'),bbox_inches='tight')
