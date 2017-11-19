#from pytriqs.gf.local import GfImFreq, iOmega_n, inverse
from pytriqs.gf.local import *
from pytriqs.archive import HDFArchive
from matplotlib import pyplot as plt
#from pytriqs.plot.mpl_interface import oplot
import numpy
#
# Create the Matsubara-frequency Green's function and initialize it
#g = GfImFreq(indices = [1], beta = 50, n_points = 1000, name = "$G_\mathrm{imp}$")
#g << inverse( iOmega_n + 0.5 )
#
#from pytriqs.plot.mpl_interface import oplot
#oplot(g, '-o',  x_window  = (0,10))

ar = HDFArchive('bethe.out.h5', 'r')

Sigma_iw = ar['dmft_out']['Sigma_iw']['up']
beta = Sigma_iw.beta
#print "beta ", beta
first_index = Sigma_iw.mesh.first_index()
last_index = Sigma_iw.mesh.last_index()
num_positive_freq = last_index + 1

x = numpy.array([(2*n+1)*numpy.pi/beta for n in range(num_positive_freq)])
y = Sigma_iw[0,0].data[-first_index:,0,0].imag

#print(x.shape)
#print(y.shape)

for i in range(len(x)):
    print x[i], y[i]

plt.plot(x**0.5, -y, marker='x')
plt.xlim([0, 5])
plt.ylim([0, 3])
plt.savefig("plot.pdf")
#plt.show()

