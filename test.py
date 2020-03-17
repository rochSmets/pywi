
import numpy as np
import pyfftw
import matplotlib.pyplot as plt


dim = [80, 60]
domain = [[-4, 4], [-2, 2]]
cut = [2.0, 0.1]

inn = pyfftw.empty_aligned(dim, dtype='complex128')
out = pyfftw.empty_aligned(dim, dtype='complex128')

fftobj = pyfftw.FFTW(inn, out, axes=(0, 1), direction='FFTW_FORWARD', flags=('FFTW_MEASURE', ), threads = 8)

x = np.linspace(domain[0][0], domain[0][1], dim[0])
y = np.linspace(domain[1][0], domain[1][1], dim[1])

def heavy(x, cut) :
   return 0.5*(np.sign(x+cut)-np.sign(x-cut))

xp = heavy(x, cut[0])
yp = heavy(y, cut[1])

yv, xv = np.meshgrid(yp, xp)

r = np.multiply(xv, yv)
i = np.zeros(dim)

#r = np.random.randn(np.prod(dim)).reshape(dim)
#i = np.random.randn(np.prod(dim)).reshape(dim)

inn = r + 1j*i

fftobj.update_arrays(inn, out)

pyfftw.FFTW.execute(fftobj)

inp = np.absolute(inn)
mod = np.absolute(out)


def reorder(data) :

   n0 = data.shape[0]
   n1 = data.shape[1]

   new = np.empty(data.shape)

   new[0:(n0-1)/2, 0:(n1-1)/2] = data[n0/2+1:  , n1/2+1: ]
   new[0:(n0-1)/2, (n1-1)/2: ] = data[n0/2+1:  , 0:n1/2+1]
   new[(n0-1)/2: , 0:(n1-1)/2] = data[0:n0/2+1:, n1/2+1: ]
   new[(n0-1)/2: , (n1-1)/2: ] = data[0:n0/2+1:, 0:n1/2+1]

   return new

out = reorder(mod)


plt.imshow(inp, interpolation = 'nearest', origin='lower')

plt.show()

