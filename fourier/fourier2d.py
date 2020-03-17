
import numpy as np
import pyfftw


class Fourier2D() :

   """
   """

   def __init__(self,
                data,
                axis,
                field,
                display) :

      """

      """

      # rank = data.ndim
      dim = data.shape

      # arrays are empty and aligned (needed for fftw)
      inn = pyfftw.empty_aligned(dim, dtype = 'complex128', n = 16)
      out = pyfftw.empty_aligned(dim, dtype = 'complex128', n = 16)

      # this is the handle on the object that contains the plan
      fftobj = pyfftw.FFTW(inn, out, axes=(0, 1), direction='FFTW_FORWARD', flags=('FFTW_MEASURE', ), threads = 8)

      # np arrays have to be copied with a[...]=b[...] to guarantee
      # the alignment of the arrays. if not, if data have been parsed
      # with a given stride, inn will be aligned, but with a stride
      # different from the initial one ; in such a case, pyfftw refuse
      # to update the array...
      if 'float' in out.dtype.name :
          inn.real[...] = data[...]
          inn.imag = np.zeros(dim, dtype = float)
      elif 'complex' in data.dtype.name :
          inn[...] = data[...]
      else :
          raise NotImplementedError("the data type should be float or complex")

      outdat = np.empty(dim, dtype = float)

      # update the arrays with the asked values
      fftobj.update_arrays(inn, out)

      # then compute the fft
      pyfftw.FFTW.execute(fftobj)

      # reorder to get m=0 component at the center
      [axis, data] = self._reorder(axis, out)

      # reshape the output data to get real, imag, ... & full domain, half of it...
      [self.axis, self.data] = self._shape(axis, data, field = field, display = display)



   def _reorder(self,
                axis,
                data) :

      """
      this method for a 2d array put nyquist at the left (& right) & m=0
      at the middle
      """

      n0 = data.shape[0]
      n1 = data.shape[1]

      newx = np.empty(data.shape[0], dtype = float)
      newy = np.empty(data.shape[1], dtype = float)

      step = np.array([2.0*np.pi/axis[0][-1], 2.0*np.pi/axis[1][-1]])

      newx[0:(n0-1)/2] = np.linspace((1-n0)/2, -1, (n0-1)/2)*step[0]
      newx[(n0-1)/2:] = np.linspace(0, n0/2, n0/2+1)*step[0]

      newy[0:(n1-1)/2] = np.linspace((1-n1)/2, -1, (n1-1)/2)*step[1]
      newy[(n1-1)/2:] = np.linspace(0, n1/2, n1/2+1)*step[1]

      newz = np.empty(data.shape, dtype = complex)

      newz[0:(n0-1)/2, 0:(n1-1)/2] = data[n0/2+1:  , n1/2+1: ]
      newz[0:(n0-1)/2, (n1-1)/2: ] = data[n0/2+1:  , 0:n1/2+1]
      newz[(n0-1)/2: , 0:(n1-1)/2] = data[0:n0/2+1:, n1/2+1: ]
      newz[(n0-1)/2: , (n1-1)/2: ] = data[0:n0/2+1:, 0:n1/2+1]

      return [[newx, newy], newz]


   def _shape(self,
              axis,
              data,
              field,
              display) :

      """
      field can be real part, imaginary part, modulus or arg
      display can be up-right, up, right or full
      """



      if field == 'real' :
          out = data.real
      elif field == 'imag' :
          out = data.imag
      elif field == 'modulus' :
          out = np.absolute(data)
      elif field == 'arg' :
          out = np.angle(data)
      else :
          raise NotImplementedError("field has to be 're', 'im', 'mod' or 'arg'")

      i0 = out.shape[0]/2
      i1 = out.shape[1]/2

      if display == 'upright' :
          axis = [axis[0][i0:], axis[1][i1:]]
          data = out[i0:,i1:]
      elif display == 'right' :
          axis = [axis[0][i0:], axis[1][:]]
          data = out[i0:,:]
      elif display == 'up' :
          axis = [axis[0][:], axis[1][i1:]]
          data = out[:,i1:]
      elif display == 'full' :
          axis = [axis[0][:], axis[1][:]]
          data = out[:,:]
      else :
          raise ValueError( "display has to be 'upright', 'right', 'up' or 'full'" )

      return [axis, data]

