

import numpy as np


# class interface :
# ------------------

#   - modelname     - gives the name of the model (pic, hybrid, etc.)
#   - display()     - display some informations about the run
#   - GetCoords()   - returns coordinate along a direction with default extent
#   - GetNbDiag()   - returns the number of dumps for fields/part. data
#   - GetMass()     - mass of 'species'
#   - GetCharge()   - charge of 'species'
#   - GetRunTime()  - returns total time of the run (even not finished)
#   - GetTimeStep() - dump step for fields, particles
#   - whereAmI()    - where is the run stored

#   - indices2coord()
#       - based on indices of the quantity, return the coordinates
#   - coord2indices()


#   - GetVe()
#   - GetVi()
#   - GetV(species)
#   - GetB(, original_grid=False)
#   - GetE(, "")
#   - GetJ()
#   - GetNe()
#   - GetNi()
#   - GetN(species)
#   - GetFlux(warning if 3D)
#   - GetP(species)
#   - GetPxx(species)
#   - GetPxy(species)
#   - GetPxz(species)
#   - GetPyy(species)
#   - GetPyz(species)
#   - GetPzz(species)
#   - GetThermalEnergy(species) (utherm)
#   - VgradV(species)
#   - divP()
#   - GetParticles(time, species, location=None)


# model dependant
# ------------------
#   - GetDebye()


# should be private :
# ------------------

#   - iter2time()
#   - time2iter()
#       - all should be accessed by time

#   - GetField(fieldname..., time)
#   - GetFieldData()
#   - divT(T)



# obsoletes ?

#   - GetFIeldFiles()
#   - GetFileTime()
#       - not valid anymore for 1 file HDF5 groups
#   - GetNbFieldFiles()

#   - GetOutOfPlaneDir() ?
#       - should always put the ignorable coordinate in Z (3rd coord)




class Run (object):

    """
    The :class:`Run` is the base class for all run objects implemented
    in pywi. It is not intended to be instanciated but defines the interface
    that all the run objects must implement. Some routines are already
    implemented in :class:`Run` (such as GetVp()) since they will be
    the same for all runs.


    In this class, should be methods that concern only basic/raw quantities
    that are available in data files or specific run properties. For instance,
    a routine returning the pressure tensor is typically one that will be
    implemented here, while one returning the reconnection rate or the degree
    of nongyrotropy will not, those are higher-level products that can be
    calculated the same way for all runs based on basics quantities

    """



    #----------------------------------------------------------
    #----------------------------------------------------------
    def __init__(self, path, runID = 'r'):
        """Exemple  :

        Creation : 2015-02-27 16:04:58.868849

        """
        self.modelname  = 'Unknown'
        self.codename   = 'Unknown code'
        self.path       = path
        self.runID      = runID
    #==========================================================






    #----------------------------------------------------------
    #----------------------------------------------------------
    def display(self):
        """
        Shows some information on screen about the run, like its path, 
        its ID, and the name of the physical model.
        """
        print('Model name       : {0}'.format(self.modelname))
        print('Run ID           : {0}'.format(self.runID))
        print('Path             : {0}'.format(self.path))
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def whereAmI(self):
        """ Exemple  : path = run.whereAmI()
        """
        return self.path
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetVe(self, time):
        """
        Returns the electron Bulk velocity at desired time

        Parameters:

            :time: is the time in inverse ion cyclotron frequency

        :Returns: numpy.ndarray, axis -1 is the component of the vector

        Exemple:

           >>> Ve = run.GetVe(0.1)

        """
        return self.GetV(time, species='electrons')
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetVi(self, time):
        """
        Returns the Ion Bulk velocity at desired time. Ions means
        that it is the average of ALL ion species

        Parameters:

            :time: is the time in inverse ion cyclotron frequency

        :Returns: numpy.ndarray, axis -1 is the component of the vector

        Exemple:

           >>> Ve = run.GetVe(0.1)

        """
        return self.GetV(time, species='ions')
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetVp(self, time):
        """
        Returns the proton Bulk velocity

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        return self.GetV(time, species='protons')
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetV(self, time, species):
        """
        Returns the Bulk velocity of a specific species

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetB(self, time, original_grid = False):
        """
        Returns the magnetic field at a specific time
        original_grid = True returns the field on the
        grid that is used in the simulation. This is not
        necessarily implemented for all models

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================




    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetE(self, time, original_grid = False):
        """
        Returns the electric field at a specific time
        original_grid = True returns the field on the
        grid that is used in the simulation. This is not
        necessarily implemented for all models

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================




    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetJ(self, time):
        """
        Returns the electric current density vector

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================






    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetNe(self, time):
        """
        Returns the electron particle density

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        return self.GetN(time, species='electrons')
    #==========================================================






    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetNi(self, time):
        """
        Returns the ion particle density, this include all
        ion species.

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        return self.GetN(time, species='ions')
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetNp(self, time):
        """
        Returns the proton particle density
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938

        """
        return self.GetN(time, species='protons')
    #==========================================================






    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetN(self, time, species):
        """
        Returns the particle density of a specific species

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetFLux(self, time):
        """
        Returns the magnetic flux function
        this raises an exception if the run is not 2D

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================






    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetP(self, time, species):
        """
        Returns the pressure tensor of a specific species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """

        Pxx = self.GetPxx(time, species)
        Pxy = self.GetPxy(time, species)
        Pxz = self.GetPxz(time, species)
        Pyy = self.GetPyy(time, species)
        Pyz = self.GetPyz(time, species)
        Pzz = self.GetPzz(time, species)

        P  = np.ndarray(dtype=self.dtype, shape=self.tfield_shape)
        P[...,0,0] = Pxx
        P[...,1,1] = Pyy
        P[...,2,2] = Pzz

        P[...,0,1] = Pxy
        P[...,0,2] = Pxz
        P[...,1,2] = Pyz

        P[...,1,0] = Pxy
        P[...,2,0] = Pxz
        P[...,2,1] = Pyz

        return P

    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPe(self, time):
        """
        Returns the electron pressure tensor
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """
        return self.GetP(time, species='electrons')
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPp(self, time):
        """
        Returns the proton pressure tensor
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """
        return self.GetP(time, species='protons')
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPxx(self, time, species):
        """
        Returns the component of the pressure tensor of a species

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPxy(self, time, species):
        """
        Returns the component of the pressure tensor of a species

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================








    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPxz(self, time, species):
        """
        Returns the component of the pressure tensor of a species

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPyy(self, time, species):
        """
        Returns the component of the pressure tensor of a species

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPyz(self, time, species):
        """
        Returns the component of the pressure tensor of a species

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPzz(self, time, species):
        """
        Returns the component of the pressure tensor of a species

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================






    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetThermalEnergy(self, time, species):
        """
        Returns the thermal energy for a specific species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """
        P = self.GetP(time, species=species)
        return 0.5 * (P[0,0,:,:,:] + P[1,1,:,:,:] + P[2,2,:,:,:])
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def VgradV(self, time, species):
        """
        Returns the (V dot Grad )(v)

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def divP(self, time, species):
        """
        Returns the divergence of the pressure tensor of a
        specific species

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def __divT(self, time, T):
        """
        Returns the divergence of the tensor T

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )
    #==========================================================





    def GetParticles(self, time, species, location=None):
        """
        Returns a group of particles from a particular species
        at some location for some time
        If no location is not specified, it returns the whole
        box

        Exemple  :

        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )







    def GetMass(self, species):
        """
        Returns the mass of particles of a particular species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )





    def GetCharge(self, species):
        """
        Returns the electric charge of particles of a particular species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938

        """
        raise NotImplementedError( "Should have implemented this" )




    def getA(self, time, comp, method) :

        if (self.ndim == 2) & (comp == 'z') :
            if method == 'fft' :
                return self.getAz_fft(time)
            elif method == 'fftw' :
                return self.getAz_fftw(time)
            elif method == 'fdm':
                return self.getAz_fdm(time)
            elif method == 'fem':
                return self.getAz_fem(time)
            else :
                raise ValueError("'method' should be in ['fft', 'fftw', 'fem', 'fdm']")
        else:
            raise ValueError('this method only works for 2d runs to calculate the z component of the vector potential')



    def getAz_fft(self, time) :

        from numpy.fft import fft2
        from numpy.fft import ifft2

        def k_(axis, sup, dim):
            N = dim[axis]
            k = np.linspace(0, N-1, N)
            kplus  =  k[:int(N/2)+1]
            kminus = -k[int((N-1)/2):0:-1]
            return np.concatenate((kplus, kminus))*2*np.pi/sup[axis]

        dim = self.ncells+1
        sup = self.domsize

        x = np.linspace(-sup[0], +sup[0], dim[0])
        y = np.linspace(-sup[1], +sup[1], dim[1])

        xv, yv = np.meshgrid(x, y, indexing = 'ij')

        bx = self.getB(time, 'x')+1j*np.zeros(yv.shape)
        by = self.getB(time, 'y')+1j*np.zeros(xv.shape)

        kx = k_(0, sup, dim)
        ky = k_(1, sup, dim)

        wy, wx = np.meshgrid(ky, kx)

        k2 = np.square(wx)+np.square(wy)
        k2[0][0] = 1.0

        BX = fft2(bx)
        BY = fft2(by)

        AZ = np.divide(1j*(wx*BY-wy*BX), k2)
        AZ[0][0] = 0.0

        return ifft2(AZ).real


    def getAz_fftw(self, time) :

       import pyfftw
       from pylab import meshgrid

        # TODO : not very nice set a fixed value... could be a env var ?
       numOfThreads = 8

       dim = self.ncells+1

       # arrays have to be aligned... for pyfftw to work
       bx = pyfftw.empty_aligned(dim, dtype = 'complex128')
       by = pyfftw.empty_aligned(dim, dtype = 'complex128')
       BX = pyfftw.empty_aligned(dim, dtype = 'complex128')
       BY = pyfftw.empty_aligned(dim, dtype = 'complex128')
       AZ = pyfftw.empty_aligned(dim, dtype = 'complex128')
       az = pyfftw.empty_aligned(dim, dtype = 'complex128')

       # 2d fftw plans
       fftForx = pyfftw.FFTW(bx, BX, axes=(0, 1), direction='FFTW_FORWARD', flags=('FFTW_MEASURE', ))
       fftFory = pyfftw.FFTW(by, BY, axes=(0, 1), direction='FFTW_FORWARD', flags=('FFTW_MEASURE', ))
       fftBack = pyfftw.FFTW(AZ, az, axes=(0, 1), direction='FFTW_BACKWARD', flags=('FFTW_MEASURE', ))

       # fill the input arrays
       bx.real[...] = self.GetB(time)[..., 0]
       by.real[...] = self.GetB(time)[..., 1]
       bx.imag = np.zeros(dim, dtype = float)
       by.imag = np.zeros(dim, dtype = float)

       # update needed before execute plan
       fftForx.update_arrays(bx, BX)
       fftFory.update_arrays(by, BY)

       pyfftw.FFTW.execute(fftForx)
       pyfftw.FFTW.execute(fftFory)

       # this intend to build the k values, first half is ascending
       # positives values, last one being the descending negative values
       axis = 0
       N = dim[axis]
       k = np.linspace(0, N-1, N)
       kplus = k[:int(N/2)+1]
       kminus =-k[int((N-1)/2):0:-1]
       kx = np.concatenate((kplus, kminus))*2*np.pi/self.domsize[axis]

       axis = 1
       N = dim[axis]
       k = np.linspace(0, N-1, N)
       kplus = k[:int(N/2)+1]
       kminus =-k[int((N-1)/2):0:-1]
       ky = np.concatenate((kplus, kminus))*2*np.pi/self.domsize[axis]

       # create 2d values to avoid nested "for" loops
       wy, wx = meshgrid(ky, kx)

       k2 = np.square(wx)+np.square(wy)
       k2[0][0] = 1.0

       AZ = np.divide(1j*(wx*BY-wy*BX), k2)

       # remove dc component while k=0 mode is not defined
       AZ[0][0] = 0.0

       fftBack.update_arrays(AZ, az)

       pyfftw.FFTW.execute(fftBack)

       _flux = az.real/(self.ncells[0]*self.ncells[1])


       return _flux



    def getAz_fdm(self, time) :

        return self.GetFlux(time)



    def GetFlux(self, time):

        #import matplotlib.pyplot as plt

        B = self.GetB(time)

        if len(B.shape) == 4:
            raise ValueError("magnetic Flux() only defined for 2D")

        shape = (B.shape[0], B.shape[1])
        flux  = np.zeros(shape=shape, dtype=B.dtype)

        dl1 = self.dl[0]
        dl2 = self.dl[1]

        n1 = shape[0]
        n2 = shape[1]

        b1 = B[...,0]
        b2 = B[...,1]


        flux[1:,0] = flux[0,0] + np.cumsum(b2[:-1,0]*dl1)
        flux[0,1:] = flux[0,0] - np.cumsum(b1[0,:-1]*dl2)


        flux = self._fast(flux, b1, b2, dl1, dl2)

        return flux



    # SEBERG (#scipy) method to calculate the 2D integral flux from B
    # http://stackoverflow.com/questions/11854522/nested-loop-with-array-indexing-in-numpy
    def _fast(self, flux, b1, b2, dl1=1., dl2=1.) :

        from scipy.ndimage import convolve

        _flux = np.zeros((flux.shape[0]+1, flux.shape[1]+1), dtype=flux.dtype)
        temp_b1 = np.zeros((b1.shape[0]+1, b1.shape[1]+1), dtype=b1.dtype)
        temp_b2 = np.zeros((b2.shape[0]+1, b2.shape[1]+1), dtype=b2.dtype)

        _flux[:-1,:-1] = flux
        convolve(_flux[:-1,:-1], [[0, 0.5], [0.5, 0]], _flux[1:,1:])

        temp_b2[1:,1:-1] = b2[:,1:]*dl1
        temp_b1[1:-1,1:] = b1[1:,:]*dl2

        conv_b = np.array([[0.0, 0.5], [0.5, 0.5]])
        convolve(temp_b2[:-1,:-1], [[0.5, 0.5], [0.5, 0.]], temp_b2[1:,1:])
        convolve(temp_b1[:-1,:-1], [[-0.5, 0.5], [0.5, 0.]], temp_b1[1:,1:])

        _flux += temp_b2
        _flux += temp_b1

        return _flux[:-1,:-1]



    def getAz_fem(self, time):

        from fipy import Grid2D, CellVariable, DiffusionTerm

        nx = self.ncells[0]+1
        ny = self.ncells[1]+1

        mesh = Grid2D(dx = self.dl[0], dy = self.dl[1], nx = nx, ny = ny)

        az_ = CellVariable(mesh = mesh, value = 0.0)

        bx_ = CellVariable(mesh = mesh)
        bx_.value = np.reshape(self.getB(time, 'x'), nx*ny, order='F')

        by_ = CellVariable(mesh = mesh)
        by_.value = np.reshape(self.getB(time, 'y'), nx*ny, order='F')

        jz_ = CellVariable(mesh = mesh)
        jz_.value = np.reshape(self.getJ(time, 'z'), nx*ny, order='F')

        az_.equation = (DiffusionTerm(coeff = 1.0) + jz_ == 0)

        # beware of the sign of the flux : always consider outward direction
        az_.faceGrad.constrain(+by_.arithmeticFaceValue * mesh._orientedFaceNormals, where = mesh.facesLeft)
        az_.faceGrad.constrain(-by_.arithmeticFaceValue * mesh._orientedFaceNormals, where = mesh.facesRight)
        az_.faceGrad.constrain(-bx_.arithmeticFaceValue * mesh._orientedFaceNormals, where = mesh.facesBottom)
        az_.faceGrad.constrain(+bx_.arithmeticFaceValue * mesh._orientedFaceNormals, where = mesh.facesTop)

        az_.equation.solve(var = az_)

        return np.reshape(az_.value, (nx, ny), order='F') - np.mean(az_)
