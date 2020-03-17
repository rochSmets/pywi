# this module codes the test particle run base

import os
import numpy as np
import scipy.ndimage as ndimage

# now for linking C to python
from numpy.ctypeslib import ndpointer
import ctypes



# for plotting routines
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D






class TPRun:
    """
    Test Particle run
    """

    #==========================================================
    #==========================================================
    def __init__(self,
                 npart,
                 charge,
                 mass,
                 run,
                 tstart,
                 tend,
                 t0,
                 r0,
                 dr0,
                 v0,
                 dv0,
                 dt,
                 loading     = '',     #randu,randn,copy
                 fieldinterp = False):

        """ constructor of the testparticle run object

        @param npart       : # of particles
        @param charge      : charge of the particles
        @param mass        : mass of the particles
        @param run         : run from which the fiels will be read
        @param tstart      : starting time of the integration
        @param tend        : end time of the integration
        @param t0          : time at which positions/velocities are given
        @param r0          : initial position
        @param dr0         : initial spatial deviation
        @param v0          : initial velocity
        @param dv0         : initial velocity deviation
        @param dt          : time step
        @param loading     : either 'randu', 'randn', or 'copy'
        @param fieldinterp : delfault(False) inteprolate fields in time or not

        if loading == 'randu' particles will be loaded randomly in a rectangle
        of size dr0
        if 'randn' is chosen, they will be loaded in a gaussian of spatial std=dr0
        if 'copy' is chosen, r0,v0 have to eb arrays of size (3,npart), dr0 and dv0
        are then disregarded


        @return: a TPRun object

        Creation : 2013-05-01 11:17:42.987369

        """
        self._npart     = npart
        self._charge    = charge
        self._mass      = mass
        self._r0        = None
        self._v0        = None

        # keep the user values for plotting etc.

        #TODO attention si loading==user
        # r0,v0 sont des tableaux et dr0,dv0 
        # ne doivent pas etre utilises !
        #- 2013-05-20 08:06:57.972466
        self._r0u       = r0
        self._v0u       = v0
        self._dr0u      = dr0
        self._dv0u      = dv0

        self.dt         = dt
        self.tstart     = tstart
        self.tend       = tend

        self.t0         = t0
        self._nt        = int((tend-tstart)/dt) + 1


        # checks the time interval is well defined 
        #----------------------------------------------------------------------
        if t0 < self.tstart or t0 > self.tend:
            print 'time (%5.3f) should be between tstart(%5.3f)\
                   and tend(%5.3f) (both included)' \
                   % (t0 ,self.tstart,self.tend)

            return None
        #----------------------------------------------------------------------
        # that is the selection time index.
        self._it0 = int((self.t0   - self.tstart)/self.dt) + 1


        # loading method
        self._loading = loading
        if loading.lower() == 'randu':
            self.load_randu(r0,dr0,v0,dv0)

        elif loading.lower() == 'randn':
            self.load_randn(r0,dr0,v0,dv0)

        elif loading.lower() == 'user':
            self._r0 = np.zeros((3,self._npart))
            self._r0[0,:] = r0[0,:]
            self._r0[1,:] = r0[1,:]
            self._v0 = np.zeros((3,self._npart))
            self._v0 = v0

        else:
            print 'Ptest : warning, no loading method specified'



        # position and velocity arrays

        self.r   = np.zeros((3, self._npart, self._nt),
                                    order='FORTRAN',
                                    dtype=np.float64)

        self.v   = np.zeros((3, self._npart, self._nt),
                                    order='FORTRAN',
                                    dtype=np.float64)


        # electric and magnetic field seen by each particle
        self.Ep   = np.zeros((3, self._npart, self._nt),
                             order = 'FORTRAN',
                             dtype=np.float64)

        self.Bp   = np.zeros((3, self._npart, self._nt),
                             order = 'FORTRAN',
                             dtype = np.float64)


        # do we interpolate fields in time ?
        self._fldinterp = fieldinterp

        # the PIC run from which he get the fields
        self._run       = run


        # if we interpolate fields in time
        # we actually need to read all the files between the time
        # interval
        if self._fldinterp == True:
            files    = run.GetFieldFiles(tstart,tend)    # here are the files
            nfiles   = len(files)

            self._E  = np.zeros((3,run.ncells[0]+1,run.ncells[1]+2,nfiles),
                                order='FORTRAN',
                                dtype=np.float32)

            self._B  = np.zeros((3,run.ncells[0]+1,run.ncells[1]+2,nfiles),
                                order='FORTRAN',
                                dtype=np.float32)

            # read the files
            for it,file in enumerate(files):

                time = run.GetFileTime(file)
                print 'reading file %s at time %5.3f' \
                        % (os.path.basename(file),time)

                self._E[:,:,:,it] = run.GetE(tend,silent='yes')
                self._B[:,:,:,it] = run.GetB(tend,silent='yes')




        # if we dont interpolate fields in time just read the last time
        else:
                self._E = run.GetE(t0, silent='yes')
                self._B = run.GetB(t0, silent='yes')

                # smooth fields may help... B is fine, E is noisy !!
                #for c in range(3):
                #    self._E[c,:,:] = ndimage.gaussian_filter(self._E[c,:,:],
                #                                            sigma=6,
                #                                            order=0)

                # checks which component is the out of plane
                if self._run.outofplane == 1:
                    Ey = self._E[1,:,:].copy()
                    Ez = self._E[2,:,:].copy()
                    self._E[1,:,:] = np.copy(Ez)
                    self._E[2,:,:] = np.copy(-Ey)

                    By = self._B[1,:,:].copy()
                    Bz = self._B[2,:,:].copy()
                    self._B[1,:,:] = np.copy(Bz)
                    self._B[2,:,:] = np.copy(-By)



        # in case we want to debug, use analytic fields
        debug = True

        if debug == False:
            self._B[0,:,:]  = 0.#1. + np.random.randn(self._B.shape[1],self._B.shape[2])*1e-3
            self._B[1,:,:]  = 0.#np.random.randn(self._B.shape[1],self._B.shape[2])*1e-3
            self._B[2,:,:]  = 1.#np.random.randn(self._B.shape[1],self._B.shape[2])*1e-3
            self._E[0,:,:]  = 0.#np.random.randn(self._B.shape[1],self._B.shape[2])*1e-2
            self._E[1,:,:]  = 1.#np.random.randn(self._B.shape[1],self._B.shape[2])*1e-2
            self._E[2,:,:]  = 0.#np.random.randn(self._B.shape[1],self._B.shape[2])*1e-2



        # linking to the C library
        pathlibdir   = os.path.dirname(os.path.realpath(__file__))
        pathlib      = os.path.join(pathlibdir,'functions.so')
        self._lib = ctypes.cdll.LoadLibrary(pathlib)
    #==========================================================







    #==========================================================
    #==========================================================
    def _velinit(self, v, dv):
        """uniform random velocities

        Creation : 2013-04-28 17:12:42.773824

        """

        # create the array
        self._v0 = np.zeros((3, self._npart))

        # x, y, z components of the velocity
        c = 0

        # loop over vx0, vy0, vz0
        for v_i,dv_i in zip(v,dv):

            v1 = v_i - dv_i/2.
            v2 = v_i + dv_i/2.

            self._v0[c,:] = np.random.random(self._npart)
            self._v0[c,:] = self._v0[c,:] * (v2 - v1) + v1
            c += 1

    #==========================================================







    #==========================================================
    #==========================================================
    def load_randu(self,
                   r,
                   dr,
                   v,
                   dv):

        """ loads particles randomly in a rectangle

        The method loads the particles in a rectangle
        of side 'dr' entered around 'r'.
        Particles are loaded uniformly in that rectangle.
        Velocities are uniformly distributed within 'dv' from 'v'

        Carefull : this method will erase any initial position
        and velocity that may have already been initialized.

        @param r    : (x,y,z) center of the rectangle
        @param dr   : (dx,dy,dz) size of the rectangle
        @param v    : (Vx,Vy,Vz) mean initial velocity
        @param dv   : (dVx,dVy,dVz) initial velocity interval

        @return: @todo

        Exemple  :

        Creation : 2013-05-01 11:21:19.671182

        """
        print 'Loading particles in a spatially uniform random distribution...'

        self._r0  = np.zeros((3, self._npart))

        c = 0. # x, y and z components

        for r_i, dr_i in zip(r, dr):

            r1            = r_i - dr_i/2.
            r2            = r_i + dr_i/2.
            self._r0[c,:] = np.random.random(self._npart)

            self._r0[c,:]  = self._r0[c,:]* (r2 - r1) + r1

            c += 1

        self._velinit(v,dv)

    #==========================================================






    #==========================================================
    #==========================================================
    def load_randn(self,
                   r,
                   dr,
                   v,
                   dv):

        """ loads particles randomly in a spatial gaussian

        The method loads the particles in a spatial gaussian
        of std 'dr' entered around 'r'.
        Particles are loaded in a normal distribution.
        Velocities are uniformly distributed within 'dv' from 'v'

        Carefull : this method will erase any initial position
        and velocity that may have already been initialized.

        @param r  : (x,y,z) center of the gaussian
        @param dr : (dx,dy,dz) standard deviation from r
        @param v  : (Vx,Vy,Vz) mean initial velocity
        @param dv : (dVx,dVy,dVz) initial velocity interval

        @return: @todo

        Exemple  :

        Creation : 2013-05-01 11:21:19.671182

        """
        print 'Loading particles in a spatially gaussian random distribution...'

        self._r0  = np.zeros((3, self._npart))

        c = 0. # x, y and z components

        for r_i, dr_i in zip(r, dr):

            r1            = r_i - dr_i/2.
            r2            = r_i + dr_i/2.
            self._r0[c,:] = np.random.random(self._npart)

            self._r0[c,:]  = self._r0[c,:]* (r2 - r1) + r1

            c += 1

        self._velinit(v,dv)


    #==========================================================




    #==========================================================
    #==========================================================
    def _validposition(self, pos):
        """checks whether the position is inside the box or not

        @return: True or False

        Exemple  : 

        Creation : 2013-05-01 14:26:02.007530

        """
        return True
    #==========================================================




    #==========================================================
    #==========================================================
    def move(self):
        """moves all the particles

        @return: @todo

        Exemple  : 

        Creation : 2013-05-01 14:28:08.724842

        """

        xr = self._run.GetCoord(axis=0)
        yr = self._run.GetCoord(axis=1)

        xmin = xr[0]
        ymin = yr[0]

        # setup the initial condition

        print ' it0  :', self._it0

        self.r[:,:,self._it0] = self._r0
        self.v[:,:,self._it0] = self._v0

        # call super-fast cython functions now haha

        # the following function will move all the particles
        # for all time steps

        self._moveall()
        self._pfields()


    #==========================================================




    #==========================================================
    #==========================================================
    def _moveall(self):
        """move all particles backward and forward as necessary

        Exemple  : 

        Creation : 2013-05-04 15:40:32.257893

        """
        func          = self._lib.moveall
        func.restype  = None
        func.argtypes = [ndpointer(ctypes.c_double),    # pos
                         ndpointer(ctypes.c_double),    # vel
                         ndpointer(ctypes.c_float),     # E
                         ndpointer(ctypes.c_float),     # B
                         ctypes.c_uint,                 # nx
                         ctypes.c_uint,                 # ny
                         ctypes.c_double,               # dx
                         ctypes.c_double,               # dy
                         ctypes.c_double,               # xmin
                         ctypes.c_double,               # ymin
                         ctypes.c_double,               # dt
                         ctypes.c_double,               # charge
                         ctypes.c_double,               # mass
                         ctypes.c_uint,                 # nt
                         ctypes.c_uint,                 # it0
                         ctypes.c_uint]                 # npart


        xc = self._run.GetCoord(axis=0)
        yc = self._run.GetCoord(axis=1)

        func(self.r,
             self.v,
             self._E,
             self._B,
             self._B[0,:,:].shape[0],
             self._B[0,:,:].shape[1],
             self._run.dl[0],
             self._run.dl[1],
             xc[0], yc[0],
             self.dt,
             self._charge,
             self._mass,
             self._nt,
             self._it0,
             self._npart)
    #==========================================================




    #==========================================================
    #==========================================================
    def _pfields(self):
        """interpolates the fields at the particle positions

        Exemple  :

        Creation : 2013-05-04 15:40:32.257893

        """
        func          = self._lib.pfields
        func.restype  = None
        func.argtypes = [ndpointer(ctypes.c_double),    # pos
                         ndpointer(ctypes.c_float),     # E
                         ndpointer(ctypes.c_float),     # B
                         ctypes.c_uint,                 # nx
                         ctypes.c_uint,                 # ny
                         ctypes.c_uint,                 # npart
                         ctypes.c_uint,                 # nt
                         ctypes.c_double,               # xmin
                         ctypes.c_double,               # ymin
                         ctypes.c_double,               # dx
                         ctypes.c_double,               # dy
                         ndpointer(ctypes.c_double),    # Ep
                         ndpointer(ctypes.c_double)]    # Bp


        xc = self._run.GetCoord(axis=0)
        yc = self._run.GetCoord(axis=1)

        func(self.r,
             self._E,
             self._B,
             self._E[0,:,:].shape[0],
             self._E[0,:,:].shape[1],
             self._npart,
             self._nt,
             xc[0],yc[0],
             self._run.dl[0],
             self._run.dl[1],
             self.Ep,
             self.Bp)

    #==========================================================





    #==========================================================
    #==========================================================
    def vnorm(self):
        """calculates the norm of the velocity of each particle

        @return: @todo

        Exemple  : 

        Creation : 2013-05-05 16:58:58.617352

        """
        v2 = self.v[0,:,:]**2 + self.v[1,:,:]**2 + self.v[2,:,:]**2
        return np.sqrt(v2)
    #==========================================================



    #==========================================================
    #==========================================================
    def kineticvariation(self):
        """calculates and return the variation of kinetic energy of each particle
        normalized at their "initial" value

        @return: @todo

        Exemple  :

        Creation : 2013-05-06 11:18:50.701666

        """

        vnorm   = self.vnorm() # get the norms
        kin     = vnorm**2/2.
        normkin = np.empty((self._npart,self._nt))

        for p in range(self._npart):
            normkin[p,:] = (kin[p,:]-kin[p,0])/kin[p,0]

        return normkin

    #==========================================================





    #==========================================================
    #==========================================================
    def kineticnormalized(self):
        """calculates and return the kinetic energy of each particle
        normalized at their "initial" value

        @return: @todo

        Exemple  :

        Creation : 2013-05-06 11:18:50.701666

        """

        vnorm   = self.vnorm() # get the norms
        kin     = vnorm**2/2.
        normkin = np.empty((self._npart,self._nt))

        for p in range(self._npart):
            normkin[p,:] = kin[p,:]/kin[p,0]

        return normkin

    #==========================================================





    #==========================================================
    #==========================================================
    def magforce(self):
        """calculates and return the magnetic force

        @return: @todo

        Exemple  : 

        Creation : 2013-05-05 17:01:57.513006

        """
        fx = self._charge * (self.v[1,:,:]*self.Bp[2,:,:] - \
                             self.v[2,:,:]*self.Bp[1,:,:])

        fy = self._charge * (self.v[2,:,:]*self.Bp[1,:,:] - \
                             self.v[1,:,:]*self.Bp[2,:,:])

        fz = self._charge * (self.v[0,:,:]*self.Bp[1,:,:] - \
                             self.v[1,:,:]*self.Bp[0,:,:])

        return (fx,fy,fz)

    #==========================================================




    #==========================================================
    #==========================================================
    def elecforce(self):
        """calculate and return the electric force

        @return: @todo

        Exemple  :

        Creation : 2013-05-05 17:04:49.261054

        """

        fx = self._charge * self.Ep[0,:,:]
        fy = self._charge * self.Ep[1,:,:]
        fz = self._charge * self.Ep[2,:,:]

        return (fx, fy, fz)
    #==========================================================




    #==========================================================
    #==========================================================
    def gettime(self):
        """returns the time

        @return: @todo

        Exemple  : 

        Creation : 2013-05-06 11:51:15.008025

        """
        return np.arange(self._nt)*self.dt +self.tstart
    #==========================================================






# ------- PLOTTING ROUTINES ----------


    #==========================================================
    #==========================================================
    def traj3D(self):
        """ plot trajectories in 3D

        @return: @todo

        Exemple  : 

        Creation : 2013-05-13 15:14:52.070841

        """
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        for p in range(self._npart):
            ax.plot(self.r[0,p,:], self.r[1,p,:],self.r[2,p,:])

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        fig.savefig('traj3D.png',dpi=300)
    #==========================================================







