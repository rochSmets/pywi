"""
.. module:: heckle
    :synopsis: class for data produced by the code Heckle
.. moduleauthor:: Nicolas Aunai <nicolas.aunai@lpp.polytechnique.fr>
"""




#import pywi.runs.run as run
import run
import h5py
import os
import numpy as np

if __name__ == '__main__':
    main()


class Heckle(run.Run):

    """
    test
    """

    #----------------------------------------------------------
    #----------------------------------------------------------
    def __init__(self, path, runID = 'r', fieldname = 'fields.h5'):

        """
        Exemple  :

        Creation : 2015-02-27 16:49:57.016701

        """
        super(Heckle, self).__init__(path, runID=runID)

        # ----- heckle private members
        self._fieldfilename   = 'fields.h5'
        self._speciefilename  = 'species.h5'
        self._timefilename    = 'time.h5'
        self._restartfilename = 'restarts.h5'
        self._electron = ['electrons', 'electron', 'e']
        self._proton   = ['protons', 'proton', 'p']
        self._alpha    = ['alphas', 'alpha', 'a']
        self._ion      = ['ions', 'ion', 'i']
        # we often have a second proton population : could we call it "alphas" ?
        # (even if backgrounds are generally not alpha particles)

        # ----- heckle public members
        self.modelname  = 'Hybrid'
        self.codename   = 'Heckle'
        self.ncells     = np.array(self._getFieldsAttrib('nbrOfCells'))
        #self.ncells     = np.array(self._getFieldsAttrib('numofcells'))
        self.domsize    = np.array(self._getFieldsAttrib('domainSize'))
        #self.domsize    = np.array(self._getFieldsAttrib('domainsize'))
        self.dl         = np.array(self._getFieldsAttrib('meshSize'))
        #self.dl         = np.array(self._getFieldsAttrib('meshsize'))
        self.hprsty     = self._getFieldsAttrib('hyperResistivity')
        #self.hprsty     = self._getFieldsAttrib('hyperresistivity')
        self.rsty       = self._getFieldsAttrib('resistivity')
        self.dumpfields = self._getFieldsAttrib('dumpFields')
        #self.dumpfields = self._getFieldsAttrib('dumpfields')

        # 1D case
        if self.ncells[1] == 1 :
            self.ncells  = np.array([self.ncells[0]])
            self.domsize = np.array([self.domsize[0]])
            self.dl      = np.array([self.dl[0]])

        # 2D case
        elif self.ncells[2] == 1 :
            self.ncells  = self.ncells[:-1]
            self.domsize = self.domsize[:-1]
            self.dl      = self.dl[:-1]

        self.ndim         = self.ncells.size
        self.sfield_shape = tuple(self.ncells+1)
        self.vfield_shape = tuple(self.ncells+1) + (3,)
        self.tfield_shape = tuple(self.ncells+1) + (3,3,)
        self.dtype        = np.float32 # what is the meaning ? only for fields.h5 ?
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def _getFieldsAttrib(self, attributeName):

        """
        Creation : 2015-02-27 19:12:14.062301
        """

        f = h5py.File(os.path.join(self.path, self._fieldfilename))
        attribute = f.attrs[attributeName]
        f.close()

        return attribute

    #==========================================================





    def getTimeGroups(self, time=None) :

        """
        get time (a single value or a list with format [tmin, tmax])
        and return a list of strings, containing all the time groupnames
        associated in "fields.h5"

        return a list with the times, and a list with the associated groups
        """

        f = h5py.File(os.path.join(self.path, self._fieldfilename))

        # build the list of all time recorded in fields.h5
        groups = f.keys()
        timesfromfile = np.empty([0])

        # first parsed string is "time", second ":" & third the time we want
        for grp in groups :
           timesfromfile = np.append(timesfromfile , float(grp.strip().split()[2]))

        if time == None:
            timesfromfile.sort()
            time=[timesfromfile[0], timesfromfile[-1]]

        if time.__len__() == 1 :
            # good times, bad times, you know i've had my share
            goodtimes = timesfromfile[np.argmin(np.fabs(np.array(time)-timesfromfile))]

            timegroups = 'time : %f' % (goodtimes)

        elif time.__len__() == 2 :
            goodtimes = timesfromfile[(timesfromfile >= time[0]) & (timesfromfile <= time[1])]

            timegroups = ['time : %f' % (tim) for tim in goodtimes]

        else :
            raise ValueError( "time has to be a list of 1 or 2 float values")

        return goodtimes, timegroups


    #----------------------------------------------------------
    #----------------------------------------------------------
    def _readFieldsDataSet(self, time, fieldname):

        """
        Creation : 2015-02-27 19:12:14.062301

           all time groups are under the "root dir" : "/"
           all dataset are in a given time group
        """

        mytimes, mygroups = self.getTimeGroups([time])

        f = h5py.File(os.path.join(self.path, self._fieldfilename))

        data = f[mygroups+'/'+fieldname][()] #data = f[mygroups+'/'+fieldname].value is deprecated

        f.close()

        # run 1D
        if self.ndim == 1 :
            mydata = data[:,0,0]

        # run 2D
        elif self.ndim == 2 :
            mydata = data[:,:,0]

        # run 3D
        elif self.ndim == 3 :
            mydata = data

        else :
            raise ValueError("what the fuck is this dim value : '%d' ?" % self.ndim)

        return mydata
    #==========================================================


    def getTimeDataSet(self, datasetname):

        f = h5py.File(os.path.join(self.path, self._timefilename))

        data = f[datasetname][()] #data = f[datasetname].value is deprecated

        f.close()

        #Bulk Energy sp = 0
        #Bulk Energy sp = 1
        #Mag. energy
        #Parallel Energy sp = 0
        #Parallel Energy sp = 1
        #Perp Energy sp = 0
        #Perp Energy sp = 1
        #Time
        #X Elec. Fluctuations
        #Y Elec. Fluctuations
        #Z Elec. Fluctuations
        #divB
        #pseudo divB

        return data



    #----------------------------------------------------------
    #----------------------------------------------------------
    def _fillVec(self, v1: object, v2: object, v3: object) -> object:

        shape =  v1.shape + (3,)
        V     = np.zeros(shape = shape, dtype = v1.dtype)
        V[...,0] = v1
        V[...,1] = v2
        V[...,2] = v3

        return V

    # private method _fillvec could be changed by _fillall
    # this one is not restricted to vector : can fill tensor of whatever rank

    def _fillall(self, mylist) :

        MyListLen = mylist.__len__()
        shape = mylist[0].shape+(MyListLen,)
        Tab = np.zeros(shape = shape, dtype = mylist[0].dtype)

        for i in range(MyListLen) :
            Tab[...,i] = mylist[i]

        return Tab
    #==========================================================




    def getB(self, time, component):

        if component == 'x' or component == 'X' or component == 0 :
            return self._readFieldsDataSet(time, 'Bx')
        elif component == 'y' or component == 'Y' or component == 1 :
            return self._readFieldsDataSet(time, 'By')
        elif component == 'z' or component == 'Z' or component == 2 :
            return self._readFieldsDataSet(time, 'Bz')
        # TODO : could component be 'all', so the method returns all components in a ndarray ?
        else :
            raise ValueError("Unknown component : should be ['x', 'y', 'z'], ['X', 'Y', 'Z'] or [0, 1, 2]")



    def getJ(self, time, component):

        if component == 'x' or component == 'X' or component == 0 :
            return self._readFieldsDataSet(time, 'Jx')
        elif component == 'y' or component == 'Y' or component == 1 :
            return self._readFieldsDataSet(time, 'Jy')
        elif component == 'z' or component == 'Z' or component == 2 :
            return self._readFieldsDataSet(time, 'Jz')
        else :
            raise ValueError("Unknown component : should be ['x', 'y', 'z'], ['X', 'Y', 'Z'] or [0, 1, 2]")



    def getE(self, time, component):

        if component == 'x' or component == 'X' or component == 0 :
            return self._readFieldsDataSet(time, 'Ex')
        elif component == 'y' or component == 'Y' or component == 1 :
            return self._readFieldsDataSet(time, 'Ey')
        elif component == 'z' or component == 'Z' or component == 2 :
            return self._readFieldsDataSet(time, 'Ez')
        else :
            raise ValueError("Unknown component : should be ['x', 'y', 'z'], ['X', 'Y', 'Z'] or [0, 1, 2]")





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetB(self, time, origin_grid = False):

        """
        Returns the magnetic field vector at desired time.

        Args:
            :param time: the time at which B is returned

        Kwargs:
            :param origin_grid (bool): [default:False] If True the\
                    magnetic components are defined on the Yee Grid.

        Returns:
            :returns B: numpy.ndarray of shape (nx,3), (nx,ny,3),(nx,ny,nz,3) \
                     for 1D, 2D and 3D runs, respectively
        """

        Bx = self._readFieldsDataSet(time, 'Bx')
        By = self._readFieldsDataSet(time, 'By')
        Bz = self._readFieldsDataSet(time, 'Bz')

        return self._fillVec(Bx, By, Bz)
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetJ(self, time):

        """
        Creation : 2015-02-27 19:12:14.062301
        """

        Jx = self._readFieldsDataSet(time, 'Jx')
        Jy = self._readFieldsDataSet(time, 'Jy')
        Jz = self._readFieldsDataSet(time, 'Jz')

        return self._fillVec(Jx, Jy, Jz)
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetE(self, time, origin_grid = False):

        """
        Creation : 2015-02-27 19:12:14.062301
        """

        Ex = self._readFieldsDataSet(time, 'Ex')
        Ey = self._readFieldsDataSet(time, 'Ey')
        Ez = self._readFieldsDataSet(time, 'Ez')

        return self._fillVec(Ex, Ey, Ez)
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetCoords(self, extent = None, axis = 0):

        """
        Creation : 2015-02-27 19:12:14.062301

           returns an array with the coordinates of grid points along
           a given axis, for a given extent
        """

        dl = self.dl[axis]
        if extent is None:
            return dl * np.arange(self.ncells[axis]+1)
        else:
            npts = (extent[1]-extent[0])/dl + 1
            return extent[0] + dl*np.arange(npts)
    #==========================================================



    def getV(self, time, specie, component):

        if (type(specie) == str and specie.lower() in self._electron) or specie == 0:
            if component == 'x' or component == 'X' or component == 0 :
                return self._readFieldsDataSet(time, 'Vx[0]')
            elif component == 'y' or component == 'Y' or component == 1 :
                return self._readFieldsDataSet(time, 'Vy[0]')
            elif component == 'z' or component == 'Z' or component == 2 :
                return self._readFieldsDataSet(time, 'Vz[0]')
            # TODO : could component be 'all', so the method returns all components in a ndarray ?
            else :
                ValueError("Unknown component : should be ['x', 'y', 'z'], or [0, 1, 2]")

        elif (type(specie) == str and specie.lower() in self._proton) or specie == 1:
            if component == 'x' or component == 'X' or component == 0 :
                return self._readFieldsDataSet(time, 'Vx[1]')
            elif component == 'y' or component == 'Y' or component == 1 :
                return self._readFieldsDataSet(time, 'Vy[1]')
            elif component == 'z' or component == 'Z' or component == 2 :
                return self._readFieldsDataSet(time, 'Vz[1]')
            else :
                ValueError("Unknown component : should be ['x', 'y', 'z'], or [0, 1, 2]")

        elif (type(specie) == str and specie.lower() in self._alpha) or specie == 2:
            if component == 'x' or component == 'X' or component == 0 :
                return self._readFieldsDataSet(time, 'Vx[2]')
            elif component == 'y' or component == 'Y' or component == 1 :
                return self._readFieldsDataSet(time, 'Vy[2]')
            elif component == 'z' or component == 'Z' or component == 2 :
                return self._readFieldsDataSet(time, 'Vz[2]')
            else :
                ValueError("Unknown component : should be ['x', 'y', 'z'], or [0, 1, 2]")

        elif (type(specie) == str and specie.lower() in self._ion):
            if component == 'x' or component == 'X' or component == 0 :
                return self._readFieldsDataSet(time, 'Vix')
            elif component == 'y' or component == 'Y' or component == 1 :
                return self._readFieldsDataSet(time, 'Viy')
            elif component == 'z' or component == 'Z' or component == 2 :
                return self._readFieldsDataSet(time, 'Viz')
            else :
                ValueError("Unknown component : should be ['x', 'y', 'z'], or [0, 1, 2]")

        else :
            raise ValueError("Unknown component : should be ['e', 'p', 'a', 'i'], or [0, 1, 2]")



    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetV(self, time, species):

        """
        Returns the Bulk velocity of a specific species
        Exemple  : Ve = run.GetV(14.2, species='electrons')
        Creation : 2015-02-27 16:15:40.856938
        """

        if (type(species) == str and species.lower() in self._electron) or species == 0:
            Vxn = 'Vx[0]'
            Vyn = 'Vy[0]'
            Vzn = 'Vz[0]'

        elif species.lower() in self._proton:
            Vxn = 'Vx[1]'
            Vyn = 'Vy[1]'
            Vzn = 'Vz[1]'

        elif species.lower() in self._ion:
            Vxn = 'Vix'
            Vyn = 'Viy'
            Vzn = 'Viz'

        else:
            raise ValueError("Unknown species '%s'" % species)

        Vx = self._readFieldsDataSet(time, Vxn)
        Vy = self._readFieldsDataSet(time, Vyn)
        Vz = self._readFieldsDataSet(time, Vzn)

        return self._fillVec(Vx, Vy, Vz)
    #==========================================================


    def getN(self, time, specie):

        if (type(specie) == str and specie.lower() in self._electron) or specie == 0:
            return self._readFieldsDataSet(time, 'n[0]')

        elif (type(specie) == str and specie.lower() in self._proton) or specie == 1:
            return self._readFieldsDataSet(time, 'n[1]')

        elif (type(specie) == str and specie.lower() in self._alpha) or specie == 2:
            return self._readFieldsDataSet(time, 'n[2]')

        elif (type(specie) == str and specie.lower() in self._ion):

            out = np.zeros_like(self._readFieldsDataSet(time, 'n[0]'))

            f = h5py.File(os.path.join(self.path, self._fieldfilename))
            times = list(f.keys())
            firstTimeGrp = '/'+times[0]
            inTheGroup = f.get(firstTimeGrp)

            for dataSet in inTheGroup.items() :
                if 'n[' in dataSet[0] and dataSet[0] != 'n[0]' :
                    print(dataSet[0])
                    out = out + self._readFieldsDataSet(time, dataSet[0])
            return out

        else :
            raise ValueError("Unknown component : should be ['e', 'p', 'a', 'i'], or [0, 1, 2]")



    # __________________________________________________________________________
    #
    def GetN(self, time, specie):

        """
        Returns the particle density of the given 'specie'
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938

        """

        if (type(specie) == str and specie.lower() in self._electron) or specie == 0:
            name = 'n[0]'

        elif (type(specie) == str and specie.lower() in self._proton) or specie == 1:
            name = 'n[1]'

        elif specie == 2:
            name = 'n[2]'

        elif specie.lower() in self._ion:
            name = 'n[0]'

        else:
            raise ValueError("Unknown specie '{0:s}'".format(specie))

        data = self._readFieldsDataSet(time, name)

        return data
    # __________________________________________________________________________





    def GetP(self, time, specie):

        """

        """

        Pxx = self.GetPxx(time, specie)
        Pxy = self.GetPxy(time, specie)
        Pxz = self.GetPxz(time, specie)
        Pyy = self.GetPyy(time, specie)
        Pyz = self.GetPyz(time, specie)
        Pzz = self.GetPzz(time, specie)

        shape = Pxx.shape

        P = np.empty(Pxx.shape+(3,3), dtype = Pxx.dtype)

        P[..., 0, 0] = Pxx
        P[..., 0, 1] = Pxy
        P[..., 0, 2] = Pxz
        P[..., 1, 0] = Pxy
        P[..., 1, 1] = Pyy
        P[..., 1, 2] = Pyz
        P[..., 2, 0] = Pxz
        P[..., 2, 1] = Pyz
        P[..., 2, 2] = Pzz

        return P







    def getP(self, time, specie, component):
        if (type(specie) == str and specie.lower() in self._electron) or specie == 0 :
            ispecie = 0
        elif (type(specie) == str and specie.lower() in self._proton) or specie == 1 :
            ispecie = 0
        elif (type(specie) == str and specie.lower() in self._alpha) or specie == 2:
            ispecie = 0
        else :
            raise ValueError("specie has to be a single existing population... as it is not an additive quantity")

        if component.lower() in ['xx', 'xy', 'xz', 'yy', 'yz', 'zz'] :
            name = 'P'+component.lower()+'['+str(ispecie)+']'
        else :
            raise ValueError("component has to be in ['xx', 'xy', 'xz', 'yy', 'yz', 'zz']")

        return self._readFieldsDataSet(time, name)












    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPxx(self, time, species):

        """
        Returns the component of the pressure tensor of a species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """

        if species.lower() in self._electron:
            name = 'Pxx[0]'

        elif species.lower() in self._proton:
            name = 'Pxx[1]'

        elif species.lower() in self._ion:
            # full Pressure tensor not additive and not written by heckle
            raise valueError("No pressure tensor for species %s" % species)

        else:
            raise ValueError("Unknown species '%s'" % species)

        data = self._readFieldsDataSet(time, name)

        return data
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPxy(self, time, species):

        """
        Returns the component of the pressure tensor of a species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """

        if species.lower() in self._electron:
            name = 'Pxy[0]'

        elif species.lower() in self._proton:
            name = 'Pxy[1]'

        elif species.lower() in self._ion:
            # full Pressure tensor not additive and not written by heckle
            raise valueError("No pressure tensor for species %s" % species)

        else:
            raise ValueError("Unknown species '%s'" % species)

        data = self._readFieldsDataSet(time, name)

        return data
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPxz(self, time, species):

        """
        Returns the component of the pressure tensor of a species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """

        if species.lower() in self._electron:
            name = 'Pxz[0]'

        elif species.lower() in self._proton:
            name = 'Pxz[1]'

        elif species.lower() in self._ion:
            # full Pressure tensor not additive and not written by heckle
            raise valueError("No pressure tensor for species %s" % species)

        else:
            raise ValueError("Unknown species '%s'" % species)

        data = self._readFieldsDataSet(time, name)

        return data
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPyy(self, time, species):

        """
        Returns the component of the pressure tensor of a species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """

        if species.lower() in self._electron:
            name = 'Pyy[0]'

        elif species.lower() in self._proton:
            name = 'Pyy[1]'

        elif species.lower() in self._ion:
            # full Pressure tensor not additive and not written by heckle
            raise valueError("No pressure tensor for species %s" % species)

        else:
            raise ValueError("Unknown species '%s'" % species)

        data = self._readFieldsDataSet(time, name)

        return data
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPyz(self, time, species):

        """
        Returns the component of the pressure tensor of a species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """

        if species.lower() in self._electron:
            name = 'Pyz[0]'

        elif species.lower() in self._proton:
            name = 'Pyz[1]'

        elif species.lower() in self._ion:
            # full Pressure tensor not additive and not written by heckle
            raise valueError("No pressure tensor for species %s" % species)

        else:
            raise ValueError("Unknown species '%s'" % species)

        data = self._readFieldsDataSet(time, name)

        return data
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetPzz(self, time, species):

        """
        Returns the component of the pressure tensor of a species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938
        """

        if species.lower() in self._electron:
            name = 'Pzz[0]'

        elif species.lower() in self._proton:
            name = 'Pzz[1]'

        elif species.lower() in self._ion:
            # full Pressure tensor not additive and not written by heckle
            raise valueError("No pressure tensor for species %s" % species)

        else:
            raise ValueError("Unknown species '%s'" % species)

        data = self._readFieldsDataSet(time, name)

        return data
    #==========================================================





    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetMass(self, species):

        """
        Returns the mass of particles of a particular species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938

        """
        # these values should be read from attribute of species.h5
        if species.lower() in self._electron:
            return 0.
        elif species.lower() in self._proton:
            return 1.
        else:
            raise ValueError("Unknown species '%s'" % species)
    #==========================================================




    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetCharge(self, species):

        """
        Returns the electric charge of particles of a particular species
        Exemple  :
        Creation : 2015-02-27 16:15:40.856938

        """
        # these values should be read from attribute of species.h5
        if species.lower() in self._electron:
            return -1.
        elif species.lower() in self._proton:
            return 1.
        else:
            raise ValueError("Unknown species '%s'" % species)
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetHyperResistivity(self):

        """
        @todo: Brief Docstring for GetHyperResistivity
        Longer description here

        @return: @todo

        Exemple  :

        Creation : 2015-03-03

        """
        return self.hprsty
    #==========================================================







    #----------------------------------------------------------
    #----------------------------------------------------------
    def GetResistivity(self):

        """
        @todo: Brief Docstring for GetHyperResistivity
        Longer description here

        @return: @todo

        Exemple  :

        Creation : 2015-03-03

        """

        return self.rsty
    #==========================================================









    #----------------------------------------------------------
    #----------------------------------------------------------
    def indices2coord(self, indices, quantity=None):

        """
        Converts an array of shape (2, N) in physical coordinates
        normalized to the ion inertial length.

        Args :
            quantity : Not used in the Heckle code since all
            quantities are defined on the same grid.

        Exemple  :

        Creation: 2015-10-08
        """

        # could more general ? indices could be a np array of whatever rank
        # (1, 2 or 3) and whatever size... even if right now, it only seems
        # to be used in recrate.py (and picgsfc)
        coords = np.ndarray(indices.shape)
        ndim   = indices.shape[0]

        for c in np.arange(ndim):
            coords[c,:] = self.dl[c]*indices[c,:]

        return coords
    #==========================================================




    def index2Coord(self,
                    index,
                    axis) :

        """
        index : 1d numpy array containing a set of indices

        axis : 0, 1 or 2 is the direction along which indices are given

        output : a 1d numpy array of same shape as index withj the coordinates
                 in physical units
        """

        index = np.array(index)
        coord = np.empty(index.shape, dtype = float)

        if axis == 0 or axis == 1 or axis == 2 :
           coord = index*self.dl[axis]
        else :
           raise ValueError("axis value is mandatory & can only be 0, 1 or 2")

        return coord


    def coord2Index(self,
                    coord,
                    axis) :

        """
        coord : 1d numpy array containing a set of coordinates

        axis : 0, 1 or 2 is the direction along which indices are given

        output : a 1d numpy array of same shape as coord with the associateds
                 indices on G1 grid
        """

        coord = np.array(coord)
        index = np.empty(coord.shape, dtype = int)

        if axis == 0 or axis == 1 or axis == 2 :
           index = int(coord/self.dl[axis])
        else :
           raise ValueError("axis value is mandatory & can only be 0, 1 or 2")

        return index


    def spawnAxis(self,
                  index,
                  coord,
                  axis) :
        """
        """

        # spawn axis from a given indices
        if index is not None :
           if index.__len__() == 2 :
              myaxis = np.linspace(index[0], index[1], index[1]-index[0]+1)
           else :
              raise ValueError("index has to be of length 2")

        # spawn axis from given coordinates
        if coord is not None :
           if coord.__len__() == 2 :
              if axis == 0 or axis == 1 or axis == 2 :
                 minbound = int(coord[0]/self.dl[axis])*self.dl[axis]
                 maxbound = int(coord[1]/self.dl[axis])*self.dl[axis]
                 myaxis = np.arange(minbound, maxbound, self.dl[axis])
              else :
                 raise ValueError("axis has to be specified")
           else :
              raise ValueError("coord has to be of length 2")

        return myaxis



    def fourierFlux(self, time) :

       import pyfftw
       from pylab import meshgrid


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



    def getAz(self, time) :

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

        #bx = self.GetB(time)[..., 0]+1j*np.zeros(yv.shape)
        #by = self.GetB(time)[..., 1]+1j*np.zeros(xv.shape)
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



    def heckleSet(self):

        f = open(os.path.join(self.path, "heckle.txt"), 'r')
        empty = f.readline()

        junk = f.readline()
        li = f.readline()
        nw = li.strip().split()
        self.ncells = [int(nw[0]), int(nw[1]), int(nw[2])]

        junk = f.readline()
        li = f.readline()
        nw = li.strip().split()
        self.domsize = [float(nw[0]), float(nw[1]), float(nw[2])]

        junk = f.readline()
        li = f.readline()
        nw = li.strip().split()
        self.bc = [int(nw[0]), int(nw[1]), int(nw[2])]

        junk = f.readline()
        li = f.readline()
        nw = li.strip().split()
        self.ts = float(nw[0])
        self.nt = int(nw[1])

        junk = f.readline()
        li = f.readline()
        nw = li.strip().split()
        self.tf = int(nw[0])
        self.tp = int(nw[1])
        self.tt = int(nw[2])

        f.close()

