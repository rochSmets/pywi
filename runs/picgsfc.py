import re
import os
import string as strm
import numpy as np
import glob
import collections
import particles

# constants




class PICGSFCrun:
    """
    PICGSFCrun is a simulation handle for a PIC simulation from Michael
    Hesse's code
    """

    Name = "PICGSFCRun"

#=============================================================================
#=============================================================================
    def __init__(self, path, runID='r', tmax=80.0):


        # These values are created
        # when the class is instantiated.
        self.modelname =  'fullPIC'
        self.codename  = 'picgsfc'
        self.runID     =  runID
        self.path      =  path



        #---------------------------------------------------------------------
        # Now we open the first file in the run directory to grab the run
        # parameters since the pic code does not have a parameter file but
        # writes all parameters within diagnostic files, this is one good way
        # to get them
        #---------------------------------------------------------------------

        # list all field files and take the first field file
        filenames    = os.listdir(self.path)
        firstfile    = (re.findall(r'fields-\d+.dat',strm.join(filenames)))[0]

        # open the file in binary mode and skip the 4B fortran header
        f   = open(os.path.join(path,firstfile), 'rb')
        f.seek(4, os.SEEK_SET)

        # Then read the parameters
        params=np.fromfile(f, np.dtype([('it','i'),('dt','f'),\
                        ('teti','f'), ('xmax','f'),\
                        ('zmax','f'), ('nx','i'),\
                        ('nz','i')]),1)

        nss = 2 # hard coded number of populations (here 1electron 1 proton)
        nx  = (params['nx'])[0]
        nz  = (params['nz'])[0]

        # ok, now we need to read wpewce etc. so that we can renormalize
        # distances etc. to ion scales

        sizeint   = 4
        sizefloat = 4
        size2ds   = nss*sizefloat*nx*nz
        size2d    = sizefloat*nx*nz

        # shift the file pointer after the field data

        # it, dt,teti,xmax,zmax,nx,nz ,
        # vxs,vys,vzs, bx,by,bz, ex,ey,ez,
        # ns, xe,ze, mass,q,time,wpewce,dfac,
        #pxx,pyy,pzz,pxy,pxz,pyz

        offset = 4 + sizeint + sizefloat*4 + \
                 2*sizeint + 3*size2ds + 6*size2d + size2ds

        f.seek( offset, os.SEEK_SET)

        x = np.fromfile(f, np.dtype('f'), nx)
        z = np.fromfile(f, np.dtype('f'), nz)

        # read those parameters
        params2 = np.fromfile(f, np.dtype([('mi','f'),
                                           ('me','f'),
                                           ('qi','f'),
                                           ('qe','f'),
                                           ('time','f8'),
                                           ('wpewce','f'),
                                           ('dfaci','f'),
                                           ('dface','f')]),1)

        # and close the file !
        f.close()


        #---------------------------------------------------------------------
        # Ok, we are ready to set the class attributes with the parameters
        # we've read and re-normalize everything to ion scales
        #---------------------------------------------------------------------

        self.teti     = (params['teti'])[0]

        self.charges  = { 'ions' : (params2['qi'])[0] , 
                          'electrons': (params2['qe'])[0]}

        self.dfac     = { 'ions' : (params2['dfaci'])[0] , 
                          'electrons' : (params2['dface'])[0]}

        self.ncells   = (nx-1, nz-1)

        self.wpewce   = (params2['wpewce'])[0]

        self.masses   = {'ions': (params2['mi'])[0], 
                         'electrons': (params2['me'])[0] }

        self.domsize  = ((params['xmax'])[0]/np.sqrt(self.masses['ions']), 
                         (params['zmax'])[0]/np.sqrt(self.masses['ions']))

        self.dl       = ((x[1]-x[0])/np.sqrt(self.masses['ions']),
                         (z[1]-z[0])/np.sqrt(self.masses['ions']))

        self.x        = x   # warning, it is in delta_e
        self.z        = z   # warning : it is in delta_e

        self.bcs      = None
        self.ts       = (params['dt'])[0]

        # to find 'tsf' we detect it within the name of 'firstfile'
        tsf = re.findall(r'\d+', firstfile)
        self.tsf = float(tsf[0])

        # now take the first particle file :
        pfiles    = re.findall(r'part-\d+.dat',strm.join(filenames))

        # no particle file found
        if (len(pfiles) == 0):
            tsp = None

        # if one or more particles files are found
        # find the first
        # detect the time
        # convert to float
        else:
            firstpfile=pfiles[0]
            tsp = re.findall(r'\d+', firstpfile)
            print tsp
            tsp = float(tsp[0])


        self.tsp = tsp


        # Now we need to set 'nt', the total number of time steps
        self.nt = np.round(tmax*self.masses['ions']*self.wpewce/self.ts)


        # set the out-of-plane direction to be 'y'
        self.outofplane = 1


        if (silent == 'No'):
            print "PICGSFCrun> teti                      :  ", self.teti
            print "PICGSFCrun> ion charge                :  ", self.charges['ions']
            print "PICGSFCrun> electron charge           :  ", self.charges['electrons']
            print "PICGSFCrun> ion dfac                  :  ", self.dfac['ions']
            print "PICGSFCrun> electron dfac             :  ", self.dfac['electrons']
            print "PICGSFCrun> ncells                    :  ", self.ncells
            print "PICGSFCrun> wpewce                    :  ", self.wpewce
            print "PICGSFCrun> domain size               :  ", self.domsize
            print "PICGSFCrun> mesh size                 :  ", self.dl
            print "PICGSFCrun> time step                 :  ", self.ts
            print "PICGSFCrun> field dump ts  (wpe^-1)   :  ", self.tsf
            print "PICGSFCrun> total nb steps (wpe^-1)   :  ", self.nt
#=============================================================================







    #==========================================================
    #==========================================================
    def model_name(self):
        """returns the type of model (here hybrid)

        This routine is required for the global API to work

        @return: either : 'hybrid', 'HallMHD', 'fullPIC'

        Exemple  :

        Creation : 2012-08-15 14:27:06.843066

        """
        return self.modelname
    #==========================================================








    # Self is used as an argument to
    # pretty much all class functions.
    # However, you do NOT need to pass
    # the argument self if you call this method
    # from a Class, because the class provides
    # the value of itself.
    def display(self):

        """

        displays the parameters of the run

        """
        print "PICGSFCrun> teti             :  ", self.teti
        print "PICGSFCrun> ion charge       :  ", self.charges['ions']
        print "PICGSFCrun> electron charge  :  ", self.charges['electrons']
        print "PICGSFCrun> ion dfac         :  ", self.dfac['ions']
        print "PICGSFCrun> electron dfac    :  ", self.dfac['electrons']
        print "PICGSFCrun> ncells           :  ", self.ncells
        print "PICGSFCrun> wpewce           :  ", self.wpewce
        print "PICGSFCrun> domain size      :  ", self.domsize
        print "PICGSFCrun> mesh size        :  ", self.dl
        print "PICGSFCrun> time step        :  ", self.ts
        print "PICGSFCrun> field dump ts    :  ", self.tsf







#===================================================================
#===================================================================
    def GetCoord(self, extent=None, axis=0):
        """

        Returns an array of coordinates along the required
        axis between extent[0] and extent[1]. If extent is None
        then the function returns the coords. for the whole simulation domain

        """
        if axis == 0:
            dl = self.dl[0]
            coord_idx = 0
        elif axis == 1:
            dl = self.dl[1]
            coord_idx = 1
        else:
            print "Invalid axis"
            return None

        if extent is None:
            return -dl + dl*np.arange(self.ncells[coord_idx]+1)
        else:
            npoints = (extent[1]-extent[0])/dl +1
            return extent[0] + dl*np.arange(npoints)

#=================================================================







#=================================================================
#=================================================================
    def GetNbDiag(self, datatype):

        """

        Returns the number of diagnostics of type 'datatype' the run is 
        supposed to produce

        Args :
            - datatype : 'fields', 'particles'

        Note : This is different from the number of diagnostics of 
               type 'datatype' that are actually in the path directory of
               this run (see GetNbFieldFiles)
        """

        if datatype == 'fields':
            return np.int(self.nt/(self.tsf/self.ts))
        elif datatype == 'particles':
            return np.int(self.nt/(self.tsp/self.ts))
        else:
            return None
#=================================================================







#=================================================================
#=================================================================
    def GetMass(self, species='ions'):

        """

        Returns the mass of the requested specie

        """
        if species in ['ions', 'ion', 'protons', 'proton'] :
            return 1.0

        elif species in ['electrons', 'electron'] :
            return self.masses['electrons']/self.masses['ions']

        else:
            print "Invalid species name"
            return None
#=================================================================






    #==========================================================
    #==========================================================
    def GetFieldFiles(self,t1,t2):
        """returns a list of the field files names for which the time
        is between t1 and t2

        @return: file list or None if t1 or t2 is not in the path dir.

        Exemple  :

        Creation : 2012-11-29 16:55:28.247903

        """
        t1e = round(t1*self.masses['ions']*self.wpewce)
        t2e = round(t2*self.masses['ions']*self.wpewce)


        file1 = os.path.join(self.path,'fields-%05d.dat' %(t1e))
        file2 = os.path.join(self.path,'fields-%05d.dat' %(t2e))


        if os.path.isfile(file1) == False:
            print 'file %s not found' %(file1)
            return None

        if os.path.isfile(file2) == False:
            print 'file %s not found' % (file2)
            return None

        files = []
        allfiles = sorted(glob.glob(os.path.join(self.path,'fields-*')))

        for f in allfiles:
            fn = re.findall('fields-\d+.dat',f)[0]
            i = int(re.findall('\d+',fn)[0])
            if i >= t1e and i <= t2e :
                files.append(f)

        return files

    #==========================================================




    #==========================================================
    #==========================================================
    def GetFileTime(self, filename):
        """returns the time at which this file was written in
        ion cyclotron periods

        Exemple  :  t = run.GetFileTime('fields-05000.dat')

        Creation : 2012-11-30 10:03:51.966027

        """
        #fn = re.findall('fields-\d+.dat$',filename)
        filename = os.path.basename(filename)
        time = int(re.findall('\d+',filename)[0])

        time_ion = time /(self.wpewce*self.masses['ions'])

        return time_ion

    #==========================================================




    #==========================================================
    #==========================================================
    def GetCharge(self, species='ions'):

        """

        Returns the charge of the requested species

        """
        if species in ['ions', 'ion', 'protons', 'proton'] :
            return self.charges['ions']

        elif species in ['electrons', 'electron'] :
            return self.charges['electrons']

        else:
            print "Invalid species name"
            return None
    #==========================================================





    #==========================================================
    #==========================================================
    def GetOutOfPlaneDir(self):

        """

        Returns the index for the out-of-plane dimension

        """

        return self.outofplane
    #==========================================================





    #==========================================================
    #==========================================================
    def GetRunTime(self):

        """

        Returns the total simulation time in inverse of ion gyrofrequency

        """


        return self.ts*self.nt/(self.masses['ions']*self.wpewce)
    #==========================================================




    #==========================================================
    #==========================================================
    def GetTimeStep(self, datatype):

        """

        Returns the time step between two diagnostic files of type 'datatype'
        in inverse ion cyclotron periods

        Args :
            datatype : 'fields', 'particles'

        """

        if datatype == 'fields':
            return self.tsf/(self.masses['ions']*self.wpewce)
        elif datatype == 'particles':
            return self.tsp/(self.masses['ions']*self.wpewce)
        else:
            return None
    #==========================================================







    #==========================================================
    def whereAmI(self):
        return self.path
    #==========================================================





    #==========================================================
    #==========================================================
    def GetNbFieldFiles(self):

        """

        returns the number of field diagnostic files
        in the path directory.

        """


        filenames    = os.listdir(self.path)
        f_filenames = re.findall(r'fields-\d+\.dat', strm.join(filenames))
        return len(f_filenames)
    #==========================================================





    #==========================================================
    #==========================================================
    def iter2time(self, it, datatype):
        if datatype == 'fields':
            return (it+1) * self.tsf/(self.wpewce*self.masses['ions'])    # +1 because the first diag file is t = dt(datatype)
        elif datatype == 'particles':
            return (it+1) * self.tsp/(self.wpewce*self.masses['ions'])
        else :
            return None
    #==========================================================






    #==========================================================
    #==========================================================
    def time2iter(self, time, datatype):

        """

        Converts a time in number of iterations


        Args :
            - time     : a given time, in inverse ion cyclotron frequency
            - datatype : the type of data, fields, particles, diag.

        Example :
            run.time2iter(80,0, 'fields')

        TODO :
            diagnostics
        """


        if datatype == 'fields':
            return int(time *self.wpewce*self.masses['ions']/ self.tsf)-1
        elif datatype == 'particles':
            return int(time *self.wpewce * self.masses['ions']/ self.tsp)-1
#       elif datatype == 'diag':
#           return int(time / self.tst)
        else :
            return None
    #==========================================================







    #==========================================================
    #==========================================================
    def indices2coord(self, indices, quantity=None):

        """

        Converts an array of shape (2, N) in physical coordinates 
        normalized to the ion inertial length.

        Args :
            quantity : default=None, is a string giving the name of the
            field we want the coordinates of. In the PICgsfc code for instance,
            the fields are defined on a Yee mesh and thus do not have the
            exact same coordinates.

            quantity can be set to 'e[x,y,z]','b[x,y,z]','moments'

        """

        coords = np.ndarray(indices.shape)

        if (quantity in ['bx', 'ez']):
            coords[0,:] = indices[0,:]*self.dl[0] - self.dl[0]
            coords[1,:] = indices[1,:]*self.dl[1] - self.dl[1]*0.5

        elif(quantity in ['by']):
            coords[0,:] = indices[0,:]*self.dl[0] - self.dl[0]*0.5
            coords[1,:] = self.dl[1] * (indices[1,:]-0.5)

        elif (quantity in ['bz', 'ex']):
            coords[0,:] = self.dl[0]*(indices[0,:]-0.5)
            coords[1,:] = self.dl[1]*(indices[1,:]-1)

        elif (quantity in ['moment', None]):
            coords[0,:] = self.dl[0]*(indices[0,:]-1)
            coords[1,:] = self.dl[1]*(indices[1,:]-1)


        return coords
    #==========================================================





    #==========================================================
    #==========================================================
    def coord2indices(self, coord):
        """Converts  a coordinate tuple into an index array of shape (2,N)
                this function returns the indices of cell corners

        Creation : 2012-08-13 10:05:22.648635

        """

        # np.rint() returns the nearest integer
        indices      = np.ndarray(coord.shape)
        indices[0,:] = np.rint(coord[0,:]/self.dl[0]+1)
        indices[1,:] = np.rint(coord[1,:]/self.dl[1]+1)

        return indices
    #==========================================================




    #==========================================================
    #==========================================================
    def GetField(self, fieldname, time, silent='No'):

        """

        returns an array containing the data for the field 'fieldname' at the
        specified 'time'


        Args :
            - fieldname : string giving the name of the quantity to read
            - time      : requested time, in inverse ion cyclotron frequency

        KEYWORDS :
            - silent : silent mode, default='No'

        """

        #get the time in wpe^-1 units
        timewpe  = time*self.masses['ions']*self.wpewce

        # build the filename
        filename = os.path.join(self.path,'fields-%05d.dat' %round(timewpe))

        # check if the file is in the directory
        if not os.path.isfile(filename):
            print "Error - file not found '" + filename+"'"
            return None

        # the file exists
        else:
            if (silent == 'No'):
                print "GetField> Retrieving ",fieldname,"..."
            data = self.GetFieldData(fieldname, filename, silent=silent)
            return data
    #==========================================================








    #==========================================================
    #==========================================================
    def GetFieldData(self, fieldname, filename, silent='No'):

        """

        Returns a 2D array containing the data for the field 'fielname'


        """


        nss       = 2   # hard coded # of populations (electrons and ions)
        sizeint   = 4
        sizefloat = 4
        size2d    = (self.ncells[0]+1)*(self.ncells[1]+1)*sizefloat
        size2ds   = size2d*nss

        fortran_offset = 4
        offset = fortran_offset + sizeint + 4*sizefloat + 2*sizeint

        fieldname = str.lower(fieldname)
        # These values are created
        # when the class is instantiated.
        fieldpos = {}
        fieldpos['bx']    = offset + 3*size2ds
        fieldpos['by']    = offset + 3*size2ds + size2d
        fieldpos['bz']    = offset + 3*size2ds + 2*size2d
        fieldpos['ex']    = offset + 3*size2ds + 3*size2d
        fieldpos['ey']    = offset + 3*size2ds + 4*size2d
        fieldpos['ez']    = offset + 3*size2ds + 5*size2d
        fieldpos['dni']   = offset + 3*size2ds + 6*size2d
        fieldpos['dne']   = offset + 3*size2ds + 6*size2d + size2d
        fieldpos['vex']   = offset + size2d
        fieldpos['vey']   = offset + nss*size2d + size2d
        fieldpos['vez']   = offset + nss*size2d*2 + size2d
        fieldpos['vex']   = offset + size2d
        fieldpos['vey']   = offset + nss*size2d + size2d
        fieldpos['vez']   = offset + nss*size2d*2 + size2d
        fieldpos['vix']   = offset
        fieldpos['viy']   = offset + nss*size2d
        fieldpos['viz']   = offset + 2*nss*size2d
        fieldpos['pexx']  = offset + 3*size2ds + 6*size2d + size2ds + (self.ncells[0]+1)*sizefloat + (self.ncells[1]+1)*sizefloat \
                    + 2*nss*sizefloat + 2*sizefloat + sizefloat + nss*sizefloat + size2d
        fieldpos['pexy']  = offset + 3*size2ds + 6*size2d + size2ds + (self.ncells[0]+1)*sizefloat + (self.ncells[1]+1)*sizefloat \
                    + 2*nss*sizefloat + 2*sizefloat + sizefloat + nss*sizefloat + 7*size2d
        fieldpos['peyx']  = offset + 3*size2ds + 6*size2d + size2ds + (self.ncells[0]+1)*sizefloat + (self.ncells[1]+1)*sizefloat \
                    + 2*nss*sizefloat + 2*sizefloat + sizefloat + nss*sizefloat + 7*size2d
        fieldpos['pexz']  = offset + 3*size2ds + 6*size2d + size2ds + (self.ncells[0]+1)*sizefloat + (self.ncells[1]+1)*sizefloat \
                    + 2*nss*sizefloat + 2*sizefloat + sizefloat + nss*sizefloat + 9*size2d
        fieldpos['pezx']  = offset + 3*size2ds + 6*size2d + size2ds + (self.ncells[0]+1)*sizefloat + (self.ncells[1]+1)*sizefloat \
                    + 2*nss*sizefloat + 2*sizefloat + sizefloat + nss*sizefloat + 9*size2d
        fieldpos['peyy']  = offset + 3*size2ds + 6*size2d + size2ds + (self.ncells[0]+1)*sizefloat + (self.ncells[1]+1)*sizefloat \
                    + 2*nss*sizefloat + 2*sizefloat + sizefloat + nss*sizefloat + 3*size2d
        fieldpos['peyz']  = offset + 3*size2ds + 6*size2d + size2ds + (self.ncells[0]+1)*sizefloat + (self.ncells[1]+1)*sizefloat \
                    + 2*nss*sizefloat + 2*sizefloat + sizefloat + nss*sizefloat + 11*size2d
        fieldpos['pezy']  = offset + 3*size2ds + 6*size2d + size2ds + (self.ncells[0]+1)*sizefloat + (self.ncells[1]+1)*sizefloat \
                    + 2*nss*sizefloat + 2*sizefloat + sizefloat + nss*sizefloat + 11*size2d
        fieldpos['pezz']  = offset + 3*size2ds + 6*size2d + size2ds + (self.ncells[0]+1)*sizefloat + (self.ncells[1]+1)*sizefloat \
                    + 2*nss*sizefloat + 2*sizefloat + sizefloat + nss*sizefloat + 5*size2d
        fieldpos['pixx']  = offset + 3*size2ds + 6*size2d + size2ds + (self.ncells[0]+1)*sizefloat + (self.ncells[1]+1)*sizefloat \
                    + 2*nss*sizefloat + 2*sizefloat + sizefloat + nss*sizefloat
        fieldpos['pixy']  = offset + 3*size2ds + 6*size2d + size2ds + (self.ncells[0]+1)*sizefloat + (self.ncells[1]+1)*sizefloat \
                    + 2*nss*sizefloat + 2*sizefloat + sizefloat + nss*sizefloat + 6*size2d
        fieldpos['piyx']  = offset + 3*size2ds + 6*size2d + size2ds + (self.ncells[0]+1)*sizefloat + (self.ncells[1]+1)*sizefloat \
                    + 2*nss*sizefloat + 2*sizefloat + sizefloat + nss*sizefloat + 6*size2d
        fieldpos['pixz']  = offset + 3*size2ds + 6*size2d + size2ds + (self.ncells[0]+1)*sizefloat + (self.ncells[1]+1)*sizefloat \
                    + 2*nss*sizefloat + 2*sizefloat + sizefloat + nss*sizefloat + 8*size2d
        fieldpos['pizx']  = offset + 3*size2ds + 6*size2d + size2ds + (self.ncells[0]+1)*sizefloat + (self.ncells[1]+1)*sizefloat \
                    + 2*nss*sizefloat + 2*sizefloat + sizefloat + nss*sizefloat + 8*size2d
        fieldpos['piyy']  = offset + 3*size2ds + 6*size2d + size2ds + (self.ncells[0]+1)*sizefloat + (self.ncells[1]+1)*sizefloat \
                    + 2*nss*sizefloat + 2*sizefloat + sizefloat + nss*sizefloat + 2*size2d
        fieldpos['piyz']  = offset + 3*size2ds + 6*size2d + size2ds + (self.ncells[0]+1)*sizefloat + (self.ncells[1]+1)*sizefloat \
                    + 2*nss*sizefloat + 2*sizefloat + sizefloat + nss*sizefloat + 10*size2d
        fieldpos['pizy']  = offset + 3*size2ds + 6*size2d + size2ds + (self.ncells[0]+1)*sizefloat + (self.ncells[1]+1)*sizefloat \
                    + 2*nss*sizefloat + 2*sizefloat + sizefloat + nss*sizefloat + 10*size2d
        fieldpos['pizz']  = offset + 3*size2ds + 6*size2d + size2ds + (self.ncells[0]+1)*sizefloat + (self.ncells[1]+1)*sizefloat \
                    + 2*nss*sizefloat + 2*sizefloat + sizefloat + nss*sizefloat +4*size2d



        if fieldname not in fieldpos:
            print "Error, '"+fieldname+"' is not a valid fieldname"
            return


        # get the data
        if silent.lower() == 'no':
            print "opening '"+filename+"'... Accessing '"+ fieldname+ "' at position ", fieldpos[fieldname]

        f   = open(filename, 'rb') # open a binary file in read mode
        f.seek(fieldpos[fieldname], os.SEEK_SET)
        data = np.fromfile(f, "float32", (self.ncells[0]+1)*(self.ncells[1]+1))
        data = data.reshape((self.ncells[0]+1, self.ncells[1]+1),order='FORTRAN')
        f.close()
        return data
    #==========================================================






    #==========================================================
    #==========================================================
    def GetVe(self, time, silent='No'):

        vex = self.GetField("vex", time, silent=silent)
        vey = self.GetField("vey", time, silent=silent)
        vez = self.GetField("vez", time, silent=silent)

        # get the electron density. We dont normalize by dfac
        # since the vex,vey,vez are the fluxes and not the bulk flows yet
        # and they are not normalized by dfac either
        dne = self.GetField('dne', time, silent=silent)


        ve = np.zeros((3,self.ncells[0]+1,self.ncells[1]+1),"float32", order='FORTRAN')

        ve[0,:,:] = vex
        ve[1,:,:] = vey
        ve[2,:,:] = vez

        # now get the mean velocity by dividing
        # by the number of particles (not density)
        ve[0,1:,1:] /= dne[1:, 1:]
        ve[1,1:,1:] /= dne[1:, 1:]
        ve[2,1:,1:] /= dne[1:, 1:]


        return ve * self.wpewce*np.sqrt(self.masses['ions'])
    #==========================================================









    #==========================================================
    #==========================================================
    def GetVi(self, time, silent='No'):

        # this is really nv and not v
        vix = self.GetField("vix", time, silent=silent)
        viy = self.GetField("viy", time, silent=silent)
        viz = self.GetField("viz", time, silent=silent)

        # get the ion number (no dfac)
        dni = self.GetField('dni', time, silent=silent)


        vi = np.zeros((3,self.ncells[0]+1,self.ncells[1]+1),"float32", order='FORTRAN')

        vi[0,:,:] = vix
        vi[1,:,:] = viy
        vi[2,:,:] = viz


        # now get the mean velocity by dividing
        # by the number of counts
        vi[0,1:,1:] /= dni[1:, 1:]
        vi[1,1:,1:] /= dni[1:, 1:]
        vi[2,1:,1:] /= dni[1:, 1:]


        return vi  * self.wpewce*np.sqrt(self.masses['ions'])
    #==========================================================








    #==========================================================
    #==========================================================
    def GetV(self, time, species='ions', silent='No'):

        if (species.lower() in ['ions', 'ion', 'protons', 'proton']):
            V = self.GetVi(time, silent=silent)

        elif(species.lower() in ['electrons', 'electron']):
            V = self.GetVe(time, silent=silent)

        else:
            print "Invalid species name"
            return None

        return V
    #==========================================================







    #==========================================================
    #==========================================================
    def GetBYee(self, time, silent='No'):

        bx = self.GetField("bx", time, silent=silent)
        by = self.GetField("by", time, silent=silent)
        bz = self.GetField("bz", time, silent=silent)

        B = np.zeros((3,self.ncells[0]+1,self.ncells[1]+1),"float32", order='FORTRAN')
        B[0,:,:] = bx
        B[1,:,:] = by
        B[2,:,:] = bz

        return B*self.wpewce
    #==========================================================








    #==========================================================
    #==========================================================
    def GetB(self, time, silent='No'):

        nx = self.ncells[0]
        nz = self.ncells[1]

        bx = self.GetField("bx", time, silent=silent)
        by = self.GetField("by", time, silent=silent)
        bz = self.GetField("bz", time, silent=silent)

        Byee = np.zeros((3,self.ncells[0]+1,self.ncells[1]+1),
                        "float32",
                         order='FORTRAN')

        B    = np.zeros((3,self.ncells[0]+1,self.ncells[1]+1),
                         "float32",
                         order='FORTRAN')

        Byee[0,:,:] = bx*self.wpewce
        Byee[1,:,:] = by*self.wpewce
        Byee[2,:,:] = bz*self.wpewce

        # B is defined on the Yee grid
        # and we want to put it at the cell's corners

        #   Bx(+), corners(o), Bz(^), By(*)

        #     +    *    +     *    +
        #
        #2    o    ^    o-----^----o-----^
        #               |          |
        #     +    *    +     *    +     *
        #               |          |
        #1    o    ^    o-----^----o-----^
        #               |          |
        #     +    *    +     *    +     *
        #               |          |
        #0    o    ^    o-----^----o-----^
        #
        #     +    *    +     *    +     *
        #
        #     o    ^    o     ^    o    o

        #     0         1          2



        BxCorner = np.zeros((nx+1, nz+1), "float32", order="FORTRAN")
        ByCorner1 = np.zeros((nx+1, nz+1), "float32", order="FORTRAN")
        ByCorner2 = np.zeros((nx+1, nz+1), "float32", order="FORTRAN")
        BzCorner = np.zeros((nx+1, nz+1), "float32", order="FORTRAN")


        BxCorner[:,1:]      = 0.5*( Byee[0,:,1:] + Byee[0,:,:-1] )
        BzCorner[1:,:]      = 0.5*( Byee[2,1:,:] + Byee[2,:-1,:] )
        ByCorner1[1:,:]     = 0.5*( Byee[1,1:,:] + Byee[1,:-1,:] )
        ByCorner2[:,1:]     = 0.5*( Byee[1,:,1:] + Byee[1,:,:-1] )
        ByCorner            = 0.5*(ByCorner1 + ByCorner2)


        B[0,:,:]   =  BxCorner
        B[1,:,:]   =  ByCorner
        B[2,:,:]   =  BzCorner

        # fills ghost points so that no troubles occurs
        # in higher level routines when dividing by B or B^2
        B[0,:,0]   = B[0,:,1]
        B[2,0,:]   = B[2,1,:]
        B[1,0,:]   = B[1,1,:]
        B[1,:,0]   = B[1,:,1]
        B[1,0,0]   = B[1,1,1]


        return B
    #==========================================================









    #==========================================================
    #==========================================================
    def GetEYee(self, time, silent='No'):

        ex = self.GetField("ex", time, silent=silent)
        ey = self.GetField("ey", time, silent=silent)
        ez = self.GetField("ez", time, silent=silent)

        E = np.zeros((3,self.ncells[0]+1,self.ncells[1]+1),"float32", order='FORTRAN')
        E[0,:,:] = ex
        E[1,:,:] = ey
        E[2,:,:] = ez

        return E * self.wpewce**2 * np.sqrt(self.masses['ions'])
    #==========================================================









    #==========================================================
    #==========================================================
    def GetE(self, time, silent='No'):

        nx = self.ncells[0]
        nz = self.ncells[1]

        ex = self.GetField("ex", time, silent=silent)
        ey = self.GetField("ey", time, silent=silent)
        ez = self.GetField("ez", time, silent=silent)

        Eyee = np.zeros((3,nx+1, nz+1),"float32", order='FORTRAN')
        E    = np.zeros((3,nx+1, nz+1),"float32", order='FORTRAN')

        Eyee[0,:,:] = ex*self.wpewce**2 * np.sqrt(self.masses['ions'])
        Eyee[1,:,:] = ey*self.wpewce**2 * np.sqrt(self.masses['ions'])
        Eyee[2,:,:] = ez*self.wpewce**2 * np.sqrt(self.masses['ions'])

        # E is defined on the Yee grid
        # Ex is at the same location as Bz
        # Ez is at the same location as Bx
        # Ey is already at the cell corners

        #   Ez(+), V(o), Ex(^), Ey(o)

        #     +    *    +     *    +
        #
        #2    o    ^    o-----^----o-----^
        #               |          |
        #     +    *    +     *    +     *
        #               |          |
        #1    o    ^    o-----^----o-----^
        #               |          |
        #     +    *    +     *    +     *
        #               |          |
        #0    o    ^    o-----^----o-----^
        #
        #     +    *    +     *    +     *
        #
        #     o    ^    o     ^    o    o

        #     0         1          2



        ExCorner = np.zeros((nx+1, nz+1), "float32", order="FORTRAN")
        EzCorner = np.zeros((nx+1, nz+1), "float32", order="FORTRAN")


        EzCorner[:,1:]      = 0.5*( Eyee[2,:,1:] + Eyee[2,:,:-1] )
        ExCorner[1:,:]      = 0.5*( Eyee[0,1:,:] + Eyee[0,:-1,:] )


        E[0,:,:]   =  ExCorner
        E[1,:,:]   =  Eyee[1,:,:]
        E[2,:,:]   =  EzCorner


        return E
    #==========================================================










    #==========================================================
    #==========================================================
    def GetJ(self, time, silent='No'):

        # these are really the currents and not the bulk flows
        # and they are not multiplied by dfac yet.
        vix = self.GetField("vix", time, silent=silent)
        viy = self.GetField("viy", time, silent=silent)
        viz = self.GetField("viz", time, silent=silent)

        vex = self.GetField("vex", time, silent=silent)
        vey = self.GetField("vey", time, silent=silent)
        vez = self.GetField("vez", time, silent=silent)

        dfe = self.dfac['electrons']
        dfi = self.dfac['ions']

        J = np.zeros((3,self.ncells[0]+1,self.ncells[1]+1),"float32", order='FORTRAN')

        J[0,:,:] = (self.charges['ions']*dfi*vix + self.charges['electrons']*dfe*vex) *self.wpewce*np.sqrt(self.masses['ions'])
        J[1,:,:] = (self.charges['ions']*dfi*viy + self.charges['electrons']*dfe*vey) *self.wpewce*np.sqrt(self.masses['ions'])
        J[2,:,:] = (self.charges['ions']*dfi*viz + self.charges['electrons']*dfe*vez) *self.wpewce*np.sqrt(self.masses['ions'])

        return J
    #==========================================================







    #==========================================================
    #==========================================================
    def GetNi(self, time, silent='No'):
        ni = self.GetField( "dni", time, silent=silent)
        ni[0,:] = ni[1,:] # ghost cells so no /0
        ni[:,0] = ni[:,1]
        return ni*self.dfac['ions']
    #==========================================================






    #==========================================================
    #==========================================================
    def GetNe(self, time, silent='No'):
        ne = self.GetField( "dne", time, silent=silent)
        ne[0,:] = ne[1,:] # ghost cells so no /0
        ne[:,0] = ne[:,1]
        return ne*self.dfac['electrons']
    #==========================================================



    #==========================================================
    #==========================================================
    def GetN(self, time, species='ions', silent='No'):

        if (species.lower() in ['ions', 'ion', 'protons', 'proton']):
            n = self.GetNe(time, silent=silent)

        elif(species.lower() in ['electrons', 'electron']):
            n = self.GetNi(time, silent=silent)

        else:
            print "Invalid species name"
            return None

        return n
    #==========================================================





    #==========================================================
    #==========================================================
    def GetDebye(self,time):
        """Returns the debye length	

        Exemple  : run.GetDebye(10.0)

        Creation : 2012-12-29 17:03:08.181182
    
        """
        ne = self.GetNe(time)
        Pe = self.GetPe(time)
        Te = 1./3.*(Pe[0,0,:,:] + Pe[1,1,:,:] + Pe[2,2,:,:])/ne
        debye = np.sqrt(Te/ne)/self.wpewce
        return debye
    
    #==========================================================
    






    # SEBERG (#scipy) method to calculate the 2D integral flux from B
    #http://stackoverflow.com/questions/11854522/nested-loop-with-array-indexing-in-numpy
    def fast(self, flux, bx, bz, dx=1., dz=1.):

        from scipy.ndimage import convolve

        _flux = np.zeros((flux.shape[0]+1, flux.shape[1]+1), dtype=flux.dtype)
        temp_bx = np.zeros((bx.shape[0]+1, bx.shape[1]+1), dtype=bx.dtype)
        temp_bz = np.zeros((bz.shape[0]+1, bz.shape[1]+1), dtype=bz.dtype)

        _flux[:-1,:-1] = flux
        convolve(_flux[:-1,:-1], [[0, 0.5], [0.5, 0]], _flux[1:,1:])

        temp_bz[1:,1:-1] = bz[:,1:]*dx
        temp_bx[1:-1,1:] = bx[1:,:]*dz

        conv_b = np.array([[0.0, 0.5], [0.5, 0.5]])
        convolve(temp_bz[:-1,:-1], [[0.5, 0.5], [0.5, 0.]], temp_bz[1:,1:])
        convolve(temp_bx[:-1,:-1], [[-0.5, 0.5], [0.5, 0.]], temp_bx[1:,1:])

        _flux += temp_bz
        _flux += temp_bx

        return _flux[:-1,:-1]



    #==========================================================
    #==========================================================
    def GetFlux(self, time, silent='No'):

        #import matplotlib.pyplot as plt

        bx = self.GetField("bx", time, silent=silent) * self.wpewce
        by = self.GetField("by", time, silent=silent) * self.wpewce
        bz = self.GetField("bz", time, silent=silent) * self.wpewce



        flux  = np.zeros((self.ncells[0]+1,self.ncells[1]+1),"float32", order='FORTRAN')

        dx = self.dl[0]
        dz = self.dl[1]

        nx = self.ncells[0]
        nz = self.ncells[1]

        j = 0


        flux[1:,0] = flux[0,0] + np.cumsum(bz[:-1,0]*dx)

        flux[0,1:] = flux[0,0] - np.cumsum(bx[0,:-1]*dz)



        flux = self.fast(flux, bx, bz, dx, dz)


        return flux
    #==========================================================








    #==========================================================
    #==========================================================
    def GetPixx(self, time, silent='No'):

        df = self.dfac['ions']

        pixx = self.GetField("pixx", time, silent=silent)
        vix  = self.GetField("vix", time, silent=silent)   #n*v no defac
        dni = self.GetField("dni", time, silent=silent)    # n no dfac

        vix[1:,1:] /= dni[1:,1:]   # now v_bulk


        pixx = self.masses['ions']*(pixx * df  -  dni*df * vix**2)*self.wpewce**2

        return pixx
    #==========================================================






    #==========================================================
    #==========================================================
    def GetPixy(self, time, silent='No'):

        df = self.dfac['ions']

        pixy = self.GetField("pixy", time, silent=silent)
        vix  = self.GetField("vix", time, silent=silent)
        viy  = self.GetField("viy", time, silent=silent)
        dni = self.GetField("dni", time, silent=silent)

        vix[1:,1:] /= dni[1:,1:]
        viy[1:,1:] /= dni[1:,1:]

        pixy = self.masses['ions']*(pixy * df  -  dni*df * vix*viy)*self.wpewce**2

        return pixy
    #==========================================================







    #==========================================================
    #==========================================================
    def GetPixz(self, time, silent='No'):

        df = self.dfac['ions']

        pixz = self.GetField("pixz", time, silent=silent)
        vix  = self.GetField("vix", time, silent=silent)
        viz  = self.GetField("viz", time, silent=silent)
        dni = self.GetField("dni", time, silent=silent)


        vix[1:,1:] /= dni[1:,1:]
        viz[1:,1:] /= dni[1:,1:]


        #pixz = self.masses['ions']*(pixz * df  -  dni * vix*viz * df**2)*self.wpewce**2
        pixz = self.masses['ions']*(pixz * df  -  dni * vix*viz * df)*self.wpewce**2


        return pixz
    #==========================================================








    #==========================================================
    #==========================================================
    def GetPiyy(self, time, silent='No'):

        df = self.dfac['ions']

        piyy = self.GetField("piyy", time, silent=silent)
        viy  = self.GetField("viy", time, silent=silent)   #n*v no defac
        dni = self.GetField("dni", time, silent=silent)    # n no dfac

        viy[1:,1:] /= dni[1:,1:]   # now v_bulk


        piyy = self.masses['ions']*(piyy * df  -  dni*df * viy**2)*self.wpewce**2

        return piyy
    #==========================================================











    #==========================================================
    #==========================================================
    def GetPiyz(self, time, silent='No'):

        df = self.dfac['ions']

        piyz = self.GetField("piyz", time, silent=silent)
        viy  = self.GetField("viy", time, silent=silent)
        viz  = self.GetField("viz", time, silent=silent)
        dni = self.GetField("dni", time, silent=silent)


        viy[1:,1:] /= dni[1:,1:]
        viz[1:,1:] /= dni[1:,1:]


        #piyz = self.masses['ions']*(piyz * df  -  dni * viy*viz * df**2)*self.wpewce**2
        piyz = self.masses['ions']*(piyz * df  -  dni * viy*viz * df)*self.wpewce**2


        return piyz
    #==========================================================







    #==========================================================
    #==========================================================
    def GetPizz(self, time, silent='No'):

        df = self.dfac['ions']

        pizz = self.GetField("pizz", time, silent=silent)
        viz  = self.GetField("viz", time, silent=silent)   #n*v no defac
        dni = self.GetField("dni", time, silent=silent)    # n no dfac

        viz[1:,1:] /= dni[1:,1:]   # now v_bulk


        pizz = self.masses['ions']*(pizz * df  -  dni*df * viz**2)*self.wpewce**2

        return pizz
    #==========================================================







    #==========================================================
    #==========================================================
    def GetPi(self, time, silent='No'):


        nx = self.ncells[0]
        nz = self.ncells[1]

        pi = np.zeros((3,3,nx+1,nz+1), "float32", order="FORTRAN")

        pi[0,0,:,:] = self.GetPixx(time, silent=silent)
        pi[0,1,:,:] = self.GetPixy(time, silent=silent)
        pi[0,2,:,:] = self.GetPixz(time, silent=silent)


        pi[1,0,:,:] = self.GetPixy(time, silent=silent)
        pi[1,1,:,:] = self.GetPiyy(time, silent=silent)
        pi[1,2,:,:] = self.GetPiyz(time, silent=silent)


        pi[2,0,:,:] = self.GetPixz(time, silent=silent)
        pi[2,1,:,:] = self.GetPiyz(time, silent=silent)
        pi[2,2,:,:] = self.GetPizz(time, silent=silent)



        pi[:,:,0,:] = pi[:,:,1,:]
        pi[:,:,:,0] = pi[:,:,:,1]
        pi[:,:,0,0] = pi[:,:,1,1]

        return pi
    #==========================================================









    #==========================================================
    #==========================================================
    def GetPexx(self, time, silent='No'):

        df = self.dfac['electrons']

        pexx = self.GetField("pexx", time, silent=silent)
        vex  = self.GetField("vex", time, silent=silent)   #n*v no defac
        dne = self.GetField("dne", time, silent=silent)    # n no dfac

        vex[1:,1:] /= dne[1:,1:]   # now v_bulk


        pexx = self.masses['electrons']*(pexx * df  -  dne*df * vex**2)*self.wpewce**2

        return pexx
    #==========================================================






    #==========================================================
    #==========================================================
    def GetPexy(self, time, silent='No'):

        df = self.dfac['electrons']

        pexy = self.GetField("pexy", time, silent=silent)
        vex  = self.GetField("vex", time, silent=silent)
        vey  = self.GetField("vey", time, silent=silent)
        dne = self.GetField("dne", time, silent=silent)

        vex[1:,1:] /= dne[1:,1:]
        vey[1:,1:] /= dne[1:,1:]

        pexy = self.masses['electrons']*(pexy * df  -  dne*df * vex*vey)*self.wpewce**2

        return pexy
    #==========================================================







    #==========================================================
    #==========================================================
    def GetPexz(self, time, silent='No'):

        df = self.dfac['electrons']

        pexz = self.GetField("pexz", time, silent=silent)
        vex  = self.GetField("vex", time, silent=silent)
        vez  = self.GetField("vez", time, silent=silent)
        dne = self.GetField("dne", time, silent=silent)


        vex[1:,1:] /= dne[1:,1:]
        vez[1:,1:] /= dne[1:,1:]


        #pexz = self.masses['electrons']*(pexz * df  -  dne * vex*vez * df**2)*self.wpewce**2
        pexz = self.masses['electrons']*(pexz * df  -  dne *df* vex*vez )*self.wpewce**2


        return pexz
    #==========================================================








    #==========================================================
    #==========================================================
    def GetPeyy(self, time, silent='No'):

        df = self.dfac['electrons']

        peyy = self.GetField("peyy", time, silent=silent)
        vey  = self.GetField("vey", time, silent=silent)   #n*v no defac
        dne = self.GetField("dne", time, silent=silent)    # n no dfac

        vey[1:,1:] /= dne[1:,1:]   # now v_bulk


        peyy = self.masses['electrons']*(peyy * df  -  dne*df * vey**2)*self.wpewce**2

        return peyy
    #==========================================================











    #==========================================================
    #==========================================================
    def GetPeyz(self, time, silent='No'):

        df = self.dfac['electrons']

        peyz = self.GetField("peyz", time, silent=silent)
        vey  = self.GetField("vey", time, silent=silent)
        vez  = self.GetField("vez", time, silent=silent)
        dne = self.GetField("dne", time, silent=silent)


        vey[1:,1:] /= dne[1:,1:]
        vez[1:,1:] /= dne[1:,1:]


        #peyz = self.masses['electrons']*(peyz * df  -  dne * vey*vez * df**2)*self.wpewce**2
        peyz = self.masses['electrons']*(peyz * df  -  dne *df * vey*vez)*self.wpewce**2


        return peyz
    #==========================================================







    #==========================================================
    #==========================================================
    def GetPezz(self, time, silent='No'):

        df = self.dfac['electrons']

        pezz = self.GetField("pezz", time, silent=silent)
        vez  = self.GetField("vez", time, silent=silent)   #n*v no defac
        dne = self.GetField("dne", time, silent=silent)    # n no dfac

        vez[1:,1:] /= dne[1:,1:]   # now v_bulk


        pezz = self.masses['electrons']*(pezz * df  -  dne*df * vez**2)*self.wpewce**2

        return pezz
    #==========================================================







    #==========================================================
    #==========================================================
    def GetPe(self, time, silent='No'):


        nx = self.ncells[0]
        nz = self.ncells[1]

        pe = np.zeros((3,3,nx+1,nz+1), "float32", order="FORTRAN")

        pe[0,0,:,:] = self.GetPexx(time, silent=silent)
        pe[0,1,:,:] = self.GetPexy(time, silent=silent)
        pe[0,2,:,:] = self.GetPexz(time, silent=silent)


        pe[1,0,:,:] = self.GetPexy(time, silent=silent)
        pe[1,1,:,:] = self.GetPeyy(time, silent=silent)
        pe[1,2,:,:] = self.GetPeyz(time, silent=silent)


        pe[2,0,:,:] = self.GetPexz(time, silent=silent)
        pe[2,1,:,:] = self.GetPeyz(time, silent=silent)
        pe[2,2,:,:] = self.GetPezz(time, silent=silent)

        pe[:,:,0,:] = pe[:,:,1,:]
        pe[:,:,:,0] = pe[:,:,:,1]
        pe[:,:,0,0] = pe[:,:,1,1]


        return pe
    #==========================================================







    #==========================================================
    #==========================================================
    def utherm(self, time, species='ions'):
        """@todo: Returns the thermal energy of the specified specie

        The thermal energy u is defined as half the trace of the pressure
        tensor : u_s = 0.5*Tr(P_s).

        @param time is the time in inverse ion pulsation
        @param species is 'ions' or 'electrons'

        @return: 1/2Trace(P_s)

        Exemple  :

        Creation : 2012-08-28 11:49:20.052955

        """

        if species == 'ions':
            P = self.GetPi(time)

        if species == 'electrons':
            P = self.GetPe(time)

        return 0.5*(P[0,0,:,:]+P[1,1,:,:]+P[2,2,:,:])
    #==========================================================





    #==========================================================
    #==========================================================
    def VGradV(self, time, species='ions', silent='No'):

        nx = self.ncells[0]
        nz = self.ncells[1]
        odx = 1./self.dl[0]
        odz = 1./self.dl[1]

        vgv  = np.zeros((3,nx+1,nz+1), "float32", order="FORTRAN")
        dxdvx = np.zeros((nx+1, nz+1), "float32", order="FORTRAN")
        dzdvx = np.zeros((nx+1, nz+1), "float32", order="FORTRAN")
        dxdvy = np.zeros((nx+1, nz+1), "float32", order="FORTRAN")
        dzdvy = np.zeros((nx+1, nz+1), "float32", order="FORTRAN")
        dxdvz = np.zeros((nx+1, nz+1), "float32", order="FORTRAN")
        dzdvz = np.zeros((nx+1, nz+1), "float32", order="FORTRAN")

        if (species.lower() in ['ions', 'ion', 'proton','protons']):
            V = self.GetVi(time, silent=silent)
        elif (species.lower() in ['electrons', 'electron']):
            V = self.GetVe(time, silent=silent)
        else:
            print "Invalid species name"
            return None

        # vgv_x = vx*dxdvx + vz*dzdvx

        dxdvx[1:-1,:] = 0.5*odx*(V[0,2:,:] - V[0,:-2,:])
        dzdvx[:,1:-1] = 0.5*odz*(V[0,:,2:] - V[0,:,:-2])


        vgv[0, :, :] = V[0,:,:] * dxdvx + V[2, :, :] * dzdvx


        # vgv_y = vx*dxvy + vz*dzvy

        dxdvy[1:-1,:] = 0.5*odx*(V[1,2:,:] - V[1,:-2,:])
        dzdvy[:,1:-1] = 0.5*odz*(V[1,:,2:] - V[1,:,:-2])

        vgv[1, :, :] = V[0, :, :] * dxdvy + V[2, :, :] * dzdvy



        # vgv_z = vx*dxvz + vz*dzvz

        dxdvz[1:-1,:] = 0.5*odx*(V[2,2:,:] - V[2,:-2,:])
        dzdvz[:,1:-1] = 0.5*odz*(V[2,:,2:] - V[2,:,:-2])

        vgv[2, :, :] = V[0, :, :] * dxdvz + V[2, :, :] * dzdvz


        return vgv
    #==========================================================













#============================================================================
#============================================================================
    def divP(self, time, species='ions',silent='No'):

        nx = self.ncells[0]
        nz = self.ncells[1]
        odx = 1./self.dl[0]
        odz = 1./self.dl[1]

        #divp  = np.zeros((3,nx+1,nz+1), "float32", order="FORTRAN")
        P     = np.zeros((3,3,nx+1,nz+1),'float32',order='F')

        if (species.lower() in ['ions','ion','protons','proton']):
            P = self.GetPi(time, silent=silent)
        elif(species.lower() in ['electrons','electron']):
            P = self.GetPe(time, silent=silent)
        else:
            print "Invalid species name"

        divp = self.divT(P)

        return divp
#============================================================================






#============================================================================
#============================================================================
    def divT(self,T):

        nx = self.ncells[0]
        nz = self.ncells[1]
        odx = 1./self.dl[0]
        odz = 1./self.dl[1]

        divT  = np.zeros((3,nx+1,nz+1), "float32", order="FORTRAN")

        dxdTxx = np.zeros((nx+1, nz+1), "float32", order="FORTRAN")
        dzdTxz = np.zeros((nx+1, nz+1), "float32", order="FORTRAN")

        dxdTxy = np.zeros((nx+1, nz+1), "float32", order="FORTRAN")
        dzdTyz = np.zeros((nx+1, nz+1), "float32", order="FORTRAN")

        dxdTxz = np.zeros((nx+1, nz+1), "float32", order="FORTRAN")
        dzdTzz = np.zeros((nx+1, nz+1), "float32", order="FORTRAN")

        Txx = T[0,0,:,:]
        Txz = T[0,2,:,:]
        Txy = T[0,1,:,:]
        Tyz = T[1,2,:,:]
        Tzz = T[2,2,:,:]

        dxdTxx[1:-1,:] = 0.5*odx*( Txx[2:,:] - Txx[:-2,:]  )
        dzdTxz[:,1:-1] = 0.5*odz*( Txz[:,2:] - Txz[:,:-2]  )

        divT[0,:,:] = dxdTxx + dzdTxz


        dxdTxy[1:-1,:] = 0.5*odx*( Txy[2:,:] - Txy[:-2,:]  )
        dzdTyz[:,1:-1] = 0.5*odz*( Tyz[:,2:] - Tyz[:,:-2]  )

        divT[1,:,:] = dxdTxy + dzdTyz



        dxdTxz[1:-1,:] = 0.5*odx*( Txz[2:,:] - Txz[:-2,:]  )
        dzdTzz[:,1:-1] = 0.5*odz*( Tzz[:,2:] - Tzz[:,:-2]  )

        divT[2,:,:] = dxdTxz + dzdTzz


        return divT
#============================================================================






    #==========================================================
    #==========================================================
    def GetParticles(self, time, species,location=None):
        """ returns the particles

        @param time is the (ion gyroperiod) time at which particles are read
        @param species is either 'electrons' or 'ions'
        @param location is the positions at which the particles are selected,
               there can be several locations, like ((xa0,ya0,xa1,ya1), ((xb0,yb0,xb1,yb1)))

        @return a list of particle groups [particles1,particles2,particles3]

        Exemple  : 
        Creation : 2013-01-11 16:44:50.458201

        """

        # checking parameters

        # check the species

        species_elec = ['electrons','electron']
        species_ions = ['protons','proton','ions','ion']

        if species.lower() not in species_elec + species_ions:
            print 'Error, unknow species "%s"' %(species)
            return None


        # check the location
        if location != None:
            if isinstance(location, collections.Iterable) == False:
                print 'Error, location should be [(x0,y0,x1,y1),]'
                return None
            else:
                x0 = []
                x1 = []
                z0 = []
                z1 = []
                for loc in location:
                    if loc[2] <= loc[0]:
                        print 'Error x1 <= x0 (%5.3f <= %5.3f)' % (loc[2],loc[0])
                        return None
                    if loc[3] <= loc[1]:
                        print 'Error z1 <= z0 (%5.3f <= %5.3f)' % (loc[3],loc[1])
                        return None

                    x0.append(loc[0])
                    x1.append(loc[2])
                    z0.append(loc[1])
                    z1.append(loc[3])


        #get the time in wpe^-1 units
        timewpe  = time*self.masses['ions']*self.wpewce

        # build the filename of the first proc file
        filename = os.path.join(self.path,'parts-%05d-p000.dat' %round(timewpe))

        if  os.path.isfile(filename) == False:
            print 'Error - no particle file at that time'
            return None


        # grab the number of particles
        f = open(filename,'rb')
        f.seek(-5*4,2) # 2 means 'from the end', -2*3 means 6 bytes backward
                       # the last byte is the size of file (fortran)
        data = np.fromfile(f, dtype=np.int32, count=5) # number of particles
        nsp  = data[0:2] # size of the particle arrays
        nsact = data[2:4]  # actual number of particles
        f.close()

        # we need to find how many processors there are
        # we proceed with the following method : 
        # first locate all parts-* files
        # then take the first in the list and find its time
        # then look at how many of file at this time there are in the list
        # that's the number of processors...


        # list all particle files
        allfiles   = glob.glob(os.path.join(self.path,'parts-*'))
        if len(allfiles) == 0:
            print 'Error, no particle file found'
            return None

        # take the first
        firstfile  = os.path.basename(allfiles[0])
        twpe       = re.findall(r'[0-9]{5}',firstfile)[0]
        allproc    = glob.glob(os.path.join(self.path,'parts-%s-*' % (twpe)))
        nbproc     = len(allproc)

         # now we want to check whether some times miss some procs
         # so lets take the first proc of all times
        alltimes = glob.glob(os.path.join(self.path,'parts-*-p000.dat'))

        #then for each time grab the time and count the procs
        # they should be equal to nbproc
        for f in alltimes:
            name = os.path.basename(f)
            twpe = re.findall(r'[0-9]{5}',name)[0]
            allproc = glob.glob(os.path.join(self.path,'parts-%s-*' % (twpe)))
            nb      = len(allproc)
            if nb != nbproc:
                print 'time %d does not have all proc files' %(int(twpe))
                return None


        # all right from now on we know how many procs there are
        # and also that all times have all proc files

        # fortran file format :
        # write(30)it,dt,teti,xmax,zmax,nx,nz,x,z,vx,vy,vz,nsp,nsact

        allproc = glob.glob(os.path.join(self.path,'parts-%05d-*' % (round(timewpe))))
        firstproc = True

        offset0 = 4+    4 + 4 + 4  +2*4+2*4

        # these lists will contains the positions of the particles
        # for each location
        xploc  = []
        zploc  = []
        vxploc = []
        vyploc = []
        vzploc = []


        for f in allproc:
            fn = os.path.basename(f)
            fp = open(f,'rb')

            print '%s is now open' % (fn)

            if species.lower() in species_ions:
                print 'reading ions...'
                # pass the stuff before the particle arrays
                offset = offset0
                fp.seek(offset,os.SEEK_SET)
                xp = np.fromfile(fp, dtype=np.float32,count=nsp[0])
                # pass electrons
                offset = offset0 + 4*nsp[0]*2
                fp.seek(offset,os.SEEK_SET)
                zp = np.fromfile(fp, dtype=np.float32,count=nsp[0])
                # pass electrons
                offset = offset0 + 4*nsp[0]*4
                fp.seek(offset,os.SEEK_SET)
                vxp = np.fromfile(fp, dtype=np.float32,count=nsp[0])
                # pass electrons
                offset = offset0 + 4*nsp[0]*6
                fp.seek(offset,os.SEEK_SET)
                vyp = np.fromfile(fp, dtype=np.float32,count=nsp[0])
                # pass electrons
                offset = offset0 + 4*nsp[0]*8
                fp.seek(offset,os.SEEK_SET)
                vzp = np.fromfile(fp, dtype=np.float32,count=nsp[0])


            elif species.lower() in species_elec:
                print "readling electrons..."
                # pass the stuff before the particle arrays
                offset = offset0 + 4*nsp[0]
                fp.seek(offset,os.SEEK_SET)
                xp = np.fromfile(fp, dtype=np.float32,count=nsp[0])
                # pass ions
                offset = offset0 + 4*nsp[0]*3
                fp.seek(offset,os.SEEK_SET)
                zp = np.fromfile(fp, dtype=np.float32,count=nsp[0])
                # pass ions
                offset = offset0 + 4*nsp[0]*5
                fp.seek(offset,os.SEEK_SET)
                vxp = np.fromfile(fp, dtype=np.float32,count=nsp[0])
                # pass ions
                offset = offset0 + 4*nsp[0]*7
                fp.seek(offset,os.SEEK_SET)
                vyp = np.fromfile(fp, dtype=np.float32,count=nsp[0])
                # pass ions
                offset = offset0 + 4*nsp[0]*9
                fp.seek(offset,os.SEEK_SET)
                vzp = np.fromfile(fp, dtype=np.float32,count=nsp[0])

            # then we have to put the coordinates in ion units
            xp /= np.sqrt(self.masses['ions'])
            zp /= np.sqrt(self.masses['ions'])

            # caution : particle position is defined in a box
            # with z=0 at the center, so let's put zp back in [0,zmax]
            zp += self.domsize[1]

            #reset to the first location
            iloc = 0

            # loop on all locations where particles are selected
            for x0i,z0i,x1i,z1i in zip(x0,z0,x1,z1):

                print "looking for particles in (%5.3f,%5.3f,%5.3f,%5.3f) in file %s "% (x0i,z0i,x1i,z1i,fn)

                # ok now we have to take those electrons which are satisfy
                # our criteria, location etc.
                if location != None:
                    idp = np.where((xp > x0i) & (xp < x1i)  \
                                   & (zp > z0i) & (zp < z1i))[0]

                print "ok found %d particles here" %(idp.size)

                # if we read the first processor
                # the arrays are not yet defined
                if firstproc == True:
                    xploc.append(xp[idp])
                    zploc.append(zp[idp])
                    vxploc.append(vxp[idp] * self.wpewce * np.sqrt(self.masses['ions']))
                    vyploc.append(vyp[idp] * self.wpewce * np.sqrt(self.masses['ions']))
                    vzploc.append(vzp[idp] * self.wpewce * np.sqrt(self.masses['ions']))


                # but if it is not the first processor, 
                # the arrays are know already so we want to concatenate
                # new data to them
                else:
                    xploc[iloc]  = np.concatenate((xploc[iloc],xp[idp]))
                    zploc[iloc]  = np.concatenate((zploc[iloc],zp[idp]))
                    vxploc[iloc] = np.concatenate((vxploc[iloc],vxp[idp] * self.wpewce * np.sqrt(self.masses['ions'])))
                    vyploc[iloc] = np.concatenate((vyploc[iloc],vyp[idp] * self.wpewce * np.sqrt(self.masses['ions'])))
                    vzploc[iloc] = np.concatenate((vzploc[iloc],vzp[idp] * self.wpewce * np.sqrt(self.masses['ions'])))

                # next location
                iloc += 1

                # that's it, we've got all particles from all locations for this
                # processor, let's close the file and loop to the next
            firstproc = False
            fp.close()


        # We've got all particles from all location and all processors
        # ok so now let's loop over the locations and group the data into
        # Particles objects

        list_particles = []

        for i in range(len(location)):

            # make a r[2,nbpart] and v[3,nbpart] to fit in Particle constructor
            r_p = np.zeros((2,xploc[i].size), dtype=np.float32)
            v_p = np.zeros((3,vxploc[i].size),dtype=np.float32)

            r_p[0,:] = xploc[i]
            r_p[1,:] = zploc[i]
            v_p[0,:] = vxploc[i]
            v_p[1,:] = vyploc[i]
            v_p[2,:] = vzploc[i]

            list_particles.append(particles.Particles(r_p,v_p,species))


        # ok we got the list now return it
        return list_particles

    #==========================================================















# Luke's version - faster than the C basicsimpementation, works
# but slower than seberg's method
# http://stackoverflow.com/questions/11854522/nested-loop-with-array-indexing-in-numpy
#   def GetFlux(self, time):
#
#       import matplotlib.pyplot as plt
#
#       bx = self.GetField("bx", time) * self.wpewce
#       by = self.GetField("by", time) * self.wpewce
#       bz = self.GetField("bz", time) * self.wpewce
#
#
#
#       flux  = np.zeros((self.ncells[0]+1,self.ncells[1]+1),"float32", order='FORTRAN')
#
#       dx = self.dl[0]
#       dz = self.dl[1]
#
#       nx = self.ncells[0]
#       nz = self.ncells[1]
#
#       j = 0
#
#
#       flux[1:,0] = flux[0,0] + np.cumsum(bz[:-1,0]*dx)
#       flux[0,1:] = flux[0,0] - np.cumsum(bx[0,:-1]*dz)
#
#       a = 0.5
#       aexp = np.arange(nz).reshape(nz, 1) - np.arange(nz).reshape(1, nz)
#       abcoeff = a**aexp
#       abcoeff[aexp<0] = 0
#       for i in np.arange(1,nx+1):
#           b = 0.5*flux[i-1, 1:] + 0.5*bz[i-1, 1:]*dx - 0.5*bx[i,:-1]*dz
#           bvals = (abcoeff * b.reshape(1, nz)).sum(axis=1)
#           n = np.arange(1, nz+1)
#           x0 = flux[i, 0]
#           flux[i, 1:] = a**n * x0 + bvals
#
#       return flux




# BASIC C ALGORITHM - WORKS BUT IS VERY SLOW
#   def GetFlux(self, time):
#
#       import matplotlib.pyplot as plt
#
#       bx = self.GetField("bx", time) * self.wpewce
#       by = self.GetField("by", time) * self.wpewce
#       bz = self.GetField("bz", time) * self.wpewce
#
#
#
#       flux  = np.zeros((self.ncells[0]+1,self.ncells[1]+1),"float32", order='FORTRAN')
#
#       dx = self.dl[0]
#       dz = self.dl[1]
#
#       nx = self.ncells[0]
#       nz = self.ncells[1]
#
#       j = 0
#
#       for i in np.arange(1, nx):
#           flux2[i,0] = flux2[i-1,0] + bz[i-1,0]*dx
#
#
#
#
#       for j in np.arange(1,nz):
#           flux2[0,j] = flux2[0,j-1] - bx[0,j-1]*dz
#
#
#
#       for i in np.arange(1,nx):
#           for j in np.arange(1,nz):
#               flux2[i,j] = 0.5*(flux2[i-1,j] + bz[i-1,j]*dx) + 0.5*(flux2[i,j-1] - bx[i,j-1]*dz)
#
#
#
#       return flux
#==============================================================================================================




# VERSION CYTHON - JAMAIS TESTE

#   def GetFlux2(self, time):
#
#       cimport numpy as np
#       cimport cython
#
#       @cython.bound_checking(False)
#       @cython.wrap_around(False) # probably wrong, but negative indexes are not used...
#       def cythonfunc(self, time):
#
#       # add cdef np.ndarray[...] to tell cython that this is an array...
#       cdef np.ndarray[ndim=2, dtype=np.float64_t] bx = self.GetField("bx", time) * self.wpewce * self.wpewce
#       cdef np.ndarray[ndim=2, dtype=np.float64_t] by = self.GetField("by", time) * self.wpewce * self.wpewce
#       cdef np.ndarray[ndim=2, dtype=np.float64_t] bz = self.GetField("bz", time) * self.wpewce * self.wpewce
#
#       # I am not sure if Cython will honor the fortran order:
#       np.ndarray[ndim=2, dtype=np.float64_t] flux = np.zeros((self.ncells[0]+1,self.ncells[1]+1),"float32", order='FORTRAN')
#
#       # add cdefs, then all math can be handled on C level.
#       cdef double dx = self.dl[0]
#       cdef double dz = self.dl[1]
#
#       # the loop must be done with known ints
#       cdef unsgined int nx = self.ncells[0]
#       cdef unsigned int nz = self.ncells[1]
#
#       cdef unsgined int j = 0
#       cdef unsgined int i = 0
#
#       for i in range(1, nx):
#       flux[i,0] = flux[i-1,0] + bz[i-1,0]*dx
#
#       for j in range(1,nz):
#       flux[0,j] = flux[0,j-1] - bx[0,j-1]*dz
#
#       for i in range(1,nx):
#       for j in range(1,nz):
#       flux[i,j] = 0.5*(flux[i-1,j] + bz[i-1,j]*dx) + 0.5*(flux[i,j-1] - bx[i,j-1]*dz)
#
#       flux[1:, 1:] = np.cumsum(bz[:-1,1:], 0)*dx - np.cumsum(bx[1:,:-1],1)*dz
#
#       print "first row"
#       print flux[0,10:20]
#       print flux2[0,10:20]
#
#       print "first column"
#       print flux[10:20,0]
#       print flux2[10:20,0]
#
#       print "middle"
#       print flux[10:20,40]
#       print flux2[10:20,40]
#
#       return flux




