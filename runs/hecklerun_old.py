
import glob
import re
import string as strm
import os
from numpy import arange,  fromfile
#from fields import GetFieldData , GetFieldNames
import numpy as np
import scipy.ndimage as ndimage
import sys


#===============================================================================
#===============================================================================
def ReadHF(path):

    """

    Read a hyb.txt file and return its content

    Args  - path : a string indicating the directory where the run is.

    Returns  - the content of the hyb.txt file

    """

    fullpath = os.path.join(path,'hyb.txt')

    f = open(fullpath, "rU")

    content = f.read()
    match = re.findall(r'[\d.]+',content)

    return match
#===============================================================================









class HeckleRun:

    """

    Class representing a Hybrid (Heckle) 2D run



    """

    Name = "HeckleRun"





    #===========================================================================
    #===========================================================================
    def __init__( self, RunID, path):
        """Builds a HeckleRun object

        Longer description here


        @param RunID @todo
        @param path @todo

        @return: @todo

        Exemple  :

        Creation : 2012-08-11 19:01:23.634991

        """

        self.modelname= 'hybrid'
        self.codename = 'heckle'
        self.runid    =  RunID
        self.path     =  path
        params        =  ReadHF(path)
        self.ncells   =  (strm.atoi(params[0]),strm.atoi(params[1]) )
        self.domsize  =  (strm.atof(params[2]),strm.atof(params[3]) )
        self.dl       =  (self.domsize[0]/float(self.ncells[0]),
                  self.domsize[1]/float(self.ncells[1]))
        self.bcs      =  (strm.atoi(params[4]),strm.atoi(params[5]) )
        self.npalloc  =  strm.atoi(params[8])
        self.ts       =  strm.atof(params[9])
        self.nt       =  strm.atoi(params[10])
        self.tsfi     =  strm.atoi(params[11])
        self.tsf      =  self.ts * self.tsfi
        self.tspi     =  strm.atoi(params[12])
        self.tsp      =  self.ts * self.tspi
        self.tsti     =  strm.atoi(params[13])
        self.tst      =  self.ts*self.tsti
        self.rsty     =  strm.atof(params[15])
        self.hprsty   =  strm.atof(params[16])
        self.nptot    =  strm.atoi(params[18])
        self.masses   =  {'ions': (strm.atof(params[20])),
                  'electrons': 0}
        self.charges  =  {'ions' : (strm.atof(params[22])),
                  'electrons':  -(strm.atof(params[22]))}
        self.outofplane = 2

    #===========================================================================




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





    #===========================================================================
    #===========================================================================
    def display(self):
        """Display the run parameters

        Longer description here


        @return: @todo

        Exemple  :

        Creation : 2012-08-11 19:03:28.883915

        """

        print "runID               : ", self.runid
        print "run path            : ", self.path
        print "# of cells           : ", self.ncells
        print "domain size            : ", self.domsize
        print "BCs                   : ", self.bcs
        print "nptot               : ", self.nptot
        print "time step             : ", self.ts
        print "field dump period   : ", self.tsf
        print "part. dump period   : ", self.tsp
        print "diag dump period    : ", self.tst
    #===========================================================================





    #===========================================================================
    #===========================================================================
    def GetOutOfPlaneDir(self):
        """Returns the index of the out of plane dimension

        @return: integer

        Exemple  : run.GetOutOfPlaneDir()

        Creation : 2012-08-11 19:41:27.040976

        """

        return self.outofplane
    #===========================================================================







    #===========================================================================
    #===========================================================================
    def GetNbDiag(self, datatype):
        """Returns the number of diagnostic files produced by a simulation

        @param datatype is either 'fields' or 'particles'

        @return: ineteger

        Exemple  : run.GetNbDiag()

        Creation : 2012-08-11 19:45:58.321810

        """

        if datatype == 'fields':
            return self.nt/self.tsfi
        elif datatype == 'particles':
            return self.nt/self.tspi
        elif datatype == 'diag':
            return self.nt/self.tsti
        else:
            return None
    #===========================================================================





    #===========================================================================
    #===========================================================================
    def GetRunTime(self):
        """Returns the total simulation time for the run

        @return: float

        Exemple  : run.GetRunTime()

        Creation : 2012-08-11 19:49:53.811341

        """

        return self.ts*self.nt
    #===========================================================================




    #===========================================================================
    #===========================================================================
    def GetTimeStep(self, datatype):
        """Returns the time step for the specified datatype

        @param datatype either 'fields' or 'particles'

        @return: float

        Exemple  : run.GetTimeStep()

        Creation : 2012-08-11 19:54:24.236749

        """

        if datatype == 'fields':
            return self.tsf
        elif datatype == 'particles':
            return self.tsp
        elif datatype == 'diag':
            return self.tst
        else:
            return None
    #===========================================================================




    #===========================================================================
    #===========================================================================
    def whereAmI(self):
        """Return the path of the run data directory

        @return: string

        Exemple  : run.whereAmI()

        Creation : 2012-08-11 19:56:44.949280

        """

        return self.path
    #===========================================================================



    #===========================================================================
    #===========================================================================
    def GetNbFieldFiles(self):
        """Returns the number of field files within the run directory

        @return: integer

        Exemple  : print run.GetNbFieldFiles()

        Creation : 2012-08-11 19:58:52.188680

        """

        filenames    = os.listdir(self.path)
        f_filenames = re.findall(r'hf\d+\.dat', strm.join(filenames))
        return len(f_filenames)
    #===========================================================================




    #==========================================================
    #==========================================================
    def GetFieldFiles(self,t1,t2):
        """returns a list of the field files names for which the time
        is between t1 and t2

        @return: file list or None if t1 or t2 is not in the path dir.

        Exemple  :

        Creation : 2012-11-29 16:55:28.247903

        """
        it1 = self.time2iter(t1,'fields')
        it2 = self.time2iter(t2,'fields')
        print it1,it2

        if it1 < 1000:
            file1 = os.path.join(self.path,'hf%03d.dat' %(it1))
        else:
            file1 = os.path.join(self.path,'hf%d.dat' %(it1))

        if it2 < 1000:
            file2 = os.path.join(self.path,'hf%03d.dat' %(it2))
        else:
            file2 = os.path.join(self.path,'hf%d.dat' %(it2))


        if os.path.isfile(file1) == False:
            print 'file %s not found' %(file1)
            return None

        if os.path.isfile(file2) == False:
            print 'file %s not found' % (file2)
            return None

        files = []
        allfiles = glob.glob(os.path.join(self.path,'hf*.dat'))

        for f in allfiles:
            fn = re.findall('hf\d+.dat',f)[0]
            i = int(re.findall('\d+',fn)[0])
            if i >= it1 and i <= it2 :
                files.append(f)

        return sorted(files)

    #==========================================================






    #==========================================================
    #==========================================================
    def GetFileTime(self, filename):
        """returns the time at which this file was written in
        ion cyclotron periods

        Exemple  :  t = run.GetFileTime('fields-05000.dat')

        Creation : 2012-11-30 10:03:51.966027

        """
        fn   = re.findall('hf\d+.dat$',filename)
        it   = int(re.findall('\d+',fn[0])[0])
        time = self.iter2time(it,'fields')

        return time

    #==========================================================









#===============================================================================
#===============================================================================
    def ncelltot(self):
        return self.ncells[0] * self.ncells[1]
#===============================================================================




    # convert iteration number into physical time
    #--------------------------------------------
#==============================================================================
#==============================================================================
    def iter2time(self, it, datatype):
        if datatype == 'fields':
            return it * self.tsf
        elif datatype == 'particles':
            return it * self.tsp
        elif datatype == 'diag':
            return it * self.tst
        else :
            return None
#==============================================================================





    # convert physical time into iter number
    #--------------------------------------------
#==============================================================================
#==============================================================================
    def time2iter(self, time, datatype):
        if datatype == 'fields':
            return int(time / self.tsf)
        elif datatype == 'particles':
            return int(time / self.tsp)
        elif datatype == 'diag':
            return int(time / self.tst)
        else :
            return None
#==============================================================================





#==============================================================================
#==============================================================================
    def indices2coord(self, indices):
        #import operator as op
        coords = np.ndarray(indices.shape)
        coords[0,:] = indices[0,:] * self.dl[0]
        coords[1,:] = indices[1,:] * self.dl[1]
        return coords
#==============================================================================







    #==========================================================
    #==========================================================
    def coord2indices(self, coord):
        """Converts a coordinate array of shape (2,N) into coordinates

        Creation : 2012-08-13 10:05:22.648635

        """

        # np.rint() returns the nearest integer
        indices = np.ndarray(coord.shape)
        indices[0,:] = np.rint(coord[0,:]/self.dl[0])
        indices[1,:] = np.rint(coord[1,:]/self.dl[1])

        return indices
    #==========================================================





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
            return dl*np.arange(self.ncells[coord_idx]+1)
        else:
            npoints = (extent[1]-extent[0])/dl +1
            return extent[0] + dl*np.arange(npoints)

#=================================================================







    # list HF files in the directory
    #--------------------------------------------
#==============================================================================
#==============================================================================
    def GetFieldFiles_old(self, time=None):
        filenames    = os.listdir(self.path)
        hf_filenames = re.findall(r'hf\d+\.dat', strm.join(filenames))

        if time is not None:
            if not isinstance(time,list):
                times = []
                times.append(time)
            else:
                times     = arange(time[0],time[1],self.tsf) # all times 
                                                             # between start 
                                                             # and end
            # all corresponding file IDs
            fileids   = [self.time2iter(t, 'fields') for t in times]  

            # all corresponding filenames
            filenames = ['hf%4d.dat' % fileid for fileid in fileids] 

            # return intersection with ALL filenames in the directory
            return list(set(filenames).intersection(set(hf_filenames))) 
        else:
            return hf_filenames
#==============================================================================








#==============================================================================
#==============================================================================
    def GetFieldNames(self):
        names = ['Bx',
        'By',
        'Bz',
        'Ax',
        'Ay',
        'Az',
        'Ex',
        'Ey',
        'Ez',
        'N',
        'Jx',
        'Jy',
        'Jz',
        'Vix',
        'Viz',
        'Viz',
        'Pexx',
        'Pexy',
        'Pexz' ,
        'Peyy',
        'Peyz',
        'Pezz',
        'Pixx',
        'Pixy',
        'Pixz',
        'Piyy',
        'Piyz',
        'Pizz',
        'Vx0',
        'Vx1',
        'Vy0',
        'Vy1',
        'Vz0',
        'Vz1',
        'dn0',
        'dn1',
        'dnr',
        'qx',
        'qy',
        'qz'
        ]

        return names
#==============================================================================







    # retrieve the field data for time 'time'
    #--------------------------------------------
#==============================================================================
#==============================================================================
    def GetField(self, fieldname, time, silent='No'):

        # get the file number
        it = self.time2iter(time, 'fields')

        # build the filename
        filename = os.path.join(self.path,'hf%03d.dat' %it)

        # check if the file is in the directory
        if not os.path.isfile(filename):
            print "Error - file not found '" + filename+"'"
            return None

        # the file exists
        else:
            data = self.GetFieldData(fieldname, filename, silent=silent)
            return data
#==============================================================================








    # convert a hf file from native heckle format to a VTK format
    #--------------------------------------------
#==============================================================================
#==============================================================================
    def FieldFile2VTK(self):

        filenames = glob.glob(os.path.join(self.whereAmI(),'hf*.dat'))#self.GetFieldFiles(times)
        filenamesvtk = [(re.search(r'hf\d+',f)).group()+'.vtk' 
                        for f in filenames]

        names = self.GetFieldNames()

        for i in range(len(filenames)): # loop over files

            # if the file exists
            if os.path.isfile(os.path.join(self.path,filenames[i])): 
                # make a VTK file
                fnvtk = os.path.join(self.path,filenamesvtk[i])
                fvtk = open(fnvtk, 'wb') #open the VTK file
                print "opening  "+ fnvtk
                n = (self.ncells[0]+1)*(self.ncells[1]+1)            # size

                fvtk.write("# vtk DataFile Version 2.0\n")
                fvtk.write("Heckle fields\n")
                fvtk.write("BINARY\n")
                fvtk.write("DATASET STRUCTURED_POINTS\n")
                fvtk.write("DIMENSIONS %d %d 1\n" % (self.ncells[0]+1, 
                                                    self.ncells[1]+1))
                fvtk.write("ORIGIN  %f %f 0\n" %(0,0))
                fvtk.write("SPACING %f %f 0\n" %(self.dl[0], self.dl[1]))
                fvtk.write("POINT_DATA %d \n" % n)

                #                fvtk.write("SCALARS Bx float\n")
                #                fvtk.write("LOOKUP_TABLE default\n")

                fn =  os.path.join(self.path,filenames[i])
                f   = open(fn, 'rb') # open a binary file in read mode
                print "opening "+fn
                for field in range(39):
                    data = np.fromfile(f, "float32", (self.ncells[0]+1)
                                                      *(self.ncells[1]+1))
                    fvtk.write("SCALARS "+names[field]+" float\n")
                    fvtk.write("LOOKUP_TABLE default\n")
                    data.byteswap(True)
                    data.tofile(fvtk)

                f.close()
                fvtk.close()
#===============================================================================





    # convert a hf file from native heckle format to a VTK format
    #--------------------------------------------
#==============================================================================
#==============================================================================
    def FieldFile2VTK3D(self):

        filenames = glob.glob(os.path.join(self.whereAmI(),'hf*.dat'))#self.GetFieldFiles(times)
        filenamesvtk = [(re.search(r'hf\d+',f)).group()+'.vtk' 
                        for f in filenames]

        names = self.GetFieldNames()

        dz = self.dl[0]
        nz = int(100./dz+1)

        for i in range(len(filenames)): # loop over files

            # if the file exists
            if os.path.isfile(os.path.join(self.path,filenames[i])): 
                # make a VTK file
                fnvtk = os.path.join(self.path,filenamesvtk[i])
                fvtk = open(fnvtk, 'wb') #open the VTK file
                print "opening  "+ fnvtk
                n = (self.ncells[0]+1)*(self.ncells[1]+1)*nz      # size

                fvtk.write("# vtk DataFile Version 2.0\n")
                fvtk.write("Heckle fields\n")
                fvtk.write("BINARY\n")
                fvtk.write("DATASET STRUCTURED_POINTS\n")
                fvtk.write("DIMENSIONS %d %d %d\n" % (self.ncells[0]+1, 
                                                    self.ncells[1]+1,nz))
                fvtk.write("ORIGIN  %f %f 0\n" %(0,0))
                fvtk.write("SPACING %f %f %f\n" %(self.dl[0], self.dl[1],dz))
                fvtk.write("POINT_DATA %d \n" % n)

                #                fvtk.write("SCALARS Bx float\n")
                #                fvtk.write("LOOKUP_TABLE default\n")

                fn =  os.path.join(self.path,filenames[i])
                f   = open(fn, 'rb') # open a binary file in read mode
                print "opening "+fn
                for field in range(3):
                    data = np.fromfile(f, "float32", (self.ncells[0]+1)
                                                      *(self.ncells[1]+1))
                    tmp = data.reshape((self.ncells[0]+1,self.ncells[1]+1))
                    data3d = np.zeros((self.ncells[0]+1,self.ncells[1]+1,nz),'float32')

                    for iz in range(nz):
                        data3d[:,:,iz] = tmp

                    fvtk.write("SCALARS "+names[field]+" float\n")
                    fvtk.write("LOOKUP_TABLE default\n")
                    data.byteswap(True)
                    data.tofile(fvtk)

                f.close()
                fvtk.close()
#===============================================================================




#===============================================================================
#===============================================================================
    def GetFieldData(self, fieldname, filename, silent='No'):


        fieldname = str.lower(fieldname)
        # These values are created
        # when the class is instantiated.
        fieldpos = {}
        fieldpos['bx']    = 0
        fieldpos['by']    = 1
        fieldpos['bz']    = 2
        fieldpos['ax']    = 3
        fieldpos['ay']    = 4
        fieldpos['az']    = 5
        fieldpos['ex']    = 6
        fieldpos['ey']    = 7
        fieldpos['ez']    = 8
        fieldpos['dn']    = 9
        fieldpos['jx']    = 10
        fieldpos['jy']    = 11
        fieldpos['jz']    = 12
        fieldpos['vix']   = 13
        fieldpos['viy']   = 14
        fieldpos['viz']   = 15
        fieldpos['pexx']  = 16
        fieldpos['pexy']  = 17
        fieldpos['peyx']  = 17
        fieldpos['pexz']  = 18
        fieldpos['pezx']  = 18
        fieldpos['peyy']  = 19
        fieldpos['peyz']  = 20
        fieldpos['pezy']  = 20
        fieldpos['pezz']  = 21
        fieldpos['pixx']  = 22
        fieldpos['pixy']  = 23
        fieldpos['piyx']  = 23
        fieldpos['pixz']  = 24
        fieldpos['pizx']  = 24
        fieldpos['piyy']  = 25
        fieldpos['piyz']  = 26
        fieldpos['piyz']  = 26
        fieldpos['pizz']  = 27
        fieldpos['vx0']   = 28
        fieldpos['vx1']   = 29
        fieldpos['vy0']   = 30
        fieldpos['vy1']   = 31
        fieldpos['vz0']   = 32
        fieldpos['vz1']   = 33
        fieldpos['dn0']   = 34
        fieldpos['dn1']   = 35
        fieldpos['dnr']   = 36
        fieldpos['qx']    = 37
        fieldpos['qy']    = 38
        fieldpos['qz']    = 39

        if fieldname not in fieldpos:
            print "Error, '"+fieldname+"' is not a valid fieldname"
            return

        # get the data
        if silent.lower() == 'no':
            print "openning '"+filename+"'... Accessing '"+ fieldname+\
            "' at position ", fieldpos[fieldname],\
            (self.ncells[0]+1)*(self.ncells[1]+1)*4*fieldpos[fieldname]


        f   = open(filename, 'rb') # open a binary file in read mode
        f.seek((self.ncells[0]+1)*(self.ncells[1]+1)*4*fieldpos[fieldname],\
                os.SEEK_SET)

        data = np.fromfile(f, "float32", (self.ncells[0]+1)*(self.ncells[1]+1))
 
        data = data.reshape((self.ncells[0]+1, 
                             self.ncells[1]+1),order='FORTRAN')
        f.close()
        return data
#===============================================================================







#===============================================================================
#===============================================================================
    def GetVe(self, time, silent='No'):
        dn = self.GetField("dn", time, silent=silent)

        jx = self.GetField("jx", time, silent=silent)
        jy = self.GetField("jy", time, silent=silent)
        jz = self.GetField("jz", time, silent=silent)

        vix = self.GetField("vix", time, silent=silent)
        viy = self.GetField("viy", time, silent=silent)
        viz = self.GetField("viz", time, silent=silent)

        vex = vix - jx/dn
        vey = viy - jy/dn
        vez = viz - jz/dn

        ve = np.zeros((3,self.ncells[0]+1,self.ncells[1]+1),
                      "float32", order='FORTRAN')
        ve[0,:,:] = vex
        ve[1,:,:] = vey
        ve[2,:,:] = vez

        return ve
#===============================================================================




#===============================================================================
#===============================================================================
    def GetVi(self, time, silent='No'):

        vix = self.GetField("vix", time, silent=silent)
        viy = self.GetField("viy", time, silent=silent)
        viz = self.GetField("viz", time, silent=silent)

        vi = np.zeros((3,self.ncells[0]+1,self.ncells[1]+1),
                       "float32", order='FORTRAN')
        vi[0,:,:] = vix
        vi[1,:,:] = viy
        vi[2,:,:] = viz

        return vi
#===============================================================================








#==============================================================================================================
#==============================================================================================================
#===============================================================================
#===============================================================================
    def GetV(self, time, specie='ions', silent='No'):

        if (specie.lower() in ['ions', 'ion', 'protons', 'proton']):
            V = self.GetVi(time, silent=silent)

        elif(specie.lower() in ['electrons', 'electron']):
            V = self.GetVe(time, silent=silent)

        else:
            print "Invalid specie name"
            return None

        return V
#===============================================================================
 #==============================================================================================================









#==============================================================================================================
#==============================================================================================================
#===============================================================================
#===============================================================================
    def GetMass(self, specie='ions'):

        """

        Returns the mass of the requested specie

        """
        if specie in ['ions', 'ion', 'protons', 'proton'] :
            return self.masses['ions']

        elif specie in ['electrons', 'electron'] :
            return self.masses['electrons']

        else:
            print "Invalid specie name"
            return None
#===============================================================================
#==============================================================================================================








#==============================================================================================================
#==============================================================================================================
#===============================================================================
#===============================================================================
    def GetCharge(self, specie='ions'):

        """

        Returns the charge of the requested specie

        """
        if specie in ['ions', 'ion', 'protons', 'proton'] :
            return self.charges['ions']

        elif specie in ['electrons', 'electron'] :
            return self.charges['electrons']

        else:
            print "Invalid specie name"
            return None
#===============================================================================
#==============================================================================================================







#==============================================================================================================
#==============================================================================================================
#===============================================================================
#===============================================================================
    def GetB(self, time, silent='No'):

        bx = self.GetField("bx", time, silent=silent)
        by = self.GetField("by", time, silent=silent)
        bz = self.GetField("bz", time, silent=silent)

        B = np.zeros((3,self.ncells[0]+1,self.ncells[1]+1),"float32", order='FORTRAN')
        B[0,:,:] = bx
        B[1,:,:] = by
        B[2,:,:] = bz

        return B
#===============================================================================
#==============================================================================================================




#==============================================================================================================
#==============================================================================================================
#===============================================================================
#===============================================================================
    def GetE(self, time, silent='No'):

        ex = self.GetField("ex", time, silent=silent)
        ey = self.GetField("ey", time, silent=silent)
        ez = self.GetField("ez", time, silent=silent)

        E = np.zeros((3,self.ncells[0]+1,self.ncells[1]+1),"float32", order='FORTRAN')
        E[0,:,:] = ex
        E[1,:,:] = ey
        E[2,:,:] = ez

        return E
#===============================================================================
#===============================================================================
#==============================================================================================================





#===============================================================================
#===============================================================================
    def GetJ(self, time, silent='No'):

        jx = self.GetField("jx", time, silent=silent)
        jy = self.GetField("jy", time, silent=silent)
        jz = self.GetField("jz", time, silent=silent)

        J = np.zeros((3,self.ncells[0]+1,self.ncells[1]+1),
                     "float32", order='FORTRAN')
        J[0,:,:] = jx
        J[1,:,:] = jy
        J[2,:,:] = jz

        return J
#===============================================================================







#===============================================================================
#===============================================================================
    def GetN(self, time, species='ions', silent='No'):

        # since it is a hybrid code, the density is the same for electrons 
        #and ions the prototype however needs the "specie" keyword for 
        #compatibilty with  higher level generic methods


        # no matter what specie is, return the only density we have in hybrid: n
        n = self.GetField( "dn", time, silent=silent)
        return n
#===============================================================================








#===============================================================================
#===============================================================================
    def GetFlux(self, time, silent='No'):
        return self.GetField("az", time, silent=silent)
#===============================================================================








#===============================================================================
#===============================================================================
    def GetPixx(self, time, silent='No'):
        return self.GetField('pixx',time, silent=silent)
#===============================================================================






#===============================================================================
#===============================================================================
    def GetPixy(self, time, silent='No'):
        return self.GetField('pixy', time, silent=silent)
#===============================================================================





#===============================================================================
#===============================================================================
    def GetPixz(self, time, silent='No'):
        return self.GetField('pixz', time, silent=silent)



#==============================================================================================================
#==============================================================================================================
    def GetPiyy(self, time, silent='No'):
        return self.GetField('piyy', time, silent=silent)
#==============================================================================================================






#==============================================================================================================
#==============================================================================================================
    def GetPiyz(self, time, silent='No'):
        return self.GetField('piyz', time, silent=silent)
#==============================================================================================================





#==============================================================================================================
#==============================================================================================================
    def GetPizz(self, time, silent='No'):
        return self.GetField('pizz', time, silent=silent)
#==============================================================================================================





#==============================================================================================================
#==============================================================================================================
    def GetPi(self, time, silent='No'):

        nx = self.ncells[0]
        ny = self.ncells[1]

        Pi = np.zeros((3,3, nx+1, ny+1), "float32", order="FORTRAN")

        pixx = self.GetPixx(time, silent=silent)
        pixy = self.GetPixy(time, silent=silent)
        pixz = self.GetPixz(time, silent=silent)

        piyy = self.GetPiyy(time, silent=silent)
        piyz = self.GetPiyz(time, silent=silent)
        pizz = self.GetPizz(time, silent=silent)


        Pi[0,0,:,:] = pixx
        Pi[0,1,:,:] = pixy
        Pi[0,2,:,:] = pixz

        Pi[1,0,:,:] = pixy
        Pi[1,1,:,:] = piyy
        Pi[1,2,:,:] = piyz

        Pi[2,0,:,:] = pixz
        Pi[2,1,:,:] = piyz
        Pi[2,2,:,:] = pizz


        return Pi
#==============================================================================================================






#===============================================================================
#===============================================================================
    def GetPexx(self, time, silent='No'):
        return self.GetField('pexx',time, silent=silent)
#===============================================================================






#===============================================================================
#===============================================================================
    def GetPexy(self, time, silent='No'):
        return self.GetField('pexy', time, silent=silent)
#===============================================================================





#===============================================================================
#===============================================================================
    def GetPexz(self, time, silent='No'):
        return self.GetField('pexz', time, silent=silent)



#==============================================================================================================
#==============================================================================================================
    def GetPeyy(self, time, silent='No'):
        return self.GetField('peyy', time, silent=silent)
#==============================================================================================================






#==============================================================================================================
#==============================================================================================================
    def GetPeyz(self, time, silent='No'):
        return self.GetField('pexz', time, silent=silent)
#==============================================================================================================





#==============================================================================================================
#==============================================================================================================
    def GetPezz(self, time, silent='No'):
        return self.GetField('pezz', time, silent=silent)
#==============================================================================================================






#==============================================================================================================
#==============================================================================================================
    def GetPe(self, time, silent='No'):

        nx = self.ncells[0]
        ny = self.ncells[1]

        Pe = np.zeros((3,3, nx+1, ny+1), "float32", order="FORTRAN")

        pexx = self.GetPexx(time, silent=silent)
        pexy = self.GetPexy(time, silent=silent)
        pexz = self.GetPexz(time, silent=silent)

        peyy = self.GetPeyy(time, silent=silent)
        peyz = self.GetPeyz(time, silent=silent)
        pezz = self.GetPezz(time, silent=silent)


        Pe[0,0,:,:] = pexx
        Pe[0,1,:,:] = pexy
        Pe[0,2,:,:] = pexz

        Pe[1,0,:,:] = pexy
        Pe[1,1,:,:] = peyy
        Pe[1,2,:,:] = peyz

        Pe[2,0,:,:] = pexz
        Pe[2,1,:,:] = peyz
        Pe[2,2,:,:] = pezz


        return Pe
#==============================================================================================================






#==============================================================================================================
#==============================================================================================================
    def GetPe_old(self, time, silent='No'):

        return self.GetField('pexx', time, silent=silent)
#==============================================================================================================








    #==========================================================
    #==========================================================
    def utherm(self, time, specie='ions'):
        """@todo: Returns the thermal energy of the specified specie

        The thermal energy u is defined as half the trace of the pressure
        tensor : u_s = 0.5*Tr(P_s).

        because the electrons are an isothermal fluid in the hybrid code
        their thermal energy is : u = 1/2 * 3*n*T_e = 3/2*P_e

        @param time is the time in inverse ion pulsation
        @param specie is 'ions' or 'electrons'

        @return: 1/2Trace(P_s)

        Exemple  : 

        Creation : 2012-08-28 11:49:20.052955

        """

        if specie == 'ions':
            Pi = self.GetPi(time)
            return 0.5*(Pi[0,0,:,:]+Pi[1,1,:,:]+Pi[2,2,:,:])

        if specie == 'electrons':
            Pe = self.GetPe(time)
            return 0.5*3.*Pe
    #==========================================================








#==============================================================================================================
#==============================================================================================================
    def GetVxB(self, time, silent='No'):

        Vi = self.GetVi(time, silent=silent)
        B  = self.GetB(time, silent=silent)

        VxB = np.zeros((3,self.ncells[0]+1,self.ncells[1]+1),"float32", order='FORTRAN')

        VxB[0,:,:] = Vi[1,:,:]*B[2,:,:] - Vi[2,:,:]*B[1,:,:]
        VxB[1,:,:] = Vi[2,:,:]*B[0,:,:] - Vi[0,:,:]*B[2,:,:]
        VxB[2,:,:] = Vi[0,:,:]*B[1,:,:] - Vi[1,:,:]*B[2,:,:]

        return VxB
#==============================================================================================================






#==============================================================================================================
#==============================================================================================================
    def GetHallTerm(self, time, silent='No'):

        J = self.GetJ(time, silent=silent)
        n = self.GetN(time, silent=silent)
        B = self.GetB(time, silent=silent)

        H = np.zeros((3,self.ncells[0]+1,self.ncells[1]+1),"float32", order='FORTRAN')

        H[0,:,:] = J[1,:,:]*B[2,:,:] - J[2,:,:]*B[1,:,:]
        H[1,:,:] = J[2,:,:]*B[0,:,:] - J[0,:,:]*B[2,:,:]
        H[2,:,:] = J[0,:,:]*B[1,:,:] - J[1,:,:]*B[2,:,:]

        H /= (n*self.charge['ions'])

        return H
#==============================================================================================================




#============================================================================
#============================================================================
    def divT(self, T):

        divT = np.zeros((3,self.ncells[0]+1,self.ncells[1]+1),\
                        "float32", order='FORTRAN')

        Txx = T[0,0,:,:]
        Txy = T[0,1,:,:]
        Tyy = T[1,1,:,:]
        Tyz = T[1,2,:,:]
        Txz = T[0,2,:,:]


        nx  = self.ncells[0]
        ny  = self.ncells[1]


        odx = 1./self.dl[0]
        ody = 1./self.dl[1]


        dxdTxx = np.zeros((nx+1, ny+1), "float32", order="FORTRAN")
        dydTxy = np.zeros((nx+1, ny+1), "float32", order="FORTRAN")

        dxdTxy = np.zeros((nx+1, ny+1), "float32", order="FORTRAN")
        dydTyy = np.zeros((nx+1, ny+1), "float32", order="FORTRAN")

        dxdTxz = np.zeros((nx+1, ny+1), "float32", order="FORTRAN")
        dydTyz = np.zeros((nx+1, ny+1), "float32", order="FORTRAN")



        dxdTxx[1:-1,:] = 0.5*odx*(Txx[2:,:] - Txx[:-2,:]) #odx*( pxx[1:,:] - pxx[:-1,:]  )
        dydTxy[:,1:-1] = 0.5*ody*(Txy[:,2:] - Txy[:,:-2])  #ody*( pxy[:,1:] - pxy[:,:-1]  )

        divT[0,:,:] = dxdTxx + dydTxy


        dxdTxy[1:-1,:] = 0.5*odx* (Txy[2:,:] - Txy[:-2,:])  #odx*( pxy[1:,:] - pxy[:-1,:]  )
        dydTyy[:,1:-1] = 0.5*ody* (Tyy[:,2:] - Tyy[:,:-2])  #ody*( pyy[:,1:] - pyy[:,:-1]  )

        divT[1,:,:] = dxdTxy + dydTyy


        dxdTxz[1:-1,:] = 0.5*odx* (Txz[2:,:] - Txz[:-2,:])  #odx*( pxz[1:,:] - pxz[:-1,:]  )
        dydTyz[:,1:-1] = 0.5*ody* (Tyz[:,2:] - Tyz[:,:-2])  #ody*( pyz[:,1:] - pyz[:,:-1]  )

        divT[2,:,:] = dxdTxz + dydTyz


        return divT
#============================================================================






#=============================================================================
#=============================================================================
    def divPi(self, time, silent='No'):

        P = self.GetPi(time, silent=silent)
        divp = self.divT(P)

        return divp
#=============================================================================





#==============================================================================================================
#==============================================================================================================
    def divPe(self, time, silent='No'):

        # for compatibility with higher level modules
        # the electron pressure force must have an out-of-plane component
        # even if for 2D runs, it is zero. Higher level modules do not know
        # that and thus expect 3 components
        dpe = np.zeros((3,self.ncells[0]+1,self.ncells[1]+1),"float32", order='FORTRAN')

        pe = self.GetField("pexx", time, silent=silent)

        odx = 1./self.dl[0]
        ody = 1./self.dl[1]

        dpe[0, 1:-1,:] = (pe[2:,:] - pe[:-2,:])*odx*0.5
        dpe[1, :,1:-1] = (pe[:,2:] - pe[:,:-2])*ody*0.5

        return dpe
#==============================================================================================================










#==============================================================================================================
#==============================================================================================================
    def divP(self, time, species='ions', silent='No'):


        if species.lower() in ['ions', 'ion', 'protons', 'proton']:
            dp = self.divPi(time, silent=silent)

        elif species.lower() in ['electrons','electron']:
            dp = self.divPe(time, silent=silent)

        else:
            print "Invalid specie name"
            return None

        return dp
#==============================================================================================================








#==============================================================================================================
#==============================================================================================================
    def VGradV(self, time, specie='ions', silent='No'):

        nx = self.ncells[0]
        ny = self.ncells[1]
        odx = 1./self.dl[0]
        ody = 1./self.dl[1]

        vgv   = np.zeros((3,nx+1,ny+1), "float32", order="FORTRAN")
        dxdvx = np.zeros((nx+1, ny+1), "float32", order="FORTRAN")
        dydvx = np.zeros((nx+1, ny+1), "float32", order="FORTRAN")
        dxdvy = np.zeros((nx+1, ny+1), "float32", order="FORTRAN")
        dydvy = np.zeros((nx+1, ny+1), "float32", order="FORTRAN")
        dxdvz = np.zeros((nx+1, ny+1), "float32", order="FORTRAN")
        dydvz = np.zeros((nx+1, ny+1), "float32", order="FORTRAN")

        if (specie.lower() in ['ions', 'ion', 'proton','protons']):
            V = self.GetVi(time, silent=silent)
        elif (specie.lower() in ['electrons', 'electron']):
            print "Electrons have no mass, this part shouldn't be reached"
            return None
        else:
            print "Invalid specie name"
            return None

        # vgv_x = vx*dxdvx + vy*dydvx

        dxdvx[1:-1,:] = 0.5*odx*(V[0,2:,:] - V[0,:-2,:])
        dydvx[:,1:-1] = 0.5*ody*(V[0,:,2:] - V[0,:,:-2])


        vgv[0, :, :] = V[0,:,:]*dxdvx + V[1, :, :]*dydvx


        # vgv_y = vx*dxvy + vy*dyvy

        dxdvy[1:-1,:] = 0.5*odx*(V[1,2:,:] - V[1,:-2,:])
        dydvy[:,1:-1] = 0.5*ody*(V[1,:,2:] - V[1,:,:-2])

        vgv[1, :, :] = V[0, :, :] * dxdvy + V[1, :, :] * dydvy



        # vgv_z = vx*dxvz + vy*dyvz

        dxdvz[1:-1,:] = 0.5*odx*(V[2,2:,:] - V[2,:-2,:])
        dydvz[:,1:-1] = 0.5*ody*(V[2,:,2:] - V[2,:,:-2])

        vgv[2, :, :] = V[0, :, :] * dxdvz + V[1, :, :] * dydvz


        return vgv
#==============================================================================================================




#==============================================================================================================
#==============================================================================================================
    def GetHyperResistivity(self, time, silent='No'):

        J = self.GetJ(time, silent=silent)

        LapJ = np.zeros((3,self.ncells[0]+1,self.ncells[1]+1),"float32",\
                        order='FORTRAN')

        dx = self.dl[0]
        dy = self.dl[1]
        dx2 = dx*dx
        dy2 = dy*dy


        Jx = J[0,:,:]
        LapJ[0,1:-1, 1:-1] = (Jx[2:,1:-1]-2*Jx[1:-1,1:-1]+Jx[:-2,1:-1])/dx2 \
                           + (Jx[1:-1,2:]-2*Jx[1:-1,1:-1]+Jx[1:-1,:-2])/dy2


        Jy = J[1,:,:]
        LapJ[1,1:-1, 1:-1] = (Jy[2:,1:-1]-2*Jy[1:-1,1:-1]+Jy[:-2,1:-1])/dx2 \
                           + (Jy[1:-1,2:]-2*Jy[1:-1,1:-1]+Jy[1:-1,:-2])/dy2

        #Jz = J[2,:,:]
        LapJ[2,1:-1, 1:-1] = (J[2,2:,1:-1]-2*J[2,1:-1,1:-1]+J[2,:-2,1:-1])/dx2 \
        +(J[2,1:-1,2:]-2*J[2,1:-1,1:-1]+J[2,1:-1,:-2])/dy2

        return -self.hprsty*LapJ
#==============================================================================================================




    #==========================================================
    #==========================================================
    def GetResistivity(self):
        """returns the resistivity for this run

        @return: scalar, float

        Exemple  : 

        Creation : 2012-08-15 14:50:28.674540

        """
        return self.rsty
    #==========================================================





#==============================================================================================================
#==============================================================================================================
    def GetAlfvenSpeed(self,time, silent='No'):
        B = self.GetB(time, silent=silent)
        n = self.GetN(time, silent=silent)
        Bmod = np.sqrt(B[0,:,:]**2 + 0*B[1,:,:]**2 + 0*B[2,:,:]**2)
        #print "shape n", np.shape(n)
        return Bmod/np.sqrt(n)
#==============================================================================================================





#==============================================================================================================
#==============================================================================================================
    def GetExB(self,time, silent='No'):
        B = self.GetB(time, silent=silent)
        E = self.GetE(time, silent=silent)

        B2 = B[0,:,:]**2 + B[1,:,:]**2 + B[2,:,:]**2

        ExB = np.zeros((3,self.ncells[0]+1,self.ncells[1]+1),"float32", order='FORTRAN')

        ExB[0,:,:] = E[1,:,:]*B[2,:,:] - E[2,:,:]*B[1,:,:]
        ExB[1,:,:] = E[2,:,:]*B[0,:,:] - E[0,:,:]*B[2,:,:]
        ExB[2,:,:] = E[0,:,:]*B[1,:,:] - E[1,:,:]*B[2,:,:]

        ExB /= B2

        return ExB
#==============================================================================================================





#==============================================================================================================
#==============================================================================================================
    def GetJdotE(self,time, silent='No'):
        J = self.GetJ(time, silent=silent)
        E = self.GetE(time, silent=silent)


        JdotE = np.zeros((self.ncells[0]+1,self.ncells[1]+1),"float32", order='FORTRAN')

        for h in range(3) :
            JdotE[:,:] +=  E[h,:,:]*J[h,:,:]

        return JdotE
#==============================================================================================================




#==============================================================================================================
#==============================================================================================================
    def GetDivB(self,time, silent='No'):
        from tools import div
        B = self.GetB(time, silent=silent)
        return div(B,self.dl)
#==============================================================================================================



#==============================================================================================================
#==============================================================================================================
    def GetEParallel(self, time, silent='No'):

        B = self.GetB(time, silent=silent)
        E = self.GetE(time, silent=silent)

        epara = np.zeros( (self.ncells[0]+1, self.ncells[1]+1), "float32", order='FORTRAN')

        for h in range(3):
            epara += B[h,:,:]*E[h,:,:]

        epara /= (B[0,:,:]*B[0,:,:]+B[1,:,:]*B[1,:,:]+B[2,:,:]*B[2,:,:])

        return epara
#==============================================================================================================




    #==========================================================
    #==========================================================
    def HasResistivity(self):
        """ tells if resistivty is on


        @return: True if eta != 0, False if eta==0

        Exemple  : run.HasResistivity()

        Creation : 2012-08-15 10:40:42.507696

        """
        if self.rsty == 0:
            return False
        else:
            return True
    #==========================================================





    #==========================================================
    #==========================================================
    def HasHyperResistivity(self):
        """tells if hyper-resistivity is on

        @return: True if on, False otherwise

        Exemple  : run.HasHyperResistivity()

        Creation : 2012-08-15 10:46:14.774664

        """
        if self.hprsty == 0:
            return False
        else:
            return True
    #==========================================================




