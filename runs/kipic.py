
# kiPIC module


import re
import os
import string as strm
import numpy as np
import glob
import collections
import particles

# constants

f64 = np.float64(0)
i32 = np.int32(0)
f64s = f64.itemsize
i32s = i32.itemsize



class KIPIC:
    """
    class for the runs performed with the code kiPIC
    """

    Name = "KIPIC"

#=============================================================================
#=============================================================================
    def __init__(self, RunID, path, silent='No', tmax=80.0):


        self.modelname  = 'fullPIC'
        self.codename   = 'kipic'
        self.runid      = RunID
        self.path       = path

        # read any of the field files and get simulation parameters
        files = glob.glob(os.path.join(path,"field*.dat"))

        f                  = open(files[0], 'rb')
        time               = np.fromfile(f, dtype=np.float64, count=1)
        self.dtsave_fields = np.fromfile(f, dtype=np.float64, count=1)[0]
        ncells_loc         = np.fromfile(f, dtype=np.int32  , count=2)
        self.ncells_glo    = np.fromfile(f, dtype=np.int32  , count=2)
        r0_loc             = np.fromfile(f, dtype=np.float64, count=2)
        rm_loc             = np.fromfile(f, dtype=np.float64, count=2)
        self.r0_glo        = np.fromfile(f, dtype=np.float64, count=2)
        self.rm_glo        = np.fromfile(f, dtype=np.float64, count=2)
        self.c             = np.fromfile(f, dtype=np.float64, count=1)[0]
        f.close()

        dl0 = (self.rm_glo[0] - self.r0_glo[0])/float(self.ncells_glo[0]-2)
        dl1 = (self.rm_glo[1] - self.r0_glo[1])/float(self.ncells_glo[1]-2)

        self.dl = (dl0,dl1)


        # now open one moments file to get some information about
        # the plasma as well

        files = glob.glob(os.path.join(path, "moments*.dat"))

        f                = open(files[0], 'rb')
        f.seek(f64s, os.SEEK_SET)
        self.dtsave_mom  = np.fromfile(f, dtype=np.float64, count=1)[0]
        f.seek( 10*f64s + 4*i32s, os.SEEK_SET)
        self.nspecies    = np.fromfile(f, dtype=np.int32, count=1)[0]


        f.close()








    #==========================================================
    #==========================================================
    def GetE(self, time, YEE=False, silent='No'):
        """ Returns the electric field at requested 'time' """

        # first find the files for that time for all processors
        # then for each processor file :
        #       read the header
        #       read the electric field
        # then concatenate the electric field in one big array

        files = glob.glob(os.path.join(self.path,"field*%011.5f.dat" % time))

        if len(files) == 0:
            print 'Error - no files found for this time (%f)' %time
            return None

        # read the header and parameters


        Eyee  = np.zeros((3,self.ncells_glo[0], self.ncells_glo[1]))

        for ifile,fn in enumerate(files):
            f           = open(fn, 'rb')
            time        = np.fromfile(f, dtype=np.float64, count=1)
            dtsave      = np.fromfile(f, dtype=np.float64, count=1)
            ncells_loc  = np.fromfile(f, dtype=np.int32  , count=2)
            ncells_glo  = np.fromfile(f, dtype=np.int32  , count=2)
            r0_loc      = np.fromfile(f, dtype=np.float64, count=2)
            rm_loc      = np.fromfile(f, dtype=np.float64, count=2)
            r0_glo      = np.fromfile(f, dtype=np.float64, count=2)
            rm_glo      = np.fromfile(f, dtype=np.float64, count=2)
            c           = np.fromfile(f, dtype=np.float64, count=1)

            # shape des tableaux locaux
            size_loc  = ncells_loc[0]*ncells_loc[1]
            shape_loc = (ncells_loc[0], ncells_loc[1])

            # on cree un tableau local
            eloc = np.zeros((3, ncells_loc[0], ncells_loc[1]))

            # dans lequel on stocke ce qu'on lit dans le fichier
            # pour ex, ey et ez
            e0          = np.fromfile(f, dtype=np.float64, count=size_loc)
            e1          = np.fromfile(f, dtype=np.float64, count=size_loc)
            e2          = np.fromfile(f, dtype=np.float64, count=size_loc)

            f.close()

            eloc[0,:,:] = e0.reshape(shape_loc)
            eloc[1,:,:] = e1.reshape(shape_loc)
            eloc[2,:,:] = e2.reshape(shape_loc)


            # maintenant il faut stocker ces tableaux locaux dans
            # un gros tableaux global. il faut calculer les indices
            # de debut et fin dans le tableau global correspondant
            # aux tableaux locaux

            # ces indices vont etre different pour chacune des composantes
            # car elles ne sont pas definies au meme point.


            # pour ez
            # selon x : premier point sur la frontiere gauche
            #           dernier point sur la frontiere droite
            # selon y   premier point sur la frontiere bottom
            #           dernier point sur la frontiere top

            # x0part    = x0glo - dx + ix0*dx
            # x1part    = x1glo - dx + ix1*dx
            # y0part    = y0glo - dy + iy0*dy
            # y1part    = y1glo - dy + iy1*dy

            i0 = np.rint((r0_loc[0] - r0_glo[0])/self.dl[0] + 1)
            i1 = np.rint((rm_loc[0] - r0_glo[0])/self.dl[0] + 1)
            j0 = np.rint((r0_loc[1] - r0_glo[1])/self.dl[1] + 1)
            j1 = np.rint((rm_loc[1] - r0_glo[1])/self.dl[1] + 1)


            # attention on prend pas les i=j=0 pour ezloc car en dehors frontiere
            Eyee[2,i0:i1+1,j0:j1+1] = eloc[2,1:,1:]


            # maintenant ex et by
            # selon x : premier point a frontiere gauche - dx/2
            #           dernier point sur frontiere droite + dx/2
            # selon y : premier et dernier points sur les frontieres bot et top

            # x0part - dx/2 = x0glo - dx + ix0*dx
            # x1part + dx/2 = x0glo - dx + ix1*dx
            # y0part        = x0glo - dx + ix0*dx
            # y1part        = x0glo - dx + ix1*dx

            i0 = np.rint((r0_loc[0] - r0_glo[0])/self.dl[0])
            i1 = np.rint((rm_loc[0] - r0_glo[0])/self.dl[0] + 1)

            Eyee[0,i0:i1+1,j0:j1+1] = eloc[0,0:,1:]


            # maintenant ey et bx
            # selon x : premier et dernier points frontieres gauches et droites
            # selon y : premier et derniers points bot-dy/2 et top+dy/2

            i0 = np.rint((r0_loc[0] - r0_glo[0])/self.dl[0] + 1)
            i1 = np.rint((rm_loc[0] - r0_glo[0])/self.dl[0] + 1)
            j0 = np.rint((r0_loc[1] - r0_glo[1])/self.dl[1])
            j1 = np.rint((rm_loc[1] - r0_glo[1])/self.dl[1] + 1)


            Eyee[1,i0:i1+1,j0:j1+1] = eloc[1,1:,0:]


        if YEE == False:
            E  = np.zeros((3,self.ncells_glo[0],self.ncells_glo[1]),
                          dtype=np.float64)

            # Eyee is defined on the Yee grid
            # E0 is at the same location as B1
            # E1 is at the same location as B0
            # E2 is already at the cell corners

            #   E1(+), E0(^), E2(o)

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


            shape = (self.ncells_glo[0], self.ncells_glo[1])
            E0  = np.zeros(shape, dtype=np.float64)
            E1  = np.zeros(shape, dtype=np.float64)

            E0[1:,:]      = 0.5*( Eyee[0,1:,:] + Eyee[0,:-1,:] )
            E1[:,1:]      = 0.5*( Eyee[1,:,1:] + Eyee[1,:,:-1] )


            E[0,:,:]   =  E0
            E[1,:,:]   =  E1
            E[2,:,:]   =  Eyee[2,:,:]

            return E

        return Eyee
    #==========================================================






    #==========================================================
    #==========================================================
    def GetB(self, time, YEE=False, silent='No'):
        """ Returns the magnetic field at requested 'time' """

        # first find the files for that time for all processors
        # then for each processor file :
        #       read the header
        #       read the magnetic field
        # then concatenate the magnetic field in one big array

        files = glob.glob(os.path.join(self.path,"field*%011.5f.dat" % time))

        if len(files) == 0:
            print 'Error - no files found for this time (%f)' %time
            return None

        # read the header and parameters

        Byee  = np.zeros((3,self.ncells_glo[0], self.ncells_glo[1]))

        for ifile,fn in enumerate(files):
            if silent.lower() == 'no':
                print fn
            f           = open(fn, 'rb')
            time        = np.fromfile(f, dtype=np.float64, count=1)
            dtsave      = np.fromfile(f, dtype=np.float64, count=1)
            ncells_loc  = np.fromfile(f, dtype=np.int32  , count=2)
            ncells_glo  = np.fromfile(f, dtype=np.int32  , count=2)
            r0_loc      = np.fromfile(f, dtype=np.float64, count=2)
            rm_loc      = np.fromfile(f, dtype=np.float64, count=2)
            r0_glo      = np.fromfile(f, dtype=np.float64, count=2)
            rm_glo      = np.fromfile(f, dtype=np.float64, count=2)
            c           = np.fromfile(f, dtype=np.float64, count=1)


            # shape des tableaux locaux
            size_loc  = ncells_loc[0]*ncells_loc[1]
            shape_loc = (ncells_loc[0], ncells_loc[1])

            f.seek(3*np.float64(1).itemsize*size_loc,os.SEEK_CUR)

            # on cree un tableau local
            bloc = np.zeros((3, ncells_loc[0], ncells_loc[1]))

            # dans lequel on stocke ce qu'on lit dans le fichier
            # pour bx, by et bz
            b0          = np.fromfile(f, dtype=np.float64, count=size_loc)
            b1          = np.fromfile(f, dtype=np.float64, count=size_loc)
            b2          = np.fromfile(f, dtype=np.float64, count=size_loc)

            f.close()

            bloc[0,:,:] = b0.reshape(shape_loc)
            bloc[1,:,:] = b1.reshape(shape_loc)
            bloc[2,:,:] = b2.reshape(shape_loc)


            # maintenant il faut stocker ces tableaux locaux dans
            # un gros tableaux global. il faut calculer les indices
            # de debut et fin dans le tableau global correspondant
            # aux tableaux locaux

            # ces indices vont etre different pour chacune des composantes
            # car elles ne sont pas definies au meme point.

            # maintenant by
            # selon x : premier point a frontiere gauche - dx/2
            #           dernier point sur frontiere droite + dx/2
            # selon y : premier et dernier points sur les frontieres bot et top

            # x0part - dx/2 = x0glo - dx + ix0*dx
            # x1part + dx/2 = x0glo - dx + ix1*dx
            # y0part        = x0glo - dx + ix0*dx
            # y1part        = x0glo - dx + ix1*dx

            i0 = np.rint((r0_loc[0] - r0_glo[0])/self.dl[0])
            i1 = np.rint((rm_loc[0] - r0_glo[0])/self.dl[0] + 1)
            j0 = np.rint((r0_loc[1] - r0_glo[1])/self.dl[1] + 1)
            j1 = np.rint((rm_loc[1] - r0_glo[1])/self.dl[1] + 1)


            Byee[1,i0:i1+1,j0:j1+1] = bloc[1,0:,1:]


            # maintenant bx
            # selon x : premier et dernier points frontieres gauches et droites
            # selon y : premier et derniers points bot-dy/2 et top+dy/2

            i0 = np.rint((r0_loc[0] - r0_glo[0])/self.dl[0] + 1)
            i1 = np.rint((rm_loc[0] - r0_glo[0])/self.dl[0] + 1)
            j0 = np.rint((r0_loc[1] - r0_glo[1])/self.dl[1])
            j1 = np.rint((rm_loc[1] - r0_glo[1])/self.dl[1] + 1)

            Byee[0,i0:i1+1,j0:j1+1] = bloc[0,1:,0:]

            # maintenant bz
            # selon x premier et dernier points left - dx/2, right + dx/2
            # selon y premier et dernier points bot  - dy/2,  top  + dz/2

            i0 = np.rint((r0_loc[0] - r0_glo[0])/self.dl[0])
            i1 = np.rint((rm_loc[0] - r0_glo[0])/self.dl[0] + 1.0)
            j0 = np.rint((r0_loc[1] - r0_glo[1])/self.dl[1])
            j1 = np.rint((rm_loc[1] - r0_glo[1])/self.dl[1] + 1.0)

            Byee[2,i0:i1+1,j0:j1+1] = bloc[2,:,:]


        if YEE == False:
            shape = (3,self.ncells_glo[0], self.ncells_glo[1])
            B  = np.zeros(shape, dtype=np.float64)

            # Byee is defined on the Yee grid

            #   B0(+), B1(^), B2(o)

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


            shape = (self.ncells_glo[0],self.ncells_glo[1])
            B0  = np.zeros(shape, dtype=np.float64)
            B1  = np.zeros(shape, dtype=np.float64)

            B1[1:,:]      = 0.5*( Byee[1,1:,:] + Byee[1,:-1,:] )
            B0[:,1:]      = 0.5*( Byee[0,:,1:] + Byee[0,:,:-1] )


            B[0,:,:]   =  B0
            B[1,:,:]   =  B1
            B[2,:,:]   =  Byee[2,:,:]

            return B

        return Byee
    #==========================================================








    #==========================================================
    #==========================================================
    def GetNi(self, time, silent='No'):
        """ Returns the ion density at requested 'time' """

        # first find the files for that time for all processors
        # then for each processor file :
        #       read the header
        #       read the density of all ion species
        # then concatenate the density in one big array

        files = glob.glob(os.path.join(self.path,"moments*%011.5f.dat" % time))

        if len(files) == 0:
            print 'Error - no files found for this time (%f)' %time
            return None

        # read the header and parameters

        Ni  = np.zeros((self.ncells_glo[0], self.ncells_glo[1]))

        for ifile,fn in enumerate(files):
            if silent.lower() == 'no':
                print fn
            f           = open(fn, 'rb')
            f.seek(2*f64s, os.SEEK_SET)
            ncells_loc  = np.fromfile(f, dtype=np.int32  , count=2)
            ncells_glo  = np.fromfile(f, dtype=np.int32  , count=2)
            r0_loc      = np.fromfile(f, dtype=np.float64, count=2)
            rm_loc      = np.fromfile(f, dtype=np.float64, count=2)
            r0_glo      = np.fromfile(f, dtype=np.float64, count=2)
            rm_glo      = np.fromfile(f, dtype=np.float64, count=2)
            nspecies    = np.fromfile(f, dtype=np.int32  , count=1)


            # shape des tableaux locaux
            size_loc  = ncells_loc[0]*ncells_loc[1]
            shape_loc = (ncells_loc[0], ncells_loc[1])

            nloc    = np.zeros(shape_loc, dtype=np.float64)
            mass_s  = 0.

            # now loop on ion species, i.e. the nspecies-1 first ones
            for ispe in np.arange(nspecies-1):
                ispei  = np.fromfile(f, dtype=np.int32, count=1)[0]
                q      = np.fromfile(f, dtype=np.float64, count=1)[0]
                m      = np.fromfile(f, dtype=np.float64, count=1)[0]
                tmp    = np.fromfile(f, dtype=np.float64, count=size_loc)
                nloc  += tmp.reshape(shape_loc)*m
                mass_s += m

            f.close()
            nloc /= mass_s

            # maintenant il faut stocker ces tableaux locaux dans
            # un gros tableaux global. il faut calculer les indices
            # de debut et fin dans le tableau global correspondant
            # aux tableaux locaux

            # selon x : premier point sur la frontiere gauche
            #           dernier point sur la frontiere droite
            # selon y   premier point sur la frontiere bottom
            #           dernier point sur la frontiere top

            # x0part    = x0glo - dx + ix0*dx
            # x1part    = x1glo - dx + ix1*dx
            # y0part    = y0glo - dy + iy0*dy
            # y1part    = y1glo - dy + iy1*dy


            i0 = np.rint((r0_loc[0] - r0_glo[0])/self.dl[0] + 1)
            i1 = np.rint((rm_loc[0] - r0_glo[0])/self.dl[0] + 1)
            j0 = np.rint((r0_loc[1] - r0_glo[1])/self.dl[1] + 1)
            j1 = np.rint((rm_loc[1] - r0_glo[1])/self.dl[1] + 1)

            # attention on prend pas les i=j=0 pour ezloc car en dehors frontiere
            Ni[i0:i1+1, j0:j1+1] = nloc[1:,1:]

        return Ni
    #==========================================================







    #==========================================================
    #==========================================================
    def GetN(self, time, species='protons', silent='No'):
        """ Returns the ion density at requested 'time' """

        # first find the files for that time for all processors
        # then for each processor file :
        #       read the header
        #       read the density of all ion species
        # then concatenate the density in one big array

        files = glob.glob(os.path.join(self.path,"moments*%011.5f.dat" % time))

        if len(files) == 0:
            print 'Error - no files found for this time (%f)' %time
            return None

        # read the header and parameters

        Ns  = np.zeros((self.ncells_glo[0], self.ncells_glo[1]))

        # to be changed later for multispecies, have a dict with species names
        # and index.

        if species.lower() == 'protons':
            ispe = 0
        else:
            ispe = 1


        for ifile,fn in enumerate(files):
            if silent.lower() == 'no':
                print fn
            f           = open(fn, 'rb')
            f.seek(2*f64s, os.SEEK_SET)
            ncells_loc  = np.fromfile(f, dtype=np.int32  , count=2)
            ncells_glo  = np.fromfile(f, dtype=np.int32  , count=2)
            r0_loc      = np.fromfile(f, dtype=np.float64, count=2)
            rm_loc      = np.fromfile(f, dtype=np.float64, count=2)
            r0_glo      = np.fromfile(f, dtype=np.float64, count=2)
            rm_glo      = np.fromfile(f, dtype=np.float64, count=2)
            nspecies    = np.fromfile(f, dtype=np.int32  , count=1)


            # shape des tableaux locaux
            size_loc  = ncells_loc[0]*ncells_loc[1]
            shape_loc = (ncells_loc[0], ncells_loc[1])

            nloc    = np.zeros(shape_loc, dtype=np.float64)
            offset = (i32s  + 2*f64s + size_loc*f64s*4)*ispe
            f.seek(offset, os.SEEK_CUR)


            ispei  = np.fromfile(f, dtype=np.int32, count=1)[0]
            q      = np.fromfile(f, dtype=np.float64, count=1)[0]
            m      = np.fromfile(f, dtype=np.float64, count=1)[0]
            tmp    = np.fromfile(f, dtype=np.float64, count=size_loc)
            nloc   = tmp.reshape(shape_loc)

            f.close()

            # maintenant il faut stocker ces tableaux locaux dans
            # un gros tableaux global. il faut calculer les indices
            # de debut et fin dans le tableau global correspondant
            # aux tableaux locaux

            # selon x : premier point sur la frontiere gauche
            #           dernier point sur la frontiere droite
            # selon y   premier point sur la frontiere bottom
            #           dernier point sur la frontiere top

            # x0part    = x0glo - dx + ix0*dx
            # x1part    = x1glo - dx + ix1*dx
            # y0part    = y0glo - dy + iy0*dy
            # y1part    = y1glo - dy + iy1*dy


            i0 = np.rint((r0_loc[0] - r0_glo[0])/self.dl[0] + 1)
            i1 = np.rint((rm_loc[0] - r0_glo[0])/self.dl[0] + 1)
            j0 = np.rint((r0_loc[1] - r0_glo[1])/self.dl[1] + 1)
            j1 = np.rint((rm_loc[1] - r0_glo[1])/self.dl[1] + 1)



            # attention on prend pas les i=j=0 pour ezloc car en dehors frontiere
            Ns[i0:i1+1, j0:j1+1] = nloc[1:,1:]

        return Ns
    #==========================================================








    #==========================================================
    #==========================================================
    def GetV(self, time, species='protons', silent='No',YEE=False):
        """ Returns the bulk velocity at requested 'time' """

        # first find the files for that time for all processors
        # then for each processor file :
        #       read the header
        #       read the density of all ion species
        # then concatenate the density in one big array

        files = glob.glob(os.path.join(self.path,"moments*%011.5f.dat" % time))

        if len(files) == 0:
            print 'Error - no files found for this time (%f)' %time
            return None

        # read the header and parameters

        Vyee  = np.zeros((3,self.ncells_glo[0], self.ncells_glo[1]))

        # to be changed later for multispecies, have a dict with species names
        # and index.

        if species.lower() == 'protons':
            ispe = 0
        else:
            print 'electrons !'
            ispe = 1


        for ifile,fn in enumerate(files):
            if silent.lower() == 'no':
                print fn
            f           = open(fn, 'rb')
            f.seek(2*f64s, os.SEEK_SET)
            ncells_loc  = np.fromfile(f, dtype=np.int32  , count=2)
            ncells_glo  = np.fromfile(f, dtype=np.int32  , count=2)
            r0_loc      = np.fromfile(f, dtype=np.float64, count=2)
            rm_loc      = np.fromfile(f, dtype=np.float64, count=2)
            r0_glo      = np.fromfile(f, dtype=np.float64, count=2)
            rm_glo      = np.fromfile(f, dtype=np.float64, count=2)
            nspecies    = np.fromfile(f, dtype=np.int32  , count=1)


            # shape des tableaux locaux
            size_loc  = ncells_loc[0]*ncells_loc[1]
            shape_loc = (3,ncells_loc[0], ncells_loc[1])

            vloc    = np.zeros(shape_loc, dtype=np.float64)

            offset = (i32s  + 2*f64s + size_loc*f64s*4)*ispe \
                    + i32s+2*f64s + size_loc*f64s

            f.seek(offset, os.SEEK_CUR)

            v0 = np.fromfile(f, dtype=np.float64, count=size_loc)
            v1 = np.fromfile(f, dtype=np.float64, count=size_loc)
            v2 = np.fromfile(f, dtype=np.float64, count=size_loc)


            f.close()

            vloc[0,:,:] = v0.reshape(shape_loc[1:])
            vloc[1,:,:] = v1.reshape(shape_loc[1:])
            vloc[2,:,:] = v2.reshape(shape_loc[1:])


            # maintenant il faut stocker ces tableaux locaux dans
            # un gros tableaux global. il faut calculer les indices
            # de debut et fin dans le tableau global correspondant
            # aux tableaux locaux

            # ces indices vont etre different pour chacune des composantes
            # car elles ne sont pas definies au meme point.


            # pour vz
            # selon x : premier point sur la frontiere gauche
            #           dernier point sur la frontiere droite
            # selon y   premier point sur la frontiere bottom
            #           dernier point sur la frontiere top

            # x0part    = x0glo - dx + ix0*dx
            # x1part    = x1glo - dx + ix1*dx
            # y0part    = y0glo - dy + iy0*dy
            # y1part    = y1glo - dy + iy1*dy

            i0 = np.rint((r0_loc[0] - r0_glo[0])/self.dl[0] + 1)
            i1 = np.rint((rm_loc[0] - r0_glo[0])/self.dl[0] + 1)
            j0 = np.rint((r0_loc[1] - r0_glo[1])/self.dl[1] + 1)
            j1 = np.rint((rm_loc[1] - r0_glo[1])/self.dl[1] + 1)


            # attention on prend pas les i=j=0 pour ezloc car en dehors frontiere
            Vyee[2,i0:i1+1,j0:j1+1] = vloc[2,1:,1:]


            # maintenant vx
            # selon x : premier point a frontiere gauche - dx/2
            #           dernier point sur frontiere droite + dx/2
            # selon y : premier et dernier points sur les frontieres bot et top

            # x0part - dx/2 = x0glo - dx + ix0*dx
            # x1part + dx/2 = x0glo - dx + ix1*dx
            # y0part        = x0glo - dx + ix0*dx
            # y1part        = x0glo - dx + ix1*dx

            i0 = np.rint((r0_loc[0] - r0_glo[0])/self.dl[0])
            i1 = np.rint((rm_loc[0] - r0_glo[0])/self.dl[0] + 1)

            Vyee[0,i0:i1+1,j0:j1+1] = vloc[0,0:,1:]


            # maintenant vy
            # selon x : premier et dernier points frontieres gauches et droites
            # selon y : premier et derniers points bot-dy/2 et top+dy/2

            i0 = np.rint((r0_loc[0] - r0_glo[0])/self.dl[0] + 1)
            i1 = np.rint((rm_loc[0] - r0_glo[0])/self.dl[0] + 1)
            j0 = np.rint((r0_loc[1] - r0_glo[1])/self.dl[1])
            j1 = np.rint((rm_loc[1] - r0_glo[1])/self.dl[1] + 1)


            Vyee[1,i0:i1+1,j0:j1+1] = vloc[1,1:,0:]


        if YEE == False:
            V = np.zeros((3,self.ncells_glo[0],self.ncells_glo[1]),
                          dtype=np.float64)

            # Vyee is defined on the Yee grid
            # V0 is at the same location as B1
            # V1 is at the same location as B0
            # V2 is already at the cell corners

            #   V1(+), V0(^), V2(o)

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


            shape = (self.ncells_glo[0], self.ncells_glo[1])
            V0  = np.zeros(shape, dtype=np.float64)
            V1  = np.zeros(shape, dtype=np.float64)

            V0[1:,:]      = 0.5*( Vyee[0,1:,:] + Vyee[0,:-1,:] )
            V1[:,1:]      = 0.5*( Vyee[1,:,1:] + Vyee[1,:,:-1] )


            V[0,:,:]   =  V0
            V[1,:,:]   =  V1
            V[2,:,:]   =  Vyee[2,:,:]

            return V

        return Vyee
    #==========================================================










    #==========================================================
    #==========================================================
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
            return -dl + dl*np.arange(self.ncells_glo[coord_idx])
        else:
            npoints = (extent[1]-extent[0])/dl + 1
            return extent[0] + dl*np.arange(npoints)
    #==========================================================


