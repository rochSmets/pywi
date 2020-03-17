

# this module allows to manipulate particles
# TODO : - projection of the velocities in a local B basis
#        - visualization of a distribution function
#             - for which we choose the basis and axis
#             - and also stuff like arrow for mean V and B
#             - do an object Distribution that allows that


import colorp as colp
import numpy as np

#==========================================================
#==========================================================
def where(particle_list, pos):
    """gives the list index corresponding to the location pos

    @return: @todo

    Exemple  : 

    Creation : 2013-01-28 08:30:36.578276

    """
    npart = len(particle_list)
    xmin  = np.zeros(npart)
    xmax  = np.zeros(npart)
    ymin  = np.zeros(npart)
    ymax  = np.zeros(npart)

    x,y = pos

    for  i, p in enumerate(particle_list):
        xmin[i] = np.min(p._pos[0,:])
        xmax[i] = np.max(p._pos[0,:])
        ymin[i] = np.min(p._pos[1,:])
        ymax[i] = np.max(p._pos[1,:])

    id = np.where((x > xmin) & ( x < xmax) & (y < ymax) & (y > ymin))

    return id

#==========================================================





class Particles(object):
    """This class represents a group of particles"""


    #==========================================================
    #==========================================================
    def __init__(self, r, v,species):
        """@todo: to be defined """

        self._pos     = r
        self._vel     = v
        self._species = species
    #==========================================================




    #==========================================================
    #==========================================================
    def getspecies(self):
        """ returns the name of the particle species
        Exemple  : name = particles.getspecies()
        Creation : 2013-01-11 16:34:09.364493
        """
        return self._species
    #==========================================================



    #==========================================================
    #==========================================================
    def __add__(self, part):
        """concatenate two groups of particles

        Exemple  :

        Creation : 2013-05-14 09:22:21.201044

        """

        s1 = self._pos.shape
        s2 = part._pos.shape
        s = (s1[0],s1[1]+s2[1])
        rtot = np.zeros(s)

        s1 = self._vel.shape
        s2 = part._vel.shape
        s = (s1[0],s1[1]+s2[1])
        vtot = np.zeros(s)

        for c in range(2):
            rtot[c,:] = np.concatenate((self._pos[c,:],part._pos[c,:]))
        for c in range(3):
            vtot[c,:] = np.concatenate((self._vel[c,:],part._vel[c,:]))


        spe  = self._species

        return Particles(rtot,vtot,spe)
    #==========================================================




    #==========================================================
    #==========================================================
    def GetArea(self):
        """returns the rectangle surrounding all particles in space

        @return: @todo

        Exemple  : 

        Creation : 2013-05-20 14:44:54.170839

        """
        r = self._pos

        rect = ( (r[0,:].min(), r[0,:].max()),
                 (r[1,:].min(), r[1,:].max()))

        return rect

    #==========================================================



