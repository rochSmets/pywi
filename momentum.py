import numpy as np
import scipy.ndimage as ndimage
import collections
from pypic.misctools.tools import dotprod,vecprod


#------------------------------------------------------------------------------
#
#  MOMENTUM
#
#  This module allows to calculate the different terms of the momentum equation
#  of a given specie of the plasma.
#
#  Genericity : This module is generic, i.e. does not depend on a particular
#               implementation. It assumes the data layout is a uniform
#               cartesian mesh
#-------------------------------------------------------------------------------



def _fixtime(time):
    """checks whether time is an interval or a single time
       returns a list of one time if time is a scalar"""
    if isinstance(time, collections.Iterable) == False:
        time = [time]

    return time


def _smooth(v, sigma):
    """ apply a gaussian filter on the quantity"""
    for c in range(3):
        v[c,:,:] = ndimage.gaussian_filter(v[c,:,:],sigma=sigma, order=0)

    return v




#----------------------------------------------------------
#----------------------------------------------------------
def electric(run, time, species,smooth='no',sigma=None):
    """@todo: returns the electric force acting on 'species'

    @param run @todo
    @param time @todo
    @param species @todo
    @param smooth @todo
    @param sigma @todo

    @return: @todo

    Exemple  :

    Creation : 2013-01-23 14:50:27.427773

    """
    time   = _fixtime(time)
    q      = run.GetCharge(species)

    tmp = run.GetE(t[0])
    nqE    = np.zeros(shape=tmp.shape, dtype=tmp.dtype)

    for t in time:
        E       = run.GetE(t)
        N       = run.GetN(t, species)
        nqE    += N*E*q

    nqE /= len(time)

    if smooth.lower() == 'yes':
        nqE = _smooth(nqE,sigma)

    return nqE
#==========================================================







#----------------------------------------------------------
#----------------------------------------------------------
def magnetic(run, time, species,smooth='no',sigma=None):
    """@todo: returns the magnetic force acting on 'species'

    @param run @todo
    @param time @todo
    @param species @todo
    @param smooth @todo
    @param sigma @todo

    @return: @todo

    Exemple  : 

    Creation : 2013-01-23 14:50:27.427773

    """
    time   = _fixtime(time)
    q      = run.GetCharge(species)

    tmp    = run.GetE(t[0])
    nqVxB  = np.zeros(shape=tmp.shape, dtype=tmp.dtype)

    for t in time:
        B       = run.GetB(t)
        V       = run.GetV(t, species)
        N       = run.GetN(t, species)

        nqVxB += vecprod(V,B)*N*q
        #nqVxB[0,:,:] += q*N*(V[1,:,:]*B[2,:,:]-V[2,:,:]*B[1,:,:])
        #nqVxB[1,:,:] += q*N*(V[2,:,:]*B[0,:,:]-V[0,:,:]*B[2,:,:])
        #nqVxB[2,:,:] += q*N*(V[0,:,:]*B[1,:,:]-V[1,:,:]*B[0,:,:])

    nqVxB /= len(time)

    if smooth.lower() == 'yes':
        nqVxB = _smooth(nqVxB, sigma)

    return nqVxB
#==========================================================







#----------------------------------------------------------
#----------------------------------------------------------
def pressure(run, time, species,smooth='no',sigma=None):
    """@todo: returns the pressure force acting on 'species'

    @param run @todo
    @param time @todo
    @param species @todo
    @param smooth @todo
    @param sigma @todo

    @return: @todo

    Exemple  : 

    Creation : 2013-01-23 14:50:27.427773

    """
    time   = _fixtime(time)
    tmp    = run.GetE(t[0])
    divP   = np.zeros(shape=tmp.shape, dtype=tmp.dtype)

    for t in time:
        divP    += run.divP(t, species)

    divP /= len(time)

    if smooth.lower() == 'yes':
        divP = _smooth(divP, sigma)

    return -divP
#==========================================================




#----------------------------------------------------------
#----------------------------------------------------------
def inertia_steady(run, time, species,smooth='no',sigma=None):
    """@todo: returns the steady inertia acting on 'species'

    @param run @todo
    @param time @todo
    @param species @todo
    @param smooth @todo
    @param sigma @todo

    @return: @todo

    Exemple  : 

    Creation : 2013-01-23 14:50:27.427773

    """
    time     = _fixtime(time)
    mass     = run.GetMass(species)

    tmp      = run.GetE(t[0])
    inertia  = np.zeros(shape=tmp.shape, dtype=tmp.dtype)

    for t in time:
        N        = run.GetN(t, species)
        inertia += mass*run.VGradV(t, species)

    inertia /= len(time)

    if smooth.lower() == 'yes':
        inertia = _smooth(inertia, sigma)

    return inertia
#==========================================================










#----------------------------------------------------------
#----------------------------------------------------------
def totalforces(run, time, species,smooth='no',sigma=None):
    """@todo: returns the sum of all forces acting on 'species'

    @param run @todo
    @param time @todo
    @param species @todo
    @param smooth @todo
    @param sigma @todo

    @return: @todo

    Exemple  : 

    Creation : 2013-01-23 14:50:27.427773

    """
    fE = electric(run, time, species, smooth=smooth, sigma=sigma)
    fB = magnetic(run, time, species, smooth=smooth, sigma=sigma)
    fP = pressure(run, time, species, smooth=smooth, sigma=sigma)


    return fE + fB + fP
#==========================================================







