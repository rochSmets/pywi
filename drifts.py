

import numpy as np
import scipy.ndimage as ndimage
from misctools.tools import vecprod,dotprod


# this module calculates the fluid drift




#----------------------------------------------------------
#----------------------------------------------------------
def diamagnetic(B,divP,n,q,smooth='yes',sigma=2):
    """calculates the diamagnetic drift velocity

    given the magnetic field, divergence of the pressure tensor, density
    and charge, the routine calculates the 2D array of diamagnetic drift
    velocity vector. It applies a gaussian filter on 2points by default.

    @param B magnetic field 2D array
    @param divP divergence of pressure tensor 2D array
    @param n density 2D array
    @param q charge scalar
    @param smooth 'yes' (default) or 'no'
    @param sigma of the gaussian filter in # of points (default 2)

    @return: the 2D array containing the diamagnetic drift velocity vector

    Exemple  : V = diamagnetic(B,divP,n,q,smooth='yes',sigma=2)

    Creation : 2013-01-17 13:52:42.070119

    """
    B2 = B[0,:,:,:]**2 + B[1,:,:,:]**2 + B[2,:,:,:]**2

    #V = np.zeros(B.shape)

    #V[0,:,:,:] = -(divP[1,:,:,:] * B[2,:,:,:] - divP[2,:,:,:] * B[1,:,:,:])/(B2*q*n)
    #V[1,:,:,:] = -(divP[2,:,:,:] * B[0,:,:,:] - divP[0,:,:,:] * B[2,:,:,:])/(B2*q*n)
    #V[2,:,:,:] = -(divP[0,:,:,:] * B[1,:,:,:] - divP[1,:,:,:] * B[0,:,:,:])/(B2*q*n)

    V = - vecprod(divP, B)/(B2*q*n)

    if smooth.lower() == 'yes':
        for c in range(3):
            V[c,:,:] = ndimage.gaussian_filter(V[c,:,:],\
                            sigma=sigma, order=0)


    return V
#==========================================================






#==========================================================
#==========================================================
def ExB(E,B):
    """calculates the E x B drift velocity (2D array)

    @param E (2D array) electric field
    @param B (2D array) magnetic field

    @return: the ExB/B^2 velocity vector (2D array)

    Exemple  :

    Creation : 2013-01-17 13:57:09.900622

    """

   # B2 = B[0,:,:]**2 + B[1,:,:]**2 + B[2,:,:]**2
    B2 = dotprod(B,B)

    #ExB = np.zeros(B.shape)

    #ExB[0,:,:] = (E[1,:,:]*B[2,:,:] - E[2,:,:]*B[1,:,:])/B2
    #ExB[1,:,:] = (E[2,:,:]*B[0,:,:] - E[0,:,:]*B[2,:,:])/B2
    #ExB[2,:,:] = (E[0,:,:]*B[1,:,:] - E[1,:,:]*B[0,:,:])/B2

    return vecprod(E,B)/B2
#==========================================================








#==========================================================
#==========================================================
def Vperp(V,B):
    """returns the component of the velocity V that is perp to B

    @param V velocity vector (3,:,:)
    @param B magnetic field vector (3,:,:)

    @return: velocity perp to B vector (3,:,:)

    Exemple  :

    Note : the perp part of the velocity is obtained from subtracting
    the parallel part from the total velovity.

    Creation : 2013-01-17 13:58:44.295446

    """

    #Bnorm = np.sqrt(B[0,:,:]**2 + B[1,:,:]**2 + B[2,:,:]**2)
    Bnorm = np.sqrt(dotprod(B,B))

    Vperp = np.zeros(B.shape)

    #VdotB  = V[0,:,:]*B[0,:,:] + V[1,:,:]*B[1,:,:] + V[2,:,:]*B[2,:,:]
    VdotB = dotprod(V,B)

    Vpara = np.zeros(B.shape)

    Vpara[0,:,:,:] = VdotB * B[0,:,:,:]/Bnorm
    Vpara[1,:,:,:] = VdotB * B[1,:,:,:]/Bnorm
    Vpara[2,:,:,:] = VdotB * B[2,:,:,:]/Bnorm

    Vperp  = V - Vpara


    return Vperp
#==========================================================


















