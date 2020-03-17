import numpy as np
import scipy.linalg as la
import pypic.fields.fields as fields


# this module does stuff with the pressure tensor



#----------------------------------------------------------
#----------------------------------------------------------
def beta(P,B):
    """return the beta  = 1/3*Tr(P)/(B^2/2)

    @param P is the pressure tensor
    @param B is the magnetic field

    @return: @todo

    Exemple  :

    Creation : 2012-08-28 12:37:47.132074

    """
    trP = fields.trace(P)
    B2  = fields.norm(B)**2


    return 2./3.*trP/(B2)
#==========================================================






#----------------------------------------------------------
#----------------------------------------------------------
def nongyro_decomposition(P,B):
    """Decomposition of the pressure tensor in a gyrotropic
    and a non-gyrotropic part

    @param P tensor
    @param B magnetic field

    @return: returns a tuple (Png,Pg,Ppara,Pperp) where Png and Pg are the
    non-gyrotropic and gyrotropic parts of P, respectively and Ppara and Pperp
    are the parallel and perpendicular components, repsectively

    Exemple  :

    Note : the non gyrotropic part ot the tensor gives the same
    results as Michael Hesse's IDL code.

    Creation : 2012-08-24 18:23:03.290921

    """


    # have some names for components
    # makes it easier to read complicated formulae
    Bx = B[...,0]
    By = B[...,1]
    Bz = B[...,2]

    B2 = fields.norm(B)**2

    # do the same for pressure tensor components
    Pxx = P[...,0,0]
    Pxy = P[...,0,1]
    Pxz = P[...,0,2]
    Pyy = P[...,1,1]
    Pyz = P[...,1,2]
    Pzz = P[...,2,2]


    # bb^T P bb == P_in_b
    # and the first component is P_\parallel


    P_para = Bx**2*Pxx + Bz**2*Pzz + By**2*Pyy + 2*Bx*By*Pxy + 2*Bx*Bz*Pxz + \
             2*By*Bz*Pyz

    P_para /= B2 #it is bb not BB so divde by B^2


    # ok so we now have P_para, and we now that the trace is 
    # invariant through change of basis, so Trace(P) is actually
    # equal to P_para + 2*P_perp, or at least that is a definition of P_perp

    P_perp = 0.5*(Pxx+Pyy+Pzz - P_para)



    # ok so now the general form of P in B basis is : 
    # P = P_ng + P_perp*I + (P_para-P_perp)bb where bb is the tensor BB/B2
    # and I the identity.
    # so we are going to get the non gyrotropic from that formula

    P_pp = np.zeros(P.shape, 'float32',order='F')

    P_pp[...,0,0] = (P_perp + (P_para-P_perp)*Bx*Bx/B2)
    P_pp[...,0,1] = (0      + (P_para-P_perp)*Bx*By/B2)
    P_pp[...,1,0] = P_pp[...,0,1]
    P_pp[...,0,2] = (0      + (P_para-P_perp)*Bx*Bz/B2)
    P_pp[...,2,0] = P_pp[...,0,2]
    P_pp[...,1,1] = (P_perp + (P_para-P_perp)*By*By/B2)
    P_pp[...,1,2] = (0      + (P_para-P_perp)*By*Bz/B2)
    P_pp[...,2,1] = P_pp[...,1,2]
    P_pp[...,2,2] = (P_perp + (P_para-P_perp)*Bz*Bz/B2)


    Pxx_ng = Pxx - (P_perp + (P_para-P_perp)*Bx*Bx/B2)
    Pxy_ng = Pxy - (0      + (P_para-P_perp)*Bx*By/B2)
    Pxz_ng = Pxz - (0      + (P_para-P_perp)*Bx*Bz/B2)
    Pyy_ng = Pyy - (P_perp + (P_para-P_perp)*By*By/B2)
    Pyz_ng = Pyz - (0      + (P_para-P_perp)*By*Bz/B2)
    Pzz_ng = Pzz - (P_perp + (P_para-P_perp)*Bz*Bz/B2)



    P_ng = P - P_pp


    return (P_ng, P_pp,P_para, P_perp)

#==========================================================






#----------------------------------------------------------
#----------------------------------------------------------
def paraperp(P,B):
    """returns the parallel and perpendicular pressure

    @param P the pressure tensor
    @param B the magnetic field

    Exemple  : para,perp = paraperp(P,B)

    Creation : 2012-08-25 11:16:34.493568

    """


    # have some names for components
    # makes it easier to read complicated formulae
    Bx = B[...,0]
    By = B[...,1]
    Bz = B[...,2]

    B2 = fields.norm(B)**2

    # do the same for pressure tensor components
    Pxx = P[...,0,0]
    Pxy = P[...,0,1]
    Pxz = P[...,0,2]
    Pyy = P[...,1,1]
    Pyz = P[...,1,2]
    Pzz = P[...,2,2]


    # bb^T P bb == P_in_b
    # and the first component is P_\parallel


    P_para = Bx**2*Pxx + Bz**2*Pzz + By**2*Pyy + 2*Bx*By*Pxy + 2*Bx*Bz*Pxz + \
             2*By*Bz*Pyz

    P_para /= B2 #it is bb not BB so divde by B^2


    # ok so we now have P_para, and we now that the trace is 
    # invariant through change of basis, so Trace(P) is actually
    # equal to P_para + 2*P_perp, or at least that is a definition of P_perp

    P_perp = 0.5*(Pxx+Pyy+Pzz - P_para)



    return  (P_para, P_perp)
#==========================================================





#----------------------------------------------------------
#----------------------------------------------------------
def anisotropy(P, B):
    """ returns 3*abs(P_perp - P_para)/(2*P_perp+P_para)

    Exemple  :

    Creation : 2012-08-25 11:03:17.253041

    """
    P_para, P_perp = paraperp(P,B)

    return 2.*np.abs(P_perp-P_para)/(P_perp+P_para)
#==========================================================








#----------------------------------------------------------
#----------------------------------------------------------
def nongyrotropy2(P,B):
    """returns a non-gyrotropy coefficient

    @param P is the pressure tensor in x,y,z
    @param B is the magnetic field

    @return: it returns the norm of the eigenvalues of the non-gyrotropic par
    of the pressure tensor.

    Exemple  :

    Creation : 2012-08-25 11:39:18.370158

    """

    # gets the non gyrotropic part of P and the perp/para components
    P_ng, P_pp, P_para, P_perp = nongyro_decomposition(P,B)


    ev = _cubicfromPs(P_ng) # finds the eigenvalues of the pressure tensor

    # calculates the Frobenius norm of the nongyrotropic part
    ng = np.sqrt(ev[0]**2 + ev[1]**2 + ev[2]**2)

    # and normalize by the thermal energy
    ng = ng *2./fields.trace(P)

    return ng
#==========================================================





#----------------------------------------------------------
#----------------------------------------------------------
def nongyrotropy(P,B):
    """returns a non-gyrotropy coefficient

    @param P is the pressure tensor in x,y,z
    @param B is the magnetic field

    @return: it returns the norm of the eigenvalues of the non-gyrotropic par
    of the pressure tensor.

    Exemple  :

    Creation : 2012-08-25 11:39:18.370158

    """

    # gets the non gyrotropic part of P and the perp/para components
    P_ng, P_pp, P_para, P_perp = nongyro_decomposition(P,B)

    ng = 0
    for i in range(3):
        for j in range(3):
            ng += P_ng[i,j,:,:]*P_ng[i,j,:,:]

    ng = np.sqrt(ng)

    # and normalize by the thermal energy
    ng = ng *2./fields.trace(P)

    return ng
#==========================================================






# --------------PRIVATE FUNCTIONS ----------------------


# cubic : finds the three real roots of a cubic polynom
def _cubic(a,b,c,d):
    pi2 = np.pi*2
    oth = 1./3.

    a1 = b/a
    a2 = c/a
    a3 = d/a

    q = (a1**2 -3*a2)/9.
    r = (2*a1**3 - 9*a1*a2 + 27*a3)/54.

    alpha = r/q**1.5

    id = np.where(np.abs(alpha)>1)
    alpha[id[0],id[1]]  = 1

    theta = np.arccos(alpha)
    x1 = -2 * np.sqrt(q) * np.cos(theta/3.) - a1/3.
    x2 = -2.* np.sqrt(q) * np.cos((theta+ pi2)/3.)-a1/3.
    x3 = -2.*np.sqrt(q) * np.cos((theta-pi2)/3.)-a1/3.

    return (x1,x2,x3)




# this function setup the coef of the cubic polynom
# from the pressure tensor.
def _cubicfromPs(P):

    a = np.ones(P[0,0,:,:].shape)*-1

    b = P[...,0,0]+P[...,1,1]+P[...,2,2]

    c = - P[...,0,0]*P[...,1,1] - P[...,1,1]*P[...,2,2] - \
          P[...,0,0]*P[...,2,2] + P[...,0,2]**2 +         \
          P[...,1,2]**2 + P[...,1,0]**2

    d = P[...,0,1]*P[...,1,2]*P[...,0,2] + \
        P[...,0,1]*P[...,1,2]*P[...,0,2] - \
        P[...,0,2]**2*P[...,1,1] -         \
        P[...,1,2]**2*P[...,0,0] -         \
        P[...,1,0]**2*P[...,2,2]

    return _cubic(a,b,c,d)






