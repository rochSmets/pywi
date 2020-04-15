
import collections
from scipy.ndimage import gaussian_filter1d as gf
import numpy as np
import scipy.ndimage as ndimage
import fields.fields as fields


# ______________________________________________________________________________
#
# checks whether 'time' is one scalar or a list of time (floats)
# if it is a list, the Ohm's term is averaged on the times of this list
# ______________________________________________________________________________



def drift(run, time, specie, comp, terms = 'all', smooth = 'no', sigma = 1):

    """
    returns the drift term of the Ohm's law, for the given specie
    at time if time is float or averaged on the list of time if time is a list (of floats)

    @return  : np.ndarray of shape run.vfield_shape

    """


    if isinstance(time, collections.Iterable) == False:
        time = [time]
        ntot = 1
    else:
        ntot = time.size

    term = np.zeros(shape = run.sfield_shape, dtype = run.dtype)

    for t in time:
        if specie.lower() == 'ion'      or specie.lower() == 'ions'      or specie.lower == 'i' or\
           specie.lower() == 'electron' or specie.lower() == 'electrons' or specie.lower == 'e' :
            Vx = run.getV(t, specie, 'x')
            Vy = run.getV(t, specie, 'y')
            Vz = run.getV(t, specie, 'z')
        else :
            raise ValueError('specie can be "e", "electron", "electrons", "i", "ion", "ions"')

        Bx = run.getB(t, 'x')
        By = run.getB(t, 'y')
        Bz = run.getB(t, 'z')

        if comp == 'x' or comp == 'X' :
            if terms == 'all' :
                term += Vz*By-Vy*Bz
            elif terms == 'VyBz':
                term += -Vy*Bz
            elif terms == 'VzBy':
                term += Vz*By
            else :
                raise ValueError('for X component of ideal term, terms can be "all", "VyBz" or "VzBy"')

        elif comp == 'y' or comp == 'Y' :
            if terms == 'all' :
                term += Vx*Bz-Vz*Bx
            elif terms == 'VzBx':
                term += -Vz*Bx
            elif terms == 'VxBz':
                term += Vx*Bz
            else :
                raise ValueError('for Y component of ideal term, terms can be "all", "VzBx" or "VxBz"')

        elif comp == 'z' or comp == 'Z' :
            if terms == 'all' :
                term += Vy*Bx-Vx*By
            elif terms == 'VxBy':
                term += -Vz*Bx
            elif terms == 'VyBx':
                term += Vx*Bz
            else :
                raise ValueError('for Z component of ideal term, terms can be "all", "VxBy" or "VyBx"')

        else :
            raise ValueError('comp can only be "x", "y" or "z"')

    term /= ntot

    if smooth.lower() == 'yes':
        term = gf(term, sigma = sigma, order = 0)


    return term


def hall(run, time, comp, terms = 'all', smooth = 'no', sigma = 1):

    """
    returns the Hall term of the Ohm's law,
    at time if time is float or averaged on the list of time if time is a list (of floats)

    @return  : np.ndarray of shape run.vfield_shape

    """

    if isinstance(time, collections.Iterable) == False:
        time = [time]
        ntot = 1
    else:
        ntot = time.size

    term = np.zeros(shape = run.sfield_shape, dtype = run.dtype)

    for t in time:
        Jx = run.getJ(t, 'x')
        Jy = run.getJ(t, 'y')
        Jz = run.getJ(t, 'z')

        Bx = run.getB(t, 'x')
        By = run.getB(t, 'y')
        Bz = run.getB(t, 'z')

        n  = run.getN(t, 'electron')

        if comp.lower() == 'x' :
            if terms == 'all' :
                term += (Jy*Bz-Jz*By)/n
            elif terms == 'JyBz':
                term += (Jy*Bz)/n
            elif terms == 'JzBy':
                term += -(Jz*By)/n
            else :
                raise ValueError('for X component of Hall term, "terms" can be "all", "JyBz" or "JzBy"')

        elif comp.lower() == 'y' :
            if terms == 'all' :
                term += (Jz*Bx-Jx*Bz)/n
            elif terms == 'JzBx':
                term += (Jz*Bx)/n
            elif terms == 'JxBz':
                term += -(Jx*Bz)/n
            else :
                raise ValueError('for Y component of Hall term, "terms" can be "all", "JzBx" or "JxBz"')

        elif comp.lower() == 'z' :
            if terms == 'all' :
                term += (Jx*By-Jy*Bx)/n
            elif terms == 'JxBy':
                term += (Jx*By)/n
            elif terms == 'JyBx':
                term += -(Jy*Bx)/n
            else :
                raise ValueError('for Z component of Hall term, "terms" can be "all", "JxBy" or "JyBx"')
        else :
            raise ValueError('comp can only be "x", "y" or "z"')

    term /= ntot

    if smooth.lower() == 'yes':
        term = gf(term, sigma = sigma, order = 0)

    return term


def gradP(run, time, specie, comp, terms = 'all', smooth = 'no', sigma = 1):

    """
    returns the gradP term of the Ohm's law, for the given specie
    at time if time is float or averaged on the list of time if time is a list (of floats)

    @return  : np.ndarray of shape run.vfield_shape

    """

    if isinstance(time, collections.Iterable) == False:
        time = [time]
        ntot = 1
    else:
        ntot = time.size

    term = np.zeros(shape = run.sfield_shape, dtype = run.dtype)

    for t in time:
        n = run.getN(t, specie)

        Pxx = run.getP(t, specie, 'xx')
        Pxy = run.getP(t, specie, 'xy')
        Pxz = run.getP(t, specie, 'xz')
        Pyy = run.getP(t, specie, 'yy')
        Pyz = run.getP(t, specie, 'yz')
        Pzz = run.getP(t, specie, 'zz')

        if run.ndim == 1 :
            dxPxx = fields.deriv(Pxx, run.dl, order = 1, axis = 0)
            dxPxy = fields.deriv(Pxy, run.dl, order = 1, axis = 0)
            dxPxz = fields.deriv(Pxz, run.dl, order = 1, axis = 0)

            if comp.lower() == 'x' :
                if terms == 'all' :
                    term += dxPxx/n
                elif terms == 'Pxx' :
                    term += dxPxx/n
                else :
                    raise ValueError('for "x" comp in 1d run, terms can only be "all" or "Pxx"')

            elif comp.lower() == 'y' :
                if terms == 'all' :
                    term += dxPxy/n
                elif terms == 'Pxy' :
                    term += dxPxy/n
                else :
                    raise ValueError('for "y" comp in 1d run, terms can only be "all" or "Pxy"')

            elif comp.lower() == 'z' :
                if terms == 'all' :
                    term += dxPxz/n
                elif terms == 'Pxz' :
                    term += dxPxz/n
                else :
                    raise ValueError('for "y" comp in 1d run, terms can only be "all" or "Pxz"')

            else :
                raise ValueError('comp can only be "x", "y" or "z"')

        if run.ndim == 2 :
            dxPxx = fields.deriv(Pxx, run.dl, order = 1, axis = 0)
            dxPxy = fields.deriv(Pxy, run.dl, order = 1, axis = 0)
            dxPxz = fields.deriv(Pxz, run.dl, order = 1, axis = 0)

            dyPxy = fields.deriv(Pxy, run.dl, order = 1, axis = 1)
            dyPyy = fields.deriv(Pyy, run.dl, order = 1, axis = 1)
            dyPyz = fields.deriv(Pyz, run.dl, order = 1, axis = 1)

            if comp.lower() == 'x' :
                if terms == 'all' :
                    term += (dxPxx+dyPxy)/n
                elif terms == 'Pxx' :
                    term += dxPxx/n
                elif terms == 'Pxy':
                    term += dyPxy / n
                else :
                    raise ValueError('for "x" comp in 1d run, terms can only be "all", "Pxx" or "Pxy"')

            elif comp.lower() == 'y' :
                if terms == 'all' :
                    term += (dxPxy+dyPyy)/n
                elif terms == 'Pxy' :
                    term += dxPxy/n
                elif terms == 'Pyy':
                    term += dyPyy / n
                else :
                    raise ValueError('for "y" comp in 1d run, terms can only be "all", "Pxy" or "Pyy"')

            elif comp.lower() == 'z' :
                if terms == 'all' :
                    term += (dxPxz+dyPyz)/n
                elif terms == 'Pxz' :
                    term += dxPxz/n
                elif terms == 'Pyz':
                    term += dyPyz / n
                else :
                    raise ValueError('for "z" comp in 1d run, terms can only be "all", "Pxz" or "Pyz"')

            else :
                raise ValueError('comp can only be "x", "y" or "z"')

        if run.ndim == 3 :
            dxPxx = fields.deriv(Pxx, run.dl, order = 1, axis = 0)
            dxPxy = fields.deriv(Pxy, run.dl, order = 1, axis = 0)
            dxPxz = fields.deriv(Pxz, run.dl, order = 1, axis = 0)

            dyPxy = fields.deriv(Pxy, run.dl, order = 1, axis = 1)
            dyPyy = fields.deriv(Pyy, run.dl, order = 1, axis = 1)
            dyPyz = fields.deriv(Pyz, run.dl, order = 1, axis = 1)

            dzPxz = fields.deriv(Pxz, run.dl, order = 1, axis = 2)
            dzPyz = fields.deriv(Pyz, run.dl, order = 1, axis = 2)
            dzPzz = fields.deriv(Pzz, run.dl, order = 1, axis = 2)

            if comp.lower() == 'x' :
                if terms == 'all' :
                    term += (dxPxx+dyPxy+dzPxz)/n
                elif terms == 'Pxx' :
                    term += dxPxx/n
                elif terms == 'Pxy' :
                    term += dyPxy/n
                elif terms == 'Pxz' :
                    term += dzPxz/n

            elif comp.lower() == 'y' :
                if terms == 'all' :
                    term += (dxPxy+dyPyy+dzPyz)/n
                elif terms == 'Pxy' :
                    term += dxPxy/n
                elif terms == 'Pyy' :
                    term += dyPyy/n
                elif terms == 'Pyz' :
                    term += dzPyz/n

            elif comp.lower() == 'z' :
                if terms == 'all' :
                    term += (dxPxz+dyPyz+dzPzz)/n
                elif terms == 'Pxz' :
                    term += dxPxz/n
                elif terms == 'Pyz' :
                    term += dyPyz/n
                elif terms == 'Pzz' :
                    term += dzPzz/n

            else :
                raise ValueError('comp can only be "x", "y" or "z"')

    term /= ntot
    term /= run.getCharge('electrons')

    if smooth.lower() == 'yes':
        for c in range(3):
            term[...,c] = gf(term[..., c], sigma = sigma, order = 0)


    return term


def resistive(run, time, comp, smooth='no',sigma=1):

    """
    returns the resistive term of the Ohm's law,
    at time if time is float or averaged on the list of time if time is a list (of floats)

    @return  : np.ndarray of shape run.vfield_shape

    """

    if isinstance(time, collections.Iterable) == False:
        time = [time]
        ntot = 1
    else:
        ntot = time.size

    term = np.zeros(shape = run.sfield_shape, dtype = run.dtype)

    for t in time:

        term += run.getJ(t, comp)

    term *= run.getResistivity()/ntot

    if smooth.lower() == 'yes':
        term = gf(term, sigma = sigma, order = 0)

    return term


def hyper(run, time, comp, terms = 'all', smooth='no', sigma=1):

    """
    returns the hyper-resistive term of the Ohm's law,
    at time if time is float or averaged on the list of time if time is a list (of floats)

    @return  : np.ndarray of shape run.vfield_shape

    """


    if isinstance(time, collections.Iterable) == False:
        time = [time]
        ntot = 1
    else:
        ntot = time.size

    term = np.zeros(shape = run.sfield_shape, dtype = run.dtype)

    for t in time:

        if run.ndim == 1:
            if comp.lower() == 'x' :
                Jx = run.getJ(t, 'x')
                if terms == 'all' :
                    d2xJx = fields.deriv(Jx, dl, 2, 0)
                    term += d2xJx
                elif terms == 'd2xJx' :
                    d2xJx = fields.deriv(Jx, dl, 2, 0)
                    term += d2xJx
                else :
                    raise ValueError('for X component of Hyper-resistive term, "terms" can be "all" or "d2xJx"')

            elif comp.lower() == 'y' :
                Jy = run.getJ(t, 'y')
                if terms == 'all' :
                    d2xJy = fields.deriv(Jy, dl, 2, 0)
                    term += d2xJy
                elif terms == 'd2xJy' :
                    d2xJy = fields.deriv(Jy, dl, 2, 0)
                    term += d2xJy
                else :
                    raise ValueError('for Y component of Hyper-resistive term, "terms" can be "all" or "d2xJy"')

            elif comp.lower() == 'z' :
                Jz = run.getJ(t, 'z')
                if terms == 'all' :
                    d2xJz = fields.deriv(Jz, dl, 2, 0)
                    term += d2xJz
                elif terms == 'd2xJz' :
                    d2xJz = fields.deriv(Jz, dl, 2, 0)
                    term += d2xJz
                else :
                    raise ValueError('for Z component of Hyper-resistive term, "terms" can be "all" or "d2xJz"')

            else :
                raise ValueError('comp can only be "x", "y" or "z"')

        if run.ndim == 2:
            if comp.lower() == 'x' :
                Jx = run.getJ(t, 'x')
                if terms == 'all' :
                    d2xJx = fields.deriv(Jx, dl, 2, 0)
                    d2yJx = fields.deriv(Jx, dl, 2, 1)
                    term += d2xJx+d2yJx
                elif terms == 'd2xJx' :
                    d2xJx = fields.deriv(Jx, dl, 2, 0)
                    term += d2xJx
                elif terms == 'd2yJx' :
                    d2yJx = fields.deriv(Jx, dl, 2, 1)
                    term += d2yJx
                else :
                    raise ValueError('for X component of Hyper-resistive term, "terms" can be "all", "d2xJx" or "d2yjJx"')

            elif comp.lower() == 'y' :
                Jy = run.getJ(t, 'y')
                if terms == 'all' :
                    d2xJy = fields.deriv(Jy, run.dl, 2, 0)
                    d2yJy = fields.deriv(Jy, run.dl, 2, 1)
                    term += d2xJy+d2yJy
                elif terms == 'd2xJy' :
                    d2xJy = fields.deriv(Jy, dl, 2, 0)
                    term += d2xJy
                elif terms == 'd2yJy' :
                    d2yJy = fields.deriv(Jy, dl, 2, 1)
                    term += d2yJy
                else :
                    raise ValueError('for Y component of Hyper-resistive term, "terms" can be "all", "d2xJy", "d2yjJy" or "d2zJy"')

            elif comp.lower() == 'z' :
                Jz = run.getJ(t, 'z')
                if terms == 'all' :
                    d2xJz = fields.deriv(Jz, run.dl, 2, 0)
                    d2yJz = fields.deriv(Jz, run.dl, 2, 1)
                    term += d2xJz+d2yJz
                elif terms == 'd2xJz' :
                    d2xJz = fields.deriv(Jz, dl, 2, 0)
                    term += d2xJz
                elif terms == 'd2yJz' :
                    d2yJz = fields.deriv(Jz, dl, 2, 1)
                    term += d2yJz
                else :
                    raise ValueError('for Z component of Hyper-resistive term, "terms" can be "all", "d2xJz" or "d2yjJz"')

            else :
                raise ValueError('comp can only be "x", "y" or "z"')

        if run.ndim == 3:
            if comp.lower() == 'x' :
                Jx = run.getJ(t, 'x')
                if terms == 'all' :
                    d2xJx = fields.deriv(Jx, dl, 2, 0)
                    d2yJx = fields.deriv(Jx, dl, 2, 1)
                    d2zJx = fields.deriv(Jx, dl, 2, 2)
                    term += d2xJx+d2yJx+d2zJx
                elif terms == 'd2xJx' :
                    d2xJx = fields.deriv(Jx, dl, 2, 0)
                    term += d2xJx
                elif terms == 'd2yJx' :
                    d2yJx = fields.deriv(Jx, dl, 2, 1)
                    term += d2yJx
                elif terms == 'd2zJx':
                    d2zJx = fields.deriv(Jx, dl, 2, 2)
                    term += d2zJx
                else :
                    raise ValueError('for X component of Hyper-resistive term, "terms" can be "all", "d2xJx", "d2yjJx" or "d2zJx"')


            elif comp.lower() == 'y' :
                Jy = run.getJ(t, 'y')
                if terms == 'all' :
                    d2xJy = fields.deriv(Jy, dl, 2, 0)
                    d2yJy = fields.deriv(Jy, dl, 2, 1)
                    d2zJy = fields.deriv(Jy, dl, 2, 2)
                    term += d2xJy+d2yJy+d2zJy
                elif terms == 'd2xJy' :
                    d2xJy = fields.deriv(Jy, dl, 2, 0)
                    term += d2xJy
                elif terms == 'd2yJy' :
                    d2yJy = fields.deriv(Jy, dl, 2, 1)
                    term += d2yJy
                elif terms == 'd2zJy':
                    d2zJy = fields.deriv(Jy, dl, 2, 2)
                    term += d2zJy
                else :
                    raise ValueError('for Y component of Hyper-resistive term, "terms" can be "all", "d2xJy", "d2yjJy" or "d2zJy"')

            elif comp.lower() == 'z' :
                Jz = run.getJ(t, 'z')
                if terms == 'all' :
                    d2xJz = fields.deriv(Jz, dl, 2, 0)
                    d2yJz = fields.deriv(Jz, dl, 2, 1)
                    d2zJz = fields.deriv(Jz, dl, 2, 2)
                    term += d2xJz+d2yJz+d2zJz
                elif terms == 'd2xJz' :
                    d2xJz = fields.deriv(Jz, dl, 2, 0)
                    term += d2xJz
                elif terms == 'd2yJz' :
                    d2yJz = fields.deriv(Jz, dl, 2, 1)
                    term += d2yJz
                elif terms == 'd2zJz':
                    d2zJz = fields.deriv(Jz, dl, 2, 2)
                    term += d2zJz
                else :
                    raise ValueError('for Z component of Hyper-resistive term, "terms" can be "all", "d2xJz", "d2yjJz" or "d2zJz"')

            else :
                raise ValueError('comp can only be "x", "y" or "z"')


    term *= run.getHyperResistivity()/ntot

    if smooth.lower() == 'yes':
        term = gf(term, sigma = sigma, order = 0)

    return term















def electronInertia(run, time, smooth='no', sigma=1):

    """
    returns the electron inertia term
    NOT WORKING AS IS... TODO

    @return  : numpy array with electron inertia term
    Exemple  :
    Creation : 2012-08-15 14:52:15.062321

    """

    if isinstance(time, collections.Iterable) == False:
        time = [time]
        ntot = 1
    else:
        ntot = time.size

    einertia = np.zeros(shape=run.vfield_shape, dtype=run.dtype)

    for t in time:
        Ve  = run.GetVe(time)
        vgv = fields.VgradV(Ve, run.dl)
        q   = run.GetCharge('electrons')
        print(q)
        m_e = run.GetMass('electrons')
        print(m_e)
        print(type(m_e), type(q), type(vgv))

        einertia += m_e/q * vgv


    einertia /= ntot

    if smooth.lower() == 'yes':
        for c in range(3):
            einertia[..., c] = gf(inertia[..., c], sigma=sigma, order=0)

    return einertia


def electron_pressure_gyro(run, time, smooth='no', sigma=1, silent='no'):
    """calculates the contribution of the gyrotropic electron pressure
    to the electron pressure term in Ohm's law


    @param silent 'no' or 'yes'

    @return: 2D numpy array containing the electron pressure term

    Exemple  :

    Creation : 2012-08-15 10:21:33.373750

    """

    if isinstance(time, collections.Iterable) == False:
        time = [time]
        ntot = 1
    else:
        ntot = time.size

    shape = (3,run.ncells[0]+1, run.ncells[1]+1)

    gradP = np.zeros(shape,'float32',order='F')

    for t in time:
        Pe    = run.GetPe(t, silent=silent)
        B     = run.GetB(t, silent=silent)
        P_ng,P_pp,P_para,P_perp = pressure.nongyro_decomposition(Pe,B)

        divPeng = run.divT(P_pp)
        n     = run.GetN(t, species='electrons', silent=silent)

        gradP += divPeng/n

    gradP /= ntot

    gradP /= run.GetCharge('electrons')

    if smooth.lower() == 'yes':
        for c in range(3):
         gradP[c,:,:] = ndimage.gaussian_filter(gradP[c,:,:],\
                            sigma=sigma, order=0)


    return gradP


def electron_pressure_nongyro(run, time, smooth='no', sigma=1):
    """calculates the contribution of the non-gyrotropic electron pressure
    to the electron pressure term in Ohm's law


    @param silent 'no' or 'yes'

    @return: 2D numpy array containing the electron pressure term

    Exemple  :

    Creation : 2012-08-15 10:21:33.373750

    """

    if isinstance(time, collections.Iterable) == False:
        time = [time]
        ntot = 1
    else:
        ntot = time.size

    gradP = np.zeros(shape=run.vfield_shape,dtype=run.dtype)

    for t in time:
        Pe    = run.GetPe(t)
        B     = run.GetB(t)
        P_ng,P_pp,P_para,P_perp = pressure.nongyro_decomposition(Pe,B)


        divPeng = run.divT(P_ng)
        n     = run.GetN(t, species='electrons', silent=silent)

        gradP += divPeng/n

    gradP /= ntot

    gradP /= run.GetCharge('electrons')

    if smooth.lower() == 'yes':
        for c in range(3):
         gradP[c,:,:] = ndimage.gaussian_filter(gradP[c,:,:],\
                            sigma=sigma, order=0)


    return gradP



