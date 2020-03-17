


import numpy as np
from scipy.ndimage import gaussian_filter1d as gf1d



# ______________________________________________________________________________
#
def deriv(scalarField, dl, order = 1, axis = 1):

    """
    return the derivative of the scalarField

    Longer description here

    @param scalarField (np.ndarray) contains the field (scalar) to derive
    @param dl is a scalar (float) containing the mesh size
    @param order (int)  order of the derivative
    @param axis (int) is the index of the axis along which to derive

    @return  : @todo
    Exemple  :
    Creation : 2015-03-02

    """



    if len(scalarField.shape) == 1:
        out = deriv1D(scalarField, order, dl)

    elif len(scalarField.shape) == 2:
        if axis > 1:
            raise ValueError("wrong axis {0:1d} for array of dim {1:1d}".format(axis, len(scalarField.shape)))
        out = deriv2D(scalarField, order, axis, dl)

    elif len(scalarField.shape) == 3:
        if axis > 2:
            raise ValueError("wrong axis {0:1d} for array of dim {1:1d}".format(axis, len(scalarField.shape)))
        out = deriv3D(scalarField, order, axis, dl)

    else:
        raise ValueError("bad shape of the scalarField")

    return out
# ______________________________________________________________________________





# ______________________________________________________________________________
#
def deriv1D(field, order, dl):

    """
    @todo: Brief Docstring for deriv1D
    Longer description here

    @param field @todo
    @param order @todo
    @param dl @todo

    @return  : @todo
    Exemple  :
    Creation : 2015-03-02

    """

    return gf1d(field, order=order, sigma=1, axis=0)/pow(dl, order)
# ______________________________________________________________________________





# ______________________________________________________________________________
#
def deriv2D(field, order, axis, dl):

    """
    @todo: Brief Docstring for deriv1d
    Longer description here
    dl is a float; namely the grid size in the axis direction

    @param field @todo
    @param order @todo
    @param dl @todo

    @return: @todo

    Exemple  :

    Creation : 2015-03-02

    """
    # sigma=1 ??????????????????????????????????????????????????????????????????

    return gf1d(field, order=order, sigma=1, axis=axis)/pow(dl[axis], order)
# ______________________________________________________________________________





# ______________________________________________________________________________
#
def deriv3D(field, order, axis, dl):

    """
    @todo: Brief Docstring for deriv1d
    Longer description here

    @param field @todo
    @param order @todo
    @param dl @todo

    @return: @todo

    Exemple  :

    Creation : 2015-03-02

    """

    return gf(field, order=order,sigma=1, axis=axis)/pow(dl, order)
# ______________________________________________________________________________





# ______________________________________________________________________________
#
def divV(vectorField, dl):

    """
    @todo: returns the divergence of a vector (1st order tensor)
    Longer description here @todo

    @param vectorField is a np.ndarray - last index is for the vector components

    @return  : a scalar, of the same type as vectorField
    Exemple  :
    Creation : 2015-03-02

    """


    ndims = len(vectorField.shape)-1

    if ndims == 1:
        out = deriv(vectorField[..., 0], dl[0], 1, 0)

    elif ndims == 2:
        out = deriv(vectorField[..., 0], dl[0], 1, 0) \
            + deriv(vectorField[..., 1], dl[1], 1, 1)

    elif ndims == 3:
        out = deriv(vectorField[..., 0], dl[0], 1, 0) \
            + deriv(vectorField[..., 1], dl[1], 1, 1) \
            + deriv(vectorField[..., 2], dl[2], 1, 2)


    return out
# ______________________________________________________________________________





# ______________________________________________________________________________
#
def curl(vectorField, dl):

    """
    @todo: Brief Docstring for curl
    Longer description here

    @param vectorField @todo
    @param *dl @todo

    @return  : @todo
    Exemple  :
    Creation : 2015-03-02

    """


    ndims = len(vectorField.shape)-1
    shape = vectorField.shape
    dtype = vectorField.dtype

    out = np.zeros(shape=shape, dtype=dtype)

    if ndims == 1:
        out[..., 0] = 0
        out[..., 1] = -deriv(vectorField[..., 2], dl, 1, 0)
        out[..., 2] =  deriv(vectorField[..., 1], dl, 1, 0)

    elif ndims == 2:
        out[..., 0] =  deriv(vectorField[..., 2], dl[1], 1, 1)

        out[..., 1] = -deriv(vectorField[..., 2], dl[0], 1, 0)

        out[..., 2] =  deriv(vectorField[..., 1], dl[0], 1, 0) \
                      -deriv(vectorField[..., 0], dl[1], 1, 1)

    elif ndims == 3:
        out[..., 0] =  deriv(vectorField[..., 2], dl[1], 1, 1) \
                      -deriv(vectorField[..., 1], dl[2], 1, 2)

        out[..., 1] = -deriv(vectorField[..., 2], dl[0], 1, 0) \
                      +deriv(vectorField[..., 0], dl[2], 1, 2)

        out[..., 2] =  deriv(vectorField[..., 1], dl[0], 1, 0) \
                      -deriv(vectorField[..., 0], dl[1], 1, 1)

    return out
# ______________________________________________________________________________





# ______________________________________________________________________________
#
def scalarLaplacian(scalarField, dl):

    """
    @todo: Brief Docstring for scalarLaplacian
    Longer description here

    @param field @todo
    @param *dl) @todo

    @return  : @todo
    Exemple  :
    Creation : 2015-03-02

    """


    ndims = len(scalarField.shape)

    if ndims == 1:
        out = deriv(scalarField, dl, 2, 0)

    elif ndims == 2:
        out = deriv(scalarField, dl[0], 2, 0) \
            + deriv(scalarField, dl[1], 2, 1) \

    elif ndims == 3:
        out = deriv(scalarField, dl[0], 2, 0) \
            + deriv(scalarField, dl[1], 2, 1) \
            + deriv(scalarField, dl[2], 2, 2)

    else:
        raise ValueError("dimension (%d) is to big" % ndims)

    return out
# ______________________________________________________________________________





# ______________________________________________________________________________
#
def vectorLaplacian(vectorField, dl):

    """
    @todo: Brief Docstring for vectorLaplacian
    Longer description here

    @param vectorField @todo
    @param *dl @todo

    @return  : @todo
    Exemple  :
    Creation : 2015-03-02

    """


    ndims = len(vectorField.shape)-1
    shape = vectorField.shape
    dtype = vectorField.dtype

    out   = np.zeros(shape=shape, dtype=dtype)

    if ndims == 1:
        out[...,0] = scalarLaplacian(vectorField[..., 0], dl)
        out[...,1] = scalarLaplacian(vectorField[..., 1], dl)
        out[...,2] = scalarLaplacian(vectorField[..., 2], dl)


    elif ndims == 2:
        out[...,0] = scalarLaplacian(vectorField[..., 0], dl)
        out[...,1] = scalarLaplacian(vectorField[..., 1], dl)
        out[...,2] = scalarLaplacian(vectorField[..., 2], dl)


    elif ndims == 3:
        out[...,0] = scalarLaplacian(vectorField[..., 0], dl)
        out[...,1] = scalarLaplacian(vectorField[..., 1], dl)
        out[...,2] = scalarLaplacian(vectorField[..., 2], dl)

    return out

# ______________________________________________________________________________





# ______________________________________________________________________________
#
def dot(V1, V2):

    """
    @todo: Brief Docstring for dot
    Longer description here

    @param V1 @todo
    @param V2 @todo

    @return: @todo

    Exemple  :

    Creation: 2015-03-03

    """


    return V1[..., 0]*V2[..., 0]+V1[..., 1]*V2[..., 1]+V1[..., 2]*V2[..., 2]
# ______________________________________________________________________________





# ______________________________________________________________________________
#
def VgradV(V, dl):

    """
    @todo: Brief Docstring for VgradV
    Longer description here

    @param V @todo
    @param dl @todo

    @return: @todo

    Exemple  :

    Creation : 2015-03-03

    """


    ndims = len(V.shape)-1
    shape = V.shape
    dtype = V.dtype

    out  = np.zeros(shape=shape, dtype=dtype)

    if ndims == 1:
        for c in np.arange(3):
            out[..., c] = V[..., 0]*deriv(V[...,c], dl[0], order=1, axis=0)

    elif ndims == 2:
        for c in np.arange(3):
            out[..., c] = V[..., 0]*deriv(V[...,c], dl[0], order=1, axis=0) \
                        + V[..., 1]*deriv(V[...,c], dl[1], order=1, axis=1)

    elif ndims == 3:
        for c in np.arange(3):
            out[..., c] = V[..., 0]*deriv(V[...,c], dl[0], order=1, axis=0) \
                        + V[..., 1]*deriv(V[...,c], dl[1], order=1, axis=1) \
                        + V[..., 2]*deriv(V[...,c], dl[2], order=1, axis=2)

    return out
# ______________________________________________________________________________





# ______________________________________________________________________________
#
def divT(T, dl):

    """
    @todo: returns the divergence of a 2nd order tensor
    Longer description here

    @param T @todo
    @param dl @todo

    @return  : @todo
    Exemple  :
    Creation : 2015-03-03

    """

    ndims = len(T.shape)-2
    divT  = np.zeros(shape=T.shape[:-1], dtype=T.dtype)

    if ndims == 1:
        divT[..., 0] = deriv(T[..., 0, 0], dl[0], order=1, axis=0)
        divT[..., 1] = deriv(T[..., 0, 1], dl[0], order=1, axis=0)
        divT[..., 2] = deriv(T[..., 0, 2], dl[0], order=1, axis=0)

    elif ndims == 2:

        divT[..., 0] = deriv(T[..., 0, 0], dl[0], order=1, axis=0) \
                     + deriv(T[..., 0, 1], dl[1], order=1, axis=1)

        divT[..., 1] = deriv(T[..., 0, 1], dl[0], order=1, axis=0) \
                     + deriv(T[..., 1, 1], dl[1], order=1, axis=1)

        divT[..., 2] = deriv(T[..., 0, 2], dl[0], order=1, axis=0) \
                     + deriv(T[..., 1, 2], dl[1], order=1, axis=1)

    elif ndims == 3:

        divT[..., 0] = deriv(T[..., 0, 0], dl[0], order=1, axis=0) \
                     + deriv(T[..., 0, 1], dl[1], order=1, axis=1) \
                     + deriv(T[..., 0, 2], dl[2], order=1, axis=2)

        divT[..., 1] = deriv(T[..., 0, 1], dl[0], order=1, axis=0) \
                     + deriv(T[..., 1, 1], dl[1], order=1, axis=1) \
                     + deriv(T[..., 1, 2], dl[2], order=1, axis=2)

        divT[..., 2] = deriv(T[..., 0, 2], dl[0], order=1, axis=0) \
                     + deriv(T[..., 1, 2], dl[1], order=1, axis=1) \
                     + deriv(T[..., 2, 2], dl[2], order=1, axis=2)
    return divT
# ______________________________________________________________________________





# ______________________________________________________________________________
#
def trace(T):

    """
    calculates the Trace of a tensor field

    @param T tensor
    @return: trace of the tensor field

    Exemple  : tr = fields.trace(Pe)

    Creation : 2015-03-05
    """


    return T[..., 0] + T[..., 1] + T[..., 2]
# ______________________________________________________________________________





# ______________________________________________________________________________
#
def norm(V):

    """
    calculates the norm of a vector field

    @param V the vector field

    @return: sqrt(Vx**2 + Vy**2 + Vz**2)

    Exemple  :

    Creation : 2015-03-05

    """


    return np.sqrt(V[..., 0]**2 + V[..., 1]**2 + V[..., 2]**2)
# ______________________________________________________________________________

