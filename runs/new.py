#
#! /usr/bin/env python
#
#
import numpy as np
import matplotlib.pyplot as plt
from fipy import Grid2D, CellVariable, DiffusionTerm, FixedFlux
#
#
#
def GetBx(B0, L, dl) :

    N = [np.int(np.round(i)) for i in np.divide(L, dl)]
    Wx = B0*np.linspace(-L[1]/2, L[1]/2, num = N[1]+1)
    Bx = np.tile(Wx, (N[0]+1, 1))

    return Bx
#
#
def GetBy(B0, L, dl) :

    N = [np.int(np.round(i)) for i in np.divide(L, dl)]
    Wy = B0*np.linspace(-L[0]/2, L[0]/2, num = N[0]+1)
    By = np.repeat(Wy, N[1]+1)
    By = np.reshape(By, [N[0]+1, N[1]+1])

    return By
#
#
def GetJz(Bx, By, dl) :

    dBxdy = np.gradient(Bx, dl[0], dl[1])[1]
    dBydx = np.gradient(By, dl[0], dl[1])[0]

    Jz = dBydx-dBxdy

    return Jz
#
#
def GetAz(Bx, By, Jz, Nc, dl) :

    mesh = Grid2D(dx = dl[0], dy = dl[1], nx = Nc[0], ny = Nc[1])

    _Az = CellVariable(mesh=mesh, value=0.0)

    _Bx = CellVariable(mesh=mesh)
    _Bx.value = np.reshape(Bx, Nc[0]*Nc[1], order = 'F')

    _By = CellVariable(mesh=mesh)
    _By.value = np.reshape(By, Nc[0]*Nc[1], order = 'F')

    _Jz = CellVariable(mesh=mesh)
    _Jz.value = np.reshape(Jz, Nc[0]*Nc[1], order = 'F')

    _Az.equation = (DiffusionTerm(coeff = 1.0)+_Jz == 0)

    #diffcoeff = CellVariable(mesh=mesh, value = 1.0)
    #diffTerm = DiffusionTerm(coeff = diffcoeff)
    #diffTerm = DiffusionTerm(coeff = 1.0)
    #diffcoeff.constrain(0., mesh.exteriorFaces)

    # beware of the sign of the flux : always consider outward direction
    _Az.faceGrad.constrain( _By.arithmeticFaceValue, where = mesh.facesLeft)
    _Az.faceGrad.constrain(-_By.arithmeticFaceValue, where = mesh.facesRight)
    _Az.faceGrad.constrain(-_Bx.arithmeticFaceValue, where = mesh.facesBottom)
    _Az.faceGrad.constrain( _Bx.arithmeticFaceValue, where = mesh.facesTop)

    #_Az.faceGrad.constrain( _By.arithmeticFaceValue[mesh.facesLeft] , where=mesh.facesLeft)
    #_Az.faceGrad.constrain(-_By.arithmeticFaceValue[mesh.facesRight] , where=mesh.facesRight)
    #_Az.faceGrad.constrain(-_Bx.arithmeticFaceValue[mesh.facesBottom] , where=mesh.facesBottom)
    #_Az.faceGrad.constrain( _Bx.arithmeticFaceValue[mesh.facesTop] , where=mesh.facesTop)

    #_Az.equation = (diffTerm+_Jz == 0)
    #_Az.equation = (DiffusionTerm(diffcoeff)+ (mesh.exteriorFaces*exteriorFlux).divergence + _Jz== 0)
    #_Az.equation = (diffTerm + _Jz == 0)

    #BCs = [FixedFlux(value= _By.getFaceValue(), faces=mesh.getFacesLeft()),
    #       FixedFlux(value=-_By.getFaceValue(), faces=mesh.getFacesRight()),
    #       FixedFlux(value=-_Bx.getFaceValue(), faces=mesh.getFacesBottom()),
    #       FixedFlux(value= _Bx.getFaceValue(), faces=mesh.getFacesTop())]

    #_Az.equation.solve(var=_Az, boundaryConditions=BCs)
    _Az.equation.solve(var=_Az)

    Az = np.reshape(_Az.value, Nc, order = 'F')

    return Az
#
#
#
Lb = [20, 10]
dl = [.1, .1]

Nc = [np.int(np.round(i))+1 for i in np.divide(Lb, dl)]

Bx = GetBx(1.0, Lb, dl)
By = GetBy(0.2, Lb, dl)
Jz = GetJz(Bx, By, dl)
Az = GetAz(Bx, By, Jz, Nc, dl)

image = plt.imshow(np.transpose(Az), origin = 'lower')
plt.colorbar(image)
plt.show()
plt.savefig('new.pdf')
