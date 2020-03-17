#!/opt/local/bin/python
# encoding: utf-8

import pypic.RunModules.heckle.hecklerun as hr
import pypic.colorplot.colorp as colp
import pypic.drifts as mdrift


def main():
    r  = hr.Heckle('/Users/nicolasaunai/Desktop/')
    P  = r.GetP(0.0, species='protons')
    Pe = r.GetPe(0.0)
    Ne = r.GetNe(0.0)
    F  = r.GetFlux(0.)
    Vi = r.GetVp(0.)
    B  = r.GetB(0.)
    J  = r.GetJ(0.)

    Vperp = mdrift.Vperp(Vi,B)

    x = r.GetCoords(axis=0)
    y = r.GetCoords(axis=1)

    colp.tofile(Ne[:,:,0],x,y,filename="test_ne.png",fieldlines=F)
    colp.tofile(J[2,:,:,0],x,y,filename="test_Jz.png",fieldlines=F)




if __name__ == '__main__':
    main()
