#!/opt/local/bin/python
# encoding: utf-8


import pypic.RunModules.heckle.hecklerun as hr
import fields.fields as fields
import pypic.colorplot.colorp as colp
import pypic.ohm as ohm
import numpy as np


def main():

    r  = hr.Heckle('/Users/nicolasaunai/Documents/code/hybride/heckle/bin/')
    #r = hr.Heckle('/Users/nicolasaunai/Desktop/')
    B  = r.GetB(0.)
    J  = r.GetJ(0.)
    Np = r.GetNp(0.)
    x = r.GetCoords(axis=0)
    y = r.GetCoords(axis=1)

    db    = fields.div(B,r.dl)
    curlB = fields.curl(B, r.dl)
    LapJ  = fields.vectlaplacian(J,r.dl)

    Ve   = r.GetVe(0.)
    Pe   = r.GetPe(0.)
    Pp   = r.GetPp(0.)
    Vp   = r.GetVp(0.)
    vgv  = fields.VgradV(Ve, r.dl)
    divt = fields.divT(Pe, r.dl)

    print db.max()
    colp.tofile(db, x, y, filename='db.png')
    colp.tofile(curlB[...,2], x, y, filename='curlb.png')
    colp.tofile(J[...,2], x, y, filename='jz.png')
    colp.tofile(1e-4*LapJ[...,2], x, y, filename='hprsty.png')
    colp.tofile(Np[...], x, y, filename='Np.png')
    colp.tofile(Pp[...,0,0], x, y, filename='Pp.png')
    colp.tofile(Pe[...,0,0], x, y, filename='Pe.png')
    colp.tofile(Vp[...,0], x, y, filename='Vp.png')



    for it,t in enumerate( np.arange(11)*0.01):
        print t,it
        B  = r.GetB(t)
        Np = r.GetNp(t)
        db  = fields.div(B, r.dl)
        colp.tofile(B[...,1],x,y, filename='By_%04d.png' % it)
        colp.tofile(B[...,0],x,y, filename='Bx_%04d.png' % it)
        colp.tofile(B[...,2],x,y, filename='Bz_%04d.png' % it)
        colp.tofile(Np,x,y, filename='Np_%04d.png' % it)
        colp.tofile( db,x,y, filename='divB_%04d.png' % it)





if __name__ == '__main__':
    main()
