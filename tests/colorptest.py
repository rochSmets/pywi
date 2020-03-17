

import colorp as cp
import numpy as np
import hecklerun as hr
import momentum as mmt


def test():
    # Create some sample data
#    dx = np.linspace(0,1,20)
#    X,Y = np.meshgrid(dx,dx)
#    Z  = X**2 - Y
#    Z3 = Y
#    Z2 = X


    r = hr.HeckleRun('r','/Volumes/drobo/nico/asymmetric/062/')

    Vi = r.GetVi(50.)
    Vix= Vi[0,:,:]
    Viy = Vi[1,:,:]

    E = r.GetE(50.0)
    Ex = E[0,:,:]
    Ey = E[1,:,:]

    dtf = r.GetTimeStep('fields')


    x   = r.dl[0]*np.arange(0, r.ncells[0]+1)
    y   = r.dl[1]*np.arange(0, r.ncells[1]+1)

    Flux = r.GetFlux(50.)

    c = cp.ColorP(Vix,x, y,
                  xlabel='x',
                  ylabel='y',
                  title=r'$V_{ix}$ and $\mathbf{E}$',
                  vmin = -1,vmax=1,
                  extent=(60,85,15,35),
                  vecfield=(Ex,Ey),
                  veccolorcode=np.sqrt(Ex**2 + Ey**2),
                  fieldlines=Flux, 
                  save=True
                  )

    c.display()

#    c.show()


def main():
    print "super test"
    test()

if __name__ == '__main__':
    main()
