

import colorp as cp
import numpy as np
import hecklerun as hr
import momentum as mmt
import matplotlib.pyplot as plt


def plot_momentum(run, mmt, position):

    mag      = mmt.magnetic()
    elec     = mmt.electric()
    press    = mmt.pressure()
    acc      = mmt.acceleration()
    tot      = mmt.totalforce()

    ix, iy = run.coord2indices((position,0.))

    tot2 = mag + elec + press

    plt.figure()
    plt.plot(mag[1,ix,:], label=r'$nq\mathbf{V}\times\mathbf{B}$')
    plt.plot(elec[1,ix,:], label=r'$nq\mathbf{E}$', linestyle='dashed')
    plt.plot(press[1,ix,:], label=r'$\nabla\cdot \mathbf{P}$')
    plt.plot(acc[1,ix,:], label=r'$mn\mathbf{v}\cdot\nabla\mathbf{v}$')
    plt.plot(tot[1,ix,:], label=r'$\Sigma$ Forces')
    plt.plot(tot2[1,ix,:], label=r'$\Sigma$ Forces 2')

    plt.axis([265,390,-1,0.4])

    plt.legend()
    plt.show()




def test():
    # Create some sample data
#    dx = np.linspace(0,1,20)
#    X,Y = np.meshgrid(dx,dx)
#    Z  = X**2 - Y
#    Z3 = Y
#    Z2 = X


#    r = hr.HeckleRun('r','/Volumes/drobo/nico/asymmetric/078/')
    r = hr.HeckleRun('r','/Users/naunai/Documents/data/071/')

    Vi = r.GetVi(26.)
    Vix= Vi[0,:,:]
    Viy = Vi[1,:,:]


    immt = mmt.Momentum(r,np.arange(25,27,0.2), 'ions', smoothing='yes')

#    force = r.divP(42.5,'ions')#immt.pressure(sigma=1) #r.GetField("piyy",42.5)#immt.pressure(sigma=1)
    force  = immt.electric(sigma=0)
    force  = immt.acceleration(sigma=0)
    force  = immt.pressure(sigma=1)

    dtf = r.GetTimeStep('fields')


    x   = r.dl[0]*np.arange(0, r.ncells[0]+1)
    y   = r.dl[1]*np.arange(0, r.ncells[1]+1)


    Flux = r.GetFlux(26.)

    c = cp.ColorP(force[1,:,:],x, y,
                  xlabel='x',
                  ylabel='y',
                  title='Test Momentum on Heckle',
                  vmin = -0.9,vmax=0.9,
                  extent=(55,95,15,35) ,#extent=(20,40,5,20),
                  vecfield=(force[0,:,:],force[1,:,:]),
                  veccolorcode=np.sqrt(Vix**2 + Viy**2),
                  vecscale =1.0,
                  vecbin=6,
                  fieldlines=Flux#, 
                  #save=True
                  )

    c.display()

    plot_momentum(r, immt, 10.0)

#    c.show()


def main():
    print "super test"
    test()

if __name__ == '__main__':
    main()
