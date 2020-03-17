

#------------------------------
#

# -----------------------------
def main():
    import picgsfc as pr
    import numpy as np
    import momentum as mmt
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import hecklerun as hr
    import scipy.ndimage as ndimage


#    run = pr.PICGSFCrun('by10',\
#                            '/Volumes/drobo/nico/asymmetric/by10/data/')

    run = hr.HeckleRun('r', '/Volumes/drobo/nico/asymmetric/79')

    immt = mmt.Momentum(run, np.arange(40,42,0.5), 'ions', 
                        smoothing='yes', silent='yes')


    datac = run.GetB(21., silent='yes')
    datac = datac[1,:,:]

    plt.imshow(np.transpose(datac),     # data to show \
       cmap='jet',                     # color map\
       origin='lower',                 # where is the zero\
       interpolation='nearest',        # pixels \
       extent=[0,64,0,25.6],           # range for x and y 
       aspect='equal',                 # auto = aspect ratio given by 'extent'
       vmin=np.min(datac),vmax=np.max(datac))
#
#  ax1.set_title(r'$J_z$', fontsize=25)
#  ax1.set_xlabel(r'$x/\delta_i$', fontsize=20)
#  ax1.set_ylabel(r'$y/\delta_i$', fontsize=20)

    plt.colorbar()
    plt.savefig('acc')
#   plt.show()


    plt.close()

    # take the indice of x=35.0, don't care about iy
    ix, iy = run.coord2indices((35.0,0.0))

    print "We do a cut at ix(x=35.) = ", ix

    forces = immt.totalforce(sigma=1)
    acc    = immt.acceleration(sigma=1)


    mag   = (immt.magnetic(sigma=0))[2, ix, :]
    elec  = (immt.electric(sigma=0))[2, ix, :]
    press = (immt.pressure(sigma=[1.,0.5]))[2, ix, :]
    acc   = (immt.acceleration(sigma=[1.,0.5]))[2, ix, :]


    z = run.GetCoord(0,50, axis='y')

    print np.shape(z), np.shape(mag)


    plt.plot(z,mag, label='Magnetic force')
    plt.plot(z,elec, label='Electric force')
    plt.plot(z,press, label='Pressure force')
    plt.plot(z,acc, label='Acceleration (steady part)')

    plt.axis((8.8,16.8,-2,2))
    plt.legend()

    plt.show()
    #plt.savefig('force.png')


# ----------------------------






if __name__ == "__main__":
    main()
