# this is a test of the streamline module









import numpy as np
import hecklerun as hr


import streamlines as stl



def getlines():

    run062 = hr.HeckleRun('run062',\
                          '/Volumes/drobo/nico/asymmetric/062/')



    Vi = run062.GetVi(50.)

    dx = run062.dl[0]
    dy = run062.dl[1]

    seeds = [ (76./dx,30./dy),(77./dx,30./dy)]

    Vix = Vi[0,:,:]
    Viy = Vi[1,:,:]

    print "Calculating fieldlines..."

    st = stl.fieldlines(Vix,Viy,seeds)

    print "fieldlines calculated."

    return st

    #streamlines = stl.Streamlines(x,y, Vi[0,:,:], Vi[1,:,:], seeds)

    # Puts a scalar value on the fieldlines
    # streamlines.attach_scalar(scalarfield)


    # attach a vector to the streamline
    # streamlines.attach_vector(vectorfield)


    # streamlines.add2plot(ax)


def main():
    print "super main()"


if __name__ == "__main__":
    main()
