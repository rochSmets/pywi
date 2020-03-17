# this file is intended to test the energy module




import hecklerun as hr
import energy as nrj








def main():

    run = hr.HeckleRun('r','/Volumes/drobo/nico/asymmetric/078')

    # this will read the data
    # calculate :vsp    
    nrjbox = nrj.box(run,
                     time,
                     extent = (x0,y0,x1,y1),
                     specie='ions',
                     )

    
    print "super main"



if __name__ == '__main__':
    main()

