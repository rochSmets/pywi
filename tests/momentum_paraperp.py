


# this file tests the para and perp component of the acceleration and 
# pressure force.



import picgsfc as pr
import numpy as np
import momentum as mmt
import colorp as cp







def main():

    run = pr.PICGSFCrun('r','/Volumes/drobo/nico/asymmetric/by10')

    immt = mmt.Momentum(run, np.arange(50,53.5,0.5),'ions')








if __name__ == '__main__':
    main()



