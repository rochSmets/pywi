

#------------------------------
# 
    
# -----------------------------
def main():
    import hecklerun as hr
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    import recrate as rr

    run078 = hr.HeckleRun('078', '/Volumes/drobo/nico/asymmetric/078/')
    run076 = hr.HeckleRun('076', '/Volumes/drobo/nico/asymmetric/076/')

    rate078 = rr.RecRate(run078,0,42.5, silent='no')    
    rate076 = rr.RecRate(run076, 0, 40.0, silent='yes')



    plt.plot(rate078.t,rate078.rate)
    plt.plot(rate078.t[1:]-0.5*run078.GetTimeStep('fields'),rate078.ratef)
    
#   plt.axis((0,39,-10,-8))
    plt.axis((0,39,-0.07,0.06))
    plt.show()
    
# ----------------------------






if __name__ == "__main__":
    main()
    
