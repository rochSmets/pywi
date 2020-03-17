

import numpy as np
import os
from pickle import dump, load
import copy
import glob
import math
from pywi.tools.misc.tools import smooth
from pywi.tools.reconnection import xpoint as xtrack






#==================================================================
#==================================================================
def reload(path):

    """

    loads a 'RecRate' object from a specified path

    Args :
        - path

    Example :
        - rate = load('/Users/me/Documents/here/')

    """


    filename = os.path.join(path,'rate.dat')

    with open(filename, 'rb') as f:
        rate = load(f)
        return rate
#==================================================================







#==================================================================
#==================================================================

def fromrun(run, times, path=None, av = 0.1, silent='No'):

    """
    returns a reconnection rate object for the specified run during
    the specified time interval.


    Genericity :
        This routine should work with all codes. (see prerequisites)

    Args :
        - run    : run object
        - times  : times for which to calculate the rate (ion units)

    Keywords :

        - path : a user custom default path for the RecRate object. If None is specified, the
             default path will be the return of run.WhereAmI()

        - silent : (default 'No'). If yes, prevent the routine from being in verbose mode.


    Prerequisite :
        To be generic, the run object ('run') should provide the routine with
        the following methods and attributes :

        - run.GetRunTime()
        - run.wherAmI()
        - run.GetField()
        - run.indices2coord()
        - run.dl
    """

#       from scipy.ndimage.filters import gaussian_filter

    # safety check
    # t1 and t2 must be in a valid interval
    if (silent.lower() == 'no'):
        print("RecRate> Checking interval validity...")



    ntot = times.size


    # the 'rate' array and xpoint array are allocated
    # for the whole simulation. This does not require
    # much space and is useful for reloading/completing
    # unfinished rate calculations

    rate   = np.zeros(ntot)
    xpoint = np.zeros((2,ntot))
    flux   = np.zeros(ntot)


    if path is None:
        path = run.whereAmI()


#        for i in range(it1, it2+1, 1):
    for i,t in enumerate(times):

        # reconnection electric field and magnetic flux
        erec   = run.GetE(t)
        f      = run.GetFlux(t)

        # calculate the X point position and put it in real coords.
        xpointij    = xtrack.GetXPoint(f)
        xpointij_array = np.zeros((2,1))
        xpointij_array[0,0] = xpointij[0]
        xpointij_array[1,0] = xpointij[1]
        #xpoint[:,i] = run.indices2coord(xpointij)
        xpoint_tmp = run.indices2coord(xpointij_array)
        xpoint[:,i] = xpoint_tmp[:,0]


        # the average area
        avarea = (int(xpointij[0]-av/run.dl[0]),\
                  int(xpointij[0]+av/run.dl[0])+2,\
                  int(xpointij[1]-av/run.dl[1]),\
                  int(xpointij[1]+av/run.dl[1])+2)


        # the reconection rate as the out of plane electric field
        #averaged around the X point.
        rate[i] = -np.mean(erec[avarea[0]:avarea[1],\
                                avarea[2]:avarea[3],2])

        # same average for the magnetic flux
        flux[i] = np.mean(f[avarea[0]:avarea[1],\
                               avarea[2]:avarea[3]])


        if (silent.lower() == 'no'):
             print(" RecRate> t =  ",t,\
                   "  rate = ",    rate[i], \
                   "  xpoint = ",  xpoint[:,i])

    # the rate as d(flux)/dt
    dflux = np.diff(flux)
    ratef = smooth(dflux/(times[1:]-times[:-1]))


    return RecRate(times,rate,ratef,flux,av,xpoint,path)

#============================================








class RecRate:
    """

    This class represents o "reconnection rate object".

    Genericity :
        This class should work with all codes.

    """

    #==================================================================
    #==================================================================
    def GetTimeInterval(self):
        """

        returns a tuple containing the time interval over which the
        reconnection rate is calculated

        """
        return (self.t[0], self.t[-1])
    #==================================================================






    #==================================================================
    #==================================================================
    def dump(self,path=None):

        """

        saves the 'RecRate' object at a specified 'path' using the
        Pickle mdule. If no path is specified the the default path
        is used.

        Args :
            -

        Keywords :
            - path : string indicating where the RecRate object is saved

        Example :
            rate.dump() # saves the rate in the default dir.
            rate.dump('/Users/me/Documents/here/')
        """

        if path == None:
            path = self.path

        filename = os.path.join(path,'rate.dat')

        with open(filename, 'wb') as f:
            dump(self,f)
    #==================================================================




    #==========================================================
    #==========================================================
    def __init__(self,time,rate,ratef,flux,av,xpoint,path):
        """builds a RecRate object

        Exemple  : 

        Creation : 2012-11-30 16:03:33.357645

        """
        self.t1     = time[0]
        self.t2     = time[-1]
        self.t      = time
        self.path   = path
        self.xpoint = xpoint
        self.rate   = rate
        self.ratef  = ratef
#       if np.mean(flux-flux[0]) < 0:
#           self.flux = -flux
#       else:
        self.flux = flux - flux[0]
        self.av     = av
    #==========================================================








    #==========================================================
    #==========================================================
    def check(self, save=False):
        """Checks if the rec. electric field is equal to the deriv. of the 
        reconnected flux as it should be. Plot is on screen or save to a file

        @return: Nothing

        Exemple  : rate.check(save=True)

        Creation : 2012-08-14 10:55:28.557358

        """
        import matplotlib.pyplot as plt
        import os


        hdt = 0.5*(self.t[1]-self.t[0])

        plt.plot(self.t, self.rate, label='Reconnection electric field')
        plt.plot(self.t[1:]-hdt, -self.ratef, \
                label=r'$\partial \Phi/\partial t$')

        tmax=np.max(self.t)
        tmin=0.
        #plt.axis((tmin,tmax,-0.05,0.08))

        plt.xlabel(r'$t/\Omega_{ci}^{-1}$')
        plt.ylabel(r'$E_{r}, \partial\Phi/\partial t$')
        plt.title('Reconnection electric field and derivative of reconnected flux')
        plt.legend()

        if save is True:
            plt.savefig(os.path.join(self.path,'checkrate.png'))
        else:
            plt.show()

        plt.close()
    #==========================================================






    #==========================================================
    #==========================================================
    def shift(self, delay, inplace=False):
        """shifts the time array of the RecRate object


        Sometime the difference between two reconnection rate can be 
        understood as a simple time shift. For instance, using two different
        magnetic perturbations may result in the same rate but the one with
        the stronger perturbation is just faster sooner.

        This function allows a better comparison between rates by shifting
        the time array of a RecRate object by a given 'delay'.

        The time shift can be done inplace (inplace=True), meaning that the
        object will be modified, or can result in a new object having the 
        same properties as the original but a time shifted.


        @param delay is the time shift
        @param inplace : False, a new object is returned, True, time is modified
        in the original object.

        @return: Depending on 'inplace', the original object or a new one with
                 a time shifter by 'delay'

        Exemple  : new_rate = rate.shift(-10, inplace=False)
                   rate.shift(-12.4,inplace=True) #modifies the object

        Creation : 2012-08-17 10:05:56.037314

        """
        if inplace is False:
            tmp_rate = copy.deepcopy(self)

        else :
            tmp = self

        tmp_rate.t  += delay
        tmp_rate.t1 += delay
        tmp_rate.t2 += delay


        # it is not enough just to shift the time, since we want 
        # the flux at 0 to be 0 we also have to substract the value
        # flux[delay] to the reconnected flux

        id = np.where(np.abs(tmp_rate.t-np.abs(delay)) == \
                      np.min(np.abs(tmp_rate.t-np.abs(delay))))

       # print "The ID is", id
       # print tmp_rate.flux
       # tmp_rate.flux -= np.mean(tmp_rate.flux[id])
       # print tmp_rate.flux

        return tmp_rate
    #==========================================================







    def __add__(self, rateobj):

        if rateobj.path != self.path:
            print ('Paths are differents')
            return None

        if math.fabs(rateobj.t[-1] - self.t[0]) < 1e-5:
            first  = rateobj
            second = self
        elif math.fabs(self.t[-1]- rateobj.t[0]) < 1e-5:
            first  = self
            second = rateobj
        else:
            print ('not contiguous rates')
            return None

        path   = first.path
        av     = first.av
        t      = np.concatenate((first.t      ,second.t[1:]))
        rate   = np.concatenate((first.rate   ,second.rate[1:]))
        #ratef  = np.concatenate((first.ratef  ,second.ratef[1:]))
        flux   = np.concatenate((first.flux   ,second.flux[1:]))
        x      = np.concatenate((first.xpoint[0,:],second.xpoint[0,1:]))
        y      = np.concatenate((first.xpoint[1,:],second.xpoint[1,1:]))
        xpoint = np.zeros((2,y.size),"float32")
        print("tailles : ", first.xpoint[1,:].size,second.xpoint[1,:].size,y.size)
        xpoint[0,:] = x
        xpoint[1,:] = y

        # the rate as d(flux)/dt
        dflux = np.diff(flux)
        ratef = smooth(dflux/(t[1:]-t[:-1]))


        return RecRate(t,rate,ratef,flux,av,xpoint,path)








# -----------------------------------------------------------------------------
#
#
#     END OF THE METHODS, BEGINING OF MODULE FUNCTIONS
#
#
# -----------------------------------------------------------------------------








#==========================================================
#==========================================================
def exists(path):
    """returns True if a rate exists at the specified path

    @param path a string indicating where to test whether a rate object exists

    @return: True if a RecRate exists

    Exemple  : if exists('/test/here'): print("a rate exists here.")

    Creation : 2012-08-13 17:24:03.039137

    """
    import os

    if os.path.exists(path) == False:
        print("Invalid path")
        return False


    files=os.listdir(path)
    if ('rate.dat' in files) :
        return True
    else:
        return False

#==========================================================










# this function is not a member of the classe
#============================================
#============================================
def plot(*rates, **kwargs):

    """

    plot several reconnection rates, flux etc.

    *rate : RecRate objects
    **kwargs :

            - type :
                - (default) 'rate' plots the rate as the electric field VS time
                - 'flux' plots the reconnected flux as a function of time
                - 'ratef' plots the rec. flux derivative as a function of time
                - 'ratevsflux' plots the flux derivative vs reconnected flux

            - range :
                - [t0,t1, vmin, vmax]
                - (default) [0,tmax, -0.1, 0.1]

            - title :
                - title of the plot
                - (default) ''

            - xlabel :
                - string for x label
                - (default) $t/\Omega_{ci}^{-1}$ for types "rate","ratef","flux"
                  or $\Phi(t)$ for type 'ratevsflux'

            - ylabel :
                - (default) if type 'rate' $E_{rec}$
                - (default) if type 'ratef' $\partial \Phi/\partial t$
                - (default) if type 'flux' $\Phi(t)$
                - (default) if type 'ratevsflux' $\partial \Phi/\partial t$

            - labels :
                - list of labels for the N RecRate objects passed in *rates

            - filename :
                setting filename to (e.g.) 'rate.png' will save the figure
                to a file instead of displaying it to the screen


    """

    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax  = fig.add_subplot(111)


    # should it plot the reconnection elec. field
    # or the derivative of the magnetic flux
    # or the magnetic flux 
    if 'type' in kwargs:
        type = kwargs['type']
    else:
        type = 'rate'


    # plot the reconnection electric field
    if type.lower() == 'rate':
        print("RecRate > plotting reconnection electric field vs time")
        for rate  in rates:
            if np.mean(rate.rate) > 0:
                data = rate.rate
            else:
                data = -rate.rate
            ax.plot(rate.t, data)

    # plot the reconnected flux
    elif type.lower() == 'flux':
        print("RecRate > plotting reconnected flux vs time")
        for flux in rates:

            # search t=0
            id = np.where(np.abs(flux.t) == \
                      np.min(np.abs(flux.t)))

            #if np.mean(flux.flux-flux.flux[id]) < 0:
#           if np.mean(flux.flux-flux.flux[0]) < 0:
#               data = -flux.flux
#           else:
#               data = flux.flux
            #plt.plot(flux.t, data - data[id])
            ax.plot(flux.t, flux.flux)#data - data[0])


    # plot the derivative of the reconnected flux
    elif type.lower() == 'ratef':
        print("RecRate > plotting reconnection rate vs time (flux derivative)")
        for ratef in rates:
            hdt = 0.5*(ratef.t[1]-ratef.t[0])

            if np.mean(ratef.ratef) > 0:
                data = ratef.ratef
            else:
                data = -ratef.ratef
            ax.plot(ratef.t[1:]-hdt, data)


    elif type.lower() == 'ratevsflux':
        print("RecRate > plotting reconnection rate (flux derivative) vs reconnected flux")
        for rate in rates:
            id = np.where(np.abs(rate.t) == np.min(np.abs(rate.t)))

            if np.mean(rate.ratef) > 0:
                data = rate.ratef
            else:
                data = -rate.ratef
            if np.mean(rate.flux-rate.flux[id]) < 0:
                dataf = -rate.flux
            else:
                dataf = rate.flux
            ax.plot(dataf[1:]-dataf[1],data)


    # setup the range
    if 'range' in kwargs:
        range = kwargs['range']
    else:
        range = (0,np.max(rates[0].t), -0.1, 0.1)

    # setup the title
    if 'title' in kwargs:
        title = kwargs['title']
    else:
        title=''

    # setup the X label
    if 'xlabel' in kwargs:
        ax.set_xlabel(kwargs['xlabel'])
    else:
        if type.lower() in ['rate','ratef','flux']:
            ax.set_xlabel(r'$t/\Omega_{ci}^{-1}$')
        elif type.lower() in ['ratevsflux']:
            ax.set_xlabel(r'$\Phi(t)$')


    # setup the Y label
    # this depend on the type of plot chosen
    if 'ylabel' in kwargs:
        ax.set_ylabel(kwargs['ylabel'])
    else:
        if type == 'rate':
            ax.set_ylabel(r'$E_{rec}$')
        elif type == 'ratef':
            ax.set_ylabel(r'$\partial \Phi/\partial t$')
        elif type == 'flux':
            ax.set_ylabel(r'$\Phi(t)$')
        elif type.lower() == ['ratevsflux']:
            ax.set_ylabel(r'$\partial \Phi/\partial t$')

    # setup legend labels
    if 'labels' in kwargs:
        ax.legend(kwargs['labels'])

    ax.set_title(title)
    ax.axis(range)

    fig.tight_layout()

    # should we display or save to file?
    if 'filename' in kwargs:
        print(kwargs['filename'])
        fig.savefig(kwargs['filename'])
    else:
        fig.show()

    # close the window
    #fig.close()





