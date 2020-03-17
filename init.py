#!/opt/local/bin/python
# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt






class CondInit(object):


    def __init__(self, domsize,ncells):
        self.L      = domsize
        self.ncells = ncells
        self.resize(domsize,ncells)


    def resize(self,domsize,ncells):
        dl0 = domsize[0]/float(ncells[0]-2.)
        dl1 = domsize[1]/float(ncells[1]-2.)
        self.dl = (dl0,dl1)
        self.domsize = domsize
        self.ncells  = ncells
        self.x0 = -self.dl[0] + np.arange(ncells[0])*self.dl[0]
        self.x1 = -self.dl[1] + np.arange(ncells[1])*self.dl[1]


    def debye(self):
        return np.sqrt(self.Te()/self.density())/self.c

    def wpe(self):
        return np.sqrt(1./self.memi)/self.c

    def Te(self):
        return np.ones(self.ncells[1])


    def Ti(self):
        return np.ones(self.ncells[1])


    def B(self):
        return np.ones(self.ncells[1])


    def density(self):
        return np.ones(self.ncells[1])


    def dlmax(self):
        return 1.5*np.min(self.debye())

    def cfl_light(self):
        return 1./self.c *1./np.sqrt(1./self.dl[0]**2  + 1./self.dl[1]**2)

    def __str__(self):
        st =  'INITIALISATION\n'
        st += '--------------\n'

        if min(self.dl) < self.dlmax():
            st += 'mesh size (%f,%f) can be increased\n' % self.dl

        if min(self.dl) > self.dlmax():
            st += 'WARNING : mesh size (%f,%f) is larger than the recommended one\n' % self.dl

        st += 'Maximum Mesh Size : 1.5*min(Debye) = 1.5*%f = %f\n' % (np.min(self.debye()),self.dlmax())
        st += 'Plasma time step should be smaller than 0.1min(wpe^-1) = %f\n' % (0.1*np.min(self.wpe()))
        st += 'Field time step should be smaller than CFL_Light = %f, recommended 0.3CFL = %f\n' % (self.cfl_light(), 0.3*self.cfl_light())
        return st



    def plot(self):
        fig = plt.figure()
        ax  = fig.add_subplot(111)
        ax.plot(self.x1, self.B())
        ax.plot(self.x1, self.density())
        ax.plot(self.x1, self.Ti())
        ax.plot(self.x1, self.Te())
        ax.set_ylim( (-2., 2.) )
        fig.savefig('init.eps')







class Thermal(CondInit):


    #==========================================================
    #==========================================================
    def __init__(self,teti=0.2,me=0.01,c=25.,n=1.,ttot=0.5,
                 domsize=(10.,10.),ncells=(1000, 1000)):


        super(Thermal, self).__init__(domsize,ncells)
        self.memi    = me
        self.domsize = domsize
        self.ncells  = ncells
        self.teti    = teti
        self.c       = c
        self.n       = n
        self.ttot    = ttot
    #==========================================================

    def density(self):
        return np.ones(self.ncells[1])

    def Te(self):
        return self.teti*self.ttot/(self.teti+1.) * np.ones(self.ncells[1])


    def Ti(self):
        return self.Te()*self.ttot/self.teti

    def B(self):
        return np.ones(self.ncells[1])






class DoubleHarris(CondInit):


    #==========================================================
    #==========================================================
    def __init__(self,teti=0.2,nb=0.2,me=0.01,c=25.,l=1.,ttot=0.5,
                 domsize=(10.,10.),ncells=(1000, 1000)):


        super(DoubleHarris, self).__init__(domsize,ncells)
        self.memi    = me
        self.domsize = domsize
        self.ncells  = ncells
        self.teti    = teti
        self.c       = c
        self.l       = l
        self.nb      = nb
        self.ttot    = ttot
    #==========================================================




    def density(self):
        return 1./np.cosh((self.x1-self.x1.max()*0.25)/self.l)**2\
             + 1./np.cosh((self.x1-self.x1.max()*0.75)/self.l)**2\
             + self.nb


    def Te(self):
        return self.teti*self.ttot/(1.+self.teti) + np.zeros(self.x1.size)


    def Ti(self):
        return self.ttot/(1.+self.teti) + np.zeros(self.x1.size)


    def B(self):
        return  np.tanh((self.x1-self.x1.max()*0.25)/self.l)\
              - np.tanh((self.x1-self.x1.max()*0.75)/self.l)









class Couche1(CondInit):


    #==========================================================
    #==========================================================
    def __init__(self,te=0.2,beta=1.,me=0.01,c=25.,l=1.,
                 domsize=(10.,10.),ncells=(1000, 1000)):


        super(Couche1, self).__init__(domsize,ncells)
        self.memi    = me
        self.domsize = domsize
        self.ncells  = ncells
        self.te      = te
        self.c       = c
        self.beta    = beta
        self.l       = l
    #==========================================================

    def density(self):
        return np.ones(self.ncells[1])


    def Te(self):
        return self.te*np.ones(self.ncells[1])


    def Ti(self):
        b2 = self.B()**2
        cte = (self.beta + 1.)* 1.**2/2.

        return (cte - b2/2.)/self.density() - self.Te()

    def B(self):
        x1mid = 0.5*self.x1.max()
        return np.tanh((self.x1-x1mid)/self.l)



