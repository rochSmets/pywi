# this module allows distribution function manipulations


import numpy as np
import os
import particles
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as mpatches



#==========================================================
#==========================================================
def show(particles,
         coords        = None,
         magfield      = None,
         elecfield     = None,
         Tpp           = None,
         log           = False,
         bins          = None,
         vrange        = None,
         filename      = None,
         path          = None,
         interpolation = None,
         gyrotropy     = False,
         ptest_select  = None,
         ptest_rv      = None,
         ax            = None):
    """display a distribution function

    display a distribution function of the given 'particles',
    in a given coordinate system (coords), e.g. 'VxVy'. It can
    either display it on screen or save it to a file

    @param particles @todo
    @param bins @todo
    @param coords @todo
    @param filename @todo
    @param path @todo
    @param gyrotropy

    @return: @todo

    Exemple  :

    Creation : 2013-01-15 20:50:27.792462

    """
    dist = Distribution(particles,
                        coords        = coords,
                        magfield      = magfield,
                        elecfield     = elecfield,
                        Tpp           = Tpp,
                        log           = log,
                        bins          = bins,
                        vrange        = vrange,
                        filename      = filename,
                        path          = path,
                        interpolation = interpolation,
                        gyrotropy     = gyrotropy,
                        ptest_select  = ptest_select,
                        ptest_rv      = ptest_rv)

    dist.show(ax=ax)

    return dist
#==========================================================






class Distribution(object):
    """this is a class for a distribution function"""

    #==========================================================
    #==========================================================
    def __init__(self,
                 particles           ,
                 coords        = None,
                 magfield      = None,
                 elecfield     = None,
                 Tpp           = None, 
                 log           = False,
                 bins          = None,
                 vrange        = None,
                 filename      = None,
                 path          = None,
                 interpolation = None,
                 gyrotropy     = False,
                 ptest_select  = None,
                 ptest_rv      = None):
        """ initializes a distribution function object

        @param Particles is the group of particles we want the distribution
               from
        @param bin is the binning in velocity units
        @param coords is the coordinate system, like 'VxVy'
        Exemple  :
        Creation : 2013-01-15 12:11:48.786137
        """

        self._particles  = particles
        self._filename   = filename
        self._path       = path
        self._magfield   = magfield
        self._elecfield  = elecfield
        self._Tpp        = Tpp
        self._gyrotropy  = gyrotropy

    #---------------interpolation-mode-----------------------------------------
        if interpolation == None:
            self._interp = 'nearest'
        elif interpolation.lower() in ['nearest','bilinear','bicubic']:
            self._interp = interpolation
        else:
            print "unknown interpolation mode"
            return None

    #---------------coordinate-system------------------------------------------

        accepted_coords = ['vxvy','vxvz','vyvz','v1v2','vav1','vav2']

        if coords == None:
            self._coords = 'VxVy'
            self._c = (0,1)

        elif coords.lower() not in accepted_coords:
            print "Unknown coordinate : %s" % (coords)
            return None

        else:
            self._coords = coords
            if coords.lower()   == 'vxvy':
                self._c = (0,1)
            elif coords.lower() == 'vxvz':
                self._c = (0,2)
            elif coords.lower() == 'vyvz':
                self._c = (1,2)
    #--------------------------------------------------------------------------



    #-----------binning-and-range----------------------------------------------

        if bins != None:
            self._bins      = bins
        else:
            self._bins = (60,60) # arbitrary and needs to be fixed

        self._vrange = vrange

        self._dv = ((vrange[0][1]-vrange[0][0])/self._bins[0],
                    (vrange[1][1]-vrange[1][0])/self._bins[1])

    #--------------------------------------------------------------------------



    #--------histogram---------------------------------------------------------
        self._log = log
        self._calculate()
    #--------------------------------------------------------------------------


    #--------ptest_select------------------------------------------------------
        self._ptest_select = ptest_select
    #--------------------------------------------------------------------------


    #--------trajectories------------------------------------------------------
        self._ptest_rv = ptest_rv
    #--------------------------------------------------------------------------


    #==========================================================






    #==========================================================
    #==========================================================
    def _calculate(self):
        """@todo: private function that calculates the histogram

        Exemple  : 

        Creation : 2013-01-16 11:10:51.171984

        """
        hist,x,y = np.histogram2d(self._particles._vel[self._c[0]],
                                  self._particles._vel[self._c[1]],
                                  bins=self._bins,
                                  range=self._vrange)

        self._hist      = hist
        self._x         = x
        self._y         = y

        if self._log == True:
            idno0 = np.where(self._hist >= 1)
            self._hist[idno0[0],idno0[1]] = \
                    np.log10(self._hist[idno0[0],idno0[1]])
    #==========================================================






    #==========================================================
    #==========================================================
    def set_interpolation(self, interpolation):
        """set the interpolation mode

        @param interpolation can be 'bilinear' or 'nearest' or 'bicubic'

        Exemple  : 

        Creation : 2013-01-16 11:44:55.441053

        """
        if interpolation in ['nearest','bilinear','bicubic']:
            self._interp = interpolation
    #==========================================================





    #==========================================================
    #==========================================================
    def set_bins(self, bins):
        """set the bins value

        @param bins is the value for the velocity bins, e.g. (10,20)
               Note that set_bins() recalculates the histgram with the new bins
        @return: @todo

        Exemple  : 

        Creation : 2013-01-15 14:19:20.388793

        """
        #TODO more checking of 'bins'
        self._bins= bins
        self._calculate()
    #==========================================================





    #==========================================================
    #==========================================================
    def set_coords(self, coords):
        """set the coordinate system

        @param coords is the coordinate system in which the distribution
        is plotted

        @return: @todo

        Exemple  : 

        Creation : 2013-01-15 14:20:42.685017

        """
        # TODO more security check on that
        self._coords = coords
    #==========================================================






    #==========================================================
    #==========================================================
    def show(self, ax = None):
        """show the distribution function

        Exemple  : 

        Creation : 2013-01-15 13:46:08.940382

        """


        if ax == None:
            fig   = plt.figure()
            axe   = fig.add_subplot(111)
        else:
            axe = ax
            fig = ax.get_figure()

        extent = (self._x[0],self._x[-1],self._y[0],self._y[-1])


        im = axe.imshow(np.transpose(self._hist),
                       aspect         = 'equal',
                       interpolation  = self._interp,
                       vmin           = self._hist.min(),
                       vmax           = self._hist.max(),
                       extent         = extent,
                       origin         = 'lower',
                       cmap           = 'gist_yarg')

        # if there is a magnetic field, plot an arrow
        if self._magfield != None:
            pass

        # add the colorbar
        divider = make_axes_locatable(axe)
        cax     = divider.append_axes("right", size="5%", pad=0.5)
        fig.colorbar(im, cax)


        # add axis labels
        if self._coords == None:
            xlabel = r'$V_x$'
            ylabel = r'$V_y$'
        elif self._coords.lower() == 'vxvy':
            xlabel = r'$V_x$'
            ylabel = r'$V_y$'
        elif self._coords.lower() == 'vxvz':
            xlabel = r'$V_x$'
            ylabel = r'$V_z$'
        elif self._coords.lower() == 'vyvz':
            xlabel = r'$V_y$'
            ylabel = r'$V_z$'

        axe.set_xlabel(xlabel)
        axe.set_ylabel(ylabel)


        # add (0,0) axis
        axe.hlines(0,extent[0],extent[1],linestyle='solid',color='#5ec23d')
        axe.vlines(0,extent[2],extent[3],linestyle='solid',color='#5ec23d')


        # add an arrow for the mean in-plane velocity
        v1av = np.mean(self._particles._vel[self._c[0]])
        v2av = np.mean(self._particles._vel[self._c[1]])
        mine = min([abs(L)  for L in extent])
        vnorm = 1/np.sqrt(v1av**2 + v2av**2)* 1./2.*mine
        #print v1av*vnorm,v2av*vnorm,v1av,v2av
        axe.arrow(v1av,v2av,v1av*vnorm, v2av*vnorm,
                head_width  = 0.5,
                head_length = 0.5,
                linewidth   = 3.5,
                length_includes_head = True,
                fc='#ff0000', ec='#ff0000',zorder=101)


        # add bulk velocity frame axis
        axe.hlines(v2av,extent[0],extent[1],linestyle='dashed',color='r')
        axe.vlines(v1av,extent[2],extent[3],linestyle='dashed',color='r')


        # add title saying the # of particles
        #TODO fix : .size with particles.size()
#       title = '%d %s' % (self._particles._pos[0].size,
                                            #self._particles.getspecies())

#       title = title + r', $|V_{plane}| = %5.3f$' % (np.sqrt(v1av**2 + v2av**2))

        vx = np.mean(self._particles._vel[0])
        vy = np.mean(self._particles._vel[1])
        vz = np.mean(self._particles._vel[2])
        v  = np.sqrt(vx**2 + vy**2 + vz**2)
#       title = title + r', $|V| = %5.3f$' % (v)

        # add temperatures
        T1 = np.std(self._x)
        T2 = np.std(self._y)
#       title = title + r', $T_1 = %5.3f$, $T_2 = %5.3f$' %(T1,T2)


        # add average position
        pos1 = np.mean(self._particles._pos[0,:])
        pos2 = np.mean(self._particles._pos[1,:])
#       title = title + r', $(%5.3f,%4.3f)$' %(pos1,pos2)


        # now add an arrow for the magnetic field
        if self._magfield != None:
            # find position of particles in average
            # interpolate the field at this position
            # normalize the field to the modulus of V (each com * |V|/|B|)
            # draw the arrow

            mine = min([abs(L)  for L in extent])
            bnorm = 1/np.sqrt(self._magfield[self._c[0]]**2 \
                             +self._magfield[self._c[1]]**2)* 1./2.*mine

            #print v1av*vnorm,v2av*vnorm,v1av,v2av
            axe.arrow(v1av,v2av,self._magfield[self._c[0]]*bnorm,
                    self._magfield[self._c[1]]*bnorm,
                    head_width  = 0.5,
                    head_length = 0.5,
                    linewidth   = 3.5,
                    length_includes_head=True,
                    fc='#ff00ff', ec='#ff00ff',zorder=100)

            #angle B^V
            pass

        # puts a small rectangle at the ExB velocity
        if self._elecfield != None and self._magfield != None:

            b2    = (self._magfield[0]**2 \
                  +    self._magfield[1]**2 \
                  +    self._magfield[2]**2)

            ExB = [0,0,0]
            ExB[0] = (self._elecfield[1]*self._magfield[2] - \
                      self._elecfield[2]*self._magfield[1])/b2

            ExB[1] = (self._elecfield[2]*self._magfield[0] - \
                      self._elecfield[0]*self._magfield[2])/b2

            ExB[2] = (self._elecfield[0]*self._magfield[1] - \
                      self._elecfield[1]*self._magfield[0])/b2

            # the square we create will be 4 dv by 4dv
            center = (ExB[self._c[0]],ExB[self._c[1]])

            radius = 2*(self._dv[1] + self._dv[0])/2.

            #we draw a blue rectangle
            art = mpatches.Circle(center, radius,ec="none",fc='b',zorder=200)
            axe.add_patch(art)


        # ---------------------------------------------------------------------
        # overplot les trajectoires dans l'espace des vitesses
        if self._ptest_rv != None:

            r,v = self._ptest_rv


            #TODO npart is private... shouldn't be used
            #- 2013-05-20 17:10:42.737331
            nptest      = v.shape[1]
            select_area = self._particles.GetArea()

            for ip in range(nptest):

                # first find the times where the particle
                # is in the selection area
                times = np.where((r[0,ip,:] > select_area[0][0])  &
                                 (r[0,ip,:] < select_area[0][1])  &
                                 (r[1,ip,:] > select_area[1][0])  &
                                 (r[1,ip,:] < select_area[1][1]))


                vp = v[:,ip,times[0]]

                axe.plot(vp[self._c[0],:],vp[self._c[1],:],
                        color='#01fed8', lw=1.4)

        # ---------------------------------------------------------------------


        # add a small rectangle on the selection area of the test particles 
        # if any is specified
        if self._ptest_select != None:
            lowerleft = self._ptest_select[0:2]
            width  = self._ptest_select[2]
            height = self._ptest_select[3]
            rect = mpatches.Rectangle(lowerleft,width,height,fc='none',
                                      ec='#ffbf00',zorder=1000)
            axe.add_patch(rect)

        # ---------------------------------------------------------------------


        #axe.set_title(title)


        # add the gyrotropy ellipse
        if self._Tpp != None and self._magfield != None:

            # find the angle between B and the X axis whatever it is
            Bnorm = self._magfield[0]**2 + self._magfield[1]**2 + self._magfield[2]**2
            Bnorm = np.sqrt(Bnorm)
            Bopp  = self._magfield[self._c[1]]
            sine  = Bopp/Bnorm
            alpha = np.arcsin(sine)*180./np.pi

            T_para, T_perp = self._Tpp
            minv1,minv2 = min(self._vrange[0]),min(self._vrange[1])
            minv = min((minv1,minv2))

            eh = T_perp *0.5*minv/T_para
            ew = 0.5*minv

            eh = T_perp
            ew = T_para

            ellipse = mpatches.Ellipse((v1av,v2av),
                                        3*ew,
                                        3*eh,
                                        angle=alpha,
                                        ec = '#f5ba0a',lw=2.0, fc='none')

            axe.add_patch(ellipse)


        # now check if we have to show on screen
        # or save to file
        # or just nothing if the ax comes from somewhere else

        if ax == None:
            if self._filename != None: # we save to a file
                fig.savefig(os.path.join(self._path,self._filename))
            else:
                fig.show()


    #==========================================================



