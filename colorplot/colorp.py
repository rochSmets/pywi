import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import numpy as np
import os
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable



#----------------------------------------------------------
#----------------------------------------------------------
def show(qty,
         x,y,
         vmin=None,
         vmax=None,
         extent = None,
         vecfield=None,
         veccolorcode=None,
         vecscale=1.,
         vecbin = 4,
         fieldlines=None,
         fieldlinesnb=50,
         xlabel='',
         ylabel='',
         title='',
         filename='colorplot.png',
         path = './',
         point = None,
         boxes = None,
         cmap = cm.jet,
         silent='no'):

    """ shows a color plot on screen

    this function makes a color plot

    Parameters
    ----------
    qty : ndarray
        2D scalar field of shape (n1,n2) to show
    x   : ndarray
        coordinate along axis 0, must be of size n1
    y   : ndarray
        coordinate along axis 1, must be of size n2

    Optional parameters
    -------------------
    vmin : scalar, optional
        minimum value of the color scale
        default is min(qty)
    vmax : scalar, optional
        maximum value of the color scale
        default is max(qty)
    extent : tuple or list, optional
        [x0,x1,y0,y1], default is [min(x),max(x),min(y),max(y)]


    vectfield : (ndarray, ndarray), optional
        shape of the arrays is (n1,n2)
        vectfield[0] and vectfield[1] are the x and y components of the vector
        field.
        default is no vector field plotted.
    veccolorcode : ndarray, optional
        arrows of the vector field will be colored according to the specified
        arrays. Shape must be (n1,n2).
    vecscale : scalar,optional
        scale of the arrows of the vector field
        default is 1.
    vecbin : scalar, optional
        rebin the vector field every vecbin points
        default is 4

    fieldlines : ndarray, optional
        magnetic flux array, shape == (n1,n2)
        will draw magnetic field lines as isocontours of this flux function
        works *only in 2D*
    fieldlinesnb, scalar, optional
        number of contours in the whole simulation domain
        default is 50

    xlabel : string, optional
        label on the X axis, default is ''
    ylabel : string, optional
        label on the Y axis, default is ''
    title : string, optional
        title for the plot

    filename : string, optional
        name of the file to which the plot will be save if
        the user calls returned_object.tofile()
        default is 'colorplot.png'
    path : string, optional
        path where the file will be saved if tofile() is called
        default is './'

    point, 


    Exemple
    -------


    Returns
    -------

    returns a ColorP object.


    Creation : 2013-01-19 15:50:45.075925

    """

    cplot = ColorP(qty,
                   x,y,
                   vmin=vmin,
                   vmax=vmax,
                   extent =extent,
                   vecfield=vecfield,
                   veccolorcode=veccolorcode,
                   vecscale=vecscale,
                   vecbin = vecbin,
                   fieldlines=fieldlines,
                   fieldlinesnb=fieldlinesnb,
                   xlabel=xlabel,
                   ylabel=ylabel,
                   title=title,
                   filename=filename,
                   path = path,
                   point = point,
                   boxes = boxes,
                   cmap  = cmap,
                   silent=silent)

    cplot.makefigure()
    cplot.show()

    return cplot
#==========================================================






#==========================================================
#==========================================================
def tofile(qty,
         x,y,
         vmin=None,
         vmax=None,
         extent = None,
         vecfield=None,
         veccolorcode=None,
         vecscale=1.,
         vecbin = 4,
         fieldlines=None,
         fieldlinesnb=50,
         xlabel='',
         ylabel='',
         title='',
         filename='colorplot.png',
         path = './',
         point = None,
         boxes = None,
         cmap  = cm.jet,
         silent='no'):

    """@todo: save a color plot to a file

    @param vmax @todo
    @return: @todo
    Exemple  : 
    Creation : 2013-01-19 15:50:45.075925

    """

    cplot = ColorP(qty,
                   x,y,
                   vmin=vmin,
                   vmax=vmax,
                   extent = extent,
                   vecfield=vecfield,
                   veccolorcode=veccolorcode,
                   vecscale= vecscale,
                   vecbin = vecbin,
                   fieldlines=fieldlines,
                   fieldlinesnb=fieldlinesnb,
                   xlabel=xlabel,
                   ylabel=ylabel,
                   title=title,
                   filename=filename,
                   path =path,
                   point = point,
                   boxes = boxes,
                   cmap  = cmap,
                   silent=silent)

    cplot.makefigure()
    cplot.tofile()

    return cplot
#==========================================================







class ColorP(object):
    """color plot of the quantity"""

    def __init__(self,
                 qty,
                 x,y,
                 vmin=None, vmax=None,
                 extent = None,
                 vecfield=None,
                 veccolorcode=None,
                 vecscale=1.,
                 vecbin = 4,
                 fieldlines=None,
                 fieldlinesnb=50,
                 xlabel='',
                 ylabel='',
                 title='',
                 filename='colorplot.png',
                 path = './',
                 point = None,
                 boxes = None,
                 cmap  = cm.jet,
                 silent='no'):

        """@todo: to be defined """

        # general information of the figure
        self._qty          = qty     # quantity color coded
        self._x            = x       # X and Y coordinates in physical units.
        self._y            = y
        self._title        = title
        self._xlabel       = xlabel
        self._ylabel       = ylabel
        self._cmap         = cmap    # colormap
        self._cb           = None    # colorbar
        self._fig          = None    # the Figure object
        self._im           = None    # color plot AxesImage
        self._ax           = None


        # choose the minmax as extent if user didn't specify
        if extent is None:
            self._extent = (np.min(x), np.max(x),
                           np.min(y), np.max(y))
        else:
            self._extent = extent


        # add a point
        self._point = point

        #add rectangles
        self._boxes_corners = boxes
        self._rects         = []

        # unless specified by the user
        # the minmax values of the color are
        # the minmax values of the quantity that is plotted
        if vmin is None:
            self._vmin = np.min(qty)
        else :
            self._vmin = vmin
        if vmax is None:
            self._vmax = np.max(qty)
        else :
             self._vmax         = vmax


        # vector field
        self._vecfield     = vecfield        # components of the vector field
        self._vf_colorcode = veccolorcode
        self._vf_cmap      = cm.Greys        # color map
        self._vf_title     = ''              # title of the vector field
        self._vf_cb        = None            # colorbar
        self._binning      = vecbin          # binning of the vector field
        self._vf_im        = None
        self._vf_scale     = vecscale



        # magnetic field lines
        self._fldl         = fieldlines      # magnetic flux
        self._levels       = []              # levels of the magnetic flux
        self._fldlnb       = fieldlinesnb    # number of field lines


        # saving information
        self._path        = path
        self._filename    = filename






    #==========================================================
    #==========================================================
    def makefigure(self):
        """make the figure for the color plot

        Exemple  :

        Creation : 2012-08-16 08:53:03.629552

        Note : must keep using pyplot for creating figures
        otherwise routines show() do not work although they are methods
        of the figure.Figure

        """
        self._fig = plt.figure()
        self._ax = self._fig.add_subplot(111)

        # setup the colorplot
        dx = self._x[1]-self._x[0]
        dy = self._y[1]-self._y[0]

        i1 =int((self._extent[0]-self._x[0])/dx)
        i2 =int((self._extent[1]-self._x[0])/dx)
        j1 =int((self._extent[2]-self._y[0])/dy)
        j2 =int((self._extent[3]-self._y[0])/dy)


        self._im = self._ax.imshow(np.transpose(self._qty[i1:i2,j1:j2]),
                       cmap          = self._cmap,
                       aspect        = 'equal',
                       interpolation = 'nearest',
                       vmin          = self._vmin,
                       vmax          = self._vmax,
                       extent        = self._extent,
                       origin        = 'lower')

        self._ax.set_xlabel(self._xlabel)
        self._ax.set_ylabel(self._ylabel)
        self._ax.set_title(self._title)


        # so that colorbar fits the height of the main plot
        divider = make_axes_locatable(self._ax)
        cax     = divider.append_axes("right", size="5%", pad=0.5)
        self._fig.colorbar(self._im, cax)

        self._ax.set_xlim(self._extent[0], self._extent[1])
        self._ax.set_ylim(self._extent[2], self._extent[3])

        # setup the vector field
        if self._vecfield is not None:
            self._vectorfield()

        # setting up the flux
        if self._fldl is not None:
            self.fieldlines()

        # point
        if self._point is not None:
            art = mpatches.Circle(self._point, radius=0.1)
            self._ax.add_patch(art)

        # boxes
        if len(self._rects) != 0: # re-create old boxes
            self._restore_boxes()

        if self._boxes_corners != None:  # new boxes added
            self.add_boxes(self._boxes_corners)
            self._boxes_corners = None   #later calls to Makefigure won't add these boxes again

        self._fig.tight_layout()

#==========================================================




    #==========================================================
    #==========================================================
    def show(self):
        """@todo: show the figure on screen

        @return: @todo

        Exemple  : 

        Creation : 2013-01-19 20:16:53.952742

        """

        self._fig.show()
    #==========================================================





    #==========================================================
    #==========================================================
    def tofile(self):
        """@todo: save the figure to a file

        @return: @todo

        Exemple  : 

        Creation : 2013-01-24 11:22:02.034870

        """

        self._fig.savefig(os.path.join(self._path,self._filename),
                          dpi = 400, bbox_inches='tight')
        #plt.close(self._fig) # so that next figure is not an overplot
    #==========================================================




    #==========================================================
    #==========================================================
    def display(self):
        """display the object parameters


        Exemple  : 

        Creation : 2012-08-16 08:28:35.738418

        """

        print("qty       : ", self.qty.shape)
        print("extnent   : ", self.extent)
        print("title     : ", self.title)
        print("labels    : ", self.xlabel, self.ylabel)


        print("")
        print("magnetic field lines")
        print("levels    : ", self.levels)

        print("")
        print("saving info")
        print("filename  : ", self.filename)
        print("path      : ", self.path)


    #==========================================================




    #==========================================================
    #==========================================================
    def add_boxes(self, corners):
        """add a rectangular box on top of the colorplot

        @param corners (x0,y0,x1,y1)

        @return: @todo

        Exemple  : cplot.add_box((x0,y0,x1,y1))

        Creation : 2013-01-24 11:41:22.299075

        """
        for box in corners:
            art = mpatches.Rectangle((box[0],box[1]),
                                     box[2]-box[0], box[3]-box[1],
                                     ec='k',fc='none')
            self._ax.add_patch(art)
            self._rects.append(box)
    #==========================================================



    #==========================================================
    #==========================================================
    def _restore_boxes(self):
        """ recreate rectangles

        @return: @todo

        Exemple  : 

        Creation : 2013-01-24 12:39:44.039518

        """
        for box in self._rects:
            art = mpatches.Rectangle((box[0],box[1]),
                                     box[2]-box[0], box[3]-box[1],
                                     ec='k',fc='none')
            self._ax.add_patch(art)

    #==========================================================




    #==========================================================
    #==========================================================
    def extent(self, *args):
        """set the extent fo the plot in physical units. Returns the
        current extent. If no argument is passed returns the current extent.

        @param extent exemple : (x0,y0,x1,y1)

        @return: nothing

        Exemple  :  self.extent()
                    self.extent((12,24,30,35))

        Creation : 2012-08-15 17:55:24.079682

        """
        if (len(args) !=0 ):
            self.extent = args[0]

        return self.extent
    #==========================================================





    #==========================================================
    #==========================================================
    def title(self, *args):
        """ set and/or returns the current title of the color plot

        @param args : string

        @return: current title of the color plot

        Exemple  :  self.title()
                    self.title('new title')

        Creation : 2012-08-15 17:58:52.079867

        """
        if (len(args) != 0):
            self.title = args[0]

        return self.title
    #==========================================================





    #==========================================================
    def xlabel(self, *args):
        """ set and/or returns the current xlabel of the color plot

        @param args : string

        @return: current xlabel of the color plot

        Exemple  :  self.xlabel()
                    self.xlabel('new xlabel')

        Creation : 2012-08-15 17:58:52.079867

        """
        if (len(args) != 0):
            self.xlabel = args[0]

        return self.xlabel

    #==========================================================







    #==========================================================
    def ylabel(self, *args):
        """ set and/or returns the current ylabel of the color plot

        @param args : string

        @return: current ylabel of the color plot

        Exemple  :  self.ylabel()
                    self.ylabel('new ylabel')

        Creation : 2012-08-15 17:58:52.079867

        """
        if (len(args) != 0):
            self.ylabel = args[0]

        return self.ylabel

    #==========================================================






    #==========================================================
    def colormap(self, *args):
        """ set and/or returns the current cmap of the color plot

        @param args : string

        @return: current cmap of the color plot

        Exemple  :  self.cmap()
                    self.cmap('new cmap')

        Creation : 2012-08-15 17:58:52.079867

        """
        if (len(args) != 0):
            self.cmap = args[0]

        return self.cmap

    #==========================================================






    #==========================================================
    def _vectorfield(self):
        """attache vector field attached to
        the plot.

        Calling sequences :

        @return: @todo

        Exemple  :

        Creation : 2012-08-15 17:44:41.917733

        """

        Xn = self._x[::self._binning]
        Yn = self._y[::self._binning]


        U = self._vecfield[0]
        V = self._vecfield[1]

        Un = U[::self._binning, ::self._binning]
        Vn = V[::self._binning, ::self._binning]



        if self._vf_colorcode is not None:

            veccol = self._vf_colorcode[::self._binning,::self._binning]

            self.vf_im = self.ax.quiver(Xn, Yn,
                                    np.transpose(Un),
                                    np.transpose(Vn),
                                    np.transpose(veccol), 
                                    angles='xy',
                                    units='xy',
                                    scale=self._vf_scale,
                                    scale_units='xy',
                                    pivot='middle',
                                    cmap=self._vf_cmap)

#            self.fig.colorbar(self.vf_im)

        else:
            self.vf_im  = self._ax.quiver(Xn, Yn,
                                    np.transpose(Un),
                                    np.transpose(Vn),
                                    units='xy',
                                    angles='xy',
                                    scale=self._vf_scale,
                                    scale_units='xy',
                                    pivot='middle',
                                    cmap=self._vf_cmap)


    #==========================================================






    #==========================================================
    #==========================================================
    def fieldlines(self):
        """ display the field lines

        Exemple  : 

        Creation : 2012-08-16 12:54:58.854679

        """

        Fmin = np.min(self._fldl)
        Fmax = np.max(self._fldl)

        self._levels = np.arange(Fmin,Fmax,(Fmax-Fmin)/self._fldlnb)

        matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
        self._ax.contour(self._x, self._y, np.transpose(self._fldl),
                        colors='k', cmap=None, levels=self._levels,
                        extent=self._extent)


    #==========================================================



