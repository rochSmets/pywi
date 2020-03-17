

import numpy as np
import pypic.picgsfc as pr
import pypic.colorp as colp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import matplotlib.patches as mpatches

# for coloring the lines
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize


class trajplot(colp.ColorP):

    #==========================================================
    #==========================================================
    def __init__(self,
                 qty,
                 x,y,
                 tprun,
                 color=None,
                 vmin=None,vmax=None,
                 extent = None,
                 vecfield = None,
                 veccolorcode = None,
                 vecscale=1.,
                 vecbin = 4,
                 fieldlines=None,
                 fieldlinesnb=50,
                 xlabel='',
                 ylabel='',
                 title='',
                 filename='colorplot.png',
                 path='./',
                 point=None,
                 boxes=None,
                 silent='no'):

        """constructor
        Exemple  : 

        Creation : 2013-05-03 13:55:49.356526

        """
        super(trajplot,self).__init__(qty,
                                      x,y,
                                      vmin=vmin,
                                      vmax=vmax,
                                      extent=extent,
                                      vecfield=vecfield,
                                      veccolorcode=veccolorcode,
                                      vecscale=vecscale,
                                      vecbin=vecbin,
                                      fieldlines=fieldlines,
                                      fieldlinesnb=fieldlinesnb,
                                      xlabel=xlabel,
                                      ylabel=ylabel,
                                      title=title,
                                      filename=filename,
                                      path=path,
                                      point=point,
                                      boxes=boxes,
                                      silent=silent)



        self._tprun   = tprun
        self._color   = color
        self.makefigure()
        self.tofile()

    #==========================================================


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

        i1 = (self._extent[0]-self._x[0])/dx
        i2 = (self._extent[1]-self._x[0])/dx
        j1 = (self._extent[2]-self._y[0])/dy
        j2 = (self._extent[3]-self._y[0])/dy


        self._im = self._ax.imshow(np.transpose(self._qty[i1:i2,j1:j2]),
                       cmap =self._cmap,
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


        if len(self._tprun.r.shape) == 3:
            npart = self._tprun.r.shape[1]

            for p in range(npart):

                print 'adding trajectory %d...' % p

                x = self._tprun.r[0,p,::100]
                y = self._tprun.r[1,p,::100]

                points = np.array([x, y]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)

                mycolor = self._color[p,::100]

                vmin =  mycolor.min()
                vmax =  mycolor.max()


                lc = LineCollection(segments,
                                    cmap=plt.get_cmap('Greys_r'),
                                    norm=Normalize(vmin,vmax))

                lc.set_array(mycolor) # for now the color is y

                lc.set_linewidth(1)
                self._ax.add_collection(lc)

                # add an arrow at the end of the trajectory
                # take the in-plane velocity
                v1 = self._tprun.v[0,p,-1]
                v2 = self._tprun.v[1,p,-1]
                vnorm = np.sqrt(v1**2 + v2**2)

                arrow = mpatches.Arrow(x[-1],y[-1],
                                       0.5*v1/vnorm,0.5*v2/vnorm,
                                       width=0.1,color='k')

                self._ax.add_patch(arrow)
        else:
            x = self._tprun.r[0,:]
            y = self._tprun.r[1,:]

            points = np.array([x, y]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

            color = self._color

            vmin = color.min()
            vmax = color.max()

            lc = LineCollection(segments,
                                cmap=plt.get_cmap('gray'),
                                norm=plt.Normalize(vmin=vmin,vmax=vmax))

            lc.set_array(self._color) # for now the color is y

            lc.set_linewidth(1)
            self._ax.add_collection(lc)


        # add a point for the selection time
        art = mpatches.Circle((self._tprun.r[0,0,self._tprun._it0],
                               self._tprun.r[1,0,self._tprun._it0]),
                               radius=0.1,
                               ec='#fe01a1',fc='none',zorder=200)
        self._ax.add_patch(art)



        self._fig.tight_layout()

#==========================================================



