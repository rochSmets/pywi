#
#import matplotlib as mpl
#mpl.use('Agg')
#
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator
from matplotlib import rc
import operator
#
#
class Colp :

    """

    """

    def __init__(self,
                 coloraxis = None,
                 colordata = None,
                 bounds = None,
                 colormap = None,
                 contouraxis = None,
                 contourdata = None,
                 flines = 20,
                 arrowaxis = None,
                 arrowdata = None,
                 arrownorm = 1.0,
                 labels = None,
                 ticks = None,
                 subticks = None,
                 colorbar = True,
                 text = None,
                 xytext = None,
                 figsize = [5.0, 3.0],
                 filename = None,
                 filetype = None):


        # ..
        self.coloraxis = coloraxis
        self.colordata = colordata
        self.contouraxis = contouraxis
        self.contourdata = contourdata
        self.arrowaxis = arrowaxis
        self.arrowdata = arrowdata

        self.arrownorm = arrownorm

        # ..
        if bounds is not None :
          self.bounds = bounds

        elif colordata is not None :
          self.bounds = [np.min(self.colordata), np.max(self.colordata)]

        elif contourdata is not None :
          self.bounds = [np.min(self.contourdata), np.max(self.contourdata)]

        elif arrowdata is not None :
          self.bounds = [np.min(self.arrowdata), np.max(self.arrowdata)]

        else :
           raise ValueError( "what the fuck you wanna plot my man ?")

        # ..
        if colormap == None :
          self.colormap = cm.get_cmap('jet', 64)
        else :
          self.colormap = cm.get_cmap(colormap[0], colormap[1])

        # ..
        self.contourdata = contourdata

        # ..
        if flines == None :
           self.flines = 20

        else :
           self.flines = flines

        # ..
        if arrowdata is not None :
          self.xarrowdata = arrowdata[0]
          self.yarrowdata = arrowdata[1]
          self.arrowdata = True
        else :
          self.arrowdata = False

        # ..
        self.labels = labels

        # .. set the position of ticks
        self.ticks = ticks

        # .. set the # of subticks
        self.subticks = subticks

        # .. set the colorbar if wanted
        self.colorbar = colorbar

        # .. set the text
        self.text = text
        self.xytext = xytext

        # ..
        self.figsize = figsize

        # ..
        if filetype is not None :
          self.filetype = filetype.lower()
        else :
          self.filetype = None

        # ..
        if filename is not None:
            self.filename = filename
        else:
            self.filename = 'Colp'

        # ..
        self.Draw()


    def Draw(self):

        """
        draw the plot
        """

        # .. use latex fonts
        rc('text', usetex = True)
        rc('font', size=12)
        rc('axes', labelsize='larger')
        rc('mathtext', default='regular')

        # .. just 1 plot
        fig, ax = plt.subplots(num = 0, figsize = self.figsize, dpi = 100)

        # .. draw the color image
        if self.colordata is not None :
            im = ax.imshow(np.transpose(self.colordata),
                           aspect = 'auto',
                           interpolation = 'nearest',
                           cmap = self.colormap,
                           origin = 'lower',
                           extent = [self.coloraxis[0][0], self.coloraxis[0][-1],
                                     self.coloraxis[1][0], self.coloraxis[1][-1]],
                           vmin = self.bounds[0],
                           vmax = self.bounds[1])

        # .. draw the isocontour
        if self.contourdata is not None:
          co = ax.contour(self.contouraxis[0],
                          self.contouraxis[1],
                          np.transpose(self.contourdata),
                          self.flines,
                          colors = ('k',),
                          origin = 'lower',
                          extent = [self.contouraxis[0][0], self.contouraxis[0][-1],
                                    self.contouraxis[1][0], self.contouraxis[1][-1]],
                          linestyles = 'solid',
                          linewidths = 1)

        # .. draw the arrow field
        if self.arrowdata is True :
            ar = ax.quiver(self.arrowaxis[0],
                           self.arrowaxis[1],
                           np.transpose(self.xarrowdata)*self.arrownorm,
                           np.transpose(self.yarrowdata)*self.arrownorm)

        # .. set the ticks, the associated labels & the # of subticks
        if self.ticks is not None :
            ax.xaxis.set_major_locator(MultipleLocator(self.ticks[0]))
            ax.yaxis.set_major_locator(MultipleLocator(self.ticks[1]))

        if self.subticks is not None :
            ax.xaxis.set_minor_locator(MultipleLocator(self.subticks[0]))
            ax.yaxis.set_minor_locator(MultipleLocator(self.subticks[1]))

        # .. set the limits of the plot
        if self.coloraxis is not None :
           ax.set_xlim([self.coloraxis[0][0], self.coloraxis[0][-1]])
           ax.set_ylim([self.coloraxis[1][0], self.coloraxis[1][-1]])

        elif self.contouraxis is not None :
           ax.set_xlim([self.contouraxis[0][0], self.contouraxis[0][-1]])
           ax.set_ylim([self.contouraxis[1][0], self.contouraxis[1][-1]])

        elif self.arrowaxis is not None :
           ax.set_xlim([self.arrowaxis[0][0], self.arrowaxis[0][-1]])
           ax.set_ylim([self.arrowaxis[1][0], self.arrowaxis[1][-1]])

        else :
           raise ValueError( "what the fuck you wanna plot my man ?")

        # .. write the labels
        if self.labels is not None :
            ax.set_xlabel(self.labels[0])
            ax.set_ylabel(self.labels[1])

        if self.colorbar == True:
            # .. 3 ticks for the colorbar
            ticks = np.linspace(self.bounds[0], self.bounds[1], num = 3)

            # .. set the colorbar
            cbar = fig.colorbar(im, ticks = ticks, pad = 0.03, aspect = 40)

        # .. write the text
        if self.text is not None :
            if self.xytext is not None :
                npl = self.text.__len__()
                for i in range(npl) :
                    ax.text(self.xytext[i][0],
                            self.xytext[i][1],
                            self.text[i])
            else :
                raise ValueError("position of text has to be given")

        # .. draw the plot
        if self.filetype is not None :
            fig.savefig(self.filename+'.'+self.filetype.lower(), dpi = 200,
                        bbox_inches = 'tight')

        else : plt.show()

        # .. clear the figure
        fig.clf()

