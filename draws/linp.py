#
#! /usr/bin/env python
#
#
#import matplotlib as mpl
#mpl.use('Agg')
#
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import rc
#
#
class Linp :

    """

    """

    def __init__(self,
                 axis,
                 data,
                 shade = None,
                 diff = None,
                 drawstyle = None,
                 linlog = None,
                 linestyle = None,
                 linecolor = None,
                 linewidth = None,
                 ticks = None,
                 subticks = None,
                 strticks = None,
                 extent = None,
                 marker = None,
                 legend = None,
                 text = None,
                 xytext = None,
                 labels = None,
                 figsize = [5, 4],
                 filetype = None,
                 filename = None):

        # .. set the x data
        self.axis = axis

        # .. set the y data
        self.data = data

        # .. set the position of shaded rectangle(s)
        self.shade = shade

        # .. set the curved needed for diff
        self.diff = diff

        # .. set the position of ticks
        self.ticks = ticks

        # .. set the # of subticks
        self.subticks = subticks

        # .. set strticks
        self.strticks = strticks

        # .. set extent
        if extent == None:
            self.extent = [[self.axis[0].min(), self.axis[0].max()],\
                           [self.data[0].min(), self.data[0].max()]]
        elif extent.__len__() == 2:
            if extent[0] == None:
                extentx = [self.axis[0].min(), self.axis[0].max()]
            else:
                extentx = extent[0]

            if extent[1] == None:
                extenty = [self.data[0].min(), self.data[0].max()]
            else:
                extenty = extent[1]
            self.extent = [extentx, extenty]

        else:
            raise ValueError("extent should be a list of 2 lists of 2 elements")

        # .. set the linestyle
        if linestyle != None:
            self.linestyle = linestyle
        else:
            self.linestyle = self.axis.__len__()*["-"]

        # .. set the linecolor
        if linecolor != None:
            self.linecolor = linecolor
        else:
            self.linecolor = self.axis.__len__()*["k"]

        # .. set the linewidth
        if linewidth != None:
            self.linewidth = linewidth
        else:
            self.linewidth = npl = self.axis.__len__()*[1.0]

        # .. set the markers
        if marker != None:
            self.marker = marker
        else:
            self.marker = self.axis.__len__()*[None]

        # .. set the linlog
        if linlog != None:
            self.linlog = linlog
        else:
            self.linlog = "linlin"

        # .. set the drawstyle
        if drawstyle != None:
            self.drawstyle = drawstyle
        else:
            self.drawstyle = self.axis.__len__()*['default']

        # .. set the legend
        self.legend = legend

        # .. set the text
        self.text = text
        self.xytext = xytext

        # .. set the labels
        self.labels = labels

        # .. set the figsize
        self.figsize = figsize

        # .. which filetype
        if filetype is not None:
            self.filetype = filetype.lower()
        else:
            self.filetype = None

        # ..
        if filename != None:
            self.filename = filename
        else:
            self.filename = 'Linp'

        # .. draw the plot...
        self.Draw()


    def Draw(self):

        """
        """

        # .. use latex fonts
        rc('text', usetex = True)

        # .. just 1 plot
        fig, ax = plt.subplots(num = 0, figsize = self.figsize, dpi = 100)

        # .. plot the data
        npl = self.axis.__len__()
        for i in range(npl) :
            if self.linestyle[i] != 'bars' :
                ax.plot(self.axis[i],
                        self.data[i],
                        drawstyle = self.drawstyle[i],
                        linestyle = self.linestyle[i],
                        color = self.linecolor[i],
                        marker = self.marker[i],
                        linewidth = self.linewidth[i])
            elif self.linestyle[i] == 'bars' :
                ax.errorbar(self.axis[i],
                            0.5*self.data[i],
                            fmt='none',
                            xerr = 0,
                            yerr = 0.5*self.data[i],
                            drawstyle = self.drawstyle[i],
                            linestyle = self.linestyle[i],
                            color = self.linecolor[i],
                            marker = self.marker[i],
                            linewidth = self.linewidth[i])

        if (self.linlog == 'loglin' or self.linlog == 'loglog') :
            ax.set_xscale('log')
        if (self.linlog == 'linlog' or self.linlog == 'loglog') :
            ax.set_yscale('log')

        # fill between 2 curves
        if self.diff != None :
            ax.fill_between(self.axis[0],
                            self.data[self.diff[0]],
                            self.data[self.diff[1]],
                            where = self.data[self.diff[1]] >= self.data[self.diff[0]],
                            facecolor = '0.4',
                            interpolate = True)
            ax.fill_between(self.axis[0], self.data[self.diff[0]],
                            self.data[self.diff[1]],
                            where = self.data[self.diff[1]] <= self.data[self.diff[0]],
                            facecolor = '0.8',
                            interpolate = True)

        # .. draw shaded rectangle(s)
        if self.shade != None :
            for i in range(self.shade.__len__()):
                plt.axvspan(self.shade[i][0],
                            self.shade[i][1],
                            facecolor = "0.8",
                            alpha = 1.0,
                            linewidth = 0.0)

        # .. set the x ticks & subticks
        if self.ticks != None:
            ax.xaxis.set_major_locator(MultipleLocator(self.ticks[0]))
        if self.strticks != None:
            ax.xaxis.set_major_formatter(FormatStrFormatter(self.strticks[0]))
        if self.subticks != None:
            ax.xaxis.set_minor_locator(MultipleLocator(self.subticks[0]))

        # .. set the x ticks & subticks
        if self.ticks != None:
            ax.yaxis.set_major_locator(MultipleLocator(self.ticks[1]))
        if self.strticks != None:
            ax.yaxis.set_major_formatter(FormatStrFormatter(self.strticks[1]))
        if self.subticks != None:
            ax.yaxis.set_minor_locator(MultipleLocator(self.subticks[1]))

        # .. set the limits of the plot
        ax.set_xlim([self.extent[0][0], self.extent[0][1]])
        ax.set_ylim([self.extent[1][0], self.extent[1][1]])

        # .. set a legend if wanted
        if self.legend != None : ax.legend(self.legend)

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

        # .. write the labels
        if self.labels != None:
            ax.set_xlabel(self.labels[0])
            ax.set_ylabel(self.labels[1])

        # .. draw the plot
        if self.filetype is not None:
            fig.savefig(self.filename + '.'+self.filetype.lower(), dpi=200,
            bbox_inches = 'tight')
        else : plt.show()

        # .. clear the figure
        fig.clf()

