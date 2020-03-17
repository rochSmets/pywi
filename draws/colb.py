
N = 12
N2 = int(N/2)
from matplotlib import cm
import numpy as np
from matplotlib.colors import ListedColormap
bot = cm.get_cmap('Greys_r', N2)
top = cm.get_cmap('Purples', N2)

newcolors = np.vstack((bot(np.linspace(0, 1, N2)),
                       top(np.linspace(0, 1, N2))))
colormap = [ListedColormap(newcolors, name='zobi'), N]

