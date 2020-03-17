import numpy




#----------------------------------------------------------
#----------------------------------------------------------
def vecprod(V1,V2):
    """
    Returns the vector product V1 x V2
    Exemple  : V3 = vecprod(V1,V2)
    Creation : 2015-02-28 11:28:50.186083

    """
    V = np.zeros(shape=V1.shape, dtype=V1.dtype)
    V[0,:,:,:] = V1[1,:,:,:]*V2[2,:,:,:] - V2[1,:,:,:]*V1[2,:,:,:]
    V[1,:,:,:] = V1[2,:,:,:]*V2[0,:,:,:] - V2[0,:,:,:]*V1[2,:,:,:]
    V[2,:,:,:] = V1[0,:,:,:]*V2[1,:,:,:] - V2[1,:,:,:]*V1[0,:,:,:]
    return V
#==========================================================






#----------------------------------------------------------
#----------------------------------------------------------
def dotprod(V1,V2):
    """
    Returns the dot product V1 dot V2
    Exemple  : scalar = dotprod(V1,V2)
    Creation : 2015-02-28 11:28:50.186083

    """
    s = V1[0,:,:,:] * V2[0,:,:,:]\
            + V1[1,:,:,:] * V2[1,:,:,:] \
            + V1[2,:,:,:] * V2[2,:,:,:]
    return s
#==========================================================











#function found here : http://www.scipy.org/Cookbook/SignalSmooth
# note : corrected the size of the output, +1 added otherwise the size of the output
# does not match the size of the input (+1)
def smooth(x,window_len=11,window='hanning'):


    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
    x: the input signal
    window_len: the dimension of the smoothing window; should be an odd integer
    window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
    flat window will produce a moving average smoothing.

    output:
    the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError( "smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError( "Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError( "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')
    #return y

    return y[(int(window_len/2)-1):-(int(window_len/2)+1)]







