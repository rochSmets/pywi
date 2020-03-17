

# this module allows one to track the position of the X point
# using the magnetic flux function (out-of-plane component of the
# vector potential). This module works for 2D reconnection runs.




def GetXPoint(flux):
    """
    Calculate the position of the dominant X point from the magnetic flux fonction

    The function returns the position of the X point by looking for the
    saddle point of the magnetic flux function. This method is valid for
    2D reconnection. The definition of the saddle point itself depend on the
    coordinate system used.

    Genericity :
        This function should work with all codes (PIC,Hybrid etc.) as long
        as the magnetic flux is defined on a uniform cartesian mesh


    Args :
        flux - 2D array containing the magnetic flux function

    Returns :
        A tuple of two integers containing the coordinates of the
        X point in terms of position in the original array.

        (i,j) tuple

    TODO :
        The position (i,j) returned is the discrete saddle point
        and not the actual, continuous saddle point, which thus
        might be better estimated by a bilinear interpolation.

    """

    if (flux[1,4] > flux[1,3]):

        flux_max_i   = flux.max(axis=1)      # find the max values along y axis for all x
        flux_max_pos = flux.argmax(axis=1)   # and locate their position along y for all x

        flux_min_pos = flux_max_i.argmin()   # locate the position of the min

        return (flux_min_pos, flux_max_pos[flux_min_pos])

    else:
        flux_min_i   = flux.min(axis=1)      # find the min values along y axis for all x
        flux_min_pos = flux.argmin(axis=1)   # and locate their position along y for all x

        flux_max_pos = flux_min_i.argmax()   # locate the position of the max

        return (flux_max_pos, flux_min_pos[flux_max_pos])



