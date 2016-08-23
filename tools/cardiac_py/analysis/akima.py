import numpy as np


def create(x, y, xi=None, axis=-1, out=None):
    """
    Parameters
    ----------
    x : array like
        1D array of monotonically increasing real values.
    y : array like
        N-D array of real values. y's length along the interpolation
        axis must be equal to the length of x.
    x_new : array like
        New independent variables.
    axis : int
        Specifies axis of y along which to interpolate. Interpolation
        defaults to last axis of y.
    out : array
        Optional array to receive results. Dimension at axis must equal
        length of x.
    """
    x = np.array(x, dtype=np.float64, copy=True)
    y = np.array(y, dtype=np.float64, copy=True)
    if axis != -1 or y.ndim != 1:
        raise NotImplementedError("implemented in C extension module")
    if x.ndim != 1:
        raise ValueError("x-arrays must be one dimensional")
    n = len(x)
    if n < 3:
        raise ValueError("array too small")
    if n != y.shape[axis]:
        raise ValueError("size of x-array must match data shape")
    dx = np.diff(x)
    if any(dx <= 0.0):
        raise ValueError("x-axis not valid")
    m = np.diff(y) / dx
    mm = 2.0 * m[0] - m[1]
    mmm = 2.0 * mm - m[0]
    mp = 2.0 * m[n - 2] - m[n - 3]
    mpp = 2.0 * mp - m[n - 2]
    m1 = np.concatenate(([mmm], [mm], m, [mp], [mpp]))
    dm = np.abs(np.diff(m1))
    f1 = dm[2:n + 2]
    f2 = dm[0:n]
    f12 = f1 + f2
    ids = np.nonzero(f12 > 1e-9 * np.max(f12))[0]
    b = m1[1:n + 1]
    b[ids] = (f1[ids] * m1[ids + 1] + f2[ids] * m1[ids + 2]) / f12[ids]
    c = (3.0 * m - 2.0 * b[0:n - 1] - b[1:n]) / dx
    d = (b[0:n - 1] + b[1:n] - 2.0 * m) / dx ** 2
    if xi == None:
        coefs = np.array([d, c, b[:-1], (y[:-1])])
        coefs = coefs.T
        return coefs
    else:
        bins = np.digitize(xi, x)
        bins = np.minimum(bins, n - 1) - 1
        bb = bins[0:len(xi)]
        wj = xi - x[bb]
        return (((wj * d[bb] + c[bb]) * wj + b[bb]) * wj + y[bb])


def interpolate(breaks, coefs, xi):
    """ Interpolation function for akima spline. Requires the generation of
    the coefficient matrix by create"""
    try:
        n = len(breaks)
        bins = np.digitize(xi, breaks)
        bins = np.minimum(bins, n - 1) - 1
        bb = bins[0:len(xi)]
        wj = xi - breaks[bb]
        return (((wj * coefs[bb, 0] + coefs[bb, 1]) * wj + coefs[bb, 2]) * wj +
        coefs[bb, 3])
    except ValueError:
        print 'xi: ', xi
        print 'xi min: ', xi.min()
        print 'xi max: ', xi.max()
        print 'xi shape: ', xi.shape
        print 'breaks min: ', breaks.min()
        print 'breaks max: ', breaks.max()
        print 'breaks shape: ', breaks.shape
        raise


def a_extreme(breaks, coefs):
    """ Wrapper function for extremes function, it returns the single
    extreme for each quantity """
    [max_val, max_val_loc, min_val, min_val_loc, grad_max, grad_max_loc,
    grad_min, grad_min_loc] = extremes(breaks, coefs)
    if len(max_val) == 0:
        max_v = np.float64(0.0)
        max_v_l = np.float64(0.0)
    else:
        max_v = np.amax(max_val)
        max_v_l = max_val_loc[np.argmax(max_val)]
    if len(min_val) == 0:
        min_v = np.float64(0.0)
        min_v_l = np.float64(0.0)
    else:
        min_v = np.amin(min_val)
        min_v_l = min_val_loc[np.argmin(min_val)]
    if len(grad_max) == 0:
        grad_mx = np.float64(0.0)
        grad_mx_l = np.float64(0.0)
    else:
        grad_mx = np.amax(grad_max)
        grad_mx_l = grad_max_loc[np.argmax(grad_max)]
    if len(grad_min) == 0:
        grad_mn = np.float64(0.0)
        grad_mn_l = np.float64(0.0)
    else:
        grad_mn = np.amin(grad_min)
        grad_mn_l = grad_min_loc[np.argmin(grad_min)]
    return max_v, max_v_l, min_v, min_v_l, grad_mx, grad_mx_l, grad_mn,\
        grad_mn_l


def return_derivative(breaks, coefs):
    """ Returns the derivative at each of the original data points
    for the spline """
    [n_points] = breaks.shape
    derivative = np.zeros(n_points)
    derivative[:-1] = coefs[:, 2]
    dx = breaks[-1] - breaks[-2]
    derivative[-1] = (coefs[-1, 2] + 2 * coefs[-1, 1] * dx + 3 * coefs[-1, 0] *
    (dx ** 2))
    return derivative


def extremes(breaks, coefs):
    """ Returns the local extreme value for the spline both in terms of
    amplitude and gradient """
    [n_cells] = breaks.shape
    #print 'coefs shape: ', coefs.shape
    n_cells -= 1
    min_val = np.array([0.0])
    min_val_loc = np.array([0.0])
    max_val = np.array([0.0])
    #print 'maxval shape: ', max_val.shape
    max_val_loc = np.array([0.0])
    grad_max = np.array([0.0])
    grad_max_loc = np.array([0.0])
    grad_min = np.array([0.0])
    grad_min_loc = np.array([0.0])
    min_x = breaks[:-1]
    max_x = breaks[1:]

    # Simplified derivative
    a = 3 * coefs[:, 0]
    b = 2 * coefs[:, 1] - 6 * min_x * coefs[:, 0]
    c = coefs[:, 2] + 3 * coefs[:, 0] * (min_x ** 2) - 2 * coefs[:, 1] * min_x
    for ii in range(n_cells):
        min_x = breaks[ii]
        max_x = breaks[(ii + 1)]
        a = 3 * coefs[ii, 0]
        b = 2 * coefs[ii, 1] - 6 * min_x * coefs[ii, 0]
        c = coefs[ii, 2] + 3 * coefs[ii, 0] * (min_x ** 2) - 2 * coefs[ii, 1] \
         * min_x
        root_1 = (-b + np.sqrt((b ** 2) - 4 * a * c)) / (2 * a)
        type(root_1)
        root_2 = (-b - np.sqrt((b ** 2) - 4 * a * c)) / (2 * a)
        if np.isfinite(root_1) and (min_x <= root_1) and (root_1 <= max_x):
            d2 = 2 * a * root_1 + b
            if d2 > 0:
                min_val_loc = np.concatenate((min_val_loc, np.array([root_1])))
                min_val = np.concatenate((min_val, np.array(
                interpolate(breaks, coefs, np.array([root_1])))))
            elif d2 < 0:
                max_val_loc = np.concatenate((max_val_loc, np.array([root_1])))
                max_val = np.concatenate((max_val,
                interpolate(breaks, coefs, np.array([root_1]))))
            else:
                pass
        if np.isfinite(root_2) and (min_x <= root_2) and (root_2 <= max_x):
            d2 = 2 * a * root_2 + b
            if d2 > 0:
                min_val_loc = np.concatenate((min_val_loc, np.array([root_2])))
                min_val_t = interpolate(breaks, coefs, np.array([root_2]))
                min_val = np.concatenate((min_val, min_val_t))
            elif d2 < 0:
                max_val_loc = np.concatenate((max_val_loc, np.array([root_2])))
                max_val = np.concatenate((max_val,
                interpolate(breaks, coefs, np.array([root_2]))))
            else:
                pass
        # max min in grad
        gmml = -b / (2 * a)
        grad = a * (gmml ** 2) + b * gmml + c
        if (min_x < gmml) and (gmml < max_x):
            # have grad max min
            if (a < 0) and (grad > 0):
                grad_max = np.concatenate((grad_max, np.array([grad])))
                grad_max_loc = np.concatenate((grad_max_loc, np.array([gmml])))
            elif (a > 0) and (grad < 0):
                grad_min = np.concatenate((grad_min, np.array([grad])))
                grad_min_loc = np.concatenate((grad_min_loc, np.array([gmml])))
            else:
                pass
    #print max_val
    max_v = max_val[1:]
    max_v_l = max_val_loc[1:]
    min_v = min_val[1:]
    min_v_l = min_val_loc[1:]
    grad_mx = grad_max[1:]
    grad_mx_l = grad_max_loc[1:]
    grad_mn = grad_min[1:]
    grad_mn_l = grad_min_loc[1:]
    # at this time we are happy that this works make wrapper function
    #a_extreme_akima
    return max_v, max_v_l, min_v, min_v_l, grad_mx, grad_mx_l, grad_mn,\
        grad_mn_l
