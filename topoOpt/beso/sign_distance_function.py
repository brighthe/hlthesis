import numpy as np

def dcircle(p, cxy=[0, 0], r=1):
    x = p[..., 0]
    y = p[..., 1]
    return np.sqrt((x - cxy[0])**2 + (y - cxy[1])**2) - r

def drectangle(p, box):
    """
    Compute the signed distance from points to the rectangle.

    Parameters:
    - p (array-like): The point(s) at which to compute the signed distance.
    - box (list): The rectangle boundaries [xmin, xmax, ymin, ymax].

    Returns:
    - float or ndarray: The signed distance(s) from the point(s) to the rectangle.
    """
    x = p[..., 0]
    y = p[..., 1]
    d = dmin(y - box[2], box[3] - y)
    d = dmin(d, x - box[0])
    d = dmin(d, box[1] - x)
    return -d

def test(p, box):
    x = p[..., 0]
    y = p[..., 1]
    pi = np.pi
    d = np.sin(6*pi*2*x/(box[1]-box[0]))*np.cos(6*pi*2*y/(box[1]-box[0])) + 0.5
    return -d

def ddiff(d0, d1):
    return np.maximum(d0, -d1)

def dmin(*args):
    d = np.array(args)
    return np.min(d, axis=0)

def dmax(*args):
    d = np.array(args)
    return np.max(d, axis=0)