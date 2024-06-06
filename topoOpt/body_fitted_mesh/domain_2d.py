import numpy as np
from fealpy.geometry.domain import Domain
from typing import List, Union, Callable, Tuple, Any
from sign_distance_function import drectangle, dcircle, ddiff, dlevelset

class RectangleDomain(Domain):
    def __init__(self, domain: List[float] = [0, 1, 0, 1], \
                       hmin: float = 0.1, hmax: float = None, \
                       fh: Callable[[Any], float] = None) -> None:
        """
        Initialize a RectangleDomain object.

        Parameters:
        - domain (list): The domain boundaries [xmin, xmax, ymin, ymax].
        - hmin (float): The minimum element size.
        - hmax (float): The maximum element size.
        - fh (function): The sizing function.

        """
        super().__init__(hmin=hmin, hmax=hmax, GD=2)
        if fh is not None:
            self.fh = fh

        self.domain = domain

        # Define the bounding box of the domain.
        mx = (domain[1] - domain[0])/10
        my = (domain[3] - domain[2])/10
        self.box = [domain[0]-mx, domain[1]+mx, domain[2]-my, domain[3]+my]

        # Define the vertices and curves of the domain.
        vertices = np.array([
            (domain[0], domain[2]),
            (domain[1], domain[2]),
            (domain[1], domain[3]),
            (domain[0], domain[3]),
            ], dtype=np.float64)
        
        curves = np.array([
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            ], dtype=np.float64)

        # Define the facets of the domain.    
        self.facets = {0:vertices, 1:curves}

    def __call__(self, p):
        """
        Compute the signed distance function.
        """
        return drectangle(p, self.domain)
    
    def signed_dist_function(self, p):
        """
        Compute the signed distance function.
        """
        return self(p)

    def sizing_function(self, p):
        """
        Compute the sizing function.
        """
        return self.fh(p)

    def facet(self, dim):
        """
        Get the facets of the specified dimension.
        """
        return self.facets[dim]

    def meshing_facet_0d(self):
        """
        Get the 0-dimensional facets for meshing.
        """
        return self.facets[0]


class BoxWithCircleHolesDomain(Domain):
    def __init__(self, box: List[float] = [0, 1, 0, 1], \
                    circles: List[Tuple[float, float, float]] = [(0.5, 0.5, 0.2)], \
                    hmin: float = 0.003, hmax: float = 0.01) -> None:
        """
        Initialize a BoxWithCircleHolesDomain object.

        Parameters:
        - box (list): The bounding box of the domain [xmin, xmax, ymin, ymax].
        - circles (list): The circles in the domain [(x, y, r), ...].
        - hmin (float): The minimum element size.
        - hmax (float): The maximum element size.

        """
        super().__init__(hmin=hmin, hmax=hmax, GD=2)

        self.box = box 

        # Define the vertices of the domain.
        vertices = np.array([
            (box[0], box[2]), 
            (box[1], box[2]), 
            (box[1], box[3]),
            (box[0], box[3])],dtype=np.float64)

        self.circles = []
        for val in circles:
            self.circles = [lambda p, val=val: dcircle(p, val[0:2], val[2]) for val in circles]

        def fd(p):
            d0 = drectangle(p, box)
            for circle in self.circles:
                d0 = ddiff(d0, circle(p))
            return d0

        def fh(p):
            d0 = 1e10 
            for circle in self.circles:
                d0 = np.minimum(d0, circle(p))
            h = hmin + 0.05*d0
            h[h>hmax] = hmax 
            return h

        self.fh = fh
        self.facets = {0:vertices, 1:fd}

    def __call__(self, p):
        """
        Compute the signed distance function.
        """
        return self.facets[1](p)

    def signed_dist_function(self, p):
        """
        Compute the signed distance function.
        """
        return self(p)

    def sizing_function(self, p):
        """
        Compute the sizing function.
        """
        return self.fh(p)

    def facet(self, dim):
        """
        Get the facets of the specified dimension.
        """
        return self.facets[dim]

    def meshing_facet_0d(self):
        """
        Get the 0-dimensional facets for meshing.
        """
        return self.facets[0]


class BoxWithZeroLevelSetDomain(Domain):
    def __init__(self, box, zero_level_set, hmin=0.003, hmax=0.01):
        """
        Initialize a BoxWithZeroLevelSetDomain object.

        Parameters:
        - box (list): The domain boundaries [xmin, xmax, ymin, ymax].
        - zero_level_set (array): An array of points representing the zero level set.
        - hmin (float): The minimum element size.
        - hmax (float): The maximum element size.
        """
        super().__init__(hmin=hmin, hmax=hmax, GD=2)

        self.box = box
        vertices = np.array([
            (box[0], box[2]),
            (box[1], box[2]),
            (box[1], box[3]),
            (box[0], box[3])
        ], dtype=np.float64)

        self.zero_level_set = zero_level_set

        def fd(p):
            d0 = drectangle(p, box)
            for level_set in zero_level_set:
                d0 = ddiff(d0, dlevelset(p, level_set))
            return d0

        def fh(p):
            d0 = 1e10
            for level_set in zero_level_set:
                d0 = np.minimum(d0, dlevelset(p, level_set))
            h = hmin + 0.05 * d0
            h[h > hmax] = hmax
            return h

        self.fh = fh
        self.facets = {0: vertices, 1: fd}

    def __call__(self, p):
        """
        Compute the signed distance function.
        """
        return self.facets[1](p)

    def signed_dist_function(self, p):
        """
        Compute the signed distance function.
        """
        return self(p)

    def sizing_function(self, p):
        """
        Compute the sizing function.
        """
        return self.fh(p)

    def facet(self, dim):
        """
        Get the facets of the specified dimension.
        """
        return self.facets[dim]

    def meshing_facet_0d(self):
        """
        Get the 0-dimensional facets for meshing.
        """
        return self.facets[0]