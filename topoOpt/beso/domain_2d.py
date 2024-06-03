import numpy as np
from fealpy.geometry.domain import Domain
from typing import List, Union, Callable, Tuple, Any
from sign_distance_function import drectangle, dcircle, ddiff, test

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

        Parameters:
        - p (array-like): The point(s) at which to compute the signed distance.

        Returns:
        - float or ndarray: The signed distance(s) from the point(s) to the domain.

        """
        return drectangle(p, self.domain)
    

    def signed_dist_function(self, p):
        """
        Compute the signed distance function.

        Parameters:
        - p (array-like): The point(s) at which to compute the signed distance.

        Returns:
        - float or ndarray: The signed distance(s) from the point(s) to the domain.

        """
        return self(p)

    def sizing_function(self, p):
        """
        Compute the sizing function.

        Parameters:
        - p (array-like): The point(s) at which to compute the sizing function.

        Returns:
        - float or ndarray: The sizing function value(s) at the point(s).

        """
        return self.fh(p)

    def facet(self, dim):
        """
        Get the facets of the specified dimension.

        Parameters:
        - dim (int): The dimension of the facets to retrieve.

        Returns:
        - ndarray: The facets of the specified dimension.

        """
        return self.facets[dim]

    def meshing_facet_0d(self):
        """
        Get the 0-dimensional facets for meshing.

        Returns:
        - ndarray: The 0-dimensional facets for meshing.

        """
        return self.facets[0]

    def meshing_facet_1d(self, hmin, fh=None):
        """
        Get the 1-dimensional facets for meshing.

        Parameters:
        - hmin (float): The minimum element size.
        - fh (function): The sizing function.

        """
        pass

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
            self.circles.append(lambda p: dcircle(p, val[0:2], val[2]))

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

    def __call__(self, p: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Compute the signed distance function.

        Parameters:
        - p (float or ndarray): The point(s) at which to compute the signed distance.

        Returns:
        - float or ndarray: The signed distance(s) from the point(s) to the domain.

        """
        return self.facets[1](p)


    def signed_dist_function(self, p: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Compute the signed distance function.

        Parameters:
        - p (float or ndarray): The point(s) at which to compute the signed distance.

        Returns:
        - float or ndarray: The signed distance(s) from the point(s) to the domain.

        """
        return self(p)

    def sizing_function(self, p: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Compute the sizing function.

        Parameters:
        - p (float or ndarray): The point(s) at which to compute the sizing function.

        Returns:
        - float or ndarray: The sizing function value(s) at the point(s).

        """
        return self.fh(p)

    def facet(self, dim: int) -> np.ndarray:
        """
        Get the facets of the specified dimension.

        Parameters:
        - dim (int): The dimension of the facets to retrieve.

        Returns:
        - ndarray: The facets of the specified dimension.

        """
        return self.facets[dim]

    def meshing_facet_0d(self) -> np.ndarray:
        """
        Get the 0-dimensional facets for meshing.

        Returns:
        - ndarray: The 0-dimensional facets for meshing.

        """
        return self.facets[0]

    def meshing_facet_1d(self, hmin: float, fh: Callable[[Any], float] = None) -> None:
        """
        Get the 1-dimensional facets for meshing.

        Parameters:
        - hmin (float): The minimum element size.
        - fh (function): The sizing function.

        """
        pass

class TestDomain(Domain):
    def __init__(self, domain: List[float] = [0, 1, 0, 1], \
                    hmin: float = 0.025, hmax: float = 0.2, \
                    fh: Callable[[Any], float] = None) -> None:
        """
        Initialize a TestDomain object.

        Parameters:
        - box (list): The bounding box of the domain [xmin, xmax, ymin, ymax].
        - circles (list): The circles in the domain [(x, y, r), ...].
        - hmin (float): The minimum element size.
        - hmax (float): The maximum element size.

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


    def __call__(self, p: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Compute the signed distance function.

        Parameters:
        - p (float or ndarray): The point(s) at which to compute the signed distance.

        Returns:
        - float or ndarray: The signed distance(s) from the point(s) to the domain.

        """
        return test(p, self.domain)


    def signed_dist_function(self, p):
        """
        Compute the signed distance function.

        Parameters:
        - p (float or ndarray): The point(s) at which to compute the signed distance.

        Returns:
        - float or ndarray: The signed distance(s) from the point(s) to the domain.

        """
        return self(p)

    def sizing_function(self, p: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Compute the sizing function.

        Parameters:
        - p (float or ndarray): The point(s) at which to compute the sizing function.

        Returns:
        - float or ndarray: The sizing function value(s) at the point(s).

        """
        return self.fh(p)
    
    def facet(self, dim):
        """
        Get the facets of the specified dimension.

        Parameters:
        - dim (int): The dimension of the facets to retrieve.

        Returns:
        - ndarray: The facets of the specified dimension.

        """
        return self.facets[dim]