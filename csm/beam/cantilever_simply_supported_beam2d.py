from beam_model import CantileveredSimplySupportedBeam

from fealpy.functionspace import LagrangeFESpace

from fealpy.fem import BilinearForm
from fealpy.fem import EulerBernoulliCantileverBeamStructureIntegrator
from fealpy.fem import DirichletBC

from scipy.sparse.linalg import spsolve

import matplotlib.pyplot as plt
import numpy as np
