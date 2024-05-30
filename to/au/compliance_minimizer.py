from jax import jit, value_and_grad
from mmaOptimize import optimize

class ComplianceMinimizer:
    def __init__(self, mesh, bc, material, globalvolCons, projection):
        self.mesh = mesh
        self.bc = bc
        self.material = material
        self.globalVolumeConstraint = globalvolCons
        self.projection = projection
        self.objectiveHandle = jit(value_and_grad(self.computeCompliance))

    # Code snippet 2.1
    def computeCompliance(self, rho):
        pass
        #-----------------------#
        @jit
        # Code snippet 2.2
        def materialModel(rho):
            E = self.material['Emin'] + \
                (self.material['Emax'] - self.material['Emin'])*\
                (rho + 0.01) ** self.material['penal']
            return E
        
    #-----------------------#
    def TO(self, optimizationParams, ft):
        optimize(self.mesh, optimizationParams, ft, \
                 self.objectiveHandle, self.consHandle, self.numConstraints)