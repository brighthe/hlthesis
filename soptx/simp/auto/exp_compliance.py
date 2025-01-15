import jax
import jax.numpy as jnp
import numpy as np
import time
import matplotlib.pyplot as plt

from jax import jit, value_and_grad

from utilfuncs import computeFilter, Mesher, MMA, applySensitivityFilter

nelx, nely = 6, 3
elemSize = np.array([1., 1.])
mesh = {'nelx':nelx, 'nely':nely, 'elemSize':elemSize,\
        'ndof':2*(nelx+1)*(nely+1), 'numElems':nelx*nely}

material = {'Emax':1., 'Emin':1e-3, 'nu':0.3, 'penal':3.}

filterRadius = 1.5
H, Hs = computeFilter(mesh, filterRadius)
ft = {'type':1, 'H':H, 'Hs':Hs}

example = 1
if(example == 1):
    # tip cantilever
    force = np.zeros((mesh['ndof'],1))
    dofs = np.arange(mesh['ndof'])
    fixed = dofs[0:2*(nely+1):1]
    free = jnp.setdiff1d(np.arange(mesh['ndof']),fixed)
    force[2*(nelx+1)*(nely+1)-2*nely+1, 0 ] = -1
    symXAxis = False
    symYAxis = False
elif(example == 2):
    ndof = 2*(nelx+1)*(nely+1)
    force = np.zeros((mesh['ndof'],1))
    dofs=np.arange(mesh['ndof'])
    fixed = dofs[0:2*(nely+1):1]
    free = jnp.setdiff1d(np.arange(mesh['ndof']),fixed)
    force[2*(nelx+1)*(nely+1)- (nely+1), 0 ] = -1
    symXAxis = True
    symYAxis = False
bc = {'force':force, 'fixed':fixed,'free':free,\
          'symXAxis':symXAxis, 'symYAxis':symYAxis}

globalVolumeConstraint = {'isOn':True, 'vf':0.5}

optimizationParams = {'maxIters':200,'minIters':100,'relTol':0.05}
projection = {'isOn':False, 'beta':4, 'c0':0.5}

class ComplianceMinimizer:
    def __init__(self, mesh, bc, material, globalvolCons, projection):
        self.mesh = mesh
        self.material = material
        self.bc = bc
        M = Mesher()
        self.edofMat, self.idx = M.getMeshStructure(mesh)
        self.K0 = M.getK0(self.material)
        self.globalVolumeConstraint = globalvolCons
        self.objectiveHandle = jit(value_and_grad(self.computeCompliance))
        
        self.consHandle = self.computeConstraints
        self.numConstraints = 1
        self.projection = projection

    def computeCompliance(self, rho):

        @jit
        def projectionFilter(rho):
            if(self.projection['isOn']):
                v1 = np.tanh(self.projection['c0']*self.projection['beta'])
                nm = v1 + jnp.tanh(self.projection['beta']*(rho-self.projection['c0']))
                dnm = v1 + jnp.tanh(self.projection['beta']*(1.-self.projection['c0']))
                return nm / dnm
            else:
                return rho
        
        @jit
        def materialModel(rho):
            E = self.material['Emin'] + \
                (self.material['Emax']-self.material['Emin'])*(rho) ** self.material['penal']
                
            return E
        
        @jit
        def assembleK(E):
            K_asm = jnp.zeros((self.mesh['ndof'], self.mesh['ndof']))
            K_elem = (self.K0.flatten()[np.newaxis]).T 

            K_elem = (K_elem*E).T.flatten()
            K_asm = K_asm.at[(self.idx)].add(K_elem)

            return K_asm
        
        @jit
        def solveKuf(K): 
            u_free = jax.scipy.linalg.solve(K[self.bc['free'], :][:,self.bc['free']], \
                                                self.bc['force'][self.bc['free']])
            u = jnp.zeros((self.mesh['ndof']))
            u = u.at[self.bc['free']].set(u_free.reshape(-1))

            return u
        
        rho = projectionFilter(rho)
        E = materialModel(rho)
        K = assembleK(E)
        u = solveKuf(K)
        J = jnp.dot(self.bc['force'].T, u)[0]
        
        return J

    def computeConstraints(self, rho, epoch): 
        @jit
        def computeGlobalVolumeConstraint(rho):
            g = jnp.mean(rho) / self.globalVolumeConstraint['vf'] - 1.
            
            return g
        
        c, gradc = value_and_grad(computeGlobalVolumeConstraint)(rho)
        c, gradc = c.reshape((1, 1)), gradc.reshape((1, -1))

        return c, gradc

Opt = ComplianceMinimizer(mesh, bc, material, globalVolumeConstraint, projection)
numConstraints = Opt.numConstraints
objectiveHandle = Opt.objectiveHandle
consHandle = Opt.consHandle

rho = np.ones((mesh['nelx'] * mesh['nely']))
loop = 0
change = 1.
m = numConstraints
n = mesh['numElems']
mma = MMA()
mma.setNumConstraints(numConstraints)
mma.setNumDesignVariables(n)
mma.setMinandMaxBoundsForDesignVariables(np.zeros((n,1)), np.ones((n,1)))

xval = rho[np.newaxis].T
xold1, xold2 = xval.copy(), xval.copy()
mma.registerMMAIter(xval, xold1, xold2)
mma.setLowerAndUpperAsymptotes(np.ones((n,1)), np.ones((n,1)))
mma.setScalingParams(1.0, np.zeros((m,1)), 10000*np.ones((m,1)), np.zeros((m,1)))
mma.setMoveLimit(0.2)

mmaTime = 0

t0 = time.perf_counter()

while( (change > optimizationParams['relTol']) \
        and (loop < optimizationParams['maxIters']) \
        or (loop < optimizationParams['minIters']) ):
    loop = loop + 1

    J, dJ = objectiveHandle(rho)

    vc, dvc = consHandle(rho, loop)

    dJ, dvc = applySensitivityFilter(ft, rho, dJ, dvc)

    J, dJ = J, dJ[np.newaxis].T
    tmr = time.perf_counter()
    mma.setObjectiveWithGradient(J, dJ)
    mma.setConstraintWithGradient(vc, dvc)

    xval = rho.copy()[np.newaxis].T
    mma.mmasub(xval)

    xmma, _, _ = mma.getOptimalValues()
    xold2 = xold1.copy()
    xold1 = xval.copy()
    rho = xmma.copy().flatten()
    mma.registerMMAIter(rho, xval.copy(), xold2.copy())

    mmaTime += time.perf_counter() - tmr

    status = 'Iter {:d}; J {:.2F}; vf {:.2F}'.format(loop, J, jnp.mean(rho))
    print(status)
    if(loop%10 == 0):
        plt.imshow(-np.flipud(rho.reshape((mesh['nelx'], mesh['nely'])).T), cmap='gray')
        plt.title(status)
        plt.show()

totTime = time.perf_counter() - t0

print('total time(s): ', totTime)
print('mma time(s): ', mmaTime)
print('FE time(s): ', totTime - mmaTime)