import jax
import time
import numpy as np
import jax.numpy as jnp
import matplotlib.pyplot as plt

from jax import jit, value_and_grad

from mesher import Mesher
from mma import  MMA
from util_functions import applySensitivityFilter


class ComplianceMinimizer:
    def __init__(self, mesh, bc, material, globalvolCons, projection):
        self.mesh = mesh
        self.bc = bc
        self.material = material
        self.globalVolumeConstraint = globalvolCons
        self.projection = projection

        M = Mesher()
        self.edofMat, self.idx = M.getMeshStructure(mesh)
        self.K0 = M.getK0(self.material)
        
        self.objectiveHandle = jit(value_and_grad(self.computeCompliance))

        self.consHandle = self.computeConstraints

        self.numConstraints = 1
        

    # Code snippet 2.1
    def computeCompliance(self, rho):

        @jit
        # 投影滤波器
        def projectionFilter(rho):
            if(self.projection['isOn']):
                v1 = np.tanh(self.projection['c0']*self.projection['beta'])
                nm = v1 + jnp.tanh(self.projection['beta']*(rho-self.projection['c0']))
                dnm = v1 + jnp.tanh(self.projection['beta']*(1.-self.projection['c0']))
                return nm/dnm
            else:
                return rho

        @jit
        # SIMP 材料插值模型
        def materialModel(rho):
            E = self.material['Emin'] + \
                (self.material['Emax'] - self.material['Emin']) * \
                (rho + 0.01) ** self.material['penal']
            return E
        
        @jit
        # 组装全局刚度矩阵
        def assembleK(E):
            K_asm = jnp.zeros((self.mesh['ndof'], self.mesh['ndof']))
            K_elem = (self.K0.flatten()[np.newaxis]).T 

            K_elem = (K_elem*E).T.flatten()
            K_asm = K_asm.at[(self.idx)].add(K_elem) # UPDATED
            return K_asm
        
        @jit
        # 直接法求解线性方程组
        def solveKuf(K): 
            u_free = jax.scipy.linalg.solve(K[self.bc['free'],:][:,self.bc['free']], \
                                        self.bc['force'][self.bc['free']], check_finite=False)
            u = jnp.zeros((self.mesh['ndof']))
            u = u.at[self.bc['free']].set(u_free.reshape(-1)) # UPDATED
            return u
        
        rho = projectionFilter(rho)
        E = materialModel(rho)
        K = assembleK(E)
        u = solveKuf(K)
        J = jnp.dot(self.bc['force'].T, u)[0]

        return J
        
    def computeConstraints(self, rho, epoch): 

        @jit
        # 计算体积约束
        def computeGlobalVolumeConstraint(rho):
            g = jnp.mean(rho)/self.globalVolumeConstraint['vf'] - 1.
            return g
        
        # 体积约束的值及其灵敏度
        c, gradc = value_and_grad(computeGlobalVolumeConstraint)(rho)
        c, gradc = c.reshape((1,1)), gradc.reshape((1,-1))
        return c, gradc

    def mmaOptimize(self, optimizationParams, ft):
        
        rho = np.ones((self.mesh['nelx']*self.mesh['nely']))
        loop = 0
        change = 1.
        m = self.numConstraints
        n = self.mesh['numElems']

        mma = MMA()
        mma.setNumConstraints(self.numConstraints)
        mma.setNumDesignVariables(n);
        mma.setMinandMaxBoundsForDesignVariables(np.zeros((n,1)), np.ones((n,1)))

        xval = rho[np.newaxis].T
        xold1, xold2 = xval.copy(), xval.copy()

        mma.registerMMAIter(xval, xold1, xold2)
        mma.setLowerAndUpperAsymptotes(np.ones((n,1)), np.ones((n,1)))
        mma.setScalingParams(1.0, np.zeros((m,1)), \
                            10000*np.ones((m,1)), np.zeros((m,1)))
        mma.setMoveLimit(0.2)

        mmaTime = 0

        t0 = time.perf_counter()

        while( (change > optimizationParams['relTol']) \
            and (loop < optimizationParams['maxIters'])\
            or (loop < optimizationParams['minIters'])):
            loop = loop + 1

            J, dJ = self.objectiveHandle(rho)

            vc, dvc = self.consHandle(rho, loop)

            dJ, dvc = applySensitivityFilter(ft, rho, dJ, dvc)

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
                plt.imshow(-np.flipud(rho.reshape((self.mesh['nelx'], self.mesh['nely'])).T), cmap='gray');
                plt.title(status)
                plt.show()
            
        totTime = time.perf_counter() - t0;

        print('total time(s): ', totTime);
        print('mma time(s): ', mmaTime);
        print('FE time(s): ', totTime - mmaTime);
        return rho;
