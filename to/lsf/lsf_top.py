import numpy as np

from fealpy.mesh import QuadrangleMesh

class TopLsf:

    def __init__(self, nelx: int = 60, nely: int = 30, volReq: float = 0.3, 
                stepLength: int = 3, numReinit: int = 2, topWeight: int = 2):
        '''
        初始化拓扑优化问题

        Parameters: 
        - nelx (int): Number of elements in the horizontal direction. Defaults to 60.
        - nely (int): Number of elements in the vertical direction. Defaults to 30.
        - volReq : 最终设计所需的体积分数. Defaults to 0.3.
        - stepLength (int): The CFL time step size used in each iteration of the evolution equation. Defaults to 3.
        - numReinit (int): The number of algorithm iterations before performing level set reinitialization. Defaults to 2.
        - topWeight (int): Weight of the topology derivative term in the evolution equation. Defaults to 2.
        '''

        self._nelx = nelx
        self._nely = nely
        self._volReq = volReq
        self._stepLength = stepLength
        self._numReinit = numReinit
        self._topWeight = topWeight
        node = np.array([[0, 2], [0, 1], [0, 0],
                         [1, 2], [1, 1], [1, 0],
                         [2, 2], [2, 1], [2, 0]], dtype=np.float64)
        cell = np.array([[0, 3, 4, 1],
                         [1, 4, 5, 2],
                         [3, 6, 7, 4],
                         [4, 7, 8, 5]], dtype=np.int_)
        self._mesh = QuadrangleMesh(node=node, cell=cell)

