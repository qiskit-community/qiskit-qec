# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2022
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
"""Define the bivariate bicylce code."""

from typing import List, Optional, Sequence, Tuple

import networkx as nx
import numpy as np
import rustworkx as rx
from rustworkx.visualization import graphviz_draw
from scipy import sparse

from qiskit_qec.linear.symplectic import normalizer
from qiskit_qec.operators.pauli_list import PauliList


def cs_pow(l, power=1):
        """
        Calculates and returns a power of the cyclic shift matrix of size lxl (C_l)

        Parameters
        ----------
        l     : int
            size of cyclic shift matrix.
        power : int
            Power to which Cl is raised. Defaults to 1.

        Returns
        --------
        np.array
            C_l^power

        Examples
        --------

        >>> BivariateBicycleCodeBuilder.cs_pow(3, 2) # (C_3)^2
        """

        return np.roll(np.eye(l, dtype=np.uint8), shift=power, axis=1)

def arr_to_indices(arr: np.array) -> List[List[int]]:
        """ Converts a numpy array to a list of list of indices where it is non-zero """
        return [np.where(row)[0].tolist() for row in arr]

class BBCodeVector:
    """Bivariate Bicycle code data.

    The X and Z gauge operator lists are given as lists of supports.
    There is a consistent qubit ordering of these lists so that
    we can construct gate schedules for circuits.
    """

    def __init__(
        self,
        l: int,
        m: int,
        a: Tuple[int, int],
        b: Tuple[int, int],
    ) -> None:
        """Initializes a bivariate bicycle code builder

        Args:
            l: dimension of first space of generators x and y,
            m: dimension of second space of generators x and y,
            a: First vector for additional check
            b: Second vector for additional check

        Examples:
            Example 1 (Gros code):
            >>> code = BBCodeVector(l=6, m=12, a=(3,2), b=(-1,3))
        """

        self.l = l
        self.m = m
        self.s = l*m
        self.n = 2*l*m

        self.a = a
        self.b = b

        self._create_check_matrices()

        self._x_stabilizers = arr_to_indices(self.hx)
        self._z_stabilizers = arr_to_indices(self.hz)

        self._create_logicals()

        self.create_tanner_graph()


    def _create_check_matrices(self):
        A1 = np.eye(self.s, dtype=np.uint8)
        A2 = cs_pow(self.s, self.l)
        A3 = np.kron(cs_pow(self.m, self.a[1]), cs_pow(self.l, self.a[0]))

        self.A = A1 + A2 + A3

        B1 = np.eye(self.s, dtype=np.uint8)
        B2 = np.kron(np.eye(self.m, dtype=np.uint8), cs_pow(self.l, -1))
        B3 = np.kron(cs_pow(self.m, self.b[1]), cs_pow(self.l, self.b[0]-1))
        
        self.B = B1 + B2 + B3

        self.hx = np.hstack([self.A, self.B])
        self.hz = np.hstack([self.B.T, self.A.T])
        # self.hx = np.zeros((self.s, self.n), dtype=int)
        # for i in range(self.s):
        #     # vertical connections (left data qubits)
        #     down = i
        #     up = (i+self.l) % (self.s)
        #     vec_from_down = self._coord2idx(self._coord_add(self._idx2coord(down), self.a))
        #     self.hx[i, up] = 1
        #     self.hx[i, down] = 1
        #     self.hx[i, vec_from_down] = 1
        #     # horizontal connections (right data qubits)
        #     left = (i-1) % self.l  + (i// code.l)*code.l
        #     right = i
        #     vec_from_left = self._coord2idx(self._coord_add(self._idx2coord(left), self.b))
        #     self.hx[i, left + self.s] = 1
        #     self.hx[i, right + self.s] = 1
        #     self.hx[i, vec_from_left + self.s] = 1

        # self.hz = np.zeros((self.s, self.n), dtype=int)
        # # similar

        # C1 = np.eye(self.s, dtype=int) # = 1
        # C2 = np.kron(np.eye(self.m, dtype=int), cs_pow(self.l, 1))
        # C3 = np.kron(cs_pow(self.m, -self.b[1]), cs_pow(self.l, -(self.b[0]-1)))

        # self.C = C1 + C2 + C3 # hopefully same as B.T

        # D1 = np.eye(self.s, dtype=int) # = 1
        # D2 = cs_pow(self.s, -self.l)
        # D3 = np.kron(cs_pow(self.m, -self.a[1]), cs_pow(self.l, -self.a[0]))

        # self.D = D1 + D2 + D3

        # A1 = np.eye(self.s, dtype=int) # = 1
        # A2 = cs_pow(self.s, self.l)
        # A3 = np.kron(cs_pow(self.m, self.a[1]), cs_pow(self.l, self.a[0]))

        # self.A = A1 + A2 + A3

        # B1 = np.eye(self.s, dtype=int)
        # B2 = np.kron(np.eye(self.m, dtype=int), cs_pow(self.l, -1))
        # B3 = np.kron(cs_pow(self.m, self.b[1]), cs_pow(self.l, self.b[0]-1))
        
        # self.B = B1 + B2 + B3

        # for i in range(self.s):
        #     # horizontal connections (left data qubits)
        #     left = i
        #     right = (i+1) % self.l  + (i// code.l)*code.l
        #     vec_from_right = self._coord2idx(self._coord_sub(self._idx2coord(right), self.b))
        #     self.hz[i, left] = 1
        #     self.hz[i, right] = 1
        #     self.hz[i, vec_from_right] = 1
        #     # vertical connections (right data qubits)
        #     up = i
        #     down = (i-self.l) % (self.s)
        #     vec_from_up = self._coord2idx(self._coord_sub(self._idx2coord(up), self.a))
        #     self.hz[i, up + self.s] = 1
        #     self.hz[i, down + self.s] = 1
        #     self.hz[i, vec_from_up + self.s] = 1

    def _create_logicals(self):
        full_sym = np.vstack([np.hstack([self.hx, np.zeros((self.n//2, self.n))]),
                                np.hstack([np.zeros((self.n//2, self.n)), self.hz])])
        center_, x_new, z_new = normalizer(full_sym.astype(np.bool_))

        self._logical_z = PauliList(z_new)
        self._logical_x = PauliList(x_new)

    def _idx2coord(self, idx):
        # returns the relative coordinate of an idx within group (X, Z, L or R)
        if type(idx) == list:
            idx = np.array(idx)
        return (idx % self.l, idx // self.l)
    
    def _fault_idx2coord(self, idx):
        if type(idx) == list:
            idx = np.array(idx)

        r_type = idx > self.s
        idx %= self.s

        return (idx % self.l + 0.5*r_type, idx // self.l + 0.5*r_type)
    
    def _check_index2coord(self, idx, check_type):
        if type(idx) == list:
            idx = np.array(idx)

        x_shift = check_type=='x'
        z_shift = check_type=='z'

        return (idx % self.l + 0.5*z_shift, idx // self.l + 0.5*x_shift)
    
    def _coord2idx(self, coord):
        return coord[0] + coord[1]*self.l
    
    def _coord_add(self, coord, v):
        return (coord[0] + v[0]) % self.l, (coord[1] + v[1]) % self.m
    
    def _coord_sub(self, coord, v):
        return (coord[0] - v[0]) % self.l, (coord[1] - v[1]) % self.m
    
    @property
    def x_stabilizers(self) -> PauliList:
        return self._x_stabilizers
    
    @property
    def z_stabilizers(self) -> PauliList:
        return self._z_stabilizers
    
    @property
    def x_gauges(self) -> PauliList:
        return self.x_stabilizers
    
    @property
    def z_gauges(self) -> PauliList:
        return self.z_stabilizers
    
    @property
    def logical_z(self) -> List[List[int]]:
        return BBCode.arr_to_indices(self._logical_z.matrix[:, self.n:])
    
    @property
    def logical_x(self) -> List[List[int]]:
        return BBCode.arr_to_indices(self._logical_x.matrix[:, :self.n])
    
    def get_syndrome(self, physical_error: np.ndarray, base: str='z') -> np.ndarray:
        if base == 'z':
            return self.hz @ physical_error % 2
        elif base == 'x':
            return self.hx @ physical_error % 2
        else:
            raise ValueError(f'"base" must be one of {{"x", "y"}}, not {base}.')
        
    def get_logical_error(self, physical_error: np.ndarray, base: str='z') -> np.ndarray:
        if base == 'z':
            return self._logical_z.matrix[:, self.n:] @ physical_error % 2
        elif base == 'x':
            return self._logical_x.matrix[:, :self.n] @ physical_error % 2
        else:
            raise ValueError(f'"base" must be one of {{"x", "y"}}, not {base}.')
    
    @property
    def k(self):
        raise NotImplementedError()

    @property
    def d(self):
        raise NotImplementedError()
    
    class CodeConverter:
        def __init__(self, code_a: "BBCodeVector", code_b: "BBCodeVector", graph_matcher: nx.isomorphism.GraphMatcher) -> None:
            self.code_a = code_a
            self.code_b = code_b
            self.graph_matcher = graph_matcher
            num_checks = code_a.s
            num_faults = code_a.n
            self._check_forward_exchange = np.empty(num_checks, dtype=int)
            self._check_backward_exchange = np.empty(num_checks, dtype=int)
            self._fault_forward_exchange = np.empty(num_faults, dtype=int)
            self._fault_backward_exchange = np.empty(num_faults, dtype=int)
            for key in graph_matcher.mapping:
                val = graph_matcher.mapping[key]
                if key < num_checks:
                    self._check_forward_exchange[key] = val
                    self._check_backward_exchange[val] = key
                else:
                    self._fault_forward_exchange[key-num_checks] = val - num_checks
                    self._fault_backward_exchange[key-num_checks] = val - num_checks

        def faults_forwards(self, indices_a: List[int]) -> List[int]:
            indices_b = self._fault_backward_exchange[indices_a]
            return indices_b

        def faults_backwards(self, indices_b: List[int]) -> List[int]:
            indices_a = self._fault_backward_exchange[indices_b]
            return indices_a
        
        def checks_forwards(self, indices_a: List[int]) -> List[int]:
            indices_b = self._check_forward_exchange[indices_a]
            return indices_b
        
        def checks_backwards(self, indices_b: List[int]) -> List[int]:
            indices_a = self._check_forward_exchange[indices_b]
            return indices_a
        
    def z_fault_graph(self):
        sparse_repr = sparse.csr_matrix(self.hz)
        return nx.algorithms.bipartite.from_biadjacency_matrix(sparse_repr)

    def check_equivalence(self, other: "BBCodeVector"):
        s_a = sparse.csr_matrix(self.hz)
        s_b = sparse.csr_matrix(other.hz)
        Ga = nx.algorithms.bipartite.from_biadjacency_matrix(s_a)
        Gb = nx.algorithms.bipartite.from_biadjacency_matrix(s_b)
        GM = nx.isomorphism.GraphMatcher(Ga, Gb)
        if not GM.is_isomorphic():
            return None
        
        return BBCodeVector.CodeConverter(code_a=self, code_b=other, graph_matcher=GM)

    @property
    def tanner_graph(self):
        return self._tanner_graph
    
    @property
    def tanner_graph_X(self):
        nodes = self.tanner_graph.filter_nodes(lambda node: node['subtype'] != 'Z')
        return self.tanner_graph.subgraph(nodes)
    
    @property
    def tanner_graph_Z(self):
        nodes = self.tanner_graph.filter_nodes(lambda node: node['subtype'] != 'X')
        return self.tanner_graph.subgraph(nodes)

    def __str__(self) -> str:
        """Formatted string."""
        return f"(l={self.l},m={self.m}) bivariate bicycle code"
        #return f"[[{self.n}, {self.k}, {self.d}]] heavy-hexagon compass code"

    def __repr__(self) -> str:
        """String representation."""
        val = str(self)
        val += f"\nx_gauges = {self.x_gauges}"
        val += f"\nz_gauges = {self.z_gauges}"
        val += f"\nx_stabilizers = {self.x_stabilizers}"
        val += f"\nz_stabilizers = {self.z_stabilizers}"
        val += f"\nlogical_x = {self.logical_x}"
        val += f"\nlogical_z = {self.logical_z}"
        return val
    
    def create_tanner_graph(self):
        """
        Creates the tanner graph of the code. Manually creates nodes and edges to have more flexibility and additional parameters,
        instead of something like rx.from_adjacency_matrix(np.hstack([self.hx, self.hz])).

        """
        tanner = rx.PyGraph()

        offset = self.l*self.m

        l_nodes = [{'index': i,
                    'subindex': i,
                    'type': 'data',
                    'subtype': 'L',
                    'node_attr': {'color': 'blue', 'fillcolor': 'blue', 'style': 'filled', 'shape': 'circle', 'label': f'L_{i}'}}
                    for i in range(self.s)]
        r_nodes = [{'index': i + offset, 
                    'subindex': i,
                    'type': 'data',
                    'subtype': 'R',
                    'node_attr': {'color': 'orange', 'fillcolor': 'orange', 'style': 'filled', 'shape': 'circle', 'label': f'R_{i}'}}
                    for i in range(self.s)]
        x_nodes = [{'index': i + 2*offset, 
                    'subindex': i,
                    'type': 'check',
                    'subtype': 'X',
                    'node_attr': {'color': 'red', 'fillcolor': 'red', 'style': 'filled', 'shape': 'square', 'label': f'X_{i}'}}
                    for i in range(self.s)]
        z_nodes = [{'index': i + 3*offset, 
                    'subindex': i,
                    'type': 'check',
                    'subtype': 'Z',
                    'node_attr': {'color': 'green', 'fillcolor': 'green', 'style': 'filled', 'shape': 'square', 'label': f'Z_{i}'}}
                    for i in range(self.s)]

        tanner.add_nodes_from(l_nodes)
        tanner.add_nodes_from(r_nodes)
        tanner.add_nodes_from(x_nodes)
        tanner.add_nodes_from(z_nodes)

        for c,q in zip(*np.where(self.A)): # between X and L
            tanner.add_edge(c + 2*self.s, q, None)
        for c,q in zip(*np.where(self.B)): # between X and R
            tanner.add_edge(c + 2*self.s, q + self.s, None)
        for c,q in zip(*np.where(self.B.T)): # between Z and L
            tanner.add_edge(c + 3*self.s, q, None)
        for c,q in zip(*np.where(self.A.T)): # between Z and R
            tanner.add_edge(c + 3*self.s, q + self.s, None)

        self._tanner_graph = tanner

    def draw_tanner(self):
        """
        Draws the tanner graph using rustworkx.visualization.graphviz_draw. Graphviz must be installed.
        """
        return graphviz_draw(self.tanner_graph, node_attr_fn=lambda node: node['node_attr'])

class BBCode:
    """Bivariate Bicycle code data.

    The X and Z gauge operator lists are given as lists of supports.
    There is a consistent qubit ordering of these lists so that
    we can construct gate schedules for circuits.
    """

    def __init__(
        self,
        l: int,
        m: int,
        p1: Optional[Sequence[int]] = None,
        p2: Optional[Sequence[int]] = None
    ) -> None:
        """Initializes a bivariate bicycle code builder

        If p1 is not specified, then A has to be specified manually by set_A.
        If p2 is not specified, then B has to be specified manually by set_B.

        Args:
            l: dimension of first space of generators x and y,
            m: dimension of second space of generators x and y,
            p1 (optional): p1 = (a,b,c) => A = x^a + y^b + y^c. Defaults to None.
            p2 (optional): p2 = (d,e,f) => B = y^d + x^e + x^f. Defaults to None.

        Examples:
            Example 1:
            >>> code = BivariateBicycleCodeBuilder(l=6, m=6, p1=(3,1,2), p2=(3,1,2)).build()

            Example 2:
            >>> builder = BivariateBicycleCodeBuilder(l=6, m=6)
            >>> builder.set_A(builder.gen_x(3), builder.gen_y(1), builder.gen_y(2))
            >>> builder.set_B(builder.gen_y(3), builder.gen_x(1), builder.gen_x(2))
            >>> code = builder.build()

            Examples 1 and 2 produce the same code. The example 1 syntax is easier when creating codes in the standard
            xyy, yxx form. The example 2 syntax allows for more flexible code building with non-standard forms.

        """

        self.l = l
        self.m = m
        self.s = l*m
        self.n = 2*l*m

        if p1 is not None:
            a, b, c = p1
            self.set_A(self.gen_x(a), self.gen_y(b), self.gen_y(c))
        else: self.A = None

        if p2 is not None:
            d, e, f = p2
            self.set_B(self.gen_y(d), self.gen_x(e), self.gen_x(f))
        else: self.B = None

        self._hx = None
        self._hz = None
        self._x_stabilizers = None
        self._z_stabilizers = None
        self._logical_z = None
        self._logical_x = None
        self._tanner_graph = None

    @property
    def hx(self) -> np.array:
        if not self.is_defined():
            raise AttributeError(f'A or B undefined, first set them via set_A and set_B')
        if self._hx is None: self._hx = np.hstack([self.A, self.B])
        return self._hx

    @property
    def hz(self) -> np.array:
        if not self.is_defined():
            raise AttributeError(f'A or B undefined, first set them via set_A and set_B')
        if self._hz is None: self._hz = np.hstack([self.B.T, self.A.T])
        return self._hz
    
    @property
    def x_stabilizers(self) -> PauliList:
        return BBCode.arr_to_indices(self.hx)
        # if self._x_stabilizers is None: self._x_stabilizers = PauliList(np.hstack([self.hx, np.zeros((self.n//2, self.n))]))
        # return self._x_stabilizers
    
    @property
    def z_stabilizers(self) -> PauliList:
        return BBCode.arr_to_indices(self.hz)
        # if self._z_stabilizers is None: self._z_stabilizers = PauliList(np.hstack([np.zeros((self.n//2, self.n)), self.hz]))
        # return self._z_stabilizers
    
    @property
    def x_gauges(self) -> PauliList:
        return self.x_stabilizers
    
    @property
    def z_gauges(self) -> PauliList:
        return self.z_stabilizers
    
    @property
    def logical_z(self):
        if self._logical_z is None:
            full_sym = np.vstack([np.hstack([self.hx, np.zeros((self.n//2, self.n))]),
                                  np.hstack([np.zeros((self.n//2, self.n)), self.hz])])
            center_, x_new, z_new = normalizer(full_sym.astype(np.bool_))
            self._logical_z = PauliList(z_new)
            self._logical_x = PauliList(x_new)

        #return self._logical_z
        return BBCode.arr_to_indices(self._logical_z.matrix[:, self.n:])
    
    @property
    def logical_x(self):
        if self._logical_x is None:
            full_sym = np.vstack([np.hstack([self.hx, np.zeros((self.n//2, self.n))]),
                                  np.hstack([np.zeros((self.n//2, self.n)), self.hz])])
            center_, x_new, z_new = normalizer(full_sym.astype(np.bool_))
            self._logical_z = PauliList(z_new)
            self._logical_x = PauliList(x_new)
        
        #return self._logical_x
        return BBCode.arr_to_indices(self._logical_x.matrix[:, :self.n])
    
    def get_syndrome(self, physical_error: np.ndarray, base: str='z') -> np.ndarray:
        if base == 'z':
            return self.hz @ physical_error % 2
        elif base == 'x':
            return self.hx @ physical_error % 2
        else:
            raise ValueError(f'"base" must be one of {{"x", "y"}}, not {base}.')
        
    def get_logical_error(self, physical_error: np.ndarray, base: str='z') -> np.ndarray:
        if base == 'z':
            if self._logical_z is None:
                full_sym = np.vstack([np.hstack([self.hx, np.zeros((self.n//2, self.n))]),
                                    np.hstack([np.zeros((self.n//2, self.n)), self.hz])])
                center_, x_new, z_new = normalizer(full_sym.astype(np.bool_))
                self._logical_z = PauliList(z_new)
                self._logical_x = PauliList(x_new)
            return self._logical_z.matrix[:, self.n:] @ physical_error % 2
        elif base == 'x':
            if self._logical_x is None:
                full_sym = np.vstack([np.hstack([self.hx, np.zeros((self.n//2, self.n))]),
                                    np.hstack([np.zeros((self.n//2, self.n)), self.hz])])
                center_, x_new, z_new = normalizer(full_sym.astype(np.bool_))
                self._logical_z = PauliList(z_new)
                self._logical_x = PauliList(x_new)
            return self._logical_x.matrix[:, :self.n] @ physical_error % 2
        else:
            raise ValueError(f'"base" must be one of {{"x", "y"}}, not {base}.')

    

    @property
    def x_boundary(self):
        raise NotImplementedError()
    
    @property
    def z_boundary(self):
        raise NotImplementedError()
    
    @property
    def k(self):
        raise NotImplementedError()

    @property
    def d(self):
        raise NotImplementedError()
    
    @property
    def tanner_graph(self):
        if self._tanner_graph is None:
            self.create_tanner_graph()
        return self._tanner_graph
    
    @property
    def tanner_graph_X(self):
        nodes = self.tanner_graph.filter_nodes(lambda node: node['subtype'] != 'Z')
        return self.tanner_graph.subgraph(nodes)
    
    @property
    def tanner_graph_Z(self):
        nodes = self.tanner_graph.filter_nodes(lambda node: node['subtype'] != 'X')
        return self.tanner_graph.subgraph(nodes)

    def __str__(self) -> str:
        """Formatted string."""
        return f"(l={self.l},m={self.m}) bivariate bicycle code"
        #return f"[[{self.n}, {self.k}, {self.d}]] heavy-hexagon compass code"

    def __repr__(self) -> str:
        """String representation."""
        val = str(self)
        val += f"\nx_gauges = {self.x_gauges}"
        val += f"\nz_gauges = {self.z_gauges}"
        val += f"\nx_stabilizers = {self.x_stabilizers}"
        val += f"\nz_stabilizers = {self.z_stabilizers}"
        val += f"\nlogical_x = {self.logical_x}"
        val += f"\nlogical_z = {self.logical_z}"
        val += f"\nx_boundary = {self.x_boundary}"
        val += f"\nz_boundary = {self.z_boundary}"
        return val
    
    def is_defined(self):
        return self.A is not None and self.B is not None
    
    def set_A(self, A1, A2, A3) -> None:
        """
        Set the matrix A as the sum of A1, A2 and A3 modulo 2.

        Parameters
        ----------
        A1 : np.array
            2d np.array of dimension lm x lm, power of generator x or y.
        A2 : np.array
            2d np.array of dimension lm x lm, power of generator x or y.
        A3 : np.array
            2d np.array of dimension lm x lm, power of generator x or y.

        Examples
        --------

        >>> builder.set_A(builder.gen_x(3), builder.gen_y(1), builder.gen_y(2))
        Sets A = x^3 + y + y^2
        """

        self.A = ((A1 + A2 + A3) % 2).astype(np.uint8)

    def set_B(self, B1, B2, B3) -> None:
        """
        Set the matrix B as the sum of B1, B2 and B3 modulo 2.

        Parameters
        ----------
        B1 : np.array
            2d np.array of dimension lm x lm, power of generator x or y.
        B2 : np.array
            2d np.array of dimension lm x lm, power of generator x or y.
        B3 : np.array
            2d np.array of dimension lm x lm, power of generator x or y.

        Examples
        --------

        >>> builder.set_B(builder.gen_y(3), builder.gen_x(1), builder.gen_x(2))
        Sets B = y^3 + x + x^2
        """

        self.B = ((B1 + B2 + B3) % 2).astype(np.uint8)

    def gen_x(self, power: int = 1):
        """
        Calculates and returns x^power

        Parameters
        ----------
        power : int
            Power to which x is raised. Defaults to 1.

        Returns
        --------
        np.array
            x^power

        Examples
        --------

        >>> builder.gen_x(3) # x^3
        >>> builder.gen_x(1) # x
        >>> builder.gen_x() # x
        """

        return np.kron(BBCode.cs_pow(self.l, power=power), np.eye(self.m)).astype(np.uint8)
    
    def gen_y(self, power=1):
        """
        Calculates and returns y^power

        Parameters
        ----------
        power : int
            Power to which y is raised. Defaults to 1.

        Returns
        --------
        np.array
            y^power

        Examples
        --------

        >>> builder.gen_y(3) # y^3
        >>> builder.gen_y(1) # y
        >>> builder.gen_y() # y
        """

        return np.kron(np.eye(self.l), BBCode.cs_pow(self.m, power=power)).astype(np.uint8)
    
    def create_tanner_graph(self):
        """
        Creates the tanner graph of the code. Manually creates nodes and edges to have more flexibility and additional parameters,
        instead of something like rx.from_adjacency_matrix(np.hstack([self.hx, self.hz])).

        """
        tanner = rx.PyGraph()

        offset = self.l*self.m

        l_nodes = [{'index': i,
                    'subindex': i,
                    'type': 'data',
                    'subtype': 'L',
                    'node_attr': {'color': 'blue', 'fillcolor': 'blue', 'style': 'filled', 'shape': 'circle', 'label': f'L_{i}'}}
                    for i in range(self.s)]
        r_nodes = [{'index': i + offset, 
                    'subindex': i,
                    'type': 'data',
                    'subtype': 'R',
                    'node_attr': {'color': 'orange', 'fillcolor': 'orange', 'style': 'filled', 'shape': 'circle', 'label': f'R_{i}'}}
                    for i in range(self.s)]
        x_nodes = [{'index': i + 2*offset, 
                    'subindex': i,
                    'type': 'check',
                    'subtype': 'X',
                    'node_attr': {'color': 'red', 'fillcolor': 'red', 'style': 'filled', 'shape': 'square', 'label': f'X_{i}'}}
                    for i in range(self.s)]
        z_nodes = [{'index': i + 3*offset, 
                    'subindex': i,
                    'type': 'check',
                    'subtype': 'Z',
                    'node_attr': {'color': 'green', 'fillcolor': 'green', 'style': 'filled', 'shape': 'square', 'label': f'Z_{i}'}}
                    for i in range(self.s)]

        tanner.add_nodes_from(l_nodes)
        tanner.add_nodes_from(r_nodes)
        tanner.add_nodes_from(x_nodes)
        tanner.add_nodes_from(z_nodes)

        for c,q in zip(*np.where(self.A)): # between X and L
            tanner.add_edge(c + 2*self.s, q, None)
        for c,q in zip(*np.where(self.B)): # between X and R
            tanner.add_edge(c + 2*self.s, q + self.s, None)
        for c,q in zip(*np.where(self.B.T)): # between Z and L
            tanner.add_edge(c + 3*self.s, q, None)
        for c,q in zip(*np.where(self.A.T)): # between Z and R
            tanner.add_edge(c + 3*self.s, q + self.s, None)

        self._tanner_graph = tanner

    def draw_tanner(self):
        """
        Draws the tanner graph using rustworkx.visualization.graphviz_draw. Graphviz must be installed.
        """
        return graphviz_draw(self.tanner_graph, node_attr_fn=lambda node: node['node_attr'])

    @staticmethod
    def cs_pow(l, power=1):
        """
        Calculates and returns a power of the cyclic shift matrix of size lxl (C_l)

        Parameters
        ----------
        l     : int
            size of cyclic shift matrix.
        power : int
            Power to which Cl is raised. Defaults to 1.

        Returns
        --------
        np.array
            C_l^power

        Examples
        --------

        >>> BivariateBicycleCodeBuilder.cs_pow(3, 2) # (C_3)^2
        """

        return np.roll(np.eye(l), shift=power, axis=1).astype(np.uint8)
    
    @staticmethod
    def arr_to_indices(arr: np.array) -> List[List[int]]:
        """ Converts a numpy array to a list of list of indices where it is non-zero """
        return [np.where(row)[0].tolist() for row in arr]