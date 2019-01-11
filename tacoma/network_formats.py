# -*- coding: utf-8 -*-
"""
This module provides additional temporal network classes only available in the Python implementation.
"""

from copy import deepcopy
import numpy as np

from scipy.sparse import csc_matrix
from scipy.sparse import coo_matrix
import scipy

from _tacoma import edge_changes as ec
from _tacoma import edge_lists as el
from _tacoma import edge_lists_with_histograms as el_h
from _tacoma import edge_changes_with_histograms as ec_h

import _tacoma as _tc


class sparse_adjacency_matrices():
    """
    Construct a temporal network as a list of sparse adjacency matrices.

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        An instance of a temporal network.
    sparse_generator : instance of ``scipy.sparse`` matrix, default : ``scipy.sparse.csc_matrix``
        The sparse matrix class with which to construct all adjacency matrices.
    dtype : numpy.dtype, default=float
        Data type of the matrix entries.
    """



    def __init__(self,temporal_network,sparse_generator=csc_matrix,dtype=float):


        if type(temporal_network) in [ec, ec_h]:
            result = _tc.convert_edge_changes(temporal_network, verbose=False)
        elif type(temporal_network) in [el, el_h]:
            pass
        else:
            raise ValueError('Please provide either an instance of `tacoma.edge_changes` or `tacoma.edge_lists`')

        self.t = np.array(temporal_network.t)
        self.N = deepcopy(temporal_network.N)
        self.int_to_node = deepcopy(temporal_network.int_to_node)
        self.notes = deepcopy(temporal_network.notes)
        self.time_unit = deepcopy(temporal_network.time_unit)
        self.tmax = deepcopy(temporal_network.tmax)
        self.sparse_generator = sparse_generator
        self.dtype = dtype

        self._generate_adjacency_matrices(temporal_network.edges)

    def _generate_adjacency_matrices(self,edges):

        self.adjacency_matrices = []

        for this_el in iter(edges):
            coords = np.array(this_el,dtype=int)
            i, j = coords[:,0], coords[:,1]
            iall = np.concatenate((i,j))
            jall = np.concatenate((j,i)) 
            data = np.ones_like(iall,dtype=float)

            A = self.sparse_generator((data, (iall, jall)), shape=(self.N, self.N))

            self.adjacency_matrices.append(A)



class adjacency_matrices():
    """
    Construct a temporal network as a list of adjacency matrices (``numpy.array``).

    Parameters
    ----------
    temporal_network : :class:`_tacoma.edge_changes`, :class:`_tacoma.edge_lists`, :class:`_tacoma.edge_changes_with_histograms`, or :class:`_tacoma.edge_lists_with_histograms`
        An instance of a temporal network.
    dtype : numpy.dtype, default=float
        Data type of the matrix entries.
    """

    def __init__(self,temporal_network,dtype=float):

        if type(temporal_network) in [ec, ec_h]:
            result = _tc.convert_edge_changes(temporal_network, verbose=False)
        elif type(temporal_network) in [el, el_h]:
            pass
        else:
            raise ValueError('Please provide either an instance of `tacoma.edge_changes` or `tacoma.edge_lists`')

        self.t = np.array(temporal_network.t)
        self.N = deepcopy(temporal_network.N)
        self.int_to_node = deepcopy(temporal_network.int_to_node)
        self.notes = deepcopy(temporal_network.notes)
        self.time_unit = deepcopy(temporal_network.time_unit)
        self.tmax = deepcopy(temporal_network.tmax)
        self.dtype = dtype

        self._generate_adjacency_matrices(temporal_network.edges)

    def _generate_adjacency_matrices(self,edges):

        self.adjacency_matrices = []

        for this_el in iter(edges):
            coords = np.array(this_el,dtype=int)
            i, j = coords[:,0], coords[:,1]
            iall = np.concatenate((i,j))
            jall = np.concatenate((j,i)) 
            data = np.ones_like(iall,dtype=float)

            A = coo_matrix((data, (iall, jall)), shape=(self.N, self.N))

            self.adjacency_matrices.append(A.toarray())

if __name__ == "__main__":

    import tacoma as tc


    L = tc.edge_lists()

    L.N = 3
    L.t = [0.0,1.0,2.0]
    L.tmax = 3.0
    L.edges = [ 
                [
                  (0,1)
                ],
                [
                  (1,2), (0,2)
                ],
                [
                  (0,1)
                ],
               ]

    A = adjacency_matrices(L)

    for a in A.adjacency_matrices:
        print(a)

    A = sparse_adjacency_matrices(L)

    for a in A.adjacency_matrices:
        print(a.toarray())

