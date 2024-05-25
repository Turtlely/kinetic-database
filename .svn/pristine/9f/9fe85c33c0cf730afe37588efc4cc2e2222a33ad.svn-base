# module to make the graph with networkx 
#!/usr/bin/env python

import warnings
try:
    import networkx as nx
    import matplotlib.pyplot as plt
except ImportError:
    warnings.warn('failed to load networkx, using default method', ImportWarning)
from kdb.common import Kdb
from kdb.config import *
from kdb.aselite import elements, NeighborList_helper, NeighborList
import numpy
import time, copy


class GraphMaker(Kdb):
    def __init__(self):
        pass

    def graph_reactant(self, r, nf):
        """ Given a list of selected atoms from kdbinsert, make a
        connectivity graph so that it can stored and used to be queried
        according to neighbor fudge
        
        Parameters
        ----------
        a: list
           List of atom id's in the reactant to be used as node numbers
           (e.g. [0, 1, 2, 3, ...]

        r: ASE Atoms Object
           The reactant structure to query and make graph for

        nf: float
            Neighbor fudge (don't think it's used here. Check and remove)

        """ #XXX: update description

        a = range(len(r))
        r_G = nx.Graph()
        #cutoff = NeighborList_helper.natural_cutoffs(r, mult=(1+nf))
        #nl = NeighborList(cutoff, skin=0.0)
        cutoff = NeighborList_helper.natural_cutoffs(r)
        nl = NeighborList(cutoff)
        nl.update(r)
        for atom in a:
             if USE_GROUP_SIMILARITY: 
                 r_G.add_node(atom,symbol = elements[r.get_chemical_symbols()[atom]]["elemental_type"]) #this reads elemental type if similarity is on in config
             else:
                 r_G.add_node(atom,symbol = elements[r.get_chemical_symbols()[atom]]["symbol"])
        for atom in a:
            indices, offsets = nl.get_neighbors(atom)

            for neighbor in indices:
                if (atom == neighbor):
                    continue
                else:
                    #r_G.add_edge(atom ,neighbor, distance=Kdb().atomAtomPbcDistance(r, atom, neighbor))
                    r_G.add_edge(atom, neighbor)
        return r_G


    def graph_reactant_2(self, r, nf):
        """ makes graph without NeighborList
        """
        a = range(len(r))
        r_G = nx.Graph()
        cutoff = NeighborList_helper.natural_cutoffs(r)
        for atom in a:
            if USE_GROUP_SIMILARITY:
                r_G.add_node(atom,symbol = elements[r.get_chemical_symbols()[atom]]["elemental_type"]) #this reads elemental type if similarity is on in config
            else:
                r_G.add_node(atom,symbol = elements[r.get_chemical_symbols()[atom]]["symbol"])
        for atom1 in a:
            # getting neighbors
            for atom2 in range(atom1+1, len(a)):
                r1 = cutoff[atom1]
                r2 = cutoff[atom2]
                d = Kdb.atomAtomPbcDistance(self, r, atom1, atom2)
                if d < (r1 + r2) * (1 + nf):
                    r_G.add_edge(atom1, atom2)
        return r_G

    
    def graph_kdbentry(self, r, rr, nf):
        """ Given a list of positions for a kdbentry, make a connectivity 
        graph so that it can be stored and used to compare to the queried 
        reactant according to neigbor fudge, then removes all unconnected
        nodes from the connectivity graph and calculates the pbc
        distances between an unconnected node and all connected nodes,
        for all unconnected nodes

        Parameters  
        ----------
        a: list
           List of atom id's in kdbentry to be used as node numbers (e.g.
           [0, 1, 2, 3, ...]

        r: ASE Atoms Object
           Graph is made based on r

        rr: ASE Atoms Object
            The reactant structure to query

        """

        a = range(len(r))
        r_G = nx.Graph()
        r.cell = rr.cell
        r.set_pbc((True, True, True)) # NOTE: is this always true?
        #cutoff = NeighborList_helper.natural_cutoffs(r, mult=(1+nf))
        #nl = NeighborList(cutoff, skin=0.0, self_interaction=False)
        cutoff = NeighborList_helper.natural_cutoffs(r)
        nl = NeighborList(cutoff, self_interaction=False)
        nl.update(r)

        for atom in a:
            indices, offsets = nl.get_neighbors(atom)
            if USE_GROUP_SIMILARITY:
                r_G.add_node(atom, symbol=elements[r.get_chemical_symbols()[atom]]["elemental_type"])
            else:
                r_G.add_node(atom, symbol=elements[r.get_chemical_symbols()[atom]]["symbol"])
            for neighbor in indices:
                r_G.add_edge(atom, neighbor)

        return r_G

