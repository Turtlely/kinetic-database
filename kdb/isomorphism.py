# Modified version of NetworkX's Isomorphism GraphMatcher code

# Notable modifications:
#    * User can specify specific mappings prior to testing for
#      isomorphism using "partialMapping" argument in __init__() or 
#      GraphMatcher().partialMapping, and the isomorphism test will
#      output only the mappings that are consistent with the user's
#      partial mapping
#    * GraphChecker class (not yet completely finished)  #TODO

# This work was originally coded by Christopher Ellison
# as part of the Computational Mechanics Python (CMPy) project.
# James P. Crutchfield, principal investigator.
# Complexity Sciences Center and Physics Department, UC Davis.

import sys, numpy, math
import networkx as nx
from kdb.common import Kdb

__all__ = ["GraphMatcher", "GraphChecker"]


class GraphMatcher:
    """Implementation of VF2 algorithm for matching undirected graphs.

    Suitable for Graph and MultiGraph instances.
    """

    def __init__(self, G1, G2, G1_Atoms, G2_Atoms, mobileAtoms, partialMapping=None, refnodes=None, swap_mapping=False, dc=0.3):
        """Initialize GraphMatcher.

        Parameters
        ----------
        G1,G2: NetworkX Graph or MultiGraph instances.
           The two graphs to check for isomorphism or monomorphism.

        partialMapping: dict
           Partial mapping provided by the user, where keys are nodes
           in G2 and values are nodes in G1 

        """  # TODO: update description
        self.G1 = G1
        self.G2 = G2
        self.G1_Atoms = G1_Atoms
        self.G2_Atoms = G2_Atoms
        self.G2_Atoms.cell = self.G1_Atoms.cell
        self.mobile = mobileAtoms
        self.G1_nodes = set(G1.nodes())
        self.G2_nodes = set(G2.nodes())
        L1 = len(self.G1)
        L2 = len(self.G2)

        # Set recursion limit.
        self.old_recursion_limit = sys.getrecursionlimit()
        expected_max_recursion_level = L2
        if self.old_recursion_limit < 1.5 * expected_max_recursion_level:
            # Give some breathing room.
            sys.setrecursionlimit(int(1.5 * expected_max_recursion_level))

        # Declare that we will be searching for a graph-graph isomorphism.
        self.test = "graph"

        # Partial mapping already provided by the user
        if partialMapping is not None:
            self.partialMapping = partialMapping
        else:
            self.partialMapping = dict()


        # calculate distances
        kdb = Kdb()
        for i in range(L2):
            for j in range(i, L2):
                kdb.atomAtomDistance(self.G2_Atoms, i, j)
        for i in range(L1):
            for j in range(i, L1):
                kdb.atomAtomPermDistance(self.G1_Atoms, i, j)

        # Ordering of G2 nodes
        self.G2_node_order = self.nodeOrder_plusplus()

        self.dc = dc
        self.swap_mapping = swap_mapping

        # Initialize state
        self.initialize()

    def nodeOrder_plusplus(self):
        G2_node_order = dict.fromkeys(self.partialMapping, -1)  # what is returned from this function
        remaining_nodes = set(self.G2.nodes) - set(self.partialMapping)

        node_rarity = dict()  # dict that returns rarity of a node given the node label
        for G1_node in self.G1.nodes:
            label = self.G1.nodes[G1_node]['symbol']  # TODO: if there are more node attributes, edit this
            if label in node_rarity:
                node_rarity[label] += 1
            else:
                node_rarity[label] = 1
        for finished_node in G2_node_order:
            label = self.G2.nodes[finished_node]['symbol']
            node_rarity[label] -= 1

        node_degrees = self.G2.degree
        order_number = 0  # to be incremented by 1 after a node is added to G2_node_order
        previous_order_number = -1
        process_partialMapping = True if self.partialMapping else False
        master_node = None

        while remaining_nodes:
            if process_partialMapping:
                master_node = '...a_random_node_name_that_hopefully_is_not_already_in_G2...'
                self.G2.add_edges_from([(master_node, node) for node in self.partialMapping])

                root_node = master_node
            else:
                # get the rarest node with the largest degree
                root_node = self.argmax(remaining_nodes, node_rarity, reverse=True, use_label=True, label_graph=self.G2)
                root_node = max(root_node, key=node_degrees)
                G2_node_order[root_node] = order_number
                order_number += 1
                remaining_nodes.remove(root_node)
                node_rarity[self.G2.nodes[root_node]['symbol']] -= 1

            G2_bfsTree = nx.bfs_edges(self.G2, source=root_node)
            parent_nodes = {root_node}

            # generating and processing set of all nodes at a given depth of G2_bfsTree
            # taking advantage of the fact that G2_bfsTree is a generator that yields edges in depth order
            node_at_depth = set()
            for edge in G2_bfsTree:
                if edge[0] in parent_nodes:
                    node_at_depth.add(edge[1])
                else:
                    # the following arguments except 'node_degrees' will be modified inside the function
                    parent_nodes, order_number = self.process_bfsTree(node_at_depth, parent_nodes, remaining_nodes,
                                                                      G2_node_order, node_degrees, node_rarity,
                                                                      order_number, process_partialMapping)
                    node_at_depth = {edge[1]}
                    if process_partialMapping:
                        process_partialMapping = False
            # for the very last depth
            if node_at_depth:  #XXX: isn't this condition always satisfied?
                parent_nodes, order_number = self.process_bfsTree(node_at_depth, parent_nodes, remaining_nodes,
                                                                  G2_node_order, node_degrees, node_rarity,
                                                                  order_number, process_partialMapping)
            if process_partialMapping:
                process_partialMapping = False

            previous_order_number = order_number

        if master_node is not None:
            self.G2.remove_node(master_node)
        return G2_node_order

    def argmax(self, nodes_to_consider, key, reverse=False, use_label=False, label_graph=None):
        if len(nodes_to_consider) == 1:
            return nodes_to_consider.copy()
        first = True
        for node in nodes_to_consider:
            if use_label:
                val = key[label_graph.nodes[node]['symbol']]
            else:
                val = key[node]
            if reverse:
                val = -val
            if first or val > maxval:
                maxval = val
                to_return = {node}
                first = False
            elif val == maxval:
                to_return.add(node)
        return to_return

    def process_bfsTree(self, node_at_depth, parent_nodes, remaining_nodes, G2_node_order, node_degrees, node_rarity, order_number, proc_part_map):
        if proc_part_map:
            return node_at_depth, order_number  # here, node_at_depth should be set of G2_nodes in self.partialMapping
        new_parent_nodes = node_at_depth.copy()

        # connectivity to nodes already in G2_node_order
        connM = dict()
        for node in node_at_depth:
            neighbors = parent_nodes.intersection(self.G2.adj[node])
            connM[node] = len(neighbors)

        while node_at_depth:
            selected_node = self.argmax(node_at_depth, connM)
            selected_node = self.argmax(selected_node, node_degrees)
            selected_node = self.argmax(selected_node, node_rarity, reverse=True, use_label=True, label_graph=self.G2).pop()
            remaining_nodes.remove(selected_node)
            node_at_depth.remove(selected_node)
            G2_node_order[selected_node] = order_number
            order_number += 1

            # refresh node_rarity and connM
            node_rarity[self.G2.nodes[selected_node]['symbol']] -= 1
            for node in node_at_depth:
                if node in self.G2.adj[selected_node]:
                    connM[node] += 1
        return new_parent_nodes, order_number

    def reset_recursion_limit(self):
        """Restores the recursion limit."""
        # TODO:
        # Currently, we use recursion and set the recursion level higher.
        # It would be nice to restore the level, but because the
        # (Di)GraphMatcher classes make use of cyclic references, garbage
        # collection will never happen when we define __del__() to
        # restore the recursion level. The result is a memory leak.
        # So for now, we do not automatically restore the recursion level,
        # and instead provide a method to do this manually. Eventually,
        # we should turn this into a non-recursive implementation.
        sys.setrecursionlimit(self.old_recursion_limit)

    def candidate_pairs_iter(self):
        """Iterator over candidate pairs of nodes in G1 and G2."""

        # All computations are done using the current state!
        G1_nodes = self.G1_nodes
        G2_nodes = self.G2_nodes
        min_key = self.G2_node_order.__getitem__

        # First we compute the inout-terminal sets.
        T1_inout = [node for node in self.inout_1 if node not in self.core_1]
        T2_inout = [node for node in self.inout_2 if node not in self.core_2]

        # If T1_inout and T2_inout are both nonempty.
        # P(s) = T1_inout x {min T2_inout}
        if T1_inout and T2_inout:
            node_2 = min(T2_inout, key=min_key)
            for node_1 in T1_inout:
                yield node_1, node_2

        else:
            # If T1_inout and T2_inout were both empty....
            # P(s) = (N_1 - M_1) x {min (N_2 - M_2)}
            # if not (T1_inout or T2_inout):  # as suggested by  [2], incorrect
            if 1:  # as inferred from [1], correct
                # First we determine the candidate node for G2
                other_node = min(G2_nodes - set(self.core_2), key=min_key)
                for node in self.G1:
                    if node not in self.core_1:
                        yield node, other_node
        # For all other cases, we don't have any candidate pairs.

    def initialize(self):
        """Reinitializes the state of the algorithm.

        This method should be redefined if using something other than GMState.
        If only subclassing GraphMatcher, a redefinition is not necessary.

        """

        # core_1[n] contains the index of the node paired with n, which is m,
        #           provided n is in the mapping.
        # core_2[m] contains the index of the node paired with m, which is n,
        #           provided m is in the mapping.
        self.core_1 = {}
        self.core_2 = {}

        # See the paper for definitions of M_x and T_x^{y}

        # inout_1[n]  is non-zero if n is in M_1 or in T_1^{inout}
        # inout_2[m]  is non-zero if m is in M_2 or in T_2^{inout}
        #
        # The value stored is the depth of the SSR tree when the node became
        # part of the corresponding set.
        self.inout_1 = {}
        self.inout_2 = {}
        # Practically, these sets simply store the nodes in the subgraph.

        self.state = GMState(self)

        # Provide a convenient way to access the isomorphism mapping.
        self.mapping = self.core_1.copy()


    def is_isomorphic(self):
        """Returns True if G1 and G2 are isomorphic graphs."""

        # Let's do two very quick checks!
        # QUESTION: Should we call faster_graph_could_be_isomorphic(G1,G2)?
        # For now, I just copy the code.

        # Check global properties
        if self.G1.order() != self.G2.order():
            return False

        # Check local properties
        d1 = sorted(d for n, d in self.G1.degree())
        d2 = sorted(d for n, d in self.G2.degree())
        if d1 != d2:
            return False

        try:
            x = next(self.isomorphisms_iter())
            return True
        except StopIteration:
            return False


    def isomorphisms_iter(self):
        """Generator over isomorphisms between G1 and G2."""
        # Declare that we are looking for a graph-graph isomorphism.
        self.test = "graph"
        self.initialize()
        yield from self.match()


    def match(self):
        """Extends the isomorphism mapping.

        This function is called recursively to determine if a complete
        isomorphism can be found between G1 and G2.  It cleans up the class
        variables after each recursive call. If an isomorphism is found,
        we yield the mapping.

        """
        if len(self.core_1) == len(self.G2):
            # Save the final mapping, otherwise garbage collection deletes it.
            if self.swap_mapping:
                self.mapping = self.core_2.copy()
            else:
                self.mapping = self.core_1.copy()
            # The mapping is complete.
            yield self.mapping
 
        else:
            for G1_node, G2_node in self.candidate_pairs_iter():
                if self.feasibility(G1_node, G2_node):
                    # Recursive call, adding the feasible state.
                    newstate = self.state.__class__(self, G1_node, G2_node)
                    yield from self.match()

                    # restore data structures
                    newstate.restore()

    def subgraph_is_isomorphic(self):
        """Returns True if a subgraph of G1 is isomorphic to G2."""
        try:
            x = next(self.subgraph_isomorphisms_iter())
            return True
        except StopIteration:
            return False


    def subgraph_is_monomorphic(self):
        """Returns True if a subgraph of G1 is monomorphic to G2."""
        try:
            x = next(self.subgraph_monomorphisms_iter())
            return True
        except StopIteration:
            return False

    #    subgraph_is_isomorphic.__doc__ += "\n" + subgraph.replace('\n','\n'+indent)

    def subgraph_isomorphisms_iter(self):
        """Generator over isomorphisms between a subgraph of G1 and G2."""
        # Declare that we are looking for graph-subgraph isomorphism.
        self.test = "subgraph"
        self.initialize()
        yield from self.match()


    def subgraph_monomorphisms_iter(self):
        """Generator over monomorphisms between a subgraph of G1 and G2."""
        # Declare that we are looking for graph-subgraph monomorphism.
        self.test = "mono"
        self.initialize()
        yield from self.match()


    def feasibility(self, G1_node, G2_node):
        """Returns True if adding (G1_node, G2_node) is syntactically feasible.

        This function returns True if it is adding the candidate pair
        to the current partial isomorphism/monomorphism mapping is allowable.
        The addition is allowable if the inclusion of the candidate pair does
        not make it impossible for an isomorphism/monomorphism to be found.

        Uses cutting rules from VF2++ (SJ)
        """

        # The VF2 algorithm was designed to work with graphs having, at most,
        # one edge connecting any two nodes.  This is not the case when
        # dealing with an MultiGraphs.
        #
        # Basically, when we test the look-ahead rules R_neighbor, we will
        # make sure that the number of edges are checked. We also add
        # a R_self check to verify that the number of selfloops is acceptable.
        #
        # Users might be comparing Graph instances with MultiGraph instances.
        # So the generic GraphMatcher class must work with MultiGraphs.
        # Care must be taken since the value in the innermost dictionary is a
        # singlet for Graph instances.  For MultiGraphs, the value in the
        # innermost dictionary is a list.

        ###
        # Test at each step to get a return value as soon as possible.
        ###

        # Look ahead 0

        # R_self

        # The number of selfloops for G1_node must equal the number of
        # self-loops for G2_node. Without this check, we would fail on
        # R_neighbor at the next recursion level. But it is good to prune the
        # search tree now.

        # EDIT(SJ): combining syntatic and semantic feasibilities into a single
        #           function. All semantic checks are followed by "#SEMANTIC"
        
        #SEMANTIC
        if self.G1.nodes[G1_node]['symbol'] != self.G2.nodes[G2_node]['symbol']:
            return False

        #SEMANTIC
        # mobile atom must not be mapped to a constrained atom
        if self.G1_Atoms.constraints:  # has 1 or more constrained atoms(s)
            if G2_node in self.mobile:
                if G1_node in self.G1_Atoms.constraints[0].index:
                    return False

        #SEMANTIC
        for G1_other_node, G2_other_node in self.core_1.items():
            d2 = self.G2_Atoms.distances[(G2_node, G2_other_node)]
            d1s = self.G1_Atoms.permDistances[(G1_node, G1_other_node)]
            passed = False
            for d1 in d1s:
                if abs(d1 - d2) < 2*self.dc:
                    passed = True
                    break
            if not passed:
                return False

        if self.test == "mono":
            if self.G1.number_of_edges(G1_node, G1_node) < self.G2.number_of_edges(
                G2_node, G2_node
            ):
                return False
        else:
            if self.G1.number_of_edges(G1_node, G1_node) != self.G2.number_of_edges(
                G2_node, G2_node
            ):
                return False

        # R_neighbor

        # For each neighbor n' of n in the partial mapping, the corresponding
        # node m' is a neighbor of m, and vice versa. Also, the number of
        # edges must be equal.
        if self.test != "mono":
            for neighbor in self.G1[G1_node]:
                if neighbor in self.core_1:
                    if not (self.core_1[neighbor] in self.G2[G2_node]):
                        return False
                    elif self.G1.number_of_edges(
                        neighbor, G1_node
                    ) != self.G2.number_of_edges(self.core_1[neighbor], G2_node):
                        return False

        for neighbor in self.G2[G2_node]:
            if neighbor in self.core_2:
                if not (self.core_2[neighbor] in self.G1[G1_node]):
                    return False
                elif self.test == "mono":
                    if self.G1.number_of_edges(
                        self.core_2[neighbor], G1_node
                    ) < self.G2.number_of_edges(neighbor, G2_node):
                        return False
                else:
                    if self.G1.number_of_edges(
                        self.core_2[neighbor], G1_node
                    ) != self.G2.number_of_edges(neighbor, G2_node):
                        return False

        if self.test != "mono":
            # Look ahead 1

            # R_terminout
            # The number of neighbors of n in T_1^{inout} is equal to the
            # number of neighbors of m that are in T_2^{inout}, and vice versa.
            # EDIT (SJ): checks by node attribute
            neighbor_symbols = {}  # TODO: if there are node labels other than symbol, edit this to check for all
            for neighbor in self.G1[G1_node]:
                if (neighbor in self.inout_1) and (neighbor not in self.core_1):
                    #SEMANTIC
                    neighbor_symbol = self.G1.nodes[neighbor]['symbol']
                    if neighbor_symbol in neighbor_symbols:
                        neighbor_symbols[neighbor_symbol] += 1
                    else:
                        neighbor_symbols[neighbor_symbol] = 1
                    
            for neighbor in self.G2[G2_node]:
                if (neighbor in self.inout_2) and (neighbor not in self.core_2):
                    #SEMANTIC
                    neighbor_symbol = self.G2.nodes[neighbor]['symbol']
                    try:
                        neighbor_symbols[neighbor_symbol] -= 1
                        if neighbor_symbols[neighbor_symbol] < 0:
                            return False
                    except KeyError:
                        # there is a neighbor of G2_node in inout_2 whose node attribute is not present for any neighbor of G1_node in inout _1
                        return False
            if self.test == "graph":
                if sum(neighbor_symbols.values()):
                    return False

            # Look ahead 2

            # R_new

            # The number of neighbors of n that are neither in the core_1 nor
            # T_1^{inout} is equal to the number of neighbors of m
            # that are neither in core_2 nor T_2^{inout}.
            # EDIT (SJ): checks by node attribute
            neighbor_symbols = {}  # TODO: if there are node labels other than symbol, edit this to check for all
            for neighbor in self.G1[G1_node]:
                if neighbor not in self.inout_1:
                    #SEMANTIC
                    neighbor_symbol = self.G1.nodes[neighbor]['symbol']
                    if neighbor_symbol in neighbor_symbols:
                        neighbor_symbols[neighbor_symbol] += 1
                    else:
                        neighbor_symbols[neighbor_symbol] = 1

            for neighbor in self.G2[G2_node]:
                if neighbor not in self.inout_2:
                    #SEMANTIC
                    neighbor_symbol = self.G2.nodes[neighbor]['symbol']
                    try:
                        neighbor_symbols[neighbor_symbol] -= 1
                        if neighbor_symbols[neighbor_symbol] < 0:
                            return False
                    except KeyError:
                        return False

            if self.test == "graph":
                if sum(neighbor_symbols.values()):
                    return False

        # Otherwise, this node pair is syntactically feasible!
        return True


class GMState:
    """Internal representation of state for the GraphMatcher class.

    This class is used internally by the GraphMatcher class.  It is used
    only to store state specific data. There will be at most G2.order() of
    these objects in memory at a time, due to the depth-first search
    strategy employed by the VF2 algorithm.
    """

    def __init__(self, GM, G1_node=None, G2_node=None):
        """Initializes GMState object.

        Pass in the GraphMatcher to which this GMState belongs and the
        new node pair that will be added to the GraphMatcher's current
        isomorphism mapping.
        """
        self.GM = GM

        # Initialize the last stored node pair.
        self.G1_node = None
        self.G2_node = None
        self.depth = len(GM.core_1)

        if G1_node is None or G2_node is None:
            # Then we reset the class variables or set it to partial mapping
            if not GM.partialMapping:
                GM.core_1 = {}
                GM.core_2 = {}
                GM.inout_1 = {}
                GM.inout_2 = {}
            else:
                self.initialize_from_partialMapping()

        # Watch out! G1_node == 0 should evaluate to True.
        if G1_node is not None and G2_node is not None:
            # Add the node pair to the isomorphism mapping.
            GM.core_1[G1_node] = G2_node
            GM.core_2[G2_node] = G1_node

            # Store the node that was added last.
            self.G1_node = G1_node
            self.G2_node = G2_node

            # Now we must update the other two vectors.
            # We will add only if it is not in there already!
            self.depth = len(GM.core_1)

            # First we add the new nodes...
            if G1_node not in GM.inout_1:
                GM.inout_1[G1_node] = self.depth
            if G2_node not in GM.inout_2:
                GM.inout_2[G2_node] = self.depth

            # Now we add every other node...

            # Updates for T_1^{inout}
            new_nodes = set()
            for node in GM.core_1:
                new_nodes.update(
                    [neighbor for neighbor in GM.G1[node] if neighbor not in GM.core_1]
                )
            for node in new_nodes:
                if node not in GM.inout_1:
                    GM.inout_1[node] = self.depth

            # Updates for T_2^{inout}
            new_nodes = set()
            for node in GM.core_2:
                new_nodes.update(
                    [neighbor for neighbor in GM.G2[node] if neighbor not in GM.core_2]
                )
            for node in new_nodes:
                if node not in GM.inout_2:
                    GM.inout_2[node] = self.depth
            

    def initialize_from_partialMapping(self):
        """ Initializes GM.core_1, GM.core_2, GM.inout_1, and GM.inout_2
            from user's partial mapping

        """
        self.GM.core_2 = self.GM.partialMapping
        for node2, node1 in self.GM.partialMapping.items():
            self.GM.core_1[node1] = node2
            self.GM.inout_1.update({node1 : -1})
            for node in self.GM.G1[node1]:
                self.GM.inout_1.update({node : -1})
            self.GM.inout_2.update({node2 : -1})
            for node in self.GM.G2[node2]:
                self.GM.inout_2.update({node : -1})


    def restore(self):
        """Deletes the GMState object and restores the class variables."""
        # First we remove the node that was added from the core vectors.
        # Watch out! G1_node == 0 should evaluate to True.
        if self.G1_node is not None and self.G2_node is not None:
            del self.GM.core_1[self.G1_node]
            del self.GM.core_2[self.G2_node]

        # Now we revert the other two vectors.
        # Thus, we delete all entries which have this depth level.
        for vector in (self.GM.inout_1, self.GM.inout_2):
            for node in list(vector.keys()):
                if vector[node] == self.depth:
                    del vector[node]

