#!/usr/bin/env python

from math import fmod
import numpy
from kdb.aselite import elements
from kdb.aselite import FixAtoms
from kdb.aselite import write_con, write_vasp
from kdb.common import Kdb
from kdb.kdbquery import KdbQuery
import kdb.fmodules as fmodules
from kdb.graph import GraphMaker
import networkx as nx

__all__ = ["KdbInsert"]

force_insert = False  # set to True to insert even with errors

class KdbInsert(Kdb):
    def __init__(self):
        pass

    def coordination_numbers(self, p, nf):
        nl = []
        for a in range(len(p)):
            nl.append([])
            for b in range(len(p)):
                if b != a:
                    dist = numpy.linalg.norm(p.get_positions()[a] - p.get_positions()[b])
                    if dist < (elements[p.get_chemical_symbols()[a]]["radius"] + 
                               elements[p.get_chemical_symbols()[b]]["radius"]) * (1.0 + nf):
                        nl[a].append(b)
        return [len(l) for l in nl]

    def getMappings(self, a, b, nf, dc, mappings = None):
        """ A recursive depth-first search for a complete set of mappings from atoms
            in configuration a to atoms in configuration b. Do not use the mappings
            argument, this is only used internally for recursion.

            Returns None if no mapping was found, or a dictionary mapping atom 
            indices a to atom indices b.

            Note: If a and b are mirror images, this function will still return a 
            mapping from a to b, even though it may not be possible to align them 
            through translation and rotation. """
        # If this is the top-level user call, create and loop through top-level
        # mappings.
        if mappings == None:
            # Find the least common coordination number in b.
            bCoordinations = self.coordination_numbers(b, nf)
            bCoordinationsCounts = {}
            for coordination in bCoordinations:
                if coordination in bCoordinationsCounts:
                    bCoordinationsCounts[coordination] += 1
                else:
                    bCoordinationsCounts[coordination] = 1
            bLeastCommonCoordination = list(bCoordinationsCounts.keys())[0]
            for coordination in list(bCoordinationsCounts.keys()):
                if bCoordinationsCounts[coordination] < bCoordinationsCounts[bLeastCommonCoordination]:
                    bLeastCommonCoordination = coordination
            # Find one atom in a with the least common coordination number in b. 
            # If it does not exist, return None.
            aCoordinations = self.coordination_numbers(a, nf)
            try:
                aAtom = aCoordinations.index(bLeastCommonCoordination)
            except ValueError:
                return None
            # Create a mapping from the atom chosen from a to each of the atoms with
            # the least common coordination number in b, and recurse.
            for i in range(len(bCoordinations)):
                if bCoordinations[i] == bLeastCommonCoordination:
                    # Make sure the element types are the same.
                    if a.get_chemical_symbols()[aAtom] != b.get_chemical_symbols()[i]:
                        continue
                    mappings = self.getMappings(a, b, nf, dc, {aAtom:i})
                    # If the result is not none, then we found a successful mapping.
                    if mappings is not None:
                        return mappings
            # There were no mappings.
            return None

        # This is a recursed invocation of this function.
        else:
            # Find an atom from a that has not yet been mapped.
            unmappedA = 0
            while unmappedA < len(a):
                if unmappedA not in list(mappings.keys()):
                    break
                unmappedA += 1
            # Calculate the distances from unmappedA to all mapped a atoms.
            distances = {}
            for i in list(mappings.keys()):
                distances[i] = self.atomAtomDistance(a, unmappedA, i)

            # Loop over each unmapped b atom. Compare the distances between it and 
            # the mapped b atoms to the corresponding distances between unmappedA 
            # and the mapped atoms. If everything is similar, create a new mapping
            # and recurse.
            for bAtom in range(len(b)):
                if bAtom not in list(mappings.values()):
                    for aAtom in distances:
                        # Break if type check fails.
                        if b.get_chemical_symbols()[bAtom] != a.get_chemical_symbols()[unmappedA]:
                            break
                        # Break if distance check fails
                        bDist = self.atomAtomDistance(b, bAtom, mappings[aAtom])
                        if abs(distances[aAtom] - bDist) > 0.2:#> dc original +100
                            break
                    else:
                        # All distances were good, so create a new mapping.
                        newMappings = mappings.copy()
                        newMappings[unmappedA] = bAtom
                        # If this is now a complete mapping from a to b, return it.
                        if len(newMappings) == len(a):
                            return newMappings
                        # Otherwise, recurse.
                        newMappings = self.getMappings(a, b, nf, dc, newMappings)
                        # Pass any successful mapping up the recursion chain.
                        if newMappings is not None:
                            return newMappings
            # There were no mappings.
            return None


    def stripUnselectedAtoms(self, atoms, selected):
        """ Removes any atoms from atoms that are not in selected and returns a new
        structure and a mapping from atoms in the old structure to atoms in the new 
        structure. """
        src = atoms.copy()
        dest = atoms.copy()
        while len(dest) > 0:
            dest.pop()
        mapping = {}
        index = 0
        constraints = []
        for i in selected:
            mapping[i] = index
            index += 1
            if i in src.constraints: 
                constraints.append(index)
            dest.append(src[i])
        dest.set_constraint(FixAtoms(constraints))
        return dest, mapping


    def getProcessMobileAtoms(self, r, s, p, mac):
        """ Returns a list of atom indices that move more than mac 
        between reactant and saddle, saddle and product, or 
        reactant and product. If no atoms move more than mac, returns
        the atom that moves the most. """
        mobileAtoms = []
        ibox = numpy.linalg.inv(s.get_cell())
        reactant2saddle = self.per_atom_norm(s.positions - r.positions, s.get_cell(), ibox)
        product2saddle = self.per_atom_norm(s.positions - p.positions, s.get_cell(), ibox)
        reactant2product = self.per_atom_norm(p.positions - r.positions, s.get_cell(), ibox)
        for i in range(len(s)):
            if max(reactant2saddle[i], product2saddle[i], reactant2product[i]) > mac:
                mobileAtoms.append(i)
        if len(mobileAtoms) == 0:
            mobileAtoms.append(list(reactant2product).index(max(reactant2product)))
        return mobileAtoms


    def raise_insert_error(self, err_type, subtype=None):
        types = ["process too large", "mobile too large", "not enough atoms"]
        if err_type == types[0]:
            print("KDB insert error: Process is too large for cell.")
            if subtype == "same neighbor":
                print("\t* same neighbor atom found in multiple PBC boxes")
        elif err_type == types[1]:
            print("KDB insert error: Mobile atoms span a space too large for the cell. Possible causes:")
            print("\t* too many mobile atoms")
            print("\t* movements larger than half the cell.")
        elif err_type == types[2]:
            print("KDB insert error: Not enough atoms in the process. Try increasing neighbor fudge (nf).")


    #function will be overridden in remote/local classes
    def insert_into_db(self, **args):
        print("function not yet overloaded")


    def insert(self, reactant, saddle, product, mode=None, forward_bar=0, reverse_bar=0, nf=0.2, dc=0.3, mac=0.7, kdbname='kdb.db'):

        # keep copy of original
        original_reactant = reactant.copy()
        original_saddle = saddle.copy()
        original_product = product.copy()
        original_mode = mode.copy() if mode is not None else None

        # start with mobile atoms
        ibox = numpy.linalg.inv(saddle.cell)
        mobileAtoms = self.getProcessMobileAtoms(reactant, saddle, product, mac)

        # clump mobile atoms at saddle state, and also find their clumped positions in the saddle & product states
        movementSR = fmodules.kdb.pbcs(reactant.positions[mobileAtoms] - saddle.positions[mobileAtoms], saddle.cell, ibox)
        movementSP = fmodules.kdb.pbcs(product.positions[mobileAtoms] - saddle.positions[mobileAtoms], saddle.cell, ibox)
        saddle.positions = fmodules.kdb.blind_clump(saddle.positions, mobileAtoms, saddle.cell, ibox)  # clump mobile atoms only
        reactant.positions[mobileAtoms] = saddle.positions[mobileAtoms] + movementSR
        product.positions[mobileAtoms] = saddle.positions[mobileAtoms] + movementSP

        # size check
        dimM_R = fmodules.kdbinsert.atoms_dimensions(reactant.positions[mobileAtoms], ibox)
        dimM_S = fmodules.kdbinsert.atoms_dimensions(saddle.positions[mobileAtoms], ibox)
        dimM_P = fmodules.kdbinsert.atoms_dimensions(product.positions[mobileAtoms], ibox)
        if max(dimM_R) >= 0.5 or max(dimM_S) >= 0.5 or max(dimM_P) >= 0.5:  # span more than half the cell length in any direction
            self.raise_insert_error("mobile too large")
            if not force_insert:
                return 1  # 1 is an error code

        # find neighbors and place them at correct location
        neighbors = []
        for atom in range(len(reactant)):
            if atom in mobileAtoms:
                continue
            r2 = elements[reactant.get_chemical_symbols()[atom]]["radius"]  # candidate radius
            already_is_neighbor = False

            for mobileAtom in mobileAtoms:
                r1 = elements[reactant.get_chemical_symbols()[mobileAtom]]["radius"]
                cutoff = (r1 + r2) * (1.0 + nf)

                if self.atomAtomPbcDistance(saddle, mobileAtom, atom) < cutoff:  # is a neighbor in saddle
                    v = self.atomAtomPbcVector(saddle, mobileAtom, atom)
                    temp = saddle.positions[mobileAtom] + v
                    if already_is_neighbor:
                        if fmodules.kdb.atom_atom_distance(temp - saddle.positions[atom]) > dc:
                            self.raise_insert_error("process too large", subtype="same neighbor")
                            if not force_insert:
                                return 1
                    else:
                        saddle.positions[atom] = saddle.positions[mobileAtom] + v
                        already_is_neighbor = True
                        neighbors.append(atom)

                if self.atomAtomPbcDistance(reactant, mobileAtom, atom) < cutoff:
                    v = self.atomAtomPbcVector(reactant, mobileAtom, atom)
                    movement = fmodules.kdb.pbc(saddle.positions[atom] - reactant.positions[atom], saddle.cell, ibox)
                    temp = reactant.positions[mobileAtom] + v + movement
                    if already_is_neighbor:
                        if fmodules.kdb.atom_atom_distance(temp - saddle.positions[atom]) > dc:
                            self.raise_insert_error("process too large", subtype="same neighbor")
                            if not force_insert:
                                return 1
                    else:
                        saddle.positions[atom] = temp
                        already_is_neighbor = True
                        neighbors.append(atom)

                if self.atomAtomPbcDistance(product, mobileAtom, atom) < cutoff:
                    v = self.atomAtomPbcVector(product, mobileAtom, atom)
                    movement = fmodules.kdb.pbc(saddle.positions[atom] - product.positions[atom], saddle.cell, ibox)
                    temp = product.positions[mobileAtom] + v + movement
                    if already_is_neighbor:
                        if fmodules.kdb.atom_atom_distance(temp - saddle.positions[atom]) > dc:
                            self.raise_insert_error("process too large", subtype="same neighbor")
                            if not force_insert:
                                return 1
                    else:
                        saddle.positions[atom] = temp
                        already_is_neighbor = True
                        neighbors.append(atom)
                
        # remove all other atoms
        selectedAtoms = mobileAtoms + neighbors
        if len(selectedAtoms) < 2:
            self.raise_insert_error("not enough atoms")
            if not force_insert:
                return 1
        strippedS, mapping = self.stripUnselectedAtoms(saddle, selectedAtoms)
        strippedSelectedAtoms = [mapping[i] for i in selectedAtoms]
        strippedR = strippedS.copy()
        strippedR.positions[strippedSelectedAtoms] = strippedS.positions[strippedSelectedAtoms] + \
                                                     fmodules.kdb.pbcs(reactant.positions[selectedAtoms] - saddle.positions[selectedAtoms],
                                                        saddle.cell, ibox)
        strippedP = strippedS.copy()
        strippedP.positions[strippedSelectedAtoms] = strippedS.positions[strippedSelectedAtoms] + \
                                                     fmodules.kdb.pbcs(product.positions[selectedAtoms] - saddle.positions[selectedAtoms],
                                                        saddle.cell, ibox)

        # size check 2
        dim_R = fmodules.kdbinsert.atoms_dimensions(strippedR.positions, ibox)
        dim_S = fmodules.kdbinsert.atoms_dimensions(strippedS.positions, ibox)
        dim_P = fmodules.kdbinsert.atoms_dimensions(strippedP.positions, ibox)
        if max(dim_R) >= 1.0 or max(dim_S) >= 1.0 or max(dim_P) >= 1.0:  # span more than full cell length in any direction
            self.raise_insert_error("process too large")
            if not force_insert:
                return 1  # 1 is an error code


        # shift saddle's centroid to origin
        coc = numpy.mean(strippedS.positions, axis=0)
        strippedS.positions -= coc
        strippedR.positions -= coc
        strippedP.positions -= coc

        # give stripped structures a huge box. TODO: is there a better way to remove PBC?
        strippedS.cell = numpy.identity(3) * 1024
        strippedR.cell = strippedS.cell.copy()
        strippedP.cell = strippedS.cell.copy()

        # generate clump order
        clumpOrderR = fmodules.kdb.find_clump_order(strippedR.positions)
        clumpOrderS = fmodules.kdb.find_clump_order(strippedS.positions)
        clumpOrderP = fmodules.kdb.find_clump_order(strippedP.positions)

        # convert clump order to store into database
        clumpOrderR = self.convertClumpOrder(clumpOrderR)
        clumpOrderS = self.convertClumpOrder(clumpOrderS)
        clumpOrderP = self.convertClumpOrder(clumpOrderP)

        # get mobile_list
        mob_list = [ mapping[atom] for atom in mobileAtoms ]

        # Update the mode.
        if mode is not None:
            newMode = numpy.zeros((len(selectedAtoms), 3))
            newMode[strippedSelectedAtoms] = mode[selectedAtoms]
        else:
            newMode = None

        arg_dict = {'or': original_reactant, 'os': original_saddle, 'op': original_product, 'om': original_mode,
                    'r': strippedR, 's': strippedS, 'p': strippedP, 'm': newMode, 'ma': mob_list,
                    'kdbname': kdbname, 'nf': nf, 'dc': dc, 'mac': mac, 'b_f': forward_bar, 'b_r': reverse_bar,
                    'clump_order_r': clumpOrderR, 'clump_order_s': clumpOrderS, 'clump_order_p': clumpOrderP}

        # function is overloaded in either local_insert.py or remote_insert.py
        return self.insert_into_db(**arg_dict)
