#!/usr/bin/env python
import multiprocessing as mp
import os, sys, numpy, glob, shutil, math, copy, json, warnings, time
from optparse import OptionParser
import time
from os import listdir
from os.path import isfile, join
from kdb.common import Kdb
from kdb.config import *
from kdb.aselite import elements, elements_groups, write_vasp, write_con, NeighborList_helper
import kdb.fmodules as fmodules
from functools import partial
try:
    from kdb.graph import GraphMaker
    import networkx as nx
    from kdb.isomorphism import GraphMatcher
except ImportError:
    warnings.warn('failed to import Networkx, will use default method', ImportWarning)

__all__ = ["KdbQuery"]

debug = False  # set to True to print out debugging print statements

class KdbQuery(Kdb):
    def __init__(self):
        self.return_dict = {}
    def isDistance(self, pbcvector, target, box, dc):
        for x in [-1, 0, 1]:
            for y in [-1, 0, 1]:
                for z in [-1, 0, 1]:
                    temp = pbcvector.copy()
                    temp += x * box[0]
                    temp += y * box[1]
                    temp += z * box[2]
                    if abs(numpy.linalg.norm(temp) - target) < dc:
                        return True
        return False
    def isDistance2(self, pbcvector, target, box, dc):
        ref = pbcvector.copy()
        for x in [-1, 0, 1]:
            ref1 = ref + x * box[0]
            for y in [-1, 0, 1]:
                ref2 = ref1 + y * box[1]
                for z in [-1, 0, 1]:
                    ref3 = ref2 + z * box[2]
                    if abs(numpy.linalg.norm(ref3) - target) < dc:
                        return True
        return False
    def isDistance3(self, pbcvector, target, box, dc):
        ref = pbcvector.copy()

        xRange = [-1, 0, 1]
        if numpy.linalg.norm(box[0])/2.0 > target:
           xRange = [0]
        yRange = [-1, 0, 1]
        if numpy.linalg.norm(box[1])/2.0 > target:
           yRange = [0]
        zRange = [-1, 0, 1]
        if numpy.linalg.norm(box[2])/2.0 > target:
           zRange = [0]

        for x in xRange:
            ref1 = ref + x * box[0]
            for y in yRange:
                ref2 = ref1 + y * box[1]
                for z in zRange:
                    ref3 = ref2 + z * box[2]
                    if abs(numpy.linalg.norm(ref3) - target) < dc:
                        return True
        return False
    def centroid(self, a, which=None):
        if which == None:
            which = list(range(len(a)))
        c = numpy.array([0.0, 0.0, 0.0])
        for i in which:
            c += a.positions[i]
        c /= len(which)
        return c

    def query_db(self, **args):
        print("function not yet overloaded")

    # Note this function gets overloaded when interacting with remote DB
    def output_query(self, outputdir, numMatches, suggestion, sugproduct, entry, modeTemp=None, processBarrier=None, score=None):
        write_vasp(outputdir + "/SADDLE_%d_%d" % (entry['id'], numMatches), suggestion)
        write_vasp(outputdir + "/PRODUCT_%d_%d" % (entry['id'], numMatches), sugproduct)
        if modeTemp is not None:
            self.save_mode(outputdir + "/MODE_%d_%d" % (entry['id'], numMatches), modeTemp)
        if processBarrier is not None:
            if (processBarrier[0] !=0 or processBarrier[1] !=0):
                complete_name = os.path.join(outputdir+"/Barrier_%d_%d" % (entry['id'], numMatches)+".txt")
                #print ("complete_name: ",complete_name)
                with open(complete_name,'w') as f:
                    f.write(str(processBarrier))
                f.close()
        # os.system("touch %s/.done_%d" % (outputdir, numMatches))

    def check_overlap(self, kdbAtom_pos, kdbAtom_radius, sug, atoms_to_check, ibox):
        for a in atoms_to_check:
            v = fmodules.kdb.pbc(kdbAtom_pos - sug.positions[a], sug.cell, ibox)
            cutoff = min(sug.radii[a], kdbAtom_radius)
            if numpy.dot(v, v) <= cutoff * cutoff:
                return True
        return False

    def check_duplicates_1(self, newAtoms_sad, newAtoms_prod, uniques, dc, ibox):
        # uniques is a list of Atoms objects
        for unique in uniques:
            pan_sad = self.per_atom_norm(unique[0].positions - newAtoms_sad.positions, newAtoms_sad.cell, ibox)
            pan_prod = self.per_atom_norm(unique[1].positions - newAtoms_prod.positions, newAtoms_prod.cell, ibox)
            if max(pan_sad) <= dc and max(pan_prod) <= dc:
                return True  # is a duplicate
        uniques.append([newAtoms_sad.copy(), newAtoms_prod.copy()])  # is not a duplicate
        return False

    def check_duplicates_2(self, newAtoms_sad, newAtoms_prod, mobile, uniques, dc, ibox):
        # uniques is a list of tuples, whose 1st element is the positions of an Atoms object and 2nd element is an ordered list of mobile atoms from the correponding Atoms object
        # mobile is a list of mapped mobile atoms of newAtoms
        mobile.sort()
        dc2 = dc * dc
        for unique in uniques:
            if mobile == unique[2]:
                diff = fmodules.kdb.pbcs(unique[0].positions[mobile] - newAtoms_sad.positions[mobile], newAtoms_sad.cell, ibox)
                norm2_sad = numpy.sum(numpy.square(diff), axis=1)
                diff = fmodules.kdb.pbcs(unique[1].positions[mobile] - newAtoms_prod.positions[mobile], newAtoms_prod.cell, ibox)
                norm2_prod = numpy.sum(numpy.square(diff), axis=1)
                if max(norm2_sad) <= dc2 and max(norm2_prod) <= dc2:
                    pan_sad = self.per_atom_norm(unique[0].positions - newAtoms_sad.positions, newAtoms_sad.cell, ibox)
                    pan_prod = self.per_atom_norm(unique[1].positions - newAtoms_prod.positions, newAtoms_prod.cell, ibox)
                    if max(pan_sad) <= dc and max(pan_prod) <= dc:
                        return True  # is a duplicate
        uniques.append([newAtoms_sad.copy(), newAtoms_prod.copy(), mobile])
        return False

    def check_duplicates_3(self, newAtoms_sad, newAtoms_prod, mapping_i, mobile, uniques, dc, ibox):
        # uniques is a list of tuples, whose 1st element is the positions of an Atoms object and 2nd element is an ordered list of mobile atoms from the correponding Atoms object
        # mobile is a list of mapped mobile atoms of newAtoms
        mapping_i = sorted(list(mapping_i.values()))
        mobile.sort()
        dc2 = dc * dc
        for unique in uniques:
            all_distances = True
            if mobile == unique[1]:
                if mapping_i == unique[0]:
                    diff = fmodules.kdb.pbcs(unique[2].positions[mobile] - newAtoms_sad.positions[mobile], newAtoms_sad.cell, ibox)
                    norm2_sad = numpy.sum(numpy.square(diff), axis=1)
                    diff = fmodules.kdb.pbcs(unique[3].positions[mobile] - newAtoms_prod.positions[mobile], newAtoms_prod.cell, ibox)
                    norm2_prod = numpy.sum(numpy.square(diff), axis=1)
                    if max(norm2_sad) <= dc2 and max(norm2_prod) <= dc2:
                        return True
        uniques.append([mapping_i, mobile, newAtoms_sad.copy(), newAtoms_prod.copy()])
        return False

    def query_kdb_entry(self, name, reactantNeighbors, reactantNameCount, Reactant_Graph, reactant, outputdir, nf, dc, nodupes, kdbname, GRAPH, entry):
        start = time.time()
        output = ["Process {}: Name: %s ID: %s" % (name, entry['id'])]
        final_data = []
        numMatches = 0
        mirrored = "not mirrored"
        if entry["mirror"]:
            mirrored = "mirrored"

        # Load the minimum
        kdbmin = copy.deepcopy(entry['minimum'])

        # Make sure the reactant has at least as many atoms of each type as the kdb configuration
        passedNameCount = True
        kdbNameCount = self.similar_atom_nameCount(kdbmin,USE_GROUP_SIMILARITY) # Same function call for kdbNameCount
        for name in kdbNameCount:
            if name not in reactantNameCount:
                passedNameCount = False
                break
            if kdbNameCount[name] > reactantNameCount[name]:
                passedNameCount = False
                break
        if not passedNameCount:
            output.append("%10s  name count fail" % entry['id'])
            return ("\n".join(output), self.return_dict, [entry['id'], "name count failed", str(time.time() - start)[0:5]], 0)
        else:
            pass
            #print ("Passed Name Count")
            # Check if using graph method for mappings
        mappings = []

        kdbmobile = copy.deepcopy(entry['mobile'])
        kdbmobile_radii = dict([ (m, elements[kdbmin.get_chemical_symbols()[m]]['radius']) for m in kdbmobile ])

        clump_order = self.convertClumpOrder(entry['clump_order'])

        if GRAPH:
            Selected_Graph = GraphMaker.graph_kdbentry(self, kdbmin, reactant, nf)
            GM = GraphMatcher(Reactant_Graph, Selected_Graph, reactant, kdbmin, kdbmobile, dc=dc, swap_mapping=True)

            final_data = [entry['id'], kdbmobile, Reactant_Graph.number_of_nodes(), Reactant_Graph.number_of_edges(), Selected_Graph.number_of_nodes(), Selected_Graph.number_of_edges(), nx.number_connected_components(Selected_Graph)] if debug else None

            for mapping in GM.subgraph_isomorphisms_iter():
                mappings.append(mapping)

        else:
                # For each mobile atom in kdbmin, create a list of neighboring element
                # types and the count of each type
                kdbNeighbors = {}
                for i in kdbmobile:
                    r1 = elements[kdbmin.get_chemical_symbols()[i]]["radius"]
                    kdbNeighbors[i] = {}
                    for j in range(len(kdbmin)):
                        if j == i:
                            continue
                        r2 = elements[kdbmin.get_chemical_symbols()[j]]["radius"]
                        d = numpy.linalg.norm(kdbmin.positions[i] - kdbmin.positions[j])
                        if d > (r1 + r2) * (1 + nf):
                            continue
                        if kdbmin.get_chemical_symbols()[j] not in kdbNeighbors[i]:
                            kdbNeighbors[i][kdbmin.get_chemical_symbols()[j]] = 0
                        kdbNeighbors[i][kdbmin.get_chemical_symbols()[j]] += 1
                # Finds kdbNeighbors using the elements covalent radius, NOT some pre determined Group_Radii
                kdbNeighbors = Kdb().Neighbor_Dict_similarity(kdbNeighbors,USE_GROUP_SIMILARITY) # Will do nothing if similarity off

                kdbUnmapped = list(range(len(kdbmin))) # Keep track of the kdb atoms that have been mapped
                # Create the initial mappings
                mappings = None
                db_a = kdbmobile[0] # This will be the selected mobile atom
                for m in kdbmobile:
                    mMappings = []
                    for freeAtom in list(reactantNeighbors.keys()):
                        for elementType in reactantNeighbors[freeAtom]:
                            if elementType not in kdbNeighbors[m]:
                                break
                            if kdbNeighbors[m][elementType] != reactantNeighbors[freeAtom][elementType]: # Current code
                    #if kdbNeighbors[m][elementType] > reactantNeighbors[freeAtom][elementType]:
                                break
                            else:
                                mMappings.append({m:freeAtom})
                    if mappings == None:
                        mappings = mMappings
                    if len(mMappings) < len(mappings):
                        mappings = mMappings
                        db_a = m
                kdbUnmapped.remove(db_a)

                ibox = numpy.linalg.inv(reactant.cell)
                d1 = numpy.linalg.norm(reactant.cell[0])
                d2 = numpy.linalg.norm(reactant.cell[1])
                d3 = numpy.linalg.norm(reactant.cell[2])
                micDist = min(d1,d2,d3)/2.0
                output.append("micDist: %s" % micDist)

                while len(kdbUnmapped) > 0 and len(mappings) > 0:

                    # Create a list of new mappings that will replace mappings at the end of this iteration
                    newMappings = []

                    # Select an unmapped atom from kdbmin
                    kdbAtom = kdbUnmapped.pop()

                    # Get the distance between kdbAtom and every other atom in the kdb configuration.
                    kdbDistances = {}
                    for i in range(len(kdbmin)):
                        kdbDistances[i] = numpy.linalg.norm(kdbmin.positions[kdbAtom] - kdbmin.positions[i])

                    # Loop over each mapping and try to place kdbAtom
                    for mapping in mappings:

                        # Loop over each atom in the reactant.
                        for reactantAtom in range(len(reactant)):
                            # Make sure it has not already been mapped.
                            if reactantAtom in list(mapping.values()):
                                continue

                            # Loop over the atoms in mapping and see if the distance
                            # between reactantAtom and mapping.values() atoms is the same
                            # within dc (DISTANCE_CUTOFF) of the distance between kdbAtom
                            # and mapping.keys() atoms.
                            for DA in list(mapping.keys()):
                                RA = mapping[DA]
                                pbcVector = self.atomAtomPbcVector(reactant, RA, reactantAtom)
                                if kdbDistances[DA] > micDist and PBC_MAPPING_CHECK:
                                    if not self.isDistance3(pbcVector, kdbDistances[DA], reactant.cell, dc):
                                        break
                                else:
                                    if abs(kdbDistances[DA] - self.atomAtomPbcDistance(reactant, RA, reactantAtom)) > dc:
                                        break
                            else:
                                newMapping = mapping.copy()
                                newMapping[kdbAtom] = reactantAtom
                                newMappings.append(newMapping)
                    mappings = newMappings


        debug and final_data.append(len(mappings))

        uniques = []  # will store list of unique suggestions
        dc2 = dc * dc

        # Load the mode and barriers for the process
        mode = copy.deepcopy(entry['mode'])
        if 'barrier' in entry:
            proc_barriers = copy.deepcopy(entry['barrier'])

        kdb_centroid = numpy.mean(kdbmin.positions, axis=0)
        kdb_shifted = kdbmin.positions - kdb_centroid
        kdb_shifted_T = kdb_shifted.transpose()
        del kdbmin

        kdbSaddle_shifted = entry['saddle'].positions - kdb_centroid
        kdbProduct_shifted = entry['product'].positions - kdb_centroid

        # check if entry is coplanar
        H = numpy.dot(kdb_shifted_T, kdb_shifted)
        _, S, _ = numpy.linalg.svd(H)  # numpy orders singular values from max to min
        normalized_min_sv = S[2] / len(kdb_shifted)
        if normalized_min_sv < dc2:  # arbitrary cutoff that seems to work very well
            is_coplanar = True
        else:
            is_coplanar = False

        # Loop over each mapping and try to find a rotation that aligns the
        # kdb configuration with the query configuration
        # this is done in the separate function below
        num_low_scores = 0
        num_overlap_pass = 0
        for mapping in mappings:
            output.append("mapping %s" % mapping)

            # current Python implementation guarantees that keys and values are called in the same order if dict is not modified in between
            kdb_atom_indices = list(mapping.keys())
            reactant_atom_indices = list(mapping.values())

            mapped_kdbmobile = [mapping[i] for i in kdbmobile]
            mapped_clump_order = [mapping[i] for i in clump_order]

            # clump gets rid of PBC and attempts to moves mapped reactant atoms into the same spatial orientation as kdb
            #reactantrot = fmodules.kdb.clump(reactant.positions, reactant_atom_indices, reactant.cell, reactant.ibox)
            reactantrot = fmodules.kdb.ordered_clump(reactant.positions, mapped_clump_order, reactant.cell, reactant.ibox)

            # SJ: Horn algorithm using SVD ("Least-Squares Fitting of Two 3-D Point Sets" by Arun, Huang, and Blostein)
            reactant_centroid = numpy.mean(reactantrot[reactant_atom_indices], axis=0)
            reactant_shifted = reactantrot[reactant_atom_indices] - reactant_centroid
            H = numpy.dot(kdb_shifted_T[:, kdb_atom_indices], reactant_shifted)  # including the indices ensures correct order
            U, _, V = numpy.linalg.svd(H)
            for mirror_check in range(1, -1*is_coplanar - 1, -2):  # loops twice only if is_coplanar
                V[2] *= mirror_check
                Rmat = numpy.dot(U, V)  # implicitly includes mirroring if necessary

                kdb_shifted_rotated = numpy.dot(kdb_shifted[kdb_atom_indices], Rmat)
                diff = kdb_shifted_rotated - reactant_shifted
                score = max(numpy.sum(diff*diff, axis=1))
                # print('\nscore:', score)
                # print('mapping:', mapping)
                output.append("Score: %s" % score)

                if score > dc2:
                    continue

                num_low_scores += 1

                # apply transformations on product
                kdbProduct = numpy.dot(kdbProduct_shifted, Rmat)
                kdbProduct += reactant_centroid

                # create suggested product and check overlap along the way
                overlap = False
                sugproduct = reactant.copy()  #XXX: what if the unit cell vectors change?
                sugproduct.radii = reactant.radii
                try:  # only works for graph method
                    atoms_to_check = set()
                    for atom in reactant_atom_indices:
                        neighbors = Reactant_Graph[atom]
                        atoms_to_check = atoms_to_check.union(neighbors)
                    atoms_to_check = atoms_to_check.difference(reactant_atom_indices)
                except TypeError:
                    atoms_to_check = set(range(len(sugproduct))).difference(reactant_atom_indices)
                if sugproduct.constraints:  # has one or more constrained atom(s)
                    for m in mapping:
                        if mapping[m] not in sugproduct.constraints[0].index:
                            if m in kdbmobile:
                                if self.check_overlap(kdbProduct[m], kdbmobile_radii[m], sugproduct, atoms_to_check, reactant.ibox):
                                    overlap = True
                                    break
                            sugproduct.positions[mapping[m]] = kdbProduct[m]
                else:  # no constraints
                    for m in mapping:
                        if m in kdbmobile:
                            if self.check_overlap(kdbProduct[m], kdbmobile_radii[m], sugproduct, atoms_to_check, reactant.ibox):
                                overlap = True
                                break
                        sugproduct.positions[mapping[m]] = kdbProduct[m]
                if overlap:
                    continue
                num_overlap_pass += 1

                # apply transformations on saddle
                kdbSaddle = numpy.dot(kdbSaddle_shifted, Rmat)
                kdbSaddle += reactant_centroid

                # Create the suggestion
                suggestion = reactant.copy()  #XXX: what if the unit cell vectors change?
                if suggestion.constraints:  # has one or more constrained atom(s)
                    for m in mapping:
                        if mapping[m] not in suggestion.constraints[0].index:
                            suggestion.positions[mapping[m]] = kdbSaddle[m]
                else:  # no constraints
                    suggestion.positions[reactant_atom_indices] = kdbSaddle[kdb_atom_indices]

                # Check for duplicates
                if nodupes == 1:
                    if self.check_duplicates_1(suggestion, sugproduct, uniques, dc, reactant.ibox):
                        continue
                elif nodupes == 2:
                    if self.check_duplicates_2(suggestion, sugproduct, mapped_kdbmobile, uniques, dc, reactant.ibox):
                        continue
                elif nodupes == 3:
                    if self.check_duplicates_3(suggestion, sugproduct, mapping, mapped_kdbmobile, uniques, dc, reactant.ibox):
                        continue

                #  Map the mode
                if mode is not None:
                    modeTemp = reactant.positions * 0.0
                    #for m in mapping:  # old way
                        #modeTemp[mapping[m]] = mode[m]
                    modeTemp[reactant_atom_indices] = mode[kdb_atom_indices]
                    try:
                        modeTemp /= numpy.linalg.norm(modeTemp)  # XXX: is this correct? (SJ)
                    except FloatingPointError:
                        mode = None

                # Perform the mode transformation
                if mode is not None:
                    modeTemp = numpy.dot(modeTemp, Rmat)

                # Rebox
                if REBOX_SUGGESTIONS:  #TODO: still need to fix
                    adjustment = 0.5 * (suggestion.cell[0] + suggestion.cell[1] + suggestion.cell[2])
                    temp = suggestion.positions - adjustment
                    suggestion.positions = fmodules.kdb.pbcs(temp, suggestion.cell, reactant.ibox) + adjustment
                    temp = sugproduct.positions - adjustment
                    sugproduct.positions = fmodules.kdb.pbcs(temp, suggestion.cell, reactant.ibox) + adjustment

                # Write suggestion
                if mode is not None:
                    if KEEP_BARRIERS:
                        self.output_query(outputdir, numMatches, suggestion, sugproduct, entry, modeTemp, processBarrier=proc_barriers, score=score)
                    else:
                        self.output_query(outputdir, numMatches, suggestion, sugproduct, entry, modeTemp, score=score)
                else:
                    if KEEP_BARRIERS:
                        self.output_query(outputdir, numMatches, suggestion, sugproduct, entry, processBarrier=proc_barriers, score=score)
                    else:
                        self.output_query(outputdir, numMatches, suggestion, sugproduct, entry, score=score)
                numMatches += 1

        if debug:
            final_data.append(num_low_scores)
            final_data.append(num_overlap_pass)
            final_data.append(numMatches)
            final_data.append(str(time.time() - start)[0:5])

        output.append("KDB matches: %s" % numMatches)
        return("\n".join(output), self.return_dict, final_data, numMatches)


    def query(self, reactant, outputdir = "./kdbmatches", nf=0.2, dc=0.3, nodupes = 0, kdbname = 'kdb.db', custom_config=None, pipe=None,local_query=False):
        # XXX: I think the best way forward to allow parallel processes
        # here is to make the query function return atoms objects instead
        # of writing them to file there.
        start_time = time.time()
        print()
        if custom_config is not None: # User submitted some custom config options
            print("Config:")
            print(custom_config)
            for option in custom_config:
                if option == 'DISTANCE_CUTOFF' or option == 'MOBILE_ATOM_CUTOFF' or option == 'NEIGHBOR_FUDGE':
                    continue # User cannot edit these in remote
                try:
                    exec(option + " = " + str(custom_config[option]), globals())
                except: # Invalid options arent used in code
                    pass
        else:
            print("Config: default")
        # Get the ibox to speed up pbcs
        ibox = numpy.linalg.inv(reactant.cell)
        reactant.ibox = ibox
        # Remove directory if kdbmatches is already there
        if os.path.isdir(outputdir):
            shutil.rmtree(outputdir)
        # Remake directory
        os.mkdir(outputdir)
        # Get a list of kdb entries that match the query configuration elementally
        entries, name = self.query_db(kdbname = kdbname, reactant = reactant)
        num_entries = len(entries)
        if num_entries == 0:
            print("No matching entries found. Exiting...")
            return

        # Create a list of element types and counts for the entire reactant
        # if Similarity turned on, it will convert to relvant dictionary
        reactantNameCount = self.similar_atom_nameCount(reactant,USE_GROUP_SIMILARITY)

        # initialized later as needed
        reactantNeighbors = {}
        Reactant_Graph = None
        GRAPH = False

        if USE_GRAPH:
            GRAPH = True
            try: # Try to import networkx, if it fails do not do USE_GRAPH
                from kdb.graph import GraphMaker
                import networkx as nx
                from kdb.isomorphism import GraphMatcher
            except ImportError:
                warnings.warn('failed to import Networkx, changing to default method', ImportWarning)
                GRAPH = False

        reactant.radii = NeighborList_helper.natural_cutoffs(reactant)

        if GRAPH:
            print("Using graph method")
            Reactant_Graph = GraphMaker.graph_reactant(self, reactant, nf)

        else:
            print("Using brute-force method")
            # For each nonfrozen atom in reactant, create a list of neighboring element
            # types and the count of each type.
            # TODO: this can be made N^2/2 trivially
            # TODO: this can use SAP for ortho boxes
            for i in range(len(reactant)):
                if reactant.constraints: #NK: if reactant does not have constraints
                    if i in reactant.constraints[0].index:
                        continue
                #r1 = elements[reactant.get_chemical_symbols()[i]]["radius"]
                r1 = reactant.radii[i]
                reactantNeighbors[i] = {}
                for j in range(len(reactant)):
                    if j == i:
                        continue
                    #r2 = elements[reactant.get_chemical_symbols()[j]]["radius"]
                    r2 = reactant.radii[i]
                    d = numpy.linalg.norm(fmodules.kdb.pbc(reactant.positions[i] - reactant.positions[j], reactant.cell, ibox))
                    if d > (r1 + r2) * (1 + nf):
                        continue
                    if reactant.get_chemical_symbols()[j] not in reactantNeighbors[i]:
                            reactantNeighbors[i][reactant.get_chemical_symbols()[j]] = 0
                    reactantNeighbors[i][reactant.get_chemical_symbols()[j]] += 1 
            reactantNeighbors=Kdb().Neighbor_Dict_similarity(reactantNeighbors,USE_GROUP_SIMILARITY)
        ###########################################################################
        # (Main) Loop over each kdb entry.
        ###########################################################################
        print("Querying %d entries..." % num_entries)
        cpuCount = mp.cpu_count()  # number of processors available (manually set to 1 to disable parallelization)
        with mp.Pool(cpuCount) as pool:
            dummy_dict = {}
            func = partial(self.query_kdb_entry, name, reactantNeighbors, reactantNameCount,
                           Reactant_Graph, reactant, outputdir, nf, dc , nodupes, kdbname, GRAPH,
                           )
            unique_id = 0
            matches_counter = 0
            debug and print("\nid, mobile, G1 nodes, G1 edges, G2 nodes, G2 edges, disconnected, mappings, low scores, no overlap, matches, time")
            query_timer = time.time()
            for result in pool.map(func, entries):
                debug and print(f"FINAL DATA {unique_id}:", result[2])
                unique_id += 1
                dummy_dict.update(result[1])
                matches_counter += result[3]  # increment by numMatches from each thread
            query_timer = time.time() - query_timer
        self.return_dict = dict([(key, value) for key, value in dummy_dict.items()])
        '''
        dummy_dict = {}
        unique_id = 0
        matches_counter = 0
        query_timer = time.time()
        for entry in entries:
            result = self.query_kdb_entry(name, reactantNeighbors, reactantNameCount, Reactant_Graph, reactant, outputdir, nf, dc, nodupes,
                                          kdbname, GRAPH, entry)
            debug and print(f"FINAL DATA {unique_id}:", result[2])
            unique_id += 1
            dummy_dict.update(result[1])
            matches_counter += result[3]  # increment by numMatches from each thread
        query_timer = time.time() - query_timer
        self.return_dict = dict([(key, value) for key, value in dummy_dict.items()])
        '''
        '''
        NOTE: If this part is used again, all the "kdbmatches" must be replaced with the variable "outputdir"

        # Rename files so that they appear from 0 to n
        if local_query and os.path.isdir(outputdir):
            onlyfiles = [f for f in listdir("kdbmatches") if isfile(join("kdbmatches", f))]
            SaddleFiles = []
            ProductFiles = []
            for file in onlyfiles:
                if "SADDLE" in file:
                    SaddleFiles.append(file)
                if "PRODUCT" in file:
                    ProductFiles.append(file)
             saddleCounter = 0
             productCounter = 0
             for sadFile in SaddleFiles:
                 for prodFile in ProductFiles:
                     if (sadFile.replace('SADDLE', '') == prodFile.replace('PRODUCT', '')):
                         os.rename((os.getcwd()+"/kdbmatches/" + sadFile), (os.getcwd()+"/kdbmatches/"+"SADDLE_%s" % saddleCounter))
                         os.rename((os.getcwd()+"/kdbmatches/" + prodFile), (os.getcwd()+"/kdbmatches/"+"PRODUCT_%s" % saddleCounter))
                 saddleCounter = saddleCounter + 1
             print("Number of kdb matches", saddleCounter)
             '''

        print("\nQuery Finished!")
        print("---------------------")
        print("%d matches (%.3f seconds)" % (matches_counter, (time.time() - start_time)) )
        print("Query time: %.3f seconds" % query_timer)
        print(f"Results written to {outputdir}\n")
        return (self.return_dict)
