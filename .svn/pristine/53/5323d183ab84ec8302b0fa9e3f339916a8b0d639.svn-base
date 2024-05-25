import glob
import os
import numpy
import sys
import re
import subprocess
import kdb.fmodules as fmodules
try:
    from .aselite import elements
    from .aselite import elements_groups
except ImportError:
    try:
        from .aselite import elements
        from .aselite import elements_groups
    except ImportError:
        from aselite import elements
        from aselite import elements_groups
except SystemError: #website related - www user has no kdb in its pypath
    from aselite import elements
    from aselite import elements_groups


class Kdb():
    def check_svn_version(self):
        #does not work
        svn_info = subprocess.check_output("svn info")
        print(svn_info)

    def check_version(self):
        pyversion = sys.version.split()[0].split(".")  # ex: 3.10.2 --> ['3', '10', '2']
        if int(pyversion[0]) < 3 or int(pyversion[1]) < 6:  # python 3.6 or greater is required
            return False
        return True


    def atomAtomPbcVector(self, atoms, a, b):
        if not hasattr(atoms, 'ibox'):
            atoms.ibox = numpy.linalg.inv(atoms.get_cell())
        if not hasattr(atoms, 'pbcVectors'):
            atoms.pbcVectors = {}
        if (a, b) not in atoms.pbcVectors or (b, a) not in atoms.pbcVectors:
            atoms.pbcVectors[(a, b)] = fmodules.kdb.pbc(atoms.positions[b] - atoms.positions[a], atoms.get_cell(), atoms.ibox)
            atoms.pbcVectors[(b, a)] = -atoms.pbcVectors[(a, b)]
        return atoms.pbcVectors[(a, b)]


    def atomAtomPbcDistance(self, atoms, a, b):
        if not hasattr(atoms, 'pbcDistances'):
            atoms.pbcDistances = {}
        if (a, b) not in atoms.pbcDistances or (b, a) not in atoms.pbcDistances:
            atoms.pbcDistances[(a, b)] = numpy.linalg.norm(self.atomAtomPbcVector(atoms, a, b))
            atoms.pbcDistances[(b, a)] = atoms.pbcDistances[(a, b)]
        return atoms.pbcDistances[(a, b)]


    def atomAtomDistance(self, atoms, a, b):
        if not hasattr(atoms, 'distances'):
            atoms.distances = {}
        if (a, b) not in atoms.distances or (b, a) not in atoms.distances:
            atoms.distances[(a, b)] = fmodules.kdb.atom_atom_distance(atoms.positions[a] - atoms.positions[b])
            atoms.distances[(b, a)] = atoms.distances[(a, b)]
        return atoms.distances[(a, b)]


    def atomAtomPseudoMicDistance(self, atoms, a, b):
        """
        Returns the distance according to the Minimum Image Convention
        (well, almost. See below:)
            * for speed, calculates the PBC vector first, and then permute
              in the 8 adjacent unit cells. Equivalent to the true
              micDistance for basically all unit cells
            * for VERY extremely skewed unit cells, this must be done
              iteratively until true minimum is observed. This is not
              implemented here
            
        For unit cells without any non-orthogonal angles, equivalent to
        atomAtomPbcDistance()

        Parameters
        -------------------
        atoms : Atoms object
            structure to compute MIC distance

        a : int
            first atom in atoms

        b : int
            second atom in atoms

        """

        if not hasattr(atoms, 'ibox'):
            atoms.ibox = numpy.linalg.inv(atoms.get_cell())
        if not hasattr(atoms, 'micDistances'):
            atoms.micDistances = {}

        if (a, b) not in atoms.micDistances or (b, a) not in atoms.micDistances:
            v = atoms.positions[b] - atoms.positions[a]
            dmic = fmodules.kdb.atom_atom_pseudo_mic_distance(v, atoms.cell, atoms.ibox)
            atoms.micDistances[(a, b)] = dmic
            atoms.micDistances[(b, a)] = dmic
            return dmic
        else:
            return atoms.micDistances[(a, b)]


    def ortho(self, atoms):
        """
        Returns list about the orthogonality of the atoms cell.

        atoms.ortho = [ * , * , * , * ]
        atoms.ortho[i] == 1 if ith cell basis vector is orthogonal to the other 2 (for i = 0,1,2)
        atoms.ortho[3] == 1 if reactant.ortho[i] == 1 for i = 0,1,2
        """
        # XXX: obsolete, as of revision 278, but kept here just in case
        # atoms.ortho[i] == 1 if ith cell basis vector is orthogonal to the other 2 (for i = 0,1,2)
        # atoms.ortho[3] == 1 if reactant.ortho[i] == 1 for i = 0,1,2
        cell = atoms.cell
        xx = numpy.dot(cell[0], cell[0])
        yy = numpy.dot(cell[1], cell[1])
        zz = numpy.dot(cell[2], cell[2])
        xy2 = numpy.dot(cell[0], cell[1]) ** 2
        xz2 = numpy.dot(cell[0], cell[2]) ** 2
        yz2 = numpy.dot(cell[1], cell[2]) ** 2
        xyortho = ( (2048 * xy2) < (xx * yy) )  # 2048 is arbitrary (intended to make cutoff around 89 deg)
        xzortho = ( (2048 * xz2) < (xx * zz) )
        yzortho = ( (2048 * yz2) < (yy * zz) )
        atoms.ortho = [xyortho and xzortho, xyortho and yzortho, xzortho and yzortho]
        atoms.ortho.append(all(atoms.ortho))
        return atoms.ortho
       

    def atomAtomPermDistance(self, atoms, a, b):
        """
        Returns list of distance between atoms a and b in atoms permuted in
        the 27 most-adjacent PBC cells
        """
        if not hasattr(atoms, 'ibox'):
            atoms.ibox = numpy.linalg.inv(atoms.get_cell())
        if not hasattr(atoms, 'permDistances'):
            atoms.permDistances = {}
        if (a, b) not in atoms.permDistances or (b, a) not in atoms.permDistances:
            v = atoms.positions[b] - atoms.positions[a]
            atoms.permDistances[(a, b)] = fmodules.kdb.atom_atom_perm_distance(v, atoms.cell, atoms.ibox)
            atoms.permDistances[(b, a)] = atoms.permDistances[(a, b)]
        return atoms.permDistances[(a, b)]


    def getNameList(self, atoms):
        """
        Returns a sorted list of element names.
        """
        nl = []
        for name in atoms.get_chemical_symbols():
            if name not in nl:
                nl.append(name)
        return sorted(nl)


    def nameCount(self, atoms):
        counts = {}
        for name in atoms.get_chemical_symbols():
            if not name in counts:
                counts[name] = 0
            counts[name] += 1
        return counts


    def pbc(self, r, box, ibox = None):
        """
        Applies periodic boundary conditions.
        Parameters:
            r:      the vector the boundary conditions are applied to
            box:    the box that defines the boundary conditions
            ibox:   the inverse of the box. This will be calculuated if not provided.
        """
        #if ibox == None:
        #if not hasattr(ibox, 'shape'):
        if type(ibox) != numpy.ndarray and type(ibox) != list and type(ibox) != tuple: #MJW fix
            ibox = numpy.linalg.inv(box)
        vdir = numpy.dot(r, ibox)
        vdir = (vdir % 1.0 + 1.5) % 1.0 - 0.5
        return numpy.dot(vdir, box)


    def per_atom_norm(self, v, box, ibox = None):
        '''
        Returns a length N numpy array containing per atom distance
            v:      an Nx3 numpy array
            box:    box matrix that defines the boundary conditions
            ibox:   the inverse of the box. will be calculated if not provided
        '''
        if type(ibox) != numpy.ndarray and type(ibox) != list and type(ibox) != tuple:
            ibox = numpy.linalg.inv(box)
        if numpy.squeeze(v).ndim == 1:
            pbc = fmodules.kdb.pbc
        else:
            pbc = fmodules.kdb.pbcs
        diff = pbc(v, box, ibox)
        return numpy.array([numpy.linalg.norm(d) for d in diff])


    def load_mode(self, modefilein):
        ''' 
        Reads a mode.dat file into an N by 3 numpy array
            modefilein: filename
        '''
        f = open(modefilein, 'r')
        lines = f.readlines()
        f.close()
        mode = []
        line_counter = 0
        for line in lines:
            line_counter += 1
            l = line.strip().split()
            try:
                assert len(l) == 3
            except:
                sys.exit(f"ERROR: line #{line_counter} of the mode file is not formatted correctly")
            for j in range(3):
                mode.append(float(l[j]))
        mode = numpy.array(mode)
        mode.resize(int(len(mode)/3), 3)
        return mode


    def save_mode(self, modefileout, displace_vector):
        '''
        Saves an Nx3 numpy array into a mode.dat file. 
            modefileout:     filename
            displace_vector: the mode (Nx3 numpy array)
        '''
        f = open(modefileout, 'w')
        for i in range(len(displace_vector)):
            f.write("%.3f %.3f %.3f\n" % (displace_vector[i][0],
                displace_vector[i][1], displace_vector[i][2]))

    def list_element_combinations(self, kdbdir):
        combinations = [os.path.basename(i) for i in glob.glob(os.path.join(kdbdir, "*"))]
        return combinations

    def combo_split(self, combo):
        elements = []
        for i in range(len(combo)):
            if combo[i] == combo[i].lower():
                elements[-1] += combo[i]
            else:
                elements.append(combo[i])
        return elements

    def is_symbol_subset(self, a, b):
        for symbol in a:
            if symbol not in b:
                return False
        return True

    def query_has_all(self, kdbdir, symbols):
        result = []
        combinations = [os.path.basename(i) for i in glob.glob(os.path.join(kdbdir, "*"))]
        for combo in combinations:
            elements = self.combo_split(combo)
            if not is_symbol_subset(symbols, elements):
                continue
            for N in glob.glob(os.path.join(kdbdir, combo, '*')):
                result.append(N)
        return result

    def check_kdb_combinations(self, name, kdbInputList):
        matchList = []
        for k in range(0, len(kdbInputList)):
            kdbinput_in_use = str(kdbInputList[k])
            inputList = re.findall('[A-Z][a-z]*', str(name))
            #print(kdbinput_in_use)
            count = 0
            for i in range(0,len(inputList)):
                kdbSeperationList = re.findall('[A-Z][a-z]*', str(kdbinput_in_use))
                inputList = re.findall('[A-Z][a-z]*', str(name))
                inputElement = inputList[i]
                for j in range(0,len(kdbSeperationList)):
                    if inputElement == kdbSeperationList[j]:
                        count = count + 1;
                        #print(count)
                if count == len(kdbSeperationList):
                    matchList.append(kdbinput_in_use)
        return matchList

    def string_split(self,system):
        return re.findall('[A-Z][a-z]*', str(system))
        #print string_split("LiFeO")
        #print elements['Li']['prime_number']

    def name_to_code(self,name):
        name_array = self.string_split(name)
        code =1
        for element in name_array:
            code *= elements[element]['prime_number']
        return code
        #print name_to_code("LiHO")

    def name_to_type_code(self,name):
        name_array = self.string_split(name)
        unique_type_list = []
        code = 1
        for element in name_array: #for each element in process
            type_element = elements[element]['elemental_type'] #i.e AM,AEM,NG etc
            if type_element in unique_type_list:
                continue
            else:
                unique_type_list.append(type_element)
        for type_element in unique_type_list:
            code *= elements_groups[type_element]['prime_number'] #multiply the prime num of that type
            similar_atom_code = code
        return similar_atom_code

    def name_to_type(self,name):
        type = elements[name]['elemental_type']
        return type

    def similar_atom_nameCount(self, atoms,similarity):
        if similarity:
            counts = {}
            for name in atoms.get_chemical_symbols():
                group_name = self.name_to_type(name)
                if not group_name in counts:
                    counts[group_name] = 0
                counts[group_name] += 1 
        else:
            return self.nameCount(atoms)
        return counts

    def Neighbor_Dict_similarity(self,kdbNeighbors,similarity):
        if similarity:
          similiar_dict = {}
          for key in kdbNeighbors:
            similiar_dict[key]={}
            for nkey in kdbNeighbors[key]:
                if Kdb().name_to_type(nkey) not in similiar_dict[key]:# need this to sum groups together,  1: {'Al': 5,'B': 6}=1:{'G13':11}
                    similiar_dict[key][Kdb().name_to_type(nkey)]=kdbNeighbors[key][nkey]
                else:
                    similiar_dict[key][Kdb().name_to_type(nkey)]+=kdbNeighbors[key][nkey]
        else:
          return kdbNeighbors
        return similiar_dict

    #define clump function here for both kdbquery and kdbinsert
    #the clump function based on pairs is called clump_new while the clump function with mobile atom is clump_old
    def clump_new(self, c, atoms, ibox=None):
        temp = c.copy()
        undone = atoms[:]
        done = [undone.pop(0)]
        while len(undone) > 0:
            mindist = 1e10
            for i in undone[:]:
                for j in done[:]:
                    v = self.pbc(temp.positions[i] - temp.positions[j], temp.cell, ibox)
                    dist = numpy.linalg.norm(v)
                    if dist < mindist:
                        vmin = v
                        mindist = dist
                        min_i = i
                        min_j = j
            temp.positions[min_i] = temp.positions[min_j] + vmin
            done.append(min_i)
            undone.remove(min_i)
        return temp

    def clump_old(self, c, atoms):
        temp = c.copy()
        undone = atoms[:]
        working = [undone.pop(0)] #using the first atom to remove pbc as it is the mobile atom
        while len(undone) > 0:
            if len(working) == 0:
                return
            a = working.pop()
            for i in undone[:]:
                v = self.pbc(temp.positions[i] - temp.positions[a], temp.cell)
                temp.positions[i] = temp.positions[a] + v
                working.append(i)
                undone.remove(i)
        return temp

    #old clump fuction of kdbinsert that pairs with clump_old of kdbquery
    def clump_insert_old(self, r, s, p):
        temp = r.copy()
        undone = list(range(len(temp)))
        working = [undone.pop(0)]
        while len(undone) > 0:
            if len(working) == 0:
                return 0
            a = working.pop()
            for i in undone[:]:
                vr = self.pbc(r.positions[i] - r.positions[a], r.get_cell())
                dr = numpy.linalg.norm(vr)
                vs = self.pbc(s.positions[i] - s.positions[a], s.get_cell())
                ds = numpy.linalg.norm(vs)
                vp = self.pbc(p.positions[i] - p.positions[a], p.get_cell())
                dp = numpy.linalg.norm(vp)
                temp[i].position = temp[a].position + vr
                working.append(i)
                undone.remove(i)
        return temp

    def convertClumpOrder(self, clumpOrder):
        """ Used to convert clump order from storage format (string) to insert/query format
            (list) and vice versa.
            See fortran/common.f90 for clump methods
        """ 
        if isinstance(clumpOrder, str):
            return [int(i) for i in clumpOrder.split()]
        elif isinstance(clumpOrder, list) or isinstance(clumpOrder, numpy.ndarray):
            return "".join(str(int(i)) + " " for i in clumpOrder)
