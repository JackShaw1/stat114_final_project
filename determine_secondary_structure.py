import csv
import os
from Bio.PDB import PDBParser, NeighborSearch




parser = PDBParser()

structure = parser.get_structure('struct', 'BRCA2_AF.pdb')

def determine_ss(structure, residue_index):
    ns = NeighborSearch(list(structure.get_atoms()))
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[1] == residue_index:
                    for atom in residue:
                        if atom.get_name() == 'CA':
                            counter = 0
                            neighbors = ns.search(atom.get_coord(), 6)
                            taken = []
                            for atom2 in neighbors:
                                if atom2.get_parent().get_id()[1] not in taken and atom2.get_name() == 'CA' and atom2.get_parent().get_id()[1] != atom.get_parent().get_id()[1]:
                                    counter += 1
                            print(counter)
        break


determine_ss(structure, 79)