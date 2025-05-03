import csv
import os
from Bio.PDB import PDBParser, NeighborSearch, ShrakeRupley



with open('out_with_contacts.csv', 'a', newline='') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(["GENE","ACC_ID","MOD_RSD", "NUM_POS", "NUM_NEG"])

parser = PDBParser(QUIET=True)

with open('holder_csv.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        if row[0] != 'GENE':
            if f'{row[1]}.pdb' in os.listdir('pdb_files'):
                structure = parser.get_structure('struct', f'pdb_files/{row[1]}.pdb') 
                ns = NeighborSearch(list(structure.get_atoms()))
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            if residue.get_id()[1] == int(row[2]):
                                print(residue.get_resname())
                                for atom in residue:
                                    if atom.get_name() == 'NZ':
                                        positive_counter = 0
                                        negative_counter = 0
                                        neighbors = ns.search(atom.get_coord(), 10.0)
                                        taken_res = []
                                        for atom2 in neighbors:
                                            if atom2.get_parent().get_id()[1] not in taken_res and atom2.get_parent().get_resname() in ['LYS', 'ARG']:
                                                taken_res.append(atom2.get_parent().get_id()[1])
                                                positive_counter += 1
                                        taken_res = []
                                        for atom2 in neighbors:
                                            if atom2.get_parent().get_id()[1] not in taken_res and atom2.get_parent().get_resname() in ['ASP', 'GLU']:
                                                taken_res.append(atom2.get_parent().get_id()[1])
                                                negative_counter += 1
                                        with open('out_with_contacts.csv', 'a', newline='') as outfile:
                                            writer = csv.writer(outfile)
                                            writer.writerow(row + [positive_counter, negative_counter])

