import csv
import os
from Bio.PDB import PDBParser, NeighborSearch, ShrakeRupley



with open('out_with_sasa.csv', 'a', newline='') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(["GENE","ACC_ID","MOD_RSD", "SASA"])

parser = PDBParser(QUIET=True)

with open('holder_csv.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        if row[0] != 'GENE':
            structure = parser.get_structure('struct', f'pdb_files/{row[1]}.pdb') 
            # ns = NeighborSearch(list(structure.get_atoms()))
            sr = ShrakeRupley(probe_radius=1.4, n_points=100)
            sr.compute(structure, level='A')
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.get_id()[1] == int(row[2]):
                            print(residue.get_resname())
                            for atom in residue:
                                if atom.get_name() == 'NZ':
                                    with open('out_with_sasa.csv', 'a', newline='') as outfile:
                                        writer = csv.writer(outfile)
                                        writer.writerow(row + [atom.sasa])

