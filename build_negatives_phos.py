from Bio.PDB import PDBParser, NeighborSearch
import os
import csv

with open('negatives_phos.csv', 'a', newline='') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(["GENE","ACC_ID","UNMOD_RSD"])

parser = PDBParser()

data_dict = {}
with open('holder_csv.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        if row[0] != 'GENE':
            if row[1] not in data_dict:
                data_dict[row[1]] = [row[0], []]


with open('holder_csv.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        if row[0] != 'GENE':
            data_dict[row[1]][1].append(int(row[2]))

counter = 0
if counter < 1500:
    for key, value in data_dict.items():
        structure = parser.get_structure('struct', f'pdb_files/{key}.pdb')
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() in ['SER', 'THR', 'TYR'] and residue.get_id()[1] not in value[1]:
                        counter += 1
                        with open('negatives_phos.csv', 'a', newline='') as outfile:
                            writer = csv.writer(outfile)
                            writer.writerow([value[0], key, residue.get_id()[1]])