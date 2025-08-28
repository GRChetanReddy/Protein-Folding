import os
import pandas as pd
import openpyxl
from Bio.PDB import MMCIFParser, PDBParser
from tqdm import tqdm

def parse_pdb_file(file_path, protein_id, chain_id):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(protein_id, file_path)
    atom_data = []
    
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    res_id = residue.id[1] 
                    res_name = residue.resname 
                    for atom in residue:
                        atom_name = atom.get_name()  
                        x, y, z = atom.coord
                        atom_data.append([protein_id, chain_id, res_id, res_name, atom_name, x, y, z])
    return atom_data

def parse_cif_file(file_path, protein_id, chain_id):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(protein_id, file_path)
    atom_data = []
    
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    res_id = residue.id[1]
                    res_name = residue.resname
                    for atom in residue:
                        atom_name = atom.get_name()
                        x, y, z = atom.coord
                        atom_data.append([protein_id, chain_id, res_id, res_name, atom_name, x, y, z])
    return atom_data

df = pd.read_excel("prot_id_vs_sequence.xlsx")
extracted_data = []

for index, row in tqdm(df.iterrows(), total = df.shape[0]):
    
    prot_chain = row["Protein_Chain_ID"]
    protein_id, chain_id = prot_chain[:4], prot_chain[4:]
    pdb_file_path = f"./PDB_Files/{protein_id}.pdb"
    mmcif_file_path = f"./mmCIF_Files/{protein_id}.cif"
    if os.path.exists(pdb_file_path):
        extracted_data.extend(parse_pdb_file(pdb_file_path, protein_id, chain_id))
    elif os.path.exists(mmcif_file_path):
        extracted_data.extend(parse_cif_file(mmcif_file_path, protein_id, chain_id))
    else:
        print(f"File not found for {prot_chain}")
        
columns = ["Protein ID", "Chain ID", "Residue ID", "Residue Name", "Atom Name", "X", "Y", "Z"]
df_atoms = pd.DataFrame(extracted_data, columns=columns)
df_atoms.to_csv("protein_residues.csv", index=False)
print("Data extraction complete. Saved as protein_residues.csv")

