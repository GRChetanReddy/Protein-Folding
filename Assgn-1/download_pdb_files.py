import os
import pandas as pd
import requests

protein_data = pd.read_excel('prot_id_vs_sequence.xlsx')
pdb_chain_ids = protein_data["Protein_Chain_ID"]

missing_pdb_ids = []
count = 0
with open("download_results.txt", "a+") as log_file: 
    log_file.seek(0)
    lines = log_file.read().splitlines()  
    for pdb_chain_id in pdb_chain_ids:
        pdb_id = pdb_chain_id[:4]
        pdb_path = f"./PDB_Files/{pdb_id}.pdb"
        if os.path.exists(pdb_path):
            continue
        if any(f"{pdb_id}.pdb" in line for line in lines):
            missing_pdb_ids.append(pdb_id)
            continue 
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(url)
        if response.status_code == 200:
            with open(pdb_path, "wb") as f:
                f.write(response.content)
            print(f"Download Complete: {pdb_id}.pdb")
        else:
            print(f"Download Failed: {pdb_id}.pdb")
            log_file.write(f"Download Failed: {pdb_id}.pdb\n")
        log_file.flush()

# To download the mmCIF files 
for pdb_id in missing_pdb_ids:
    cif_path = f"./mmCIF_Files/{pdb_id}.cif"
    if os.path.exists(cif_path):
        print(f"Already exists: {cif_path}")
        continue
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    response = requests.get(url)
    if response.status_code == 200:
        with open(cif_path, "wb") as f:
            f.write(response.content)
        print(f"Download Complete: {pdb_id}.cif")
    else:
        print(f"Download Failed: {pdb_id}.cif")
