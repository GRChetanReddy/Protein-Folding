import pandas as pd  
data = pd.read_csv("protein_residues.csv")
required_atoms = {"N", "CA", "C"}  
filtered_data = data[data["Atom Name"].isin(required_atoms)].copy()
filtered_data.to_csv("ncac_filtered_data.csv", index=False)