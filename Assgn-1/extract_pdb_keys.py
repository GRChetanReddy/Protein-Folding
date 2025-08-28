import pandas as pd 

with open("cullpdb_pc25.0_res0.0-2.0_len40-1000_R0.25_Xray+Nmr_d2025_02_12_chains9832.fasta", "r") as f:
    lines = f.readlines()
    sequence_dict = {}
    curr_key = ""
    for line in lines:
        if line.startswith(">"):
            curr_key = line.split(" ")[0][1:]
            sequence_dict[curr_key] = ""
        else:
            line = line.strip()
            sequence_dict[curr_key] += line

data = pd.DataFrame(list(sequence_dict.items()), columns=["Protein_Chain_ID", "Sequence"])
data.to_excel("prot_id_vs_sequence.xlsx", index=False)