import pandas as pd
import numpy as np
from tqdm import tqdm

data = pd.read_csv("ncac_filtered_data.csv")

def calculate_dihedral_angle(p1, p2, p3, p4): #Atoms A, B, C, D
    b1 = p2 - p1 # AB Vector = B - A
    b2 = p3 - p2 # BC Vector = C - B
    b3 = p4 - p3 # CD Vector = D - C
    n1 = np.cross(b1, b2) # N1 = AB x BC (plane 1)
    n1_normalized = n1/np.linalg.norm(n1) # N1 = N1/|N1|
    n2 = np.cross(b2, b3) # N2 = BC x CD
    n2_normalized = n2/np.linalg.norm(n2) # N2 = N2/|N2|
    cos_theta = np.dot(n1_normalized, n2_normalized) # x = cos(t) = N1 . N2
    sin_theta = np.dot(b2, np.cross(n1_normalized, n2_normalized))/np.linalg.norm(b2) # y = sin(t) = (BC . (N1 x N2))/|BC|
    return np.degrees(np.arctan2(sin_theta, cos_theta)) # Dihedral angle = tan-1(sin(t) cos(t))

angles = []

for i in tqdm(range(0, len(data) - 2, 3), total=(len(data) - 2) // 3): # Sliding through the center atoms C_pre-|N-CA-C|-N_next
    if (i + 2 < len(data)) and (data.loc[i, ["Protein ID", "Chain ID"]].tolist() == data.loc[i + 2, ["Protein ID", "Chain ID"]].tolist()):
        
        # Extracting coordinates of center N - phi - CA - psi - N
        N = np.array([data["X"].iloc[i], data["Y"].iloc[i], data["Z"].iloc[i]])
        CA = np.array([data["X"].iloc[i + 1], data["Y"].iloc[i + 1], data["Z"].iloc[i + 1]])
        C = np.array([data["X"].iloc[i + 2], data["Y"].iloc[i + 2], data["Z"].iloc[i + 2]])

        phi = None
        psi = None

        # Calculate Phi if C_prev exists
        if i - 1 >= 0:  # Checking if previous residue exists
            C_prev = np.array([data["X"].iloc[i - 1], data["Y"].iloc[i - 1], data["Z"].iloc[i - 1]])
            phi = calculate_dihedral_angle(C_prev, N, CA, C)

        # Calculate Psi if N_next exists
        if i + 3 < len(data): 
            N_next = np.array([data["X"].iloc[i + 3], data["Y"].iloc[i + 3], data["Z"].iloc[i + 3]])
            psi = calculate_dihedral_angle(N, CA, C, N_next)
        angles.append([data["Protein ID"].iloc[i], data["Chain ID"].iloc[i], data["Residue ID"].iloc[i], data["Residue Name"].iloc[i], phi, psi])

residue_phi_psi_angle_df = pd.DataFrame(angles, columns=["Protein ID", "Chain ID", "Residue ID", "Residue", "Phi", "Psi"])
residue_phi_psi_angle_df .to_csv("residue_wise_phi_psi_angles.csv", index=False)
