import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm

data = pd.read_csv("residue_wise_phi_psi_angles.csv")
data = data[1:len(data)-1] #Eliminating None values

#Only 20 default amino acid residues
standard_residues = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

heatmap_dir = "2D_Heatmaps"
plot3d_dir = "3D_Plots"
os.makedirs(heatmap_dir, exist_ok=True)
os.makedirs(plot3d_dir, exist_ok=True)

binsize = 3
bin_edges = np.arange(-180, 181, binsize)
num_bins = len(bin_edges) - 1

for residue_name in tqdm(standard_residues, total=len(standard_residues)):
    residue_data = data[data["Residue"] == residue_name]
    density_matrix = np.zeros((num_bins, num_bins))
    for phi, psi in zip(residue_data["Phi"], residue_data["Psi"]):
        phi_idx = int((phi + 180) // binsize)
        psi_idx = int((psi + 180) // binsize)
        density_matrix[psi_idx, phi_idx] += 1 #Increasing the count in the bins
    density_matrix = density_matrix[::-1]

    # 2D Heatmap
    plt.figure(figsize=(8, 6))
    plt.imshow(density_matrix, origin="lower", extent=[-180, 180, -180, 180], cmap="viridis_r")
    plt.colorbar(label="Count")
    plt.xlabel("Phi (°)")
    plt.ylabel("Psi (°)")
    plt.title(f"2D Ramachandran Heatmap for {residue_name}")
    plt.xticks(np.arange(-180, 181, 30))
    plt.yticks(np.arange(-180, 181, 30))
    
    heatmap_path = os.path.join(heatmap_dir, f"{residue_name}_heatmap.png")
    plt.savefig(heatmap_path, dpi=300)
    plt.close()

    #3D plot
    X, Y = np.meshgrid(bin_edges[:-1], bin_edges[:-1])
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(X, Y, density_matrix, cmap="viridis")
    ax.set_xlabel("Phi (°)")
    ax.set_ylabel("Psi (°)")
    ax.set_title(f"3D Ramachandran Density Plot for {residue_name}")
    fig.colorbar(surf, ax=ax, aspect=10, label="Count")
    ax.set_xticks(np.arange(-180, 181, 30))
    ax.set_yticks(np.arange(-180, 181, 30))

    plot3d_path = os.path.join(plot3d_dir, f"{residue_name}_3D_plot.png")
    plt.savefig(plot3d_path, dpi=300)
    plt.close()
