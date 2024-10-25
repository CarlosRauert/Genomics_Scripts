import networkx as nx
import pandas as pd
import multiprocessing as mp
import sys
import numpy as np
import math

def cycle_overlap(cycle1, cycle2):
        """Calculate percentage overlap between two cycles based on shared nodes."""
        set1, set2 = set(cycle1), set(cycle2)
        intersection = set1.intersection(set2)
        return len(intersection) / min(len(set1), len(set2))

def is_valid_cycle(new_cycle, stored_cycles, overlap_threshold=0.9):
    """
    Check if the new cycle has less than 90% overlap or should replace 
    smaller cycles with high overlap.
    """
    for i, stored_cycle in enumerate(stored_cycles):
        overlap = cycle_overlap(new_cycle, stored_cycle)
        if overlap >= overlap_threshold:
            if len(new_cycle) > len(stored_cycle):
                stored_cycles[i] = new_cycle  # Replace with larger cycle
                return False  # Reject the new cycle after replacing
            return False  # Reject the new cycle if overlap is high and it's smaller
    return True  # If no high overlap, the cycle can be added

def optimized_simple_cycles(G, output_file, overlap_threshold=0.9):
    stored_cycles = []  # Only stores cycles that pass the conditions
    # Stream cycles one by one from nx.simple_cycles (generator)
    for cycle in nx.simple_cycles(G):
        if is_valid_cycle(cycle, stored_cycles, overlap_threshold):
            print(f"Cycle found: {cycle}")  # This will print each cycle
            #  Add cycle if it's valid
            stored_cycles.append(cycle)
    # After all cycles are generated, write them to the file
    with open(output_file, 'w') as f:
        for cycle in stored_cycles:
             f.write(f"{cycle}\n")  # Write each element followed by a newline

MatrixFile_l = ["/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1KU/A1KU_adj_Matrix_12.csv",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1KW/A1KW_adj_Matrix_38.csv",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1L0/A1L0_adj_Matrix_36.csv",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1L2/A1L2_adj_Matrix_20.csv",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1L3/A1L3_adj_Matrix_48.csv",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A240/A240_adj_Matrix_15.csv",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A2IZ/A2IZ_adj_Matrix_20.csv",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A2J0/A2J0_adj_Matrix_14.csv",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A23R/A23R_adj_Matrix_55.csv",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A2J4/A2J4_adj_Matrix_13.csv",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A2QS/A2QS_adj_Matrix_38.csv",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A3LS/A3LS_adj_Matrix_29.csv",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A3LT/A3LT_adj_Matrix_125.csv",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A3LW/A3LW_adj_Matrix_8.csv",
                "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A3LY/A3LY_adj_Matrix_14.csv",]


CaseID_l = ["A1KU", "A1KW", "A1L0", "A1L2", "A1L3", "A240", "A2IZ", "A2J0", "A23R", "A2J4", "A2QS", "A3LS", "A3LT", "A3LW", "A3LY"]
CNxt_l = ["12", "38", "36", "20", "48", "15", "20", "14", "55", "13", "38", "29", "125", "8", "14"]

i=1

for i in range(len(MatrixFile_l)):
    adj_matrix_directed = pd.read_csv(MatrixFile_l[i])
    adj_matrix_directed = adj_matrix_directed.drop(["Unnamed: 0"],axis=1)
    adj_matrix_directed.index=adj_matrix_directed.columns
    CaseID = CaseID_l[i]
    CNxt = CNxt_l[i]
    filename=f"/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/{CaseID}/optimized_cycles_{CNxt}.txt"
    print("Creating NX Graph and Basis")
    G = nx.from_pandas_adjacency(adj_matrix_directed, create_using=nx.DiGraph)
    H = G.to_undirected()
    cycle_basis=nx.cycle_basis(H)
    cycle_nodes = set(node for cycle in cycle_basis for node in cycle)
    subgraph = G.subgraph(cycle_nodes)
    optimized_simple_cycles(subgraph,filename,0,7)