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
            # print(f"Cycle found: {cycle}")  # This will print each cycle
            # Add cycle if it's valid
            stored_cycles.append(cycle)
    # After all cycles are generated, write them to the file
    with open(output_file, 'w') as f:
        for cycle in stored_cycles:
             f.write(f"{cycle}\n")  # Write each element followed by a newline

CaseID_l = ["CPCT02010386T", "CPCT02010680T", "CPCT02060104T", "CPCT02060191T", "CPCT02070051T", "CPCT02070366T", "CPCT02080206T", "CPCT02080227T", "CPCT02090057T", "CPCT02340046T"]
CNxt_l = ["91", "11", "36", "11", "9", "9", "9", "10", "18", "13", "10", "91"]

i=1

for i in range(len(CaseID_l)):
    CaseID = CaseID_l[i]
    CNxt = CNxt_l[i]
    MatrixFile=f"/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/{CaseID}/{CaseID}_adj_Matrix_{CNxt}.csv"
    adj_matrix_directed = pd.read_csv(MatrixFile)
    adj_matrix_directed = adj_matrix_directed.drop(["Unnamed: 0"],axis=1)
    adj_matrix_directed.index=adj_matrix_directed.columns
    filename=f"/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/{CaseID}/optimized_cycles_{CNxt}.txt"
    print(CaseID)
    print("Creating NX Graph and Basis")
    G = nx.from_pandas_adjacency(adj_matrix_directed, create_using=nx.DiGraph)
    H = G.to_undirected()
    cycle_basis=nx.cycle_basis(H)
    cycle_nodes = set(node for cycle in cycle_basis for node in cycle)
    subgraph = G.subgraph(cycle_nodes)
    optimized_simple_cycles(subgraph,filename,0.7)

    adj_mae_directed="I love you baby"
    print(adj_mae_directed)    