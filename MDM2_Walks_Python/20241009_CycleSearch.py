import networkx as nx
import pandas as pd
import multiprocessing as mp
import sys
import numpy as np
import math

print("Reading file ", sys.argv[1])
MatrixFile = sys.argv[1]
CaseID = sys.argv[2]
CNxt = sys.argv[3]

#MatrixFile = "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/A1KU/A1KU_adj_Matrix_8.csv"
#CaseID = "A1KU"
#CNxt = 8

adj_matrix_directed = pd.read_csv(MatrixFile)
adj_matrix_directed = adj_matrix_directed.drop(["Unnamed: 0"],axis=1)
adj_matrix_directed.index=adj_matrix_directed.columns

print("Creating NX Graph and Basis")
G = nx.from_pandas_adjacency(adj_matrix_directed, create_using=nx.DiGraph)
H = G.to_undirected()
cycle_basis=nx.cycle_basis(H)
cycle_nodes = set(node for cycle in cycle_basis for node in cycle)
subgraph = G.subgraph(cycle_nodes)

simple_cycles = list(nx.simple_cycles(subgraph))

print("searching for cycles")
Minimum_Cycle_Basis = list(nx.minimum_cycle_basis(G))
Cycles= list(nx.simple_cycles(G))
print("\nElementary Cycles in the Graph:")
print(len(Cycles))

if len(Cycles)>0:
    print("found Cycle")
    #Save list to txt file
    #Open a file in write mode
    print("writing file")
    with open("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/"+CaseID+"/"+CaseID+"_"+CNxt+"_Cycles.txt", "w") as file:
        for item in Cycles:
            file.write(f"{item}\n")  # Write each element followed by a newline
else:
    print("No Cycle found")
    with open("/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/"+CaseID+"/"+CaseID+"_"+CNxt+"_Cycles.txt", "w") as file:
            file.write(f"No_Cycle_Found")