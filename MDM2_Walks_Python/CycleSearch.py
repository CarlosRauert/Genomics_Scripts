import networkx as nx
import pandas as pd
import multiprocessing as mp
import sys

print("Reading file ", sys.argv[1])
MatrixFile = sys.argv[1]
CaseID = sys.argv[2]
CNxt = sys.argv[3]

adj_matrix_directed = pd.read_csv(MatrixFile)
adj_matrix_directed = adj_matrix_directed.drop(["Unnamed: 0"],axis=1)
adj_matrix_directed.index=adj_matrix_directed.columns

print("Creating NX Graph")
G = nx.from_pandas_adjacency(adj_matrix_directed, create_using=nx.DiGraph)

print("searching for cycles")
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