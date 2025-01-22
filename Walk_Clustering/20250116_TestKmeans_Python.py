import sklearn
import pandas as pd
from sklearn.cluster import KMeans

filee=pd.read_csv("/data/cephfs-1/home/users/rauertc_c/liposarcoma-wgs/genomics/20250110_Kmed/A1KW/logpca_PC1_2.csv")
filee=filee.drop("Unnamed: 0", axis=1)

x=filee["V1"].array()
y=filee["V2"]

data = list(zip(x, y))
inertias = []

for i in range(1,11):
    kmeans = KMeans(n_clusters=i)
    kmeans.fit(data)
    inertias.append(kmeans.inertia_)