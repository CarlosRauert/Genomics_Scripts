#!/bin/sh

CaseID=$1
CNxt=$2

conda init
conda activate /data/cephfs-1/home/users/rauertc_c/work/miniconda/envs/genomics

Rscript /data/cephfs-1/home/users/rauertc_c/work/genomics/Scripts/MDM2_Walks_Python/Prepare_Adj_Python.R ${CaseID} ${CNxt}

MatrixFile=/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/${CaseID}/${CaseID}_adj_Matrix_${CNxt}.csv

chmod +x /data/cephfs-1/home/users/rauertc_c/work/genomics/Scripts/MDM2_Walks_Python/CycleSearch.py
python3 -u /data/cephfs-1/home/users/rauertc_c/work/genomics/Scripts/MDM2_Walks_Python/CycleSearch.py 2>&1 ${MatrixFile} ${CaseID} ${CNxt}

if grep -q No_Cycle_Found "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/${CaseID}/${CaseID}_${CNxt}_Cycles.txt"; then
  echo "No Cycles Were Found" # String was found
else
  echo "String Was Not Found"
  Rscript /data/cephfs-1/home/users/rauertc_c/work/genomics/Scripts/MDM2_Walks_Python/CycleSearch_MakegWalk.R ${CaseID} ${CNxt}
fi

echo "done"