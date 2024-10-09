#!/bin/sh

CaseID=$1
CNxt=$2

conda init
conda activate /data/cephfs-1/home/users/rauertc_c/work/miniforge3/envs/Genomics

# Loop until CNxt reaches or exceeds 8
while (( $(echo "$CNxt < 8" | bc -l) )); do
  # Execute the R script
  Rscript /data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/MDM2_Walks_Python/20241009_Prepare_Adj_Python.R ${CaseID} ${CNxt}

  # Define the MatrixFile path
  MatrixFile=/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/${CaseID}/${CaseID}_adj_Matrix_${CNxt}.csv

  # Make the Python script executable
  chmod +x /data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/MDM2_Walks_Python/20241009_CycleSearch.py

  # Run the Python script and capture output
  python3 -u /data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/MDM2_Walks_Python/20241009_CycleSearch.py 2>&1 ${MatrixFile} ${CaseID} ${CNxt}

  # Check if "No_Cycle_Found" is in the output file
  if grep -q No_Cycle_Found "/data/cephfs-1/home/users/rauertc_c/work/genomics/MDM2_Walks_Out/${CaseID}/${CaseID}_${CNxt}_Cycles.txt"; then
    echo "No Cycles Were Found" # String was found
  else
    echo "String Was Not Found"
    Rscript /data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/MDM2_Walks_Python/20241009_CycleSearch_MakegWalk.R ${CaseID} ${CNxt}
  fi

  # Update CNxt by multiplying it by 0.8 and rounding down to the nearest integer
  CNxt=$(echo "scale=0; $CNxt * 0.8 / 1" | bc)

  # Print current value of CNxt for debugging
  echo "Current CNxt: $CNxt"
done

echo "done"