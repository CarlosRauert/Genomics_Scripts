#!/bin/bash

conda init
conda activate genomics

a=$1
b=$2

Rscript /data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/Scripts/A1KU_Test_Debugged.R $a $b