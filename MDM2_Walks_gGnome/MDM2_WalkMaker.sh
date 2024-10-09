#!/bin/bash

conda init
conda activate genomics

a=$1

Rscript /data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/Scripts/MDM2_WalkMaker_ForShellScript.R $a