#!/bin/bash

conda init
conda activate ficture

a=$1
b=$2

Rscript /data/cephfs-1/work/groups/dubois/users/rauertc_c/genomics/Scripts/Print_Test.R $a $b