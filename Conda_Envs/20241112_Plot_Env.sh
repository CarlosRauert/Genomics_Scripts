conda ycreate -n Plot_Env r-base=4.0.2
install.packages("devtools")
export PKG_CONFIG_PATH=/data/cephfs-1/home/users/rauertc_c/work/miniforge3/envs/Plot_Env/lib/pkgconfig:$PKG_CONFIG_PATH
devtools::install_github("mskilab/gGnome")
