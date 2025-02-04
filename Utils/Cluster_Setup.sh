#!/bin/sh

tmuxname=$1

tmux new-session -A -s ${tmuxname}
# Check if $HOSTNAME contains "login"
if [[ "$HOSTNAME" == *"login"* ]]; then
  # Run the srun command if the condition is met
srun --mem=200G --ntasks=64 --time=1-00:00:00 -p medium --pty bash -i
  srun --x11 --pty bash -i 
else
  echo "The hostname does not contain 'login'."
fi
#git config --global user.name "rauertc_c"
#git clone https://github.com/CarlosRauert/Genomics_Scripts $HOME/work/Scripts_Git_Repos/Genomics_Scripts
#git clone https://github.com/CarlosRauert/Spatial_Scripts $HOME/work/Scripts_Git_Repos/Spatial_Scripts

cd $HOME/work/Scripts_Git_Repos/Genomics_Scripts
git fetch origin

cd $HOME/work/Scripts_Git_Repos/Spatial_Scripts
git fetch origin

XAUTHORITY=/data/cephfs-1/home/users/rauertc_c/.Xauthority2

export PATH="/data/cephfs-1/home/users/rauertc_c/work/miniforge3/bin:$PATH"

source /data/cephfs-1/home/users/rauertc_c/work/miniforge3/bin/activate.c~

cp -r /data/cephfs-1/home/users/rauertc_c/work/miniforge3/envs /data/cephfs-1/home/users/rauertc_c/work/envs_backup

wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -u
