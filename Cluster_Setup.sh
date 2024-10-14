#!/bin/sh

tmuxname=$1

tmux new-session -A -s ${tmuxname}
# Check if $HOSTNAME contains "login"
if [[ "$HOSTNAME" == *"login"* ]]; then
  # Run the srun command if the condition is met
  srun --mem=600G --ntasks=32 --pty --time 0-03:00:00 bash -i
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

