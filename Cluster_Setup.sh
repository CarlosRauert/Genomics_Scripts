#!/bin/sh

tmuxname=$1

tmux new-session -A -s ${tmuxname}
# Check if $HOSTNAME contains "login"
if [[ "$HOSTNAME" == *"login"* ]]; then
  # Run the srun command if the condition is met
  srun --mem=8G --ntasks=12 --pty bash -i
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

