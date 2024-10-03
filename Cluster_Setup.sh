#!/bin/sh
tmux new-session -A -s default
# Check if $HOSTNAME contains "login"
if [[ "$HOSTNAME" == *"login"* ]]; then
  # Run the srun command if the condition is met
  srun --mem=8G --ntasks=12 --pty bash -i
else
  echo "The hostname does not contain 'login'."
fi
#git config --global user.name "rauertc_c"
git clone https://github.com/CarlosRauert/Genomics_Scripts $HOME/work/Scripts_Git_Repos