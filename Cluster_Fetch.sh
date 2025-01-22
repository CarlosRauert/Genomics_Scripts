#!/bin/sh

cd $HOME/work/Scripts_Git_Repos/Genomics_Scripts
git pull origin main
git config pull.rebase false
cd $HOME/work/Scripts_Git_Repos/Spatial_Scripts
git pull origin main