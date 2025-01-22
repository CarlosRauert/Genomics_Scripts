#!/bin/bash
token=$(head -n 1 $HOME/work/Scripts_Git_Repos/20250122_Token.txt)
cd $HOME/work/Scripts_Git_Repos/Spatial_Scripts
git add .
DATE=$(date +"%Y-%m-%d")
git commit -m "Changes of ${DATE}"
git push https://${token}@github.com/CarlosRauert/Spatial_Scripts main

cd $HOME/work/Scripts_Git_Repos/Genomics_Scripts
git add .
DATE=$(date +"%Y-%m-%d")
git commit -m "Changes of ${DATE}"
git push https://${token}@github.com/CarlosRauert/Genomics_Scripts main
git add NF_Jabba/nf-jabba
