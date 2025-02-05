#!/bin/bash
token=$(head -n 1 $HOME/work/Scripts_Git_Repos/20250122_Token.txt)
cd $HOME/work/Scripts_Git_Repos/Spatial_Scripts
git add .
DATE=$(date +"%Y-%m-%d")
git commit -m "Changes of ${DATE}"
git push https://${token}@github.com/CarlosRauert/Spatial_Scripts main
git config --global credential.helper store

cd $HOME/work/Scripts_Git_Repos/Genomics_Scripts
git add .
DATE=$(date +"%Y-%m-%d")
git commit -m "Changes of ${DATE}"
git[] push https://CarlosRauert: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@github.com/CarlosRauert/Genomics_Scripts main
git add NF_Jabba/nf-jabba/*
cd NF_Jabba/nf-jabba
git add .
git commit -m "Changes of ${DATE}"
cd ../..
git add NF_Jabba/nf-jabba
git commit -m "Update submodule NF_Jabba/nf-jabba"
git push https://${token}@github.com/CarlosRauert/Genomics_Scripts main


git push ${token}@github.com/CarlosRauert/Genomics_Scripts.git main