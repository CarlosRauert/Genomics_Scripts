#!/bin/bash
token=$(head -n 1 $HOME/work/dir/token.txt)
cd $HOME/work/Scripts_Git_Repos/Spatial_Scripts
git push https://${token}@github.com/CarlosRauert/Spatial_Scripts main

cd $HOME/work/Scripts_Git_Repos/Genomics_Scripts
git push https://${token}@github.com/CarlosRauert/Genomics_Scripts main