conda create -n NF_JAB
conda activate NF_JAB
conda install nf-core
cd work 
nextflow run mskilab/nf-jabba -profile test,bih
# Get NF Jabba Github repository and remove submodule status
cd Scripts_Git_Repos/Genomics_Scripts/NF_Jabba
git clone https://CarlosRauert:$Token@github.com/mskilab-org/nf-jabba.git
git remote add submodule_origin https://github.com/mskilab-org/nf-jabba.git
git fetch submodule_origin
git merge -s ours --no-commit submodule_origin/master --allow-unrelated-histories
git rm --cached NF_Jabba/nf-jabba
cp -r NF_Jabba/nf-jabba/* NF_Jabba
git add .
git commit -m "Converted submodule to regular directory"
rm -rf NF_Jabba/nf-jabba/.git
git add .
git commit -m "Removed .git directory from submodule and converted to regular directory"
git push
# Run from cloned repo
cd NF_Jabba
nextflow run . -profile test
# Get GDC data transfer tool
conda install bioconda::gdc-client

