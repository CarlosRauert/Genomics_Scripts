conda create -n NF_JAB
conda activate NF_JAB
conda install nf-core
cd work 
nextflow run mskilab/nf-jabba -profile test,bih

cd Scripts_Git_Repos/Genomics_Scripts/NF_Jabba
git clone https://CarlosRauert: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@github.com/mskilab-org/nf-jabba.git
cd nf-jabba
nextflow run . -profile test
