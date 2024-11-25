conda create -n Liana2
conda install pip
conda install python=3.11
cd /data/cephfs-1/home/users/rauertc_c/work/miniforge3/envs/Liana/conda-meta
echo "python=3.11" > pinned
cat pinned
pip install decoupler
pip install plotnine
pip install mudata
pip install adjustText
pip install git+https://github.com/saezlab/liana-py
conda install scipy=1.13.1
conda install -c conda-forge gxx_linux-64

