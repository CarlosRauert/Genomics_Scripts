conda create -n faiss
conda install -c pytorch -c nvidia -c rapidsai -c conda-forge faiss-gpu-raft=1.9.0 pytorch pytorch-cuda=12 numpy
conda install python=3.12
pip install -U "jax[cuda12]"
pip install scib-metrics

du -d 1 -h /data/cephfs-1/work/groups/dubois/users/rauertc_c

