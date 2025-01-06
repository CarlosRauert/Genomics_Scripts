conda create -n scib python=3.12
conda activate scib
pip install scib-metrics
pip install cellcharter
pip install scvi-tools
pip install --upgrade "jax[cuda12]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
ln -s work/.cache .cache
