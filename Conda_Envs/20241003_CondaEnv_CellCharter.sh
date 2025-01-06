conda create --name CellCharter jupyter
conda activate CellCharter
conda install pytorch torchvision torchaudio cpuonly -c pytorch
pip install scvi-tools
pip install cellcharter
pip install --upgrade "jax[cuda12]"
conda install -c conda-forge faiss-gpu
conda install -c pytorch -c nvidia faiss-gpu
conda install -c conda-forge libstdcxx-ngd

conda install -c conda-forge cmake
export PATH=/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH
set(CMAKE_CUDA_ARCHITECTURES 60 61 70 75 80 86)
conda install -c conda-forge swig
conda install conda-forge::gflags
conda install -c nvidia cuda-toolkit=12.6
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

git clone https://github.com/facebookresearch/faiss.git
cd ~/work/faiss
cmake -B build -DFAISS_OPT_LEVEL=generic -DFAISS_ENABLE_PYTHON=ON -DCMAKE_CUDA_HOST_COMPILER=${CUDA_TOOLKIT_ROOT_DIR}/bin/gcc .
make -C build -j32 faiss
make -C build -j32 swigfaiss
cd build/faiss/python
python setup.py install
conda install -c pytorch -c nvidia faiss-gpu=1.9.0

conda install -c pytorch -c nvidia -c rapidsai -c conda-forge faiss-gpu-raft=1.9.0 pytorch pytorch-cuda=11 numpy