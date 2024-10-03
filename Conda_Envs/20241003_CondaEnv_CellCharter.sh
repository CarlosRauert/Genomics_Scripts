conda create --name CellCharter jupyter
conda activate CellCharter
conda install pytorch torchvision torchaudio cpuonly -c pytorch
pip install scvi-tools
pip install cellcharter
