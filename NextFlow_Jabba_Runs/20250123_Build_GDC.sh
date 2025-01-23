cd ~/work/Scripts_Git_Repos
git clone --branch feat/try-macos-11 https://github.com/NCI-GDC/gdc-client.git
cd gdc-client
cd bin

conda_init
conda create -n GDC_New python=3.8
conda activate GDC_New

mkdir -p ~/bin
ln -s /usr/bin/python3 ~/bin/python
export PATH="$HOME/bin:$PATH"
source ~/.bashrc
pip install virtualenv

./package
pip install --upgrade pip
pip install -r requirements.txt
python setup.py install


