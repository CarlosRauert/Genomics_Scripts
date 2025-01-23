conda create -n gdc-env python=3.8
conda activate gdc-env
mamba install -c conda-forge pyopenssl=18.0.0 # didn work
conda install -c conda-forge mamba
mamba install -c bioconda gdc-client=1.3
conda install -c bioconda gdc-client

unzip /data/cephfs-1/home/users/rauertc_c/work/GDC_Datatransfer/gdc-client_2.3_Ubuntu_x64-py3.8-ubuntu-20.04.zip
cd /data/cephfs-1/home/users/rauertc_c/work/GDC_Datatransfer
export PATH=$PATH:/data/cephfs-1/home/users/rauertc_c/work/GDC_Datatransfer
source ~/.bashrc

conda create -n gdc-env python=3.8
wget -c https://ftp.gnu.org/gnu/glibc/glibc-2.29.tar.gz
tar -zxvf glibc-2.29.tar.gz
cd glibc-2.29
mkdir build
cd build
../configure --prefix=/opt/glibc
../configure --prefix=/opt/glibc-2.29 CFLAGS="-Wno-error"

make
make install
unset LD_LIBRARY_PATH
unset LD_PRELOAD

mkdir -p /data/cephfs-1/home/users/rauertc_c/work/glibc_new/lib64
cp -r /usr/lib64/* /data/cephfs-1/home/users/rauertc_c/work/glibc_new/lib64/

export LD_LIBRARY_PATH=/data/cephfs-1/home/users/rauertc_c/work/2.29/build:$LD_LIBRARY_PATH
export LD_PRELOAD=/data/cephfs-1/home/users/rauertc_c/work/2.29/build/libc.so.6
export PATH=/data/cephfs-1/home/users/rauertc_c/work/2.29/build/bin:$PATH



echo $LD_LIBRARY_PATH
echo $LD_PRELOAD

token=$(cat /data/cephfs-1/home/users/rauertc_c/work/GDC_Datatransfer/gdc-user-token.2025-01-22T15_32_13.891Z.txt)

curl --remote-name --remote-header-name 'https://api.gdc.cancer.gov/data/5a993135-be43-4c41-9aa4-3fb3d1abea29'
curl https://api.gdc.cancer.gov/files/cb92f61d-041c-4424-a3e9-891b7545f351?pretty=true
curl -O -J -H "X-Auth-Token: $token" 'https://api.gdc.cancer.gov/5a993135-be43-4c41-9aa4-3fb3d1abea29'
ping doublehelix.helix.prod.aws.gel.ac
doublehelix.helix.prod.aws.gel.ac

../glibc-2.29/configure --prefix=/data/cephfs-1/home/users/rauertc_c/work/2.29/made


#include <stdio.h>
int main() {
    printf("GLIBC test successful!\n");
    return 0;
}