module load gcc
g++ LF_L4.cpp -o LF_L4
chmod +x LF_L4

module load python
module load all-pkgs
module load AstroPy
mkdir qp_file_$1

ls -l
#magical forloop
python ./MC_fNew.py
./LF_L4

tar czf qp_file_$1.tar.gz qp_file_$1

