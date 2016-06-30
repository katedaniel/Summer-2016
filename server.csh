module load gcc
g++ LF_L4.cpp -o LF_L4
chmod +x LF_L4

module load python
module load all-pkgs
module load AstroPy
mkdir qp_file_$1

ls -l
#for i in {1..300}
#do
python ./MC_fNew.py
./LF_L4
ls -l qp_file_$1
#done

tar czf qp_file_$1.tar.gz qp_file_$1/*