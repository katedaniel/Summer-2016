module load gcc
g++ LF_L4.cpp -o LF_L4
chmod +x LF_L4

module load python
module load all-pkgs
module load AstroPy
mkdir qp_file
mkdir qp_file_$1

ls -l
for i in {1..2}
do
    python ./MC_fNew.py
    ./LF_L4
done

ls -l qp_file/
mv !$ qp_file_$1/
tar -czf qp_file_$1.tar.gz qp_file_$1/
