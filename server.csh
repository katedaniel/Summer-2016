#!/bin/csh

module load python
module load all-pkgs
module load AstroPy
chmod +x LF_L4

mkdir qp_file_$1
#magical forloop
python MC_fnew.py
./LF_L4

tar czf runNfilename.tar.gz qp*

