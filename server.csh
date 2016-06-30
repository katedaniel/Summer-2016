#!/bin/csh

mkdir qp_file_$1
#magical forloop
python MC_fnew.py
./LF_L4

tar czf runNfilename.tar.gz qp*

