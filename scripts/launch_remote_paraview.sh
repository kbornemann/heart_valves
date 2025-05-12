#!/bin/bash

# Enter ParaView folder

cd /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin

# Start pvserver

srun -p amarsden pvserver --server-port=12345
