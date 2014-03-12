#!/bin/sh

# create gradients
tend grads -n 12 \
 | unu pad -min 0 -1 -max M M -b pad -v 0 -o grads.txt

# create tensors
tend helix -s 29 30 31 -ip 0.1 0.3 0.6 -mp -0.8 0.1 0.4 -r 40 -o ten.nrrd

# create synthetic DWI from tensors
tend helix -s 29 30 31 -ip 0.1 0.3 0.6 -mp -0.8 0.1 0.4 -r 50 \
 | tend anvol -a fa \
 | unu 2op x - 7000 -o b0.nrrd

tend sim -g grads.txt -r b0.nrrd -i ten.nrrd \
  -b 800 -t ushort -kvp -o dwi-D.nhdr

# Estimate the gradients from the DT-MRI data
tend estim -i dwi-D.nhdr -knownB0 false -B kvp -o ten-from-dwi.nrrd

