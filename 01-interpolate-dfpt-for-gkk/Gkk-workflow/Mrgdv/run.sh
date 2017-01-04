#!/bin/bash


MPIRUN='mpirun -n 1'
MRGDV='/path/to/abinit/build/src/98_main/mrgdv'

$MPIRUN $MRGDV < mrgdv.in &> mrgdv.out 2> stderr

