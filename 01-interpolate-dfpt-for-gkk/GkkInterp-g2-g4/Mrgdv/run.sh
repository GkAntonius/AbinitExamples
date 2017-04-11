#!/bin/bash


MPIRUN='mpirun -n 1'
MRGDV='mrgdv'

$MPIRUN $MRGDV < mrgdv.in &> mrgdv.out 2> stderr

