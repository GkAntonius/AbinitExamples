#!/bin/bash


MPIRUN='mpirun -n 1'
MRGDDB='/path/to/abinit/build/src/98_main/mrgddb'

$MPIRUN $MRGDDB < mrgddb.in &> mrgddb.out 2> stderr

