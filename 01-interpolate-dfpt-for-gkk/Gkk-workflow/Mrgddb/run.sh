#!/bin/bash


MPIRUN='mpirun -n 1'
MRGDDB='mrgddb'

$MPIRUN $MRGDDB < mrgddb.in &> mrgddb.out 2> stderr

