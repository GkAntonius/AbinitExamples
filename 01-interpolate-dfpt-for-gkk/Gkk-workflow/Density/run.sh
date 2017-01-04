#!/bin/bash


MPIRUN='mpirun -n 4'
ABINIT='/path/to/abinit/build/src/98_main/abinit'

$MPIRUN $ABINIT < calc.files &> calc.log 2> stderr

