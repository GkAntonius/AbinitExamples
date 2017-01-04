#!/bin/bash


MPIRUN='mpirun -n 4'
ABINIT='/path/to/abinit/build/src/98_main/abinit'

ln -nfs ../../Density/out_data/odat_DEN input_data/idat_DS1_DEN
ln -nfs ../../Density/out_data/odat_DEN input_data/idat_DS2_DEN
ln -nfs ../../Density/out_data/odat_DEN input_data/idat_DS3_DEN
ln -nfs ../out_data/odat_DS1_WFK input_data/idat_DS3_WFK
ln -nfs ../out_data/odat_DS1_WFK input_data/idat_DS3_WFQ

$MPIRUN $ABINIT < calc.files &> calc.log 2> stderr

