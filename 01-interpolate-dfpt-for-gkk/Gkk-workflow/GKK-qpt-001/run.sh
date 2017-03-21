#!/bin/bash


MPIRUN='mpirun -n 4'
ABINIT='abinit'

ln -nfs ../../Density/out_data/odat_DEN input_data/idat_DS1_DEN
ln -nfs ../../Density/out_data/odat_DEN input_data/idat_DS2_DEN
ln -nfs ../../Density/out_data/odat_DEN input_data/idat_DS3_DEN
ln -nfs ../out_data/odat_DS1_WFK input_data/idat_DS3_WFK
ln -nfs ../out_data/odat_DS2_WFQ input_data/idat_DS3_WFQ
ln -nfs ../../Mrgddb/out_data/odat_DDB input_data/idat_DS3_DDB
ln -nfs ../../Mrgdv/out_data/odat_DVDB input_data/idat_DS3_DVDB

$MPIRUN $ABINIT < calc.files &> calc.log 2> stderr

