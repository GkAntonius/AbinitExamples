#!/bin/bash



cd Den
bash run.sh
cd ..
cd DVSCF/WFK
bash run.sh
cd ../..
cd EPC/WFK
bash run.sh
cd ../..
cd DVSCF/qpt-0001/WFQ
bash run.sh
cd ../../..
cd DVSCF/qpt-0001/DVSCF
bash run.sh
cd ../../..
cd DVSCF/qpt-0002/WFQ
bash run.sh
cd ../../..
cd DVSCF/qpt-0002/DVSCF
bash run.sh
cd ../../..
cd DVSCF/qpt-0003/WFQ
bash run.sh
cd ../../..
cd DVSCF/qpt-0003/DVSCF
bash run.sh
cd ../../..
cd EPC/qpt-0001/WFQ
bash run.sh
cd ../../..
cd EPC/qpt-0001/EPC
bash run.sh
cd ../../..
cd EPC/qpt-0002/WFQ
bash run.sh
cd ../../..
cd EPC/qpt-0002/EPC
bash run.sh
cd ../../..
cd EPC/qpt-0003/WFQ
bash run.sh
cd ../../..
cd EPC/qpt-0003/EPC
bash run.sh
cd ../../..

