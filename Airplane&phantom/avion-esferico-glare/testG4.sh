#!/bin/bash

rm -r build/*
cd build
source ~/work/GEANT4/10.07.p04-install/bin/geant4.sh
cmake ~/work/GEANT4/10.07.p04-install/lib/Geant4-10.7.4/ ../
make
cp ../testG4b .
bash testG4b 
cd ../

