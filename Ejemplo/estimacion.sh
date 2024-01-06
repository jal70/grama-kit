#!/bin/bash
# Version G4 10.07.p04
rutaG4=~/work/GEANT4/10.07.p04-install/
lines=$(wc -l < instruccionesG4.mac)
lines=$(((lines-1-14) / 5))
rm ED_totales/*
#cp coordenadas.csv ED_totales/

for ((i=0; i<lines; i++))
do
	awk 'BEGIN{FS=","} NR=='$((i*5+2+14))', NR=='$((i*5+6+14))'{print $2}' instruccionesG4.mac
	echo "/run/initialize" > aeronaveG4/gamma.mac
	echo "/control/verbose 2" >> aeronaveG4/gamma.mac
	echo "/tracking/verbose 2" >> aeronaveG4/gamma.mac
	awk 'BEGIN{FS=","} NR=='$((i*5+2+14))', NR=='$((i*5+6+14))'{print $2}' instruccionesG4.mac  >> aeronaveG4/gamma.mac
	rm -r aeronaveG4/build/*
	cd aeronaveG4/build
	source $rutaG4""bin/geant4.sh
	cmake $rutaG4""lib/Geant4-10.7.4/ ../
	make
	./avion_esferico_glare
	primaria=$(awk 'NR==4{print $2}' gamma.mac)
	energia=$(awk 'NR==5{print $2}' gamma.mac)
	numPrim=$(awk 'NR==8{print $2}' gamma.mac)
	cd ../../
	echo "$primaria""_$energia""_$numPrim" > ED_totales/$i""_E.csv
	awk '{print $4}' aeronaveG4/build/totalED.txt >> ED_totales/$i""_E.csv
done
awk 'BEGIN{FS=","} NR=='1', NR=='14'{print $0}' instruccionesG4.mac > ED_totales/head
cd ED_totales/
paste -d' ' ../coordenadas.csv *_E.csv > totalED_t.csv
cat head totalED_t.csv > totalED.csv
rm *_E.csv head totalED_t.csv
echo "gitkeep" > .gitkeep
cd ../
rm -r aeronaveG4/build/*
echo "gitkeep" > aeronaveG4/build/.gitkeep
