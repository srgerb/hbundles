#!/bin/bash

#FA_REP_CUTOFF=2000
#PDB_LENGTH=184
SCFILE=score.sc

head -2 $SCFILE > ${SCFILE}.filtered

tail -n +3 $SCFILE | while read LINE; do
	NAME=`echo $LINE | awk '{print $NF}'`
	PDB=$NAME.pdb

	# check clash in score.sc
	#FA_REP=`echo $LINE | awk '{print $9}'`
	#if [ $(echo $FA_REP'>'${FA_REP_CUTOFF} | bc -l) -eq 1 ]; then
	#	rm $PDB
	#	continue
	#fi

	# remove pdbs that have two "TER" - ones that have broken chains
	if [ $(grep "TER" $PDB | wc -l) -gt 1 ]; then
		rm $PDB
		continue
	fi

	# length check
	#LENGTH=$(grep "CtermProteinFull" $PDB | awk '{print $1}' | awk -F'_' '{print $NF}')
	#if [[ ! "$LENGTH" -eq "$PDB_LENGTH" ]]; then
	#	rm $PDB
	#	continue
	#fi

	#echo $PDB
	echo $LINE >> ${SCFILE}.filtered
done

