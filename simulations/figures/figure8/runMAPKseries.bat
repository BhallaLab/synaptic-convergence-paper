#!/bin/bash

# First, run freq series somehow for the Branco-Hausser stim
# May have to instead scale receptor weight or # of inputs.

for ww in 5 10 15 20 25 30 35 40;
do
	echo $ww
	python allNonlin7.py --frequency 20 --seqDt 3 --seqOnT 2 --seqDx 2 --seqN 5 --xOffset 20 --runtime 30 -w $ww
done


for ff in 5 10 15 20 25 30 35 40;
do
	echo $ff
	python allNonlin7.py --frequency $ff --seqDt 3 --seqOnT 2 --seqDx 2 --seqN 5 --xOffset 20 --runtime 30
done


for dx in 1 2 3 4 5 6 8 10;
do
	echo $dx
	python allNonlin7.py --frequency 20 --seqDt 3 --seqOnT 2 --seqDx $dx --seqN 5 --xOffset 20 --runtime 30
done
	
	
