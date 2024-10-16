#!/bin/bash

# First, run freq series somehow for the Branco-Hausser stim
# May have to instead scale receptor weight or # of inputs.

for ww in 5 10 15 20 25 30 35 40;
do
	echo $ww
	python allNonlin7.py --frequency 1 --seqDt 0.01 --seqOnT 0.01 --seqDx 6 --seqN 7 --xOffset 2 --runtime 4 -w $ww
done

for dx in 1 2 3 4 5 6 8 10;
do
	echo $dx
	python allNonlin7.py --frequency 1 --seqDt 0.01 --seqOnT 0.01 --seqDx $dx --seqN 7 --xOffset 2 --runtime 4 -w 20
done
