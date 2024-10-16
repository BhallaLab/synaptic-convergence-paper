#!/bin/bash

# First, run freq series somehow for the Branco-Hausser stim
# May have to instead scale receptor weight or # of inputs.

for ww in 5 10 15 20 25 30 35 40;
do
	echo $ww
	python allNonlin7.py --frequency 1 --seqDt 0.01 --seqOnT 0.01 --seqDx 6 --seqN 7 --xOffset 2 --runtime 4 -w $ww
done
