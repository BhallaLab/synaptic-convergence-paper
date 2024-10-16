#!/bin/bash

# First, run freq series somehow for the Branco-Hausser stim
# May have to instead scale receptor weight or # of inputs.

for ww in 5 10 15 20 25 30 35 40;
do
	echo $ww
	python allNonlin7.py --frequency 10 --seqDt 1 --seqOnT 0.5 --seqDx 1 --seqN 5 --xOffset 2 --runtime 5 --groupStimulus -w $ww
done

for dx in 1 2 3 4 5;
do
	echo $dx
	python allNonlin7.py --frequency 10 --seqDt 1 --seqOnT 0.5 --seqDx $dx --seqN 5 --xOffset 2 --runtime 5 --groupStimulus
done
