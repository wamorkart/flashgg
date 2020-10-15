#!/bin/bash

training=$1
mass=$2
year=$3
nameVar=$4

python CatTrainMVA.py  -t $training -m $mass -y $year -o CatTrain_${training}_${nameVar}_${year}
cp dataset/weights/TMVAClassification_BDTG.weights.xml CatTrain_${training}_${nameVar}_${year}_TMVAClassification_BDTG.weights.xml
rm -r dataset
