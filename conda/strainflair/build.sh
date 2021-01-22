#!/bin/bash

mkdir -p $PREFIX/bin
cp scripts/*.py $PREFIX/bin
cp StrainFLAIR.sh $PREFIX/bin 
chmod +x $PREFIX/bin/*.py $PREFIX/bin/*.sh
