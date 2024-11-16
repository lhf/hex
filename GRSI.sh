#!/bin/sh

N=6

for R in `seq 0 $N`
do
	make R=$R O=hex-$R
done

