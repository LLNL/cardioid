#!/bin/bash

for exponent in 2 3 4 5 6
do
for mantissa in 5.0 2.5 1.0
do
    dt=${mantissa}e-$exponent
    echo $dt `./OneCell cellml_tt06 -50 0 2 1000 1.1 $dt 1 0 2 | ./computeActivationTime.py`
done
done
