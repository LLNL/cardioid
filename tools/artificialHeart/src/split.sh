#!/bin/sh

for n in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39
do
#    cat anatomy#0000* | awk '{if ($5=='$n') print}' > d$n
    echo $n
    cat ana* |  awk '{if ($5=='$n') print}' > d$n
done
