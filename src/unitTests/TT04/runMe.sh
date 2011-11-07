#!/bin/sh

echo "testing bb"
./TT04 bb -52 2 1 1000 3000 2e-4 50 0 > bb_endo
./TT04 bb -52 2 1 1000 3000 2e-4 50 1 > bb_mid
./TT04 bb -52 2 1 1000 3000 2e-4 50 2 > bb_epi

echo "testing mr"
./TT04 mr -52 2 1 1000 3000 2e-4 50 0 > mr_endo
./TT04 mr -52 2 1 1000 3000 2e-4 50 1 > mr_mid
./TT04 mr -52 2 1 1000 3000 2e-4 50 2 > mr_epi

./TT04 mr -52 2 2 1000 3000 2e-4 50 0 > mr_endo_2ms
./TT04 mr -52 2 2 1000 3000 2e-4 50 1 > mr_mid_2ms
./TT04 mr -52 2 2 1000 3000 2e-4 50 2 > mr_epi_2ms

echo "tesing cellml"
./TT04 cellml -52 2 1 1000 3000 2e-4 50 0 > cellml_endo
./TT04 cellml -52 2 1 1000 3000 2e-4 50 1 > cellml_mid
./TT04 cellml -52 2 1 1000 3000 2e-4 50 2 > cellml_epi
