#!/bin/sh

echo "testing bb-tt04"
./TT04 bb_tt04 -52 2 1 1000 3000 2e-2 50 0 0 > bb_tt04_endo
./TT04 bb_tt04 -52 2 1 1000 3000 2e-2 50 0 1 > bb_tt04_mid
./TT04 bb_tt04 -52 2 1 1000 3000 2e-2 50 0 2 > bb_tt04_epi

echo "testing cellml_tt04"
./TT04 cellml_tt04 -52 2 1 1000 3000 2e-2 50 0 0 > cellml_tt04_endo
./TT04 cellml_tt04 -52 2 1 1000 3000 2e-2 50 0 1 > cellml_tt04_mid
./TT04 cellml_tt04 -52 2 1 1000 3000 2e-2 50 0 2 > cellml_tt04_epi

echo "testing cellml_tt04_fe"
./TT04 cellml_tt04_fe -52 2 1 1000 3000 2e-4 50 0 0 > cellml_tt04_fe_endo
./TT04 cellml_tt04_fe -52 2 1 1000 3000 2e-4 50 0 1 > cellml_tt04_fe_mid
./TT04 cellml_tt04_fe -52 2 1 1000 3000 2e-4 50 0 2 > cellml_tt04_fe_epi

echo "testing cellml_tt06_fe"
./TT04 cellml_tt06_fe -52 2 1 1000 3000 2e-4 50 0 0 > cellml_tt06_endo
./TT04 cellml_tt06_fe -52 2 1 1000 3000 2e-4 50 0 1 > cellml_tt06_mid
./TT04 cellml_tt06_fe -52 2 1 1000 3000 2e-4 50 0 2 > cellml_tt06_epi

echo "testing tt04dev"
./TT04 tt04dev -52 2 1 1000 3000 2e-2 50 0 0 > tt04_dev_endo
./TT04 tt04dev -52 2 1 1000 3000 2e-2 50 0 1 > tt04_dev_mid
./TT04 tt04dev -52 2 1 1000 3000 2e-2 50 0 2 > tt04_dev_epi

