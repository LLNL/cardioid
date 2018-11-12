#!/bin/sh

echo "testing bb-tt04"
./OneCell bb_tt04 -52 50 1 1000 3000 2e-2 50 0 0 > bb_tt04_endo
./OneCell bb_tt04 -52 50 1 1000 3000 2e-2 50 0 1 > bb_tt04_mid
./OneCell bb_tt04 -52 50 1 1000 3000 2e-2 50 0 2 > bb_tt04_epi

echo "testing cellml_tt04"
./OneCell cellml_tt04 -52 50 1 1000 3000 2e-2 50 0 0 > cellml_tt04_endo
./OneCell cellml_tt04 -52 50 1 1000 3000 2e-2 50 0 1 > cellml_tt04_mid
./OneCell cellml_tt04 -52 50 1 1000 3000 2e-2 50 0 2 > cellml_tt04_epi

echo "testing cellml_tt04_fe"
./OneCell cellml_tt04_fe -52 50 1 1000 3000 2e-4 50 0 0 > cellml_tt04_fe_endo
./OneCell cellml_tt04_fe -52 50 1 1000 3000 2e-4 50 0 1 > cellml_tt04_fe_mid
./OneCell cellml_tt04_fe -52 50 1 1000 3000 2e-4 50 0 2 > cellml_tt04_fe_epi

echo "testing cellml_tt06_fe"
./OneCell cellml_tt06_fe -52 50 1 1000 3000 2e-4 50 0 0 > cellml_tt06_endo
./OneCell cellml_tt06_fe -52 50 1 1000 3000 2e-4 50 0 1 > cellml_tt06_mid
./OneCell cellml_tt06_fe -52 50 1 1000 3000 2e-4 50 0 2 > cellml_tt06_epi

echo "testing fhn"
./OneCell fhn -52 50 1 1000 3000 2e-2 50 0 2 > fhn

echo "testing tt04dev"
./OneCell tt04dev -52 50 1 1000 3000 2e-2 50 0 0 > tt04_dev_endo
./OneCell tt04dev -52 50 1 1000 3000 2e-2 50 0 1 > tt04_dev_mid
./OneCell tt04dev -52 50 1 1000 3000 2e-2 50 0 2 > tt04_dev_epi

