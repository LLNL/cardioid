g++ -fopenmp tt06.cc -o tt06.x -I..
icc -g -std=c++11 -fopenmp tt06.cc -o tt06.x -I.. -DSIMDOPS_ARCH_X86_AVX2 -march=native
icc -g -std=c++11 -fopenmp -O3 -mmic -S -fsource-asm -c tt06.cc -I.. -DSIMDOPS_ARCH_X86_AVX2 -march=native
gcc -std=c99 -c hpm.x86.c
icc -g -std=c++11 -fopenmp tt06.cc -o tt06.x -I.. -DSIMDOPS_ARCH_X86_AVX2 -march=native hpm.x86.o -lpapi -Wl,-rpath=/usr/tce/packages/papi/papi-5.4.3/lib
