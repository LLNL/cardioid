nvcc -arch=sm_60 kernel.cu -o repack.x
nvprof ./repack.x
