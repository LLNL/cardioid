
.SUFFIXES: .c .o .cc .cu .cpp

SRC = \
     main.cpp \
     simd_mthrd.cpp simd_mthrd_parti.cpp simd_hoist.cpp simd_baseline.cpp nosimd.cpp \
     simd_mthrd_larger.cpp simd_mthrd_30x30.cpp simd_hoist_to_mthrd.cpp simd_mthrd_30x30_overhead.cpp simd_mthrd_30disj.cpp \
     options.c

OBJS = $(SRC:.cpp=.o)

#LIB=-L/home/ckim/project/CORAL/cardioid/serial_hpm/lib/ -lhpm -Wl,-rpath=/home/ckim/project/CORAL/cardioid/serial_hpm/lib
LIB=

#CXX= /shared/opt/xlC/14.1.0.0-170623/xlC/14.1.0/bin/xlC -g
#CXX_R= /shared/opt/xlC/14.1.0.0-170623/xlC/14.1.0/bin/xlC_r -g -qsmp=omp
CXX= xlC -g
CXX_R= xlC_r -g -qsmp=omp

all: diffusion

diffusion: $(OBJS)
	$(CXX_R) $(LIB) $^ -o $@
#	$(CXX) -O3 $^ -o $@

simd_mthrd.o : simd_mthrd.cpp
#	$(CXX_R) -qaltivec -qarch=pwr9 -O5 -c $^ -qreport -qlist
	$(CXX_R) -qarch=pwr9 -O5 -c $^ -qreport -qlist

simd_mthrd_larger.o : simd_mthrd_larger.cpp
	$(CXX_R) -qaltivec -qarch=pwr9 -O5 -c $^ -qreport -qlist

simd_hoist.o : simd_hoist.cpp
	$(CXX) -qaltivec -qarch=pwr9 -O3 -qhot -c $^ -qreport -qlist

simd_hoist_to_mthrd.o : simd_hoist_to_mthrd.cpp
	$(CXX) -qaltivec -qarch=pwr9 -O3 -qhot -c $^ -qreport -qlist

simd_mthrd_30disj.o : simd_mthrd_30disj.cpp
	$(CXX_R) -qaltivec -qarch=pwr9 -O3 -qhot -c $^ -qreport -qlist

simd_mthrd_30x30.o : simd_mthrd_30x30.cpp
	$(CXX_R) -qaltivec -qarch=pwr9 -O3 -qhot -c $^ -qreport -qlist

simd_mthrd_30x30_overhead.o : simd_mthrd_30x30_overhead.cpp
	$(CXX_R) -qaltivec -qarch=pwr9 -O3 -qhot -c $^ -qreport -qlist

simd_mthrd_parti.o : simd_mthrd_parti.cpp
	$(CXX_R) -qaltivec -qarch=pwr9 -O5 -c $^ -qreport -qlist

nosimd.o : nosimd.cpp
	$(CXX) -qarch=pwr9 -O5 -c $^

simd_baseline.o:simd_baseline.cpp
	$(CXX) -qarch=pwr9 -qaltivec -O3 -qhot -c $^ -qreport -qlist


#.cc.o:
#	$(CXX) -g $(CFLAGS) $(INCLUDES) -c $<  -o $@
#
#.cpp.o:
#	$(CXX) -g $(CFLAGS) $(INCLUDES) -c $<  -o $@

#.cu.o:
#	nvcc -g -arch=sm_60 -c $<  -o $@

.c.o:
	xlc -g $(CFLAGS) $(INCLUDES) -c $<  -o $@


clean:
	rm -rf diffusion
	rm *.o

GENGETOPT ?= gengetopt

gengetopt:
	$(GENGETOPT) -i options.ggo -F options

test1:
	for x in `seq 1 21`; do export HPM_GROUP=$x; ./diffusion -x 200 -y 32 -z 32 -K simd_cpu -N 10; done
