export OMP_NUM_THREADS=40
export OMP_DISPLAY_ENV=true
export OMP_PROC_BIND=close
for x in simd_cpu simd_cpu_hoist simd_hoist_to_mthrd simd_cpu_thr simd_cpu_thr_30 simd_cpu_thr_30_oh simd_cpu_thr_lg simd_mthrd_30disj simd_cpu_thr2; do echo $x; ./diffusion -x 100 -y 242 -z 152 -K $x; done

