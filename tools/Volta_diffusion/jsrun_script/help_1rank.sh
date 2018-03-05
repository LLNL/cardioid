#!/bin/bash 
#--------------------------------------------------------------------------------
# jsrun --rs_per_host ${ppn} --np ${nmpi}  helper.sh  your.exe [args]
# You must set env variable RANKS_PER_NODE in your job script;
# optionally set BIND_SLOTS in your job script = #hwthreads per rank.
# Set USE_GOMP=yes to get proper binding for GNU and PGI OpenMP runtimes.
#--------------------------------------------------------------------------------
let ngpus=6

 cpus_per_node=160
 declare -a list0=(`seq 0 79`)
 declare -a list1=(`seq 88 167`)
 declare -a mpscpu=(80-83 168-171)
 declare -a gpulist=(0 1 2 3 4 5)

if [ -z "$PMIX_RANK" ]; then
  echo helper.sh : PMIX_RANK is not set ... exiting
fi

let world_rank=$PMIX_RANK

if [ -z "$RANKS_PER_NODE" ]; then
  if [ $world_rank = 0 ]; then
    echo helper.sh : you must set RANKS_PER_NODE ... exiting
  fi
  exit
fi

let local_size=$RANKS_PER_NODE
let local_rank=$(expr $world_rank % $local_size)

export CUDA_CACHE_PATH=/dev/shm/$USER/nvcache_$local_rank


#---------------------------------------------
# set CUDA device for each MPI rank
#---------------------------------------------
let product=$ngpus*$local_rank
let gpu=$product/$local_size
let mydevice=${gpulist[$gpu]}

#-------------------------------------------------
# assign socket and affinity mask
#-------------------------------------------------
let x2rank=2*$local_rank
let socket=$x2rank/$local_size
let ranks_per_socket=1

#---------------------------------------------
# optionally start MPS, one per socket
#---------------------------------------------
if [ "$USE_MPS" == "yes" ]; then
  if [ $local_rank = 0 ]; then
    if [ $world_rank = 0 ]; then
      echo starting mps ...
    fi  

    rm -rf /dev/shm/${USER}/mps_0
    rm -rf /dev/shm/${USER}/mps_log_0
    mkdir -p /dev/shm/${USER}/mps_0
    mkdir -p /dev/shm/${USER}/mps_log_0
    export CUDA_VISIBLE_DEVICES=0,1,2
    export CUDA_MPS_PIPE_DIRECTORY=/dev/shm/${USER}/mps_0
    export CUDA_MPS_LOG_DIRECTORY=/dev/shm/${USER}/mps_log_0
    taskset -c ${mpscpu[0]} /usr/bin/nvidia-cuda-mps-control -d
    sleep 1

    rm -rf /dev/shm/${USER}/mps_1
    rm -rf /dev/shm/${USER}/mps_log_1
    mkdir -p /dev/shm/${USER}/mps_1
    mkdir -p /dev/shm/${USER}/mps_log_1
    export CUDA_VISIBLE_DEVICES=3,4,5
    export CUDA_MPS_PIPE_DIRECTORY=/dev/shm/${USER}/mps_1
    export CUDA_MPS_LOG_DIRECTORY=/dev/shm/${USER}/mps_log_1
    taskset -c ${mpscpu[1]} /usr/bin/nvidia-cuda-mps-control -d
    sleep 1

  else
    sleep 3
  fi
  sleep 1
  printf -v myfile "/dev/shm/${USER}/mps_%d" $socket
  export CUDA_MPS_PIPE_DIRECTORY=$myfile
  unset CUDA_VISIBLE_DEVICES
else
  export CUDA_VISIBLE_DEVICES=0,1,2,3,4,5
fi

# divide available slots evenly or specify slots by env variable
if [ -z "$BIND_SLOTS" ]; then
  let cpus_per_rank=$cpus_per_node/$local_size
else
  let cpus_per_rank=$BIND_SLOTS 
fi

if [ -z "$OMP_NUM_THREADS" ]; then
  let num_threads=$cpus_per_rank
else
  let num_threads=$OMP_NUM_THREADS
fi

# BIND_STRIDE is used in OMP_PLACES ... it will be 1 if OMP_NUM_THREADS was not set
let BIND_STRIDE=$(expr $cpus_per_rank / $num_threads)
#echo BIND_STRIDE = $BIND_STRIDE

if [ $socket = 0 ]; then
  let ndx=$local_rank*$cpus_per_rank
  let start_cpu=${list0[$ndx]}
  let stop_cpu=$start_cpu+$cpus_per_rank-1
else
  let rank_in_socket=$local_rank-$ranks_per_socket
  let ndx=$rank_in_socket*$cpus_per_rank
  let start_cpu=${list1[$ndx]}
  let stop_cpu=$start_cpu+$cpus_per_rank-1
fi

#---------------------------------------------
# set OMP_PLACES or GOMP_CPU_AFFINITY
#---------------------------------------------
if [ "$USE_GOMP" == "yes" ]; then
  export GOMP_CPU_AFFINITY="$start_cpu-$stop_cpu:$BIND_STRIDE"
  unset OMP_PLACES
else
  export OMP_PLACES={$start_cpu:$num_threads:$BIND_STRIDE}
fi

#-------------------------------------------------
# set an affinity mask for each rank using taskset
#-------------------------------------------------
printf -v command "taskset -c %d-%d"  $start_cpu  $stop_cpu 
echo command = $command
echo exec = $1

executable=$1

shift

echo $command $executable "$@"
#-------------------------
# run the code
#-------------------------
$command $executable "$@"

#---------------------------------------------
# optionally stop MPS 
#---------------------------------------------
if [ "$USE_MPS" == "yes" ]; then
  if [ $local_rank = 0 ]; then
    if [ $world_rank = 0 ]; then
      echo stopping mps ...
    fi
    for ((i=0; i<2; i++))
    do  
     export CUDA_MPS_PIPE_DIRECTORY=/dev/shm/${USER}/mps_$i
     echo "quit" | /usr/bin/nvidia-cuda-mps-control
     sleep 1
     rm -rf /dev/shm/${USER}/mps_$i
     rm -rf /dev/shm/${USER}/mps_log_$i
    done
    rm -rf /dev/shm/${USER}
    unset CUDA_MPS_PIPE_DIRECTORY
  fi
fi


if [ $local_rank = 0 ]; then
 rm -rf /dev/shm/$USER/nvcache_*
fi
