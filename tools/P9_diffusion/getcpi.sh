#!/bin/bash 
#------------------------------------------------------------------
# primitive utility to extract CPI stack from hpm_job_summary files
# specifically for the Power8 hardware counter set
# Syntax: getcpi.sh hpm_job_summary.jobid
#------------------------------------------------------------------
files=`ls ${1}.*`
echo "run cycles"
for file in $files
do
  grep "(PM_RUN_CYC)" $file | awk '{printf("%.5le %s\n",$5,$6)}'
done
echo " "
echo "stall cycles, no groups completed"
for file in $files
do
  grep "(PM_CMPLU_STALL)" $file | awk '{printf("%.5le %s\n",$5,$6)}'
done
echo " "
   echo "   stall cycles due to loads or stores"
   for file in $files
   do
     grep "(PM_CMPLU_STALL_LSU)" $file | awk '{printf("   %.5le %s\n",$5,$6)}'
   done
   echo " "
   echo "   stall cycles in the floating-point pipeline"
   for file in $files
   do
     grep "(PM_CMPLU_STALL_VSU)" $file | awk '{printf("   %.5le %s\n",$5,$6)}'
   done
   echo " "
   echo "   stall cycles in the fixed-point units"
   for file in $files
   do
     grep "(PM_CMPLU_STALL_FXU)" $file | awk '{printf("   %.5le %s\n",$5,$6)}'
   done
   echo " "
   echo "   stall cycles due to instruction fetch"
   for file in $files
   do
     grep "(PM_CMPLU_STALL_BRU_CRU)" $file | awk '{printf("   %.5le %s\n",$5,$6)}'
   done
   echo " "
   echo "   stall cycles due to instruction flush"
   for file in $files
   do
     grep "(PM_CMPLU_STALL_NTCG_FLUSH)" $file | awk '{printf("   %.5le %s\n",$5,$6)}'
   done
echo " "
echo "cycles from instructions finished to group completed"
for file in $files
do
  grep "(PM_NTCG_ALL_FIN)" $file | awk '{print $5,$6}'
done
echo " "
echo "cycles stalled due to thread conflict"
for file in $files
do
  grep "(PM_CMPLU_STALL_THRD)" $file | awk '{printf("%.5le %s\n",$5,$6)}'
done
echo " "
echo "cycles stalled due to empty instruction pipeline"
for file in $files
do
  grep "(PM_GCT_NOSLOT_CYC)" $file | awk '{printf("%.5le %s\n",$5,$6)}'
done
echo " "
echo "cycles with group completions"
for file in $files
do
  grep "(PM_GRP_CMPL)" $file | awk '{printf("%.5le %s\n",$5,$6)}'
done
