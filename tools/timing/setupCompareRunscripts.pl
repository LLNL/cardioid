#! /usr/bin/perl -w
#
# create runscripts to run compareSnapshots on pairs of directories containing Cardioid snapshots
#
# written by Erik Draeger, LLNL, 8/27/2012


$bgqExe = "../../../../../bin/compareSnapshots-bgq-spi";
$pelotonExe = "../../../../../bin/compareSnapshots-peloton";

$nBGQnodes = 32;
$nPelnodes = 1;

if ($#ARGV < 0) {
   print "syntax:  setupCompareRunscripts.pl [simulation directories]\n";
}

$cnt1 = 1;
foreach $dir1 (@ARGV[0..$#ARGV]) {
   $cnt2 = 1;
   foreach $dir2 (@ARGV[0..$#ARGV]) {
      if ($dir1 ne $dir2 && -d $dir1 && -d $dir2 && $cnt1 < $cnt2) {

         $dirbase = join '',"dir",$cnt1,"-dir",$cnt2;
         system("mkdir -p $dirbase");

         $bgqbatch = join '',"sbatch.",$dirbase,".bgq";
         open BGQ, ">$dirbase/$bgqbatch";
         print BGQ "\#!/bin/bash\n";
         print BGQ "\#SBATCH --nodes=$nBGQnodes\n";
         print BGQ "\n";
         print BGQ "export OMP_NUM_THREADS=64\n";
         print BGQ "export MUSPI_NUMBATIDS=203\n";
         print BGQ "export MUSPI_NUMINJFIFOS=3\n";
         print BGQ "export MUSPI_NUMRECFIFOS=3\n\n";
         print BGQ "srun --ntasks=$nBGQnodes $bgqExe ../$dir1 ../$dir2\n";
         close BGQ;

         $pelbatch = join '',"msub.",$dirbase,".pel";
         $nPeltasks = 16*$nPelnodes;
         open PEL, ">$dirbase/$pelbatch";
         print PEL "\#!/bin/bash\n";
         print PEL "\#MSUB -l nodes=$nPelnodes\n";
         print PEL "\#MSUB -l walltime=12:00:00\n";
         print PEL "\#MSUB -A gbcq\n";
         print PEL "\n";
         print PEL "export OMP_NUM_THREADS=1\n";
         print PEL "srun -n $nPeltasks $pelotonExe ../$dir1 ../$dir2\n";
         close PEL;

         $peldebug = join '',"rundebug.",$dirbase,".pel";
         $pelout = join '',$dirbase,".out";
         open DEB, ">$dirbase/$peldebug";
         print DEB "\#!/bin/bash\n";
         print DEB "export OMP_NUM_THREADS=1\n";
         print DEB "srun -N $nPelnodes -n $nPeltasks -p pdebug $pelotonExe ../$dir1 ../$dir2 > $pelout\n";   
         close DEB;
         system("chmod u+x $dirbase/$peldebug");

      }
      $cnt2++;
   }
   $cnt1++;
}
