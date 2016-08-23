#! /usr/bin/perl -w
#
# generate Cardioid object.data input files for verification runs.  Compare output data 
# using compareSnapshots.cc (in Cardioid src directory) 
#
# written by Erik Draeger, LLNL, 8/20/2012

# global variables

$thisdir = `pwd`;  chomp $thisdir;
$makeVoidAnatomyScript = "$thisdir/tools/makeAnatomyBlockWithVoids.pl";
$bgqExe = "../../../../bin/cardioid-bgq-spi";
$pelotonExe = "../../../../bin/cardioid-peloton";
$nthreadsBGQ = 64;
$nthreadsPeloton = 4;

$nIterations = 100000;
$checkpointRate = 1000;

$weakScaling = 1;   # if set to zero, anatomy size corresponding 
                    # to $strongTaskCount will be used throughout
$strongTaskCount = 64;

foreach $anatomy ("block247", "swiss247")
{
   foreach $celltype ("random")
   {
      # these correspond to if blocks in printObject (below), not cardioid reaction types
      #foreach $reaction ("TT06RRG", "TT06RRGOpt", "TT06", "TT06Opt") 
      foreach $reaction ("TT06RRGOpt") 
      {
         foreach $fastgates (1)
         {
            foreach $rationalfns (1)
            {
               $smoothing = 1;
               foreach $setgkr (0, 1, 2)
               {
                  #foreach $ntasks (16, 32, 64)
                  foreach $ntasks (512)
                  {
                     #foreach $machine ("bgq", "peloton")
                     foreach $machine ("bgq")
                     {
                        printObject($anatomy,$celltype,$reaction,$fastgates,
                                    $rationalfns,$smoothing,$ntasks,$machine,$setgkr);
                     }
                  }
               }
            }
         }
      }
   }
}

# details of how to build object.data files for each set of parameters
sub printObject
{
   my($anatomy,$celltype,$reaction,$fastgates,$rationalfns,$smoothing,$ntasks,$machine,$setgkr) = @_;

   # skip file creation for conflicting parameter sets
   if ($reaction eq "TT06RRG" && !($fastgates == 0 && $smoothing == 0 && $rationalfns == 0)) { return; }
   if ($reaction eq "TT06" && !($fastgates == 0 && $smoothing == 0 && $rationalfns == 0)) { return; }
   if ($smoothing == 0 && $fastgates == 1) { return; }
   if ($rationalfns == 0 && $fastgates == 1) { return; }

   $nnodes = $ntasks;  # BGQ default
   if ($machine eq "peloton") { $nnodes = $ntasks/(16/$nthreadsPeloton); }
   if ($ntasks%$nnodes != 0) { $nnodes++; }
   $nthreads = $nthreadsBGQ;
   if ($machine eq "peloton") { $nthreads = $nthreadsPeloton; }

   $date = `date +%m%d%y`;  chomp $date;
   $maindir = join '','verif-runs-',$date;
   if ($weakScaling == 0) { $maindir = join '','verif-runs-strong-',$date; }
   $dirname = join '',$anatomy,'-',$celltype,'-gkr',$setgkr,'-N',$nnodes,'t',$nthreads;
   system("mkdir -p $maindir/$machine/$dirname");

# store different process grids in hashes
   $px{16} = 2;     $py{16} = 4;     $pz{16} = 2;  
   $px{32} = 4;     $py{32} = 4;     $pz{32} = 2;  
   $px{64} = 4;     $py{64} = 4;     $pz{64} = 4;  
   $px{128} = 4;    $py{128} = 8;    $pz{128} = 4;  
   $px{256} = 8;    $py{256} = 8;    $pz{256} = 4;  
   $px{512} = 8;    $py{512} = 8;    $pz{512} = 8;  
   $px{1024} = 8;   $py{1024} = 16;  $pz{1024} = 8;  
   $px{2048} = 16;  $py{2048} = 16;  $pz{2048} = 8;  
   $px{4096} = 16;  $py{4096} = 16;  $pz{4096} = 16;  
   $px{8192} = 16;  $py{8192} = 32;  $pz{8192} = 16;  
   $px{16384} = 32; $py{16384} = 32; $pz{16384} = 16;  
   $px{24576} = 32; $py{24576} = 32; $pz{24576} = 24;  
   $px{32768} = 32; $py{32768} = 32; $pz{32768} = 32;  
      
   if ($ntasks != $px{$ntasks}*$py{$ntasks}*$pz{$ntasks})
   {
      print "Process grid not defined correctly for ntasks = $ntasks:  px = $px{$ntasks}, py = $py, pz = $pz{$ntasks}\n";
      exit;
   }

   open OBJECT, ">$maindir/$machine/$dirname/object.data";

   # simulate block
   print OBJECT "simulate SIMULATE\n";
   print OBJECT "{\n";
   print OBJECT "   anatomy = $anatomy;\n";
   print OBJECT "   decomposition = grid;\n";
   print OBJECT "   diffusion = FGR;\n";
   print OBJECT "   reaction = $reaction;\n";
   print OBJECT "   stimulus = stimulus;\n";
   print OBJECT "   loop = 0;\n";
   print OBJECT "   maxLoop = $nIterations;\n";
   print OBJECT "   dt = 10 us;\n";
   print OBJECT "   time = 0;\n";
   print OBJECT "   printRate = $checkpointRate;\n";
   print OBJECT "   snapshotRate = $nIterations;\n";
   print OBJECT "   checkpointRate = $checkpointRate;\n";
   if ($reaction =~ /Opt/) 
   {
      print OBJECT "   parallelDiffusionReaction = 1;\n";
      print OBJECT "   nDiffusionCores = 2;\n";
   }
   print OBJECT "}\n\n";

   # anatomy and conductivity blocks
   if ($anatomy eq "block247")  # 247 cells per BG/Q core, 247*16 cells per MPI task
   {
      $nx = 16*$px{$ntasks}; $ny = 19*$py{$ntasks}; $nz = 13*$pz{$ntasks};
      if ($weakScaling == 0)
      {
         $nx = 16*$px{$strongTaskCount}; $ny = 19*$py{$strongTaskCount}; $nz = 13*$pz{$strongTaskCount};
      } 

      print OBJECT "$anatomy ANATOMY\n";
      print OBJECT "{\n";
      print OBJECT "   method = brick;\n";
      print OBJECT "   cellType = $celltype;\n";
      print OBJECT "   dx = 0.10;\n";
      print OBJECT "   dy = 0.10;\n";
      print OBJECT "   dz = 0.10;\n";
      print OBJECT "   nx = $nx;\n"; 
      print OBJECT "   ny = $ny;\n";
      print OBJECT "   nz = $nz;\n";
      print OBJECT "   conductivity = conductivity;\n";
      print OBJECT "}\n\n";
      print OBJECT "conductivity CONDUCTIVITY\n";
      print OBJECT "{\n";
      print OBJECT "    method = fibre;\n";
      print OBJECT "    sigmaLi = 0.0001334177;\n";
      print OBJECT "    sigmaTi = 0.0000176062;\n";
      print OBJECT "}\n\n";

   }
   elsif ($anatomy eq "swiss247")
   {
      $nvoids = 6;
      $nx = 16*$px{$ntasks}; $ny = 19*$py{$ntasks}; $nz = 13*$pz{$ntasks};
      if ($weakScaling == 0)
      {
         $nx = 16*$px{$strongTaskCount}; $ny = 19*$py{$strongTaskCount}; $nz = 13*$pz{$strongTaskCount};
      } 

      $anatdir = join '','anatomy-',$anatomy,'-',$celltype,'-n',$ntasks;
      $fullanat = join '',$maindir,'/',$anatdir;
      if (! -d $fullanat)
      {
         system("mkdir -p $fullanat");
         chdir($fullanat);
         if (!-e $makeVoidAnatomyScript)
         {
            print "$makeVoidAnatomyScript not found!\n";
            exit;
         }
         system("$makeVoidAnatomyScript $nx $ny $nz $celltype $nvoids");
         chdir("../..");
      }

      print OBJECT "$anatomy ANATOMY\n";
      print OBJECT "{\n";
      print OBJECT "   method = pio;\n";
      print OBJECT "   fileName = ..\/..\/$anatdir\/anatomy\#;\n";
      print OBJECT "   conductivity = conductivity;\n";
      print OBJECT "   dx = 0.10;\n";
      print OBJECT "   dy = 0.10;\n";
      print OBJECT "   dz = 0.10;\n";
      print OBJECT "}\n\n";
      print OBJECT "conductivity CONDUCTIVITY\n";
      print OBJECT "{\n";
      print OBJECT "    method = JHU;\n";
      print OBJECT "    sigmaLi = 0.0001334177;\n";
      print OBJECT "    sigmaTi = 0.0000176062;\n";
      print OBJECT "}\n\n";
      #print OBJECT "conductivity CONDUCTIVITY\n";
      #print OBJECT "{\n";
      #print OBJECT "   method = pio;\n";
      #print OBJECT "}\n\n";
   }

   # diffusion block
   print OBJECT "FGR DIFFUSION\n";
   print OBJECT "{\n";
   print OBJECT "   method = FGR;\n";
   print OBJECT "   diffusionScale = 714.2857143;\n";
   print OBJECT "}\n\n";

   # decomposition block
   print OBJECT "grid DECOMPOSITION \n";
   print OBJECT "{\n";
   print OBJECT "   method = grid;\n";
   print OBJECT "   nx = $px{$ntasks};\n";
   print OBJECT "   ny = $py{$ntasks};\n";
   print OBJECT "   nz = $pz{$ntasks};\n";
   print OBJECT "}\n\n";

   # reaction block
   print OBJECT "$reaction REACTION\n";
   print OBJECT "{\n";

   if ($reaction eq "TT06RRGOpt")
   {
      print OBJECT "   method = TT06Opt;\n";
      if ($rationalfns == 1) {
          print OBJECT "   tolerance = 0.0001;\n";
      }
      elsif ($rationalfns == 0) {
          print OBJECT "   tolerance = 0.0;\n";
      }
      else {
         print "Invalid choice of rationalfns = $rationalfns!\n";
         exit;
      }
      print OBJECT "   mod = $smoothing;\n";
      print OBJECT "   fastGate =$fastgates;\n"; 
      print OBJECT "   fastNonGate =$fastgates;\n";
      print OBJECT "   cellTypes = endo mid epi;\n";
      print OBJECT "}\n\n";
      if ($setgkr == 1)
      {
         print OBJECT "endo CELLTYPE { clone=endoRRG; g_Kr = 0.153;}\n";
         print OBJECT "mid CELLTYPE { clone=midRRG;  P_NaK=3.0; g_NaL=0.6; g_Kr = 0.153;}\n";
         print OBJECT "epi CELLTYPE { clone=epiRRG; g_Kr = 0.153; }\n\n";
      }
      elsif ($setgkr == 2)
      {
         print OBJECT "endo CELLTYPE { clone=endoRRG; g_Kr = 0.6;}\n";
         print OBJECT "mid CELLTYPE { clone=midRRG;  P_NaK=3.0; g_NaL=0.6; g_Kr = 0.4;}\n";
         print OBJECT "epi CELLTYPE { clone=epiRRG; g_Kr = 0.5; }\n\n";
      }
      else
      {
         print OBJECT "endo CELLTYPE { clone=endoRRG; }\n";
         print OBJECT "mid CELLTYPE { clone=midRRG;  P_NaK=3.0; g_NaL=0.6; }\n";
         print OBJECT "epi CELLTYPE { clone=epiRRG; }\n\n";
      }
   }
   elsif ($reaction eq "TT06RRG") 
   {
      print OBJECT "   method = TT06_RRG;\n";
      print OBJECT "}\n\n";
   }
   elsif ($reaction eq "TT06") 
   {
      print OBJECT "   method = TT06_CellML;\n";
      print OBJECT "   integrator = rushLarsen;\n";
      print OBJECT "}\n\n";
   }
   elsif ($reaction eq "TT06Opt")
   {
      print OBJECT "   method = TT06Opt;\n";
      if ($rationalfns == 1) {
          print OBJECT "   tolerance = 0.0001;\n";
      }
      elsif ($rationalfns == 0) {
          print OBJECT "   tolerance = 0.0;\n";
      }
      else {
         print "Invalid choice of rationalfns = $rationalfns!\n";
         exit;
      }
      print OBJECT "   mod = $smoothing;\n";
      print OBJECT "   fastGate =$fastgates;\n"; 
      print OBJECT "   fastNonGate =$fastgates;\n";
      print OBJECT "   cellTypes = endo mid epi;\n";
      print OBJECT "}\n\n";
      print OBJECT "endo CELLTYPE { clone=endoCellML; }\n";
      print OBJECT "mid CELLTYPE { clone=midCellML; }\n";
      print OBJECT "epi CELLTYPE { clone=epiCellML; }\n\n";
   }
   else
   {
      print "Reaction $reaction not defined in printObject!\n";
      exit;
   }

   # stimulus block
   if ($anatomy =~ /swiss/ || $anatomy =~ /block/)
   {
      print OBJECT "stimulus STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMax = 20;\n";
      print OBJECT "   yMax = 20;\n";
      print OBJECT "   zMax = 20;\n";
      print OBJECT "   vStim = -35.71429;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 1000;\n";
      print OBJECT "}\n\n";
   }
   else 
   {
      print "Stimulus not defined for anatomy $anatomy\n";
      exit;
   }
   close OBJECT;

   # print batch script
   if ($machine eq "bgq")
   {
      $bgqbatch = "sbatch.bgq";
      open BGQ, ">$maindir/$machine/$dirname/$bgqbatch";
      print BGQ "\#!/bin/bash\n";
      print BGQ "\#SBATCH --nodes=$nnodes\n";
      print BGQ "\n";
      print BGQ "export OMP_NUM_THREADS=$nthreadsBGQ\n";
      print BGQ "export MUSPI_NUMBATIDS=203\n";
      print BGQ "export MUSPI_NUMINJFIFOS=3\n";
      print BGQ "export MUSPI_NUMRECFIFOS=3\n\n";
      print BGQ "srun --ntasks=$ntasks $bgqExe\n";
      close BGQ;
   }
   elsif ($machine eq "peloton")
   {
      $pelbatch = "msub.pel";
      open PEL, ">$maindir/$machine/$dirname/$pelbatch";
      print PEL "\#!/bin/bash\n";
      print PEL "\#MSUB -l nodes=$nnodes\n";
      print PEL "\#MSUB -l walltime=12:00:00\n";
      print PEL "\#MSUB -A gbcq\n";
      print PEL "\n";
      print PEL "export OMP_NUM_THREADS=$nthreadsPeloton\n";
      print PEL "srun -n $ntasks $pelotonExe\n";
      close PEL;

      $peldebug = "rundebug.pel";
      open DEB, ">$maindir/$machine/$dirname/$peldebug";
      print DEB "\#!/bin/bash\n";
      print DEB "export OMP_NUM_THREADS=$nthreadsPeloton\n";
      print DEB "srun -N $nnodes -n $ntasks -p pdebug $pelotonExe object.data > slurm.out\n";   
      close DEB;
      system("chmod u+x $maindir/$machine/$dirname/$peldebug");
   }
}
