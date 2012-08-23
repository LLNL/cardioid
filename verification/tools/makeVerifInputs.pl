#! /usr/bin/perl
#
# generate Cardioid object.data input files for verification runs.  Compare output data 
# using compareSnapshots.cc (in Cardioid src directory) 
#
# written by Erik Draeger, LLNL, 8/20/2012

# global variables

$thisdir = `pwd`;  chomp $thisdir;
$makeVoidAnatomyScript = "$thisdir/tools/makeAnatomyBlockWithVoids.pl";
$bgqExe = "../../../bin/cardioid-bgq-spi";
$pelotonExe = "../../../bin/cardioid-peloton";

$nIterations = 100000;
$checkpointRate = 1000;

foreach $anatomy ("swiss247", "block247")
{
   foreach $celltype ("100", "101", "102", "random")
   {
      # these correspond to if blocks in printObject (below), not cardioid reaction types
      foreach $reaction ("TT06RRG", "TT06RRGOpt", "TT06", "TT06Opt") 
      {
         foreach $fastgates (1, 0)
         {
            foreach $rationalfns (1, 0)
            {
               foreach $ntasks (16, 32, 64)
               {
                  printObject($anatomy,$celltype,$reaction,$fastgates,$rationalfns,$ntasks);
               }
            }
         }
      }
   }
}

# details of how to build object.data files for each set of parameters
sub printObject
{
   my($anatomy,$celltype,$reaction,$fastgates,$rationalfns,$ntasks) = @_;

   # skip file creation for conflicting parameter sets
   if ($reaction eq "TT06RRG" && !($fastgates == 0 && $rationalfns == 0)) { return; }
   if ($reaction eq "TT06" && !($fastgates == 0 && $rationalfns == 0)) { return; }
   if ($rationalfns == 0 && $fastgates == 1) { return; }

   $date = `date +%m%d%y`;  chomp $date;
   $maindir = join '','verif-runs-',$date;
   $dirname = join '',$anatomy,'-',$celltype,'-',$reaction,'-','fast',$fastgates,'mod',$rationalfns,'-N',$ntasks;
   system("mkdir -p $maindir/$dirname");

   if ($ntasks == 16) { $px = 2; $py = 4; $pz = 2; }
   elsif ($ntasks == 32) { $px = 4; $py = 4; $pz = 2; }
   elsif ($ntasks == 64) { $px = 4; $py = 4; $pz = 4; }
   elsif ($ntasks == 128) { $px = 4; $py = 8; $pz = 4; }
   elsif ($ntasks == 256) { $px = 8; $py = 8; $pz = 4; }
   elsif ($ntasks == 512) { $px = 8; $py = 8; $pz = 8; }
   elsif ($ntasks == 1024) { $px = 8; $py = 16; $pz = 8; }
   elsif ($ntasks == 2048) { $px = 16; $py = 16; $pz = 8; }
   elsif ($ntasks == 4096) { $px = 16; $py = 16; $pz = 16; }
   elsif ($ntasks == 8192) { $px = 16; $py = 32; $pz = 16; }
   elsif ($ntasks == 16384) { $px = 32; $py = 32; $pz = 16; }
   elsif ($ntasks == 24576) { $px = 32; $py = 32; $pz = 24; }
   elsif ($ntasks == 32768) { $px = 32; $py = 32; $pz = 32; }
      
   if ($ntasks != $px*$py*$pz)
   {
      print "Process grid not defined correctly for ntasks = $ntasks:  px = $px, py = $py, pz = $pz\n";
      exit;
   }

   open OBJECT, ">$maindir/$dirname/object.data";

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
      #print OBJECT "   nDiffusionCores = 2;\n";
   }
   print OBJECT "}\n\n";

   # anatomy and conductivity blocks
   if ($anatomy eq "block247")  # 247 cells per BG/Q core, 247*16 cells per MPI task
   {
      $nx = 16*$px; $ny = 19*$py; $nz = 13*$pz;

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
      $nx = 16*$px; $ny = 19*$py; $nz = 13*$pz;

      $anatdir = join '','anatomy-',$anatomy,'-',$celltype,'-N',$ntasks;
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
      print OBJECT "   fileName = ..\/$anatdir\/anatomy\#;\n";
      print OBJECT "   conductivity = conductivity;\n";
      print OBJECT "   dx = 0.10;\n";
      print OBJECT "   dy = 0.10;\n";
      print OBJECT "   dz = 0.10;\n";
      print OBJECT "}\n\n";
      print OBJECT "conductivity CONDUCTIVITY\n";
      print OBJECT "{\n";
      print OBJECT "   method = pio;\n";
      print OBJECT "}\n\n";
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
   print OBJECT "   nx = $px;\n";
   print OBJECT "   ny = $py;\n";
   print OBJECT "   nz = $pz;\n";
   print OBJECT "}\n\n";

   # reaction block
   print OBJECT "$reaction REACTION\n";
   print OBJECT "{\n";

   if ($reaction eq "TT06RRGOpt")
   {
      print OBJECT "   method = TT06Opt;\n";
      print OBJECT "   tolerance = 0.0001;\n";
      print OBJECT "   mod = $rationalfns;\n";
      print OBJECT "   fastGate =$fastgates;\n"; 
      print OBJECT "   fastNonGate =$fastgates;\n";
      print OBJECT "   cellTypes = endo mid epi;\n";
      print OBJECT "}\n\n";
      print OBJECT "endo CELLTYPE { clone=endoRRG; }\n";
      print OBJECT "mid CELLTYPE { clone=midRRG; }\n";
      print OBJECT "epi CELLTYPE { clone=epiRRG; }\n\n";
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
      print OBJECT "   tolerance = 0.0001;\n";
      print OBJECT "   mod = $rationalfns;\n";
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
      print OBJECT "   duration = 2;\n";
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
   $bgqbatch = "sbatch.bgq";
   open BGQ, ">$maindir/$dirname/$bgqbatch";
   print BGQ "\#!/bin/bash\n";
   print BGQ "\#SBATCH --nodes=$ntasks\n";
   print BGQ "\n";
   print BGQ "export OMP_NUM_THREADS=64\n";
   print BGQ "export MUSPI_NUMBATIDS=203\n";
   print BGQ "export MUSPI_NUMINJFIFOS=3\n";
   print BGQ "export MUSPI_NUMRECFIFOS=3\n\n";
   print BGQ "srun --ntasks=$ntasks $bgqExe\n";
   close BGQ;

   $pelbatch = "msub.pel";
   $nnodes = $ntasks/8;
   if ($ntasks%$nnodes != 0) { $nnodes++; }
   open PEL, ">$maindir/$dirname/$pelbatch";
   print PEL "\#!/bin/bash\n";
   print PEL "\#MSUB -l nodes=$nnodes\n";
   print PEL "\#MSUB -l walltime=12:00:00\n";
   print PEL "\n";
   print PEL "export OMP_NUM_THREADS=2\n";
   print PEL "srun -n $ntasks $pelotonExe\n";
   close PEL;

   $peldebug = "rundebug.pel";
   $nnodes = $ntasks/8;
   if ($ntasks%$nnodes != 0) { $nnodes++; }
   open DEB, ">$maindir/$dirname/$peldebug";
   print DEB "\#!/bin/bash\n";
   print DEB "export OMP_NUM_THREADS=2\n";
   print DEB "srun -N $nnodes -n $ntasks -p pdebug $pelotonExe object.data > slurm.out\n";   
   close DEB;
   system("chmod u+x $maindir/$dirname/$peldebug");
}

