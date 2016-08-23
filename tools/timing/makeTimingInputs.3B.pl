#! /usr/bin/perl -w
#
# generate Cardioid object.data input files for timing runs.
#
# written by Erik Draeger, LLNL, 10/15/2012

# global variables

$AnatomyDir370M = "/p/ls1/emhm/370M/anatomy-19-Jun-2012";
$AnatomyDir3B = "/p/ls1/emhm/3B_human_0.05mm";

###DEBUG:  run out of /nfs/tmp2/emhm until Lustre comes back
#$AnatomyDir370M ="/nfs/tmp2/emhm/370M-19-Jun-2012/anatomy.370M";
#$AnatomyDir3B = "/nfs/tmp2/emhm/vgurev/human_05";



$thisdir = `pwd`;  chomp $thisdir;
$makeVoidAnatomyScript = "$thisdir/tools/makeAnatomyBlockWithVoids.pl";
$bgqSPIExe = "../../../bin/cardioid-bgq-spi";
$bgqHPMExe = "../../../bin/cardioid-bgq-spi-hpm";
$pelotonExe = "../../../bin/cardioid-peloton";
$nthreadsBGQ = 64;
$nthreadsPeloton = 4;

$useStateSensor = 0;  # add code for the new stateVariables sensor

$nIterations = 1000000;
$printRate = $nIterations/1000;
$checkpointRate = 50000000;

$weakScaling = 1;   # if set to zero, anatomy size corresponding 
                    # to $strongTaskCount will be used throughout
$strongTaskCount = 64;

$celltype = "undefined";  # placeholder

$anatomy = "370M";
#foreach $anatomy ("370M", "3B")
{
   $reaction = "TT06RRGOpt";
   $fastgates = 1;
   $rationalfns = 1;
   $smoothing = 1;
   foreach $ntasks (1024, 2048, 4096)
   #foreach $ntasks (1024, 2048, 4096, 8192, 16384, 24576, 32768, 49152, 73728, 98304)
   {
      $machine = "bgqhpm";
      #foreach $machine ("bgqspi", "bgqhpm")
      {
         $loadbal = "grid";
         #foreach $loadbal ("grid", "work")
         {
            $globalsync = -1;
            #foreach $globalsync (-1, 0, 1)
            {
               printObject($anatomy,$celltype,$reaction,$fastgates,
                           $rationalfns,$smoothing,$ntasks,$machine);
            }
         }
      }
   }
}

# details of how to build object.data files for each set of parameters
sub printObject
{
   my($anatomy,$celltype,$reaction,$fastgates,$rationalfns,$smoothing,$ntasks,$machine) = @_;

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
   $maindir = join '','timing-runs-',$date;
   if ($weakScaling == 0) { $maindir = join '','verif-runs-strong-',$date; }
   #$dirname = join '',$anatomy,'-',$reaction,'-fast1-',$machine,'-',$loadbal,'bal','-N',$nnodes;
   if ($globalsync > 0)
   {
      $dirname = join '',$anatomy,'-',$reaction,'-fast1-',$machine,'-globsync',$globalsync,'-',$loadbal,'bal','-N',$nnodes;
   }
   elsif ($globalsync == 0)
   {
      $dirname = join '',$anatomy,'-',$reaction,'-fast1-',$machine,'-globsync-once-',$loadbal,'bal','-N',$nnodes;
   }
   else
   {
      $dirname = join '',$anatomy,'-',$reaction,'-fast1-',$machine,'-',$loadbal,'bal','-N',$nnodes;
   }
   system("mkdir -p $maindir/$dirname");

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
   $px{36864} = 32; $py{36864} = 36; $pz{36864} = 32;  
   $px{49152} = 32; $py{49152} = 32; $pz{49152} = 48;  
   $px{73728} = 32; $py{73728} = 64; $pz{73728} = 36;  
   #$px{98304} = 32; $py{98304} = 64; $pz{98304} = 48;  
   $px{98304} = 48; $py{98304} = 32; $pz{98304} = 64;  
      
# store workbound balancer block sizes for each task count
   if ($anatomy eq "3B")
   {
      $wx{1024} = 92;   $wy{1024} = 88;  $wz{1024} = 398;  
      $wx{2048} = 46;  $wy{2048} = 88;  $wz{2048} = 398;  
      $wx{4096} = 46;  $wy{4096} = 44;  $wz{4096} = 398;  
      $wx{8192} = 46;  $wy{8192} = 44;  $wz{8192} = 198;  
      $wx{16384} = 46; $wy{16384} = 44; $wz{16384} = 98;  
      $wx{24576} = 38; $wy{24576} = 36; $wz{24576} = 98;  
      $wx{32768} = 46; $wy{32768} = 22; $wz{32768} = 98;  
      $wx{36864} = 30; $wy{36864} = 30; $wz{36864} = 98;  
      $wx{49152} = 26; $wy{49152} = 26; $wz{49152} = 98;  
      $wx{73728} = 21; $wy{73728} = 22; $wz{73728} = 98;  
      $wx{98304} = 48; $wy{98304} = 48; $wz{98304} = 14;  
   }
   elsif ($anatomy eq "370M")
   {
      $wx{1024} = 32;   $wy{1024} = 24;  $wz{1024} = 398;
      $wx{2048} = 32;  $wy{2048} = 24;  $wz{2048} = 198;
      $wx{4096} = 32;  $wy{4096} = 24;  $wz{4096} = 98;
      $wx{8192} = 35;  $wy{8192} = 34;  $wz{8192} = 42;
      $wx{16384} = 18; $wy{16384} = 20; $wz{16384} = 70;
      $wx{24576} = 20; $wy{24576} = 20; $wz{24576} = 42;
      $wx{32768} = 16; $wy{32768} = 18; $wz{32768} = 42;
      $wx{36864} = 20; $wy{36864} = 19; $wz{36864} = 26;
      $wx{49152} = 20; $wy{49152} = 19; $wz{49152} = 22;
      $wx{73728} = 20; $wy{73728} = 19; $wz{73728} = 14;
      $wx{98304} = 18; $wy{98304} = 16; $wz{98304} = 14;  
   }

   if ($ntasks != $px{$ntasks}*$py{$ntasks}*$pz{$ntasks})
   {
      print "Process grid not defined correctly for ntasks = $ntasks:  px = $px{$ntasks}, py = $py, pz = $pz{$ntasks}\n";
      exit;
   }

   open OBJECT, ">$maindir/$dirname/object.data";

   # simulate block
   print OBJECT "simulate SIMULATE\n";
   print OBJECT "{\n";
   print OBJECT "   anatomy = $anatomy;\n";   
   print OBJECT "   decomposition = $loadbal;\n";
   print OBJECT "   diffusion = FGR;\n";
   print OBJECT "   reaction = $reaction;\n";
   if ($anatomy eq "3B")
   {
      print OBJECT "   stimulus = box0 box1 box2 box3 box4 box5 box6 box7 box8 box10 box11 box12;
\n";
   }
   elsif ($anatomy eq "370M")
   {
      print OBJECT "   stimulus = box0 box1 box2 box3 box4 box5 box6 box7 box8 box10 box11;\n";
   }
   else {
      print OBJECT "   stimulus = stimulus;\n";
   }
   print OBJECT "   loop = 0;\n";
   print OBJECT "   maxLoop = $nIterations;\n";
   if ($anatomy eq "3B")
   {
      print OBJECT "   dt = 4 us;\n";
   }
   else 
   {
      print OBJECT "   dt = 20 us;\n";
   }
   print OBJECT "   time = 0;\n";
   print OBJECT "   checkRanges = 0;\n";
   print OBJECT "   writeTorusMap = 0;\n";
   if ($globalsync >= 0)
   {
      print OBJECT "   globalSyncRate = $globalsync;\n";
   }
   print OBJECT "   printRate = $printRate;\n";
   print OBJECT "   snapshotRate = $checkpointRate;\n";
   print OBJECT "   checkpointRate = $checkpointRate;\n";
   print OBJECT "   sensor = ecg;\n";
   if ($useStateSensor == 1)
   {
      print OBJECT "   sensor = stateVariable;\n";
   }
   if ($reaction =~ /Opt/) 
   {
      #print OBJECT "   parallelDiffusionReaction = 1;\n";
      print OBJECT "   loopType = pdr;\n";
      if ($loadbal eq "grid")
      {
         # comment this out for now, let assignCellsToTasks try to figure it out
         #print OBJECT "   nDiffusionCores = 3;\n";
         #print OBJECT "   nDiffusionCores = 4;\n";
      }
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
      print OBJECT "   method = pio;\n";
      print OBJECT "}\n\n";
   }
   elsif ($anatomy eq "370M")
   {
      print OBJECT "$anatomy ANATOMY\n";
      print OBJECT "{\n";
      print OBJECT "   method = pio;\n";
      print OBJECT "   fileName = $AnatomyDir370M\/anatomy\#;\n";
      print OBJECT "   conductivity = conductivity;\n";
      print OBJECT "   dx = 0.100966;\n";
      print OBJECT "   dy = 0.102087;\n";
      print OBJECT "   dz = 0.100805;\n";
      print OBJECT "}\n\n";
      print OBJECT "conductivity CONDUCTIVITY\n";
      print OBJECT "{\n";
      print OBJECT "   method = pio;\n";
      print OBJECT "}\n\n";
   }
   elsif ($anatomy eq "3B")
   {
      print OBJECT "$anatomy ANATOMY\n";
      print OBJECT "{\n";
      print OBJECT "   method = pio;\n";
      print OBJECT "   fileName = $AnatomyDir3B\/anatomy\#;\n";
      print OBJECT "   conductivity = conductivity;\n";
      print OBJECT "   dx = 0.0504831;\n";
      print OBJECT "   dy = 0.0500896;\n";
      print OBJECT "   dz = 0.0504024;\n";
      print OBJECT "}\n\n";
      print OBJECT "conductivity CONDUCTIVITY\n";
      print OBJECT "{\n";
      print OBJECT "   method = pio;\n";
      print OBJECT "}\n\n";
   }

   # ecg block
   print OBJECT "ecg SENSOR\n";
   print OBJECT "{\n";
   print OBJECT "   method = ECG;\n";
   print OBJECT "   evalRate =  500;\n";
   print OBJECT "   printRate = $nIterations;\n";
   print OBJECT "   stencilSize = 9;\n";
   print OBJECT "   nSensorPoints = 6;\n";
   print OBJECT "   nFiles = 3072;\n";
   print OBJECT "}\n\n";

   # diffusion block
   print OBJECT "FGR DIFFUSION\n";
   print OBJECT "{\n";
   print OBJECT "   method = FGR;\n";
   print OBJECT "   diffusionScale = 714.2857143;\n";
   print OBJECT "   variant = strip;\n";
   print OBJECT "}\n\n";

   # decomposition block
   if ($loadbal eq "grid")
   {
      print OBJECT "$loadbal DECOMPOSITION \n";
      print OBJECT "{\n";
      print OBJECT "   method = grid;\n";
      print OBJECT "   nx = $px{$ntasks};\n";
      print OBJECT "   ny = $py{$ntasks};\n";
      print OBJECT "   nz = $pz{$ntasks};\n";
      print OBJECT "   diffCost = 0.5;\n";
      print OBJECT "   nDiffCores = 2;\n";
      print OBJECT "   printTaskInfo = 1;\n";
      print OBJECT "}\n\n";
   }
   elsif ($loadbal eq "work")
   {
      print OBJECT "$loadbal DECOMPOSITION \n";
      print OBJECT "{\n";
      print OBJECT "   method = workBound;\n";
      print OBJECT "   dx = $wx{$ntasks};\n";
      print OBJECT "   dy = $wy{$ntasks};\n";
      print OBJECT "   dz = $wz{$ntasks};\n";
      print OBJECT "   nCores = 16;\n";
      print OBJECT "   nRCoresBB = 2;\n";
      print OBJECT "   alpha = 0.1;\n";
      print OBJECT "   beta = 100.0;\n";
      print OBJECT "   printStats = 0;\n";
      print OBJECT "   printTaskInfo = 0;\n";
      print OBJECT "}\n\n";      
   }
   else {
      print "loadbal = $loadbal undefined!\n";
      exit;
   }

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
      print OBJECT "   fitFile = fit.data;\n";
      print OBJECT "   fastGate =$fastgates;\n"; 
      print OBJECT "   fastNonGate =$fastgates;\n";
      print OBJECT "   cellTypes = endo mid epi;\n";
      print OBJECT "}\n\n";
      # hack for this version of the code only
      #print OBJECT "endo CELLTYPE { clone=endoRRG; }\n";
      #print OBJECT "mid CELLTYPE { clone=midRRG;  P_NaK=3.0; g_NaL=0.6; }\n";
      print OBJECT "endo CELLTYPE { clone=endoRRG; anatomyIndices=30 31;}\n";
      print OBJECT "mid CELLTYPE { clone=midRRG;  P_NaK=3.0; g_NaL=0.6;anatomyIndices=76 100 101; }\n";
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
   elsif ($anatomy eq "370M")
   {
      print OBJECT "box0 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 452;\n";
      print OBJECT "   xMax = 472;\n";
      print OBJECT "   yMin = 324;\n";
      print OBJECT "   yMax = 344;\n";
      print OBJECT "   zMin = 510;\n";
      print OBJECT "   zMax = 530;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n\n";
      print OBJECT "box1 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 476;\n";
      print OBJECT "   xMax = 496;\n";
      print OBJECT "   yMin = 315;\n";
      print OBJECT "   yMax = 335;\n";
      print OBJECT "   zMin = 605;\n";
      print OBJECT "   zMax = 625;\n";
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
      print OBJECT "box2 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 411;\n";
      print OBJECT "   xMax = 431;\n";
      print OBJECT "   yMin = 318;\n";
      print OBJECT "   yMax = 338;\n"; 
      print OBJECT "   zMin = 406;\n";
      print OBJECT "   zMax = 426;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
      print OBJECT "box3 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 556;\n";
      print OBJECT "   xMax = 576;\n";
      print OBJECT "   yMin = 424;\n";
      print OBJECT "   yMax = 444;\n"; 
      print OBJECT "   zMin = 479;\n";
      print OBJECT "   zMax = 499;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
      print OBJECT "\n";
      print OBJECT "box4 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 624;\n";
      print OBJECT "   xMax = 644;\n";
      print OBJECT "   yMin = 388;\n";
      print OBJECT "   yMax = 408;\n"; 
      print OBJECT "   zMin = 520;\n";
      print OBJECT "   zMax = 540;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
      print OBJECT "box5 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 536;\n";
      print OBJECT "   xMax = 556;\n";
      print OBJECT "   yMin = 413;\n";
      print OBJECT "   yMax = 433;\n"; 
      print OBJECT "   zMin = 367;\n";
      print OBJECT "   zMax = 387;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
      print OBJECT "\n";
      print OBJECT "box6 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 224;\n";
      print OBJECT "   xMax = 244;\n";
      print OBJECT "   yMin = 560;\n";
      print OBJECT "   yMax = 580;\n"; 
      print OBJECT "   zMin = 675;\n";
      print OBJECT "   zMax = 695;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
      print OBJECT "box7 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 232;\n";
      print OBJECT "   xMax = 252;\n";
      print OBJECT "   yMin = 453;\n";
      print OBJECT "   yMax = 473;\n"; 
      print OBJECT "   zMin = 714;\n";
      print OBJECT "   zMax = 734;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
      print OBJECT "box8 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 414;\n";
      print OBJECT "   xMax = 434;\n";
      print OBJECT "   yMin = 501;\n";
      print OBJECT "   yMax = 521;\n"; 
      print OBJECT "   zMin = 590;\n";
      print OBJECT "   zMax = 610;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
      print OBJECT "box9 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 33;\n";
      print OBJECT "   xMax = 53;\n";
      print OBJECT "   yMin = 561;\n";
      print OBJECT "   yMax = 581;\n"; 
      print OBJECT "   zMin = 579;\n";
      print OBJECT "   zMax = 599;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
      print OBJECT "\n";
      print OBJECT "box10 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 123;\n";
      print OBJECT "   xMax = 143;\n";
      print OBJECT "   yMin = 760;\n";
      print OBJECT "   yMax = 780;\n"; 
      print OBJECT "   zMin = 509;\n";
      print OBJECT "   zMax = 529;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
      print OBJECT "box11 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 904;\n";
      print OBJECT "   xMax = 944;\n";
      print OBJECT "   yMin = 657;\n";
      print OBJECT "   yMax = 698;\n"; 
      print OBJECT "   zMin = 1020;\n";
      print OBJECT "   zMax = 1060;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
   }
   elsif ($anatomy eq "3B")
   {
      print OBJECT "box0 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 904;\n";
      print OBJECT "   xMax = 944;\n";
      print OBJECT "   yMin = 657;\n";
      print OBJECT "   yMax = 698;\n";
      print OBJECT "   zMin = 1020;\n";
      print OBJECT "   zMax = 1060;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n\n";
      print OBJECT "box1 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 952;\n";
      print OBJECT "   xMax = 992;\n";
      print OBJECT "   yMin = 639;\n";
      print OBJECT "   yMax = 680;\n";
      print OBJECT "   zMin = 1210;\n";
      print OBJECT "   zMax = 1250;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n\n";
      print OBJECT "box2 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 822;\n";
      print OBJECT "   xMax = 862;\n";
      print OBJECT "   yMin = 645;\n";
      print OBJECT "   yMax = 686;\n";
      print OBJECT "   zMin = 812;\n";
      print OBJECT "   zMax = 852;\n";
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n\n";
      print OBJECT "box3 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 1112;\n";
      print OBJECT "   xMax = 1152;\n";
      print OBJECT "   yMin = 860;\n";
      print OBJECT "   yMax = 901;\n"; 
      print OBJECT "   zMin = 958;\n";
      print OBJECT "   zMax = 998;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
      print OBJECT "\n";
      print OBJECT "box4 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 1248;\n";
      print OBJECT "   xMax = 1288;\n";
      print OBJECT "   yMin = 787;\n";
      print OBJECT "   yMax = 828;\n"; 
      print OBJECT "   zMin = 1040;\n";
      print OBJECT "   zMax = 1080;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
      print OBJECT "\n";
      print OBJECT "box5 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 1072;\n";
      print OBJECT "   xMax = 1112;\n";
      print OBJECT "   yMin = 838;\n";
      print OBJECT "   yMax = 878;\n"; 
      print OBJECT "   zMin = 734;\n";
      print OBJECT "   zMax = 774;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
      print OBJECT "\n";
      print OBJECT "box6 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 448;\n";
      print OBJECT "   xMax = 488;\n";
      print OBJECT "   yMin = 1136;\n";
      print OBJECT "   yMax = 1177;\n"; 
      print OBJECT "   zMin = 1350;\n";
      print OBJECT "   zMax = 1390;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
      print OBJECT "\n";
      print OBJECT "box7 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 464;\n";
      print OBJECT "   xMax = 504;\n";
      print OBJECT "   yMin = 919;\n";
      print OBJECT "   yMax = 960;\n"; 
      print OBJECT "   zMin = 1428;\n";
      print OBJECT "   zMax = 1468;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
      print OBJECT "\n";
      print OBJECT "box8 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 828;\n";
      print OBJECT "   xMax = 868;\n";
      print OBJECT "   yMin = 1017;\n";
      print OBJECT "   yMax = 1057;\n"; 
      print OBJECT "   zMin = 1180;\n";
      print OBJECT "   zMax = 1220;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
      print OBJECT "\n";
      print OBJECT "box9 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 66;\n";
      print OBJECT "   xMax = 106;\n";
      print OBJECT "   yMin = 1138;\n";
      print OBJECT "   yMax = 1179;\n"; 
      print OBJECT "   zMin = 1158;\n";
      print OBJECT "   zMax = 1198;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
      print OBJECT "\n";
      print OBJECT "box10 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 246;\n";
      print OBJECT "   xMax = 286;\n";
      print OBJECT "   yMin = 1542;\n";
      print OBJECT "   yMax = 1583;\n"; 
      print OBJECT "   zMin = 1018;\n";
      print OBJECT "   zMax = 1058;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n\n";
      print OBJECT "box11 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 728;\n";
      print OBJECT "   xMax = 768;\n";
      print OBJECT "   yMin = 881;\n";
      print OBJECT "   yMax = 921;\n"; 
      print OBJECT "   zMin = 392;\n";
      print OBJECT "   zMax = 432;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
      print OBJECT "\n";
      print OBJECT "box12 STIMULUS\n";
      print OBJECT "{\n";
      print OBJECT "   method = box;\n";
      print OBJECT "   xMin = 776;\n";
      print OBJECT "   xMax = 816;\n";
      print OBJECT "   yMin = 13;\n";
      print OBJECT "   yMax = 53;\n"; 
      print OBJECT "   zMin = 998;\n";
      print OBJECT "   zMax = 1038;\n"; 
      print OBJECT "   vStim = -36 mV/ms;\n";
      print OBJECT "   tStart = 0;\n";
      print OBJECT "   duration = 1;\n";
      print OBJECT "   period = 2000;\n";
      print OBJECT "}\n";
   }
   else 
   {
      print "Stimulus not defined for anatomy $anatomy\n";
      exit;
   }

   if ($useStateSensor == 1)
   {
      print OBJECT "stateVariable SENSOR\n";
      print OBJECT "{\n";
      print OBJECT "   gid = 727354661;\n";
      print OBJECT "   radius = 4.0;\n";
      print OBJECT "   method = stateVariables;\n";
      print OBJECT "   fields = all;\n";
      print OBJECT "   startTime = 0.0;\n";
      print OBJECT "   endTime = 500.0;\n";
      print OBJECT "   dirname = sensorData;\n";
      print OBJECT "   printRate = 1;\n";
      print OBJECT "}\n\n";
   }

   close OBJECT;

   # print batch script
   if ($machine eq "bgqspi")
   {
      $bgqbatch = "sbatch.bgq";
      open BGQ, ">$maindir/$dirname/$bgqbatch";
      print BGQ "\#!/bin/bash\n";
      print BGQ "\#SBATCH --nodes=$nnodes\n";
      print BGQ "\n";
      print BGQ "export OMP_NUM_THREADS=$nthreadsBGQ\n";
      print BGQ "export MUSPI_NUMBATIDS=203\n";
      print BGQ "export MUSPI_NUMINJFIFOS=3\n";
      print BGQ "export MUSPI_NUMRECFIFOS=3\n";
      print BGQ "export MUSPI_NUMCLASSROUTES=3\n\n";
      print BGQ "srun --ntasks=$ntasks $bgqSPIExe\n";
      close BGQ;
   }
   # print batch script
   elsif ($machine eq "bgqhpm")
   {
      $bgqbatch = "sbatch.bgq";
      open BGQ, ">$maindir/$dirname/$bgqbatch";
      print BGQ "\#!/bin/bash\n";
      print BGQ "\#SBATCH --nodes=$nnodes\n";
      print BGQ "\n";
      print BGQ "export OMP_NUM_THREADS=$nthreadsBGQ\n";
      print BGQ "export MUSPI_NUMBATIDS=203\n";
      print BGQ "export MUSPI_NUMINJFIFOS=3\n";
      print BGQ "export MUSPI_NUMRECFIFOS=3\n";
      print BGQ "export MUSPI_NUMCLASSROUTES=3\n\n";
      print BGQ "srun --ntasks=$ntasks $bgqHPMExe\n";
      close BGQ;
   }
   elsif ($machine eq "peloton")
   {
      $pelbatch = "msub.pel";
      open PEL, ">$maindir/$dirname/$pelbatch";
      print PEL "\#!/bin/bash\n";
      print PEL "\#MSUB -l nodes=$nnodes\n";
      print PEL "\#MSUB -l walltime=12:00:00\n";
      print PEL "\#MSUB -A gbcq\n";
      print PEL "\n";
      print PEL "export OMP_NUM_THREADS=$nthreadsPeloton\n";
      print PEL "srun -n $ntasks $pelotonExe\n";
      close PEL;

      $peldebug = "rundebug.pel";
      open DEB, ">$maindir/$dirname/$peldebug";
      print DEB "\#!/bin/bash\n";
      print DEB "export OMP_NUM_THREADS=$nthreadsPeloton\n";
      print DEB "srun -N $nnodes -n $ntasks -p pdebug $pelotonExe object.data > slurm.out\n";   
      close DEB;
      system("chmod u+x $maindir/$dirname/$peldebug");
   }
}
