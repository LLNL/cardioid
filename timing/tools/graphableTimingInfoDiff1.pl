#! /usr/bin/perl -w
#
# read taskLoadInfo.data file and time.snapshot.0000xxxxxxxxx.info file (output from profileTimingInfo.pl)
# and compile timing data into graphable format.
#
# written by Erik Draeger, LLNL, 10/8/2012


if ($#ARGV < 0) {
   print "syntax:  graphableTimingInfoDiff1.pl [time.snapshot.diff1.info file(s)]\n";
}

$taskInfoFile = "taskLoadInfo.data";
if (! -e $taskInfoFile)
{
   print "Task info file $taskInfoFile not found!\n";
   exit;
}

foreach $timefile (@ARGV[0..$#ARGV]) 
{
   open LOAD, $taskInfoFile;
   open TIME, $timefile;

   $filebase_diffcalc = join '',$timefile,'.diffcalc';
   $filebase_diffstencil = join '',$timefile,'.diffstencil';
   $filebase_diffhalo = join '',$timefile,'.diffhalo';
   $filebase_diffwait = join '',$timefile,'.diffwait';

   $filename_difftot = join '',$timefile,'.difftot.dat';
   $filename_reactot = join '',$timefile,'.reactot.dat';


   open OUTD, ">$filename_difftot";
   open OUTR, ">$filename_reactot";

   # make array of filehandles for each number of threads
   #for ($tt=4; $tt<64; $tt+=4)
   #{
   #   $ttp = sprintf("%0.2i",$tt);
   #   $tfile1 = join '',$timefile,'.g1.t',$ttp,'.dat';
   #   open($fileh1[$tt],">$tfile1");

   #   $tfile2 = join '',$timefile,'.g2.t',$ttp,'.dat';
   #   open($fileh2[$tt],">$tfile2");
   #}

   $tmp = <LOAD>;  # header line
   $tmp = <TIME>;  # header line

   $linecnt = 0;
   while (($timeline = <TIME>) && ($loadline = <LOAD>))
   {
      @timedata = split ' ',$timeline;
      @loaddata = split ' ',$loadline;

      $timepe = $timedata[0];
      $loadpe = $loaddata[0];
      if ($timepe != $loadpe) 
      {
         print "File mismatch:  different pes on line $linecnt!  $timepe, $loadpe\n";
         exit;
      }

      $nthreadsDiff = $timedata[1];
      $nthreadsReac = $timedata[2];
      $timeDiff = $timedata[3];
      $timeDiffWait = $timedata[4];
      $timeDiffCalc = $timedata[5];
      $timeDiffStencil = $timedata[6];
      $timeDiffHalo = $timedata[7];
      $timeReac = $timedata[8];
      $timeReacWait = $timedata[9];

      $ncells = $loaddata[1];
      $volume = $loaddata[2];
      $density = $ncells/$volume;
      $bb_xmin = $loaddata[3];
      $bb_xmax = $loaddata[4];
      $bb_ymin = $loaddata[5];
      $bb_ymax = $loaddata[6];
      $bb_zmin = $loaddata[7];
      $bb_zmax = $loaddata[8];


      $rwork = $timeReac*$nthreadsReac;
      print OUTR "   $ncells    $rwork   $timeReac   $nthreadsReac  $volume  $density   $timepe\n";

      if (!defined($fileh1[$nthreadsReac])) {
         $ttp = sprintf("%0.2i",$nthreadsReac);
         $tfile1 = join '',$timefile,'.reac.t',$ttp,'.dat';
         open($fileh1[$nthreadsReac],">$tfile1");
      }
      print { $fileh1[$nthreadsReac] } "  $ncells    $timeReac   $rwork   $nthreadsReac  $volume  $density     $timepe\n";
      
      $dwork = $timeDiff*$nthreadsDiff;

      print OUTD "  $volume    $dwork   $timeDiff   $nthreadsDiff  $ncells  $density   $timepe\n";

      if (!defined($fileh_diffcalc[$nthreadsDiff])) {
         $ttp = sprintf("%0.2i",$nthreadsDiff);

         $tfilecalc = join '',$filebase_diffcalc,'.t',$ttp,'.dat';
         open($fileh_diffcalc[$nthreadsDiff],">$tfilecalc");
         $tfilestencil = join '',$filebase_diffstencil,'.t',$ttp,'.dat';
         open($fileh_diffstencil[$nthreadsDiff],">$tfilestencil");
         $tfilehalo = join '',$filebase_diffhalo,'.t',$ttp,'.dat';
         open($fileh_diffhalo[$nthreadsDiff],">$tfilehalo");
         $tfilewait = join '',$filebase_diffwait,'.t',$ttp,'.dat';
         open($fileh_diffwait[$nthreadsDiff],">$tfilewait");
         $tfiletot = join '',$timefile,'.diff.t',$ttp,'.dat';
         open($fileh_difftot[$nthreadsDiff],">$tfiletot");
      }
      print { $fileh_difftot[$nthreadsDiff] } "  $volume    $timeDiff  $timepe\n";
      print { $fileh_diffcalc[$nthreadsDiff] } "  $volume    $timeDiffCalc  $timepe\n";
      print { $fileh_diffstencil[$nthreadsDiff] } "  $volume    $timeDiffStencil  $timepe\n";
      print { $fileh_diffhalo[$nthreadsDiff] } "  $volume    $timeDiffHalo  $timepe\n";
      print { $fileh_diffwait[$nthreadsDiff] } "  $volume    $timeDiffWait  $timepe\n";

      $linecnt++;
   }

   close OUTD;
   close OUTR;
   close TIME;
   close LOAD;
}
