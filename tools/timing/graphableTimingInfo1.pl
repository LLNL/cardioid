#! /usr/bin/perl -w
#
# read taskLoadInfo.data file and time.snapshot.0000xxxxxxxxx.info file (output from profileTimingInfo.pl)
# and compile timing data into graphable format.
#
# written by Erik Draeger, LLNL, 10/8/2012


if ($#ARGV < 0) {
   print "syntax:  graphableTimingInfo.pl [time.snapshot.info file(s)]\n";
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

   $outfile1 = join '',$timefile,'.g1.dat';
   $outfile2 = join '',$timefile,'.g2.dat';
   open OUT1, ">$outfile1";
   open OUT2, ">$outfile2";

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

   $linecnt = 0;
   while (($timeline = <TIME>) && ($loadline = <LOAD>))
   {
      #$loadline = <LOAD>;
      @timedata = split ' ',$timeline;
      @loaddata = split ' ',$loadline;

      $timepe = $timedata[1];  chop $timepe;  # remove comma
      $loadpe = $loaddata[0];
      if ($timepe != $loadpe) 
      {
         print "File mismatch:  different pes on line $linecnt!  $timepe, $loadpe\n";
         exit;
      }

      $nthreadsDiff = $timedata[2];
      $nthreadsReac = $timedata[5];
      $timeDiff = $timedata[12];
      $timeReac = $timedata[20];

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
      print OUT1 "   $ncells    $rwork   $timeReac   $nthreadsReac  $volume  $density   $timepe\n";

      if (!defined($fileh1[$nthreadsReac])) {
         $ttp = sprintf("%0.2i",$nthreadsReac);
         $tfile1 = join '',$timefile,'.g1.t',$ttp,'.dat';
         open($fileh1[$nthreadsReac],">$tfile1");
      }
      print { $fileh1[$nthreadsReac] } "  $ncells    $timeReac   $rwork   $nthreadsReac  $volume  $density     $timepe\n";
      
      $dwork = $timeDiff*$nthreadsDiff;

      print OUT2 "  $volume    $dwork   $timeDiff   $nthreadsDiff  $ncells  $density   $timepe\n";

      if (!defined($fileh2[$nthreadsDiff])) {
         $ttp = sprintf("%0.2i",$nthreadsDiff);
         $tfile2 = join '',$timefile,'.g2.t',$ttp,'.dat';
         open($fileh2[$nthreadsDiff],">$tfile2");
      }
      print { $fileh2[$nthreadsDiff] } "  $volume    $timeDiff   $dwork   $nthreadsDiff  $ncells  $density    $timepe\n";

      $linecnt++;
   }

   close OUT1;
   close OUT2;
   close TIME;
   close LOAD;
}
