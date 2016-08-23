#! /usr/bin/perl -w
#
# read profile# files from snapshot directory, analyze timings w. diffusion details
#
# written by Erik Draeger, LLNL, 10/18/2012


if ($#ARGV < 0) {
   print "syntax:  profileTimingInfo-Diff1.pl [snapshot directories]\n";
}

foreach $dir (@ARGV[0..$#ARGV]) 
{
   if (-e $dir) 
   {
      # remove any slashes from $dir
      $dir =~ s/\///g;

      opendir(DIR,$dir) or die "Error opening $dir:  $!";
      $ii = 0;
      while (my $file = readdir(DIR)) {
         if ($file =~ /profile\#/) { $proffiles[$ii++] = $file; }
      }
      closedir(DIR);

      @sortfiles = sort @proffiles;

      $outfile = join '',"time.",$dir,".info";
      open OUT, ">$outfile";

      for ($i=0; $i<=$#sortfiles; $i++)
      {
         print "Processing $sortfiles[$i]...\n";
         
         open PFILE, "$dir/$sortfiles[$i]";
         while ($line = <PFILE>)
         {
            if ($line =~ "Performance for rank")
            {
               @data = split ' ',$line;
               $rank = $data[3];
               
               $nDiffThreads[$rank] = 0;
               $nReactThreads[$rank] = 0;
               $minDiffWait[$rank] = 1.0E+19;
               $maxDiffCalc[$rank] = 0;
               $maxDiffStencil[$rank] = 0;
               $maxDiffHalo[$rank] = 0;
               $maxDiffThread[$rank] = -1;
               $minReactWait[$rank] = 1.0E+19;
               $maxReactThread[$rank] = -1;
            }
            elsif ($line =~ "DiffusionLoop ")
            {
               $nDiffThreads[$rank]++;
               @data = split ' ',$line;
               $lastDiffTime = $data[5];
            }
            elsif ($line =~ "DiffusionCalc ")
            {
               @data = split ' ',$line;
               $lastDiffCalc = $data[5];
            }
            elsif ($line =~ "FGR_Stencil ")
            {
               @data = split ' ',$line;
               $lastDiffStencil = $data[5];
            }
            elsif ($line =~ "HaloExchange ")
            {
               @data = split ' ',$line;
               $lastDiffHalo = $data[5];
            }
            elsif ($line =~ "DiffusionWait ")
            {
               @data = split ' ',$line;
               @data2 = split ':',$data[0];
               $tid = $data2[0];
               if ($data[5] < $minDiffWait[$rank]) {
                  $minDiffWait[$rank] = $data[5];
                  $maxDiffTime[$rank] = $lastDiffTime - $data[5];
                  $maxDiffThread[$rank] = $tid;
                  $maxDiffCalc[$rank] = $lastDiffCalc;
                  $maxDiffStencil[$rank] = $lastDiffStencil;
                  $maxDiffHalo[$rank] = $lastDiffHalo;
               }
            }
            elsif ($line =~ "ReactionLoop ")
            {
               $nReactThreads[$rank]++;
               @data = split ' ',$line;
               $lastReactTime = $data[5];
            }
            elsif ($line =~ "ReactionWait ")
            {
               @data = split ' ',$line;
               @data2 = split ':',$data[0];
               $tid = $data2[0];
               if ($data[5] < $minReactWait[$rank]) {
                  $minReactWait[$rank] = $data[5];
                  $maxReactTime[$rank] = $lastReactTime - $data[5];
                  $maxReactThread[$rank] = $tid;
               }
            }
         }
         close PFILE;
      }

      # print out all timing information
      print OUT "#  Task  nThreadsD   nThreadsR  DiffTime  DiffWait  DiffCalc  DiffStencil  DiffHalo  ReactTime   ReactWait\n";
      for ($ip=0; $ip<=$#maxDiffTime; $ip++)
      {
         #print OUT "Task $ip, $nDiffThreads[$ip] diff. threads, $nReactThreads[$ip] react. threads, max diff time = $maxDiffTime[$ip] (wait = $minDiffWait[$ip]), max react time = $maxReactTime[$ip] (wait = $minReactWait[$ip])\n";
         print OUT " $ip    $nDiffThreads[$ip]   $nReactThreads[$ip]   $maxDiffTime[$ip]   $minDiffWait[$ip]    $maxDiffCalc[$ip]   $maxDiffStencil[$ip]   $maxDiffHalo[$ip]    $maxReactTime[$ip]   $minReactWait[$ip]\n";
      }

      # print out maximum and average timing information
      $maxpe = -1;
      $maxtime = -1;
      $avgtime = 0;
      $avgDifftime = 0;
      $avgReacttime = 0;
      $cnt = 0;
      for ($ip=0; $ip<$#maxDiffTime; $ip++)
      {
         $cnt++;
         if ($maxDiffTime[$ip] > $maxtime) { $maxtime = $maxDiffTime[$ip]; $maxpe = $ip; }
         if ($maxReactTime[$ip] > $maxtime) { $maxtime = $maxReactTime[$ip]; $maxpe = $ip; }

         if ($maxDiffTime[$ip] > $maxReactTime[$ip]) { $avgtime += $maxDiffTime[$ip]; }
         else { $avgtime += $maxReactTime[$ip]; }
         
         $avgDifftime += $maxDiffTime[$ip];
         $avgReacttime += $maxReactTime[$ip];

      }
      $avg = $avgtime/$cnt;
      $avgDiff = $avgDifftime/$cnt;
      $avgReact = $avgReacttime/$cnt;

      print OUT "\n";
      printf OUT "Average task time = %0.2f (avg diff = %0.2f, avg react = %0.2f)\n",$avg,$avgDiff,$avgReact;
      printf OUT "Maximum time spent on task $maxpe:  max time = %0.2f  (diff time = %0.2f, react time = %0.2f)\n",$maxtime,$maxDiffTime[$maxpe],$maxReactTime[$maxpe];
      close OUT;
  }
}
