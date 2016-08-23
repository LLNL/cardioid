#! /usr/bin/perl -w
#
# generate report of results from different compareSnapshots runs
#
# written by Erik Draeger, LLNL, 9/5/2012


if ($#ARGV < 0) {
   print "syntax:  reportCompareResults.pl [directories containing verif.Vm.dat files]\n";
}

foreach $dir (@ARGV[0..$#ARGV]) {
  if (-e $dir) 
  {
     # remove any slashes from $dir
     $dir =~ s/\///g;
     $runfile = join '',$dir,"/rundebug.",$dir,".pel";
     $vmfile = join '',$dir,"/verif.Vm.dat";
     if (-e $runfile && -e $vmfile) 
     {
        $srunline = `grep srun $runfile`;

        @data = split 'TT06',$srunline;

        @tmp1 = split ' ',$data[1];
        $rundetails1 = $tmp1[0];
        @tmp2 = split ' ',$data[2];
        $rundetails2 = $tmp2[0];

        # rundetails will look like e.g. RRG-fast0mod0rfns0-N4t4
        @run1 = split '-',$rundetails1;
        @run2 = split '-',$rundetails2;

        @vmdata = split ' ',`tail -1 $vmfile`;

        $common = '';
        $diff1 = '';
        $diff2 = '';
        for ($i=0; $i<=$#run1; $i++)
        {
           if ($run1[$i] eq $run2[$i])
           {
              if ($common eq '')
              {
                 $common = $run1[$i];
              }
              else 
              {
                 $common = join '-',$common,$run1[$i];
              }
           }
           else 
           {
              if ($diff1 eq '')
              {
                 $diff1 = $run1[$i];
                 $diff2 = $run2[$i];
              }
              else
              {
                 $diff1 = join '-',$diff1,$run1[$i];
                 $diff2 = join '-',$diff2,$run2[$i];
              }
           }
        }

        print "$dir:   $common, $diff1 vs. $diff2\n";

        for ($i=0; $i<=$#vmdata; $i++)
        {
           print "   $vmdata[$i] ";
        }
        print "\n\n";
     }
  }
}
