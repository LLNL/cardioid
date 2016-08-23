#! /usr/bin/perl -w

# plotecg_7panel.pl plots the ecgs for all 7 leads on a single file.  Each directory
# should contain:
#   - electrode_list.txt : list of electrodes
#   - object.data : simulation parameters
#   - electrode#000123456 : file containing ecg data for each electrode
#
# written by Erik Draeger, LLNL, 1/13/2014
#   -r2162 alterations and additions by Jonathan Cranford, LLNL, 9/21/2015 (altered from file /p/lscratchv/emhm/GrandChallenge/scripts/plotecg_7panel.pl)
#   	-added capability to specify which heart beat # is the first heart beat in the graph, rather than always showing the $graphNBeats final heart beats of simulation
#  	 -added extra checks for confirming $dt, $units, $printRate, and $BCL are same between regular and control directories
#   	-added various other automated error checking functionality
#  	 -added ability to print to title what beat(s) showing in graph

use POSIX;

if ($#ARGV < 0) {
   print "syntax:  plotecg_7panel.pl [directories with simulation and ecg data]\n";
   exit;
}

$interactive = 0;  # set to 1 to spawn xmgrace windows instead of plotting to file
$startBeat = 1;	   # beat number to start at in graphs
$graphNBeats = 1;  # number of BCLs to show in each graph, including current beat and going forward in time
$plotControl = 1;  # if control directory is symlinked in ecgdir, include on plot
$pout = 0;	   # set to 1 to print out values of key variables to terminal for bug detection, set to 0 to disable printing for a cleaner look

if ($pout)
{
  print "\nstartBeat=$startBeat, graphNBeats=$graphNBeats\n";
}

# input filenames
$object = "object.data";
$electrodeList = "electrode_list.txt";
$electrodeBase = "electrode\#";

# ecg electrode map
$electrodeName[0] = "V1";
$electrodeName[1] = "V2";
$electrodeName[2] = "V3";
$electrodeName[3] = "V4";
$electrodeName[4] = "V5";
$electrodeName[5] = "V6";
$electrodeName[6] = "I";
$electrodeName[7] = "I";

ECGDIR: foreach $ecgdir (@ARGV[0..$#ARGV]) {  

   $ymin = 1.E+09;
   $ymax = -1.E+09;

   # if control dir is symlinked, include on plot
   $controldir = "$ecgdir/control";
   $controlFound = 0;
   if (-e $controldir) 
   {
       $controlFound = 1;
       if ($plotControl == 1)
       {
          print "Control directory found, will include on plots.\n";
       }   
   }

   print "Processing $ecgdir...\n";
   if (-e $ecgdir)
   {
      # read dt and units parameters from object.data file in regular
      ($dt,$unit) = grabTimestep("$ecgdir/$object");
      if ($dt < 0)
      {
         print "Error reading time step from $ecgdir/$object.\n";
         next ECGDIR;
      }
      if ($unit eq "us")
      {
         $unit = "ms";
         $dt /= 1000;
      }
      elsif ($unit eq "s")
      {
         $unit = "ms";
         $dt *= 1000;
      }
      if (-e $controldir)
      {
      	# read dt and units parameters from object.data in control directory
      	($dtC,$unitC) = grabTimestep("$controldir/$object");
      	if ($dtC < 0)
      	{
          print "Error reading time step from $controldir/$object.\n";
          next ECGDIR;
        }
        if ($unitC eq "us")
        {
           $unitC = "ms";
           $dtC /= 1000;
        }
        elsif ($unitC eq "s")
        {
           $unitC = "ms";
           $dtC *= 1000;
        }
        if ($dt != $dtC)
        {
      	   print "\nERROR:  dt in regular directory ($dt) and control directory ($dtC) do not match, probably comparing apples to oranges\n";
	   print "Fatal error, exiting\n";
	   exit 1;
        }
        if ($unit ne $unitC)
        {
           print "\nERROR:  units in regular directory ($unit) and control directory ($unitC) do not match, probably comparing apples to oranges\n";
	   print "Fatal error, exiting\n";
	   exit 1;
        }
      }
      
      # read printRate parameter from object.data in regular directory
      $printRate = grabECGPrintRate("$ecgdir/$object");
      if ($printRate < 0)
      {
         print "Error reading ECG print rate from $object.\n";
         next ECGDIR;
      }
      if (-e $controldir)
      {
        # read printRate parameter from object.data in regular directory
        $printRateC = grabECGPrintRate("$controldir/$object");
        if ($printRateC < 0)
        {
           print "Error reading ECG print rate from $controldir/$object.\n";
           next ECGDIR;
        }
        if ($printRate != $printRateC)
        {
       	 print "\nERROR:  printRate in regular directory ($printRate) and control directory ($printRateC) do not match, probably comparing apples to oranges\n";
	 print "Fatal error, exiting\n";
	 exit 1;
        }
      } 

      # read BCL parameter from object.data in regular directory
      $BCL = grabMinStimPeriod("$ecgdir/$object");     # BCL should be in ms
      if ($BCL < 0)
      {
         print "Error reading BCL from $object.\n";
         next ECGDIR;
      }
      if (-e $controldir)
      {
        # read BCL parameter from object.data in control directory
        $BCL = grabMinStimPeriod("$ecgdir/$object");     # BCL should be in ms
        $BCLC = grabMinStimPeriod("$controldir/$object");     # BCL should be in ms
        if ($BCLC < 0)
        {
           print "Error reading BCL from $controldir/$object.\n";
           next ECGDIR;
        }
         if ($BCL != $BCLC)
        {
           print "\nERROR:  BCL in regular directory ($BCL) and control directory ($BCLC) do not match, probably comparing apples to oranges\n";
   	   print "Fatal error, exiting\n";
	   exit 1;
        }
      }

      $rescaleFactor = $printRate*$dt;
      if (-e $controldir)
      {
        $rescaleFactorC = $printRateC*$dtC;
      }
      
      # open electrode files, count lines for regular directory
      $electrodeListFull = "$ecgdir/$electrodeList";
      if (-e $electrodeListFull)
      {
         open ELIST,$electrodeListFull or die "Error opening $electrodeListFull: $!\n";;
         $ecnt = 0;
         while ($line=<ELIST>)
         {
            chomp $line;
            $enum = $line;
            $efiles[$ecnt] = sprintf("$ecgdir/$electrodeBase%9.9i",$enum);
            if (-e $efiles[$ecnt])
            {
               open ETRODE,$efiles[$ecnt];
               $linecnt = 0;
               while ($eline = <ETRODE>)
               {
                  chomp $eline;
                  if ($eline > $ymax) { $ymax = $eline; }
                  if ($eline < $ymin) { $ymin = $eline; }
                  $linecnt++;
               }
               close ETRODE;
               $elen[$ecnt] = $linecnt;
               $ecnt++;
            }
            else
            {
               print "Error:  file $efiles[$ecnt] not found!\n";
            }
         }
         close ELIST;
      }
      else
      {
         print "Error:  file $electrodeListFull not found!\n";
         next ECGDIR;
      }
      
      if (-e $controldir)
      {
        # open electrode files, count lines for control directory
        $ymin = 1.E+09;
        $ymax = -1.E+09;
        $electrodeListFullC = "$controldir/$electrodeList";
        if (-e $electrodeListFullC)
        {
           open ELISTC,$electrodeListFullC or die "Error opening $electrodeListFullC: $!\n";;
           $ccnt = 0;
           while ($line=<ELISTC>)
           {
              chomp $line;
              $cnum = $line;
              $cfiles[$ccnt]  = sprintf("$controldir/$electrodeBase%9.9i",$cnum);
              if (-e $cfiles[$ccnt])
              {
                 open ETRODEC,$cfiles[$ccnt];
                 $linecnt = 0;
                 while ($cline = <ETRODEC>)
                 {
                    chomp $cline;
                    if ($cline > $ymax) { $ymax = $cline; }
                    if ($cline < $ymin) { $ymin = $cline; }
                    $linecnt++;
                 }
                 close ETRODEC;
                 $clen[$ccnt] = $linecnt;
                 $ccnt++;
              }
              else
              {
                 print "Error:  file $cfiles[$ccnt] not found!\n";
              }
           }
           close ELISTC;
        }
        else
        {
           print "Error:  file $electrodeListFullC not found!\n";
           next ECGDIR;
        }
      }

      $nbeats = $elen[0]*$rescaleFactor/$BCL;
      if (-e $controldir)
      {
        $nbeatsC = $clen[0]*$rescaleFactorC/$BCLC;
      }
      if ( ($graphNBeats > ($nbeats-$startBeat+1)) && $startBeat != 1) {
      	print "\nERROR in regular directory:  \$graphNBeats+\$startBeat-1 must be <= \$nbeats!\n";
	print "Current state is \$graphNBeats=$graphNBeats, \$startBeat=$startBeat, \$nbeats=$nbeats\n";
        print "Fatal error, exiting\n";
      	exit 1;
      }
      if (-e $controldir)
      {
        if ( ($graphNBeats > ($nbeatsC-$startBeat+1)) && $startBeat != 1) {
       	  print "\nERROR in control directory:  \$graphNBeats+\$startBeat-1 must be <= \$nbeatsC!\n";
  	  print "Current state is \$graphNBeats=$graphNBeats, \$startBeat=$startBeat, \$nbeatsC=$nbeatsC\n";
          print "Fatal error, exiting\n";
      	  exit 1;
        }
      }
      if ($pout)
      {
        print "\nelen=$elen[0], rescaleFactor=$rescaleFactor, BCL=$BCL, dt=$dt, unit=$unit, printRate=$printRate, nbeats=$nbeats\n";
        if (-e $controldir)
        {
	  print "\nclen=$clen[0], rescaleFactorC=$rescaleFactorC, BCLC=$BCLC, dtC=$dtC, unitC=$unitC, printRateC=$printRateC, nbeats=$nbeats\n";
	}
      }
      @ecgsets = ($nbeats,$graphNBeats);
      plot7PanelGraceBatch(@ecgsets);


   }
   else
   {
      print "Error:  directory $ecgdir not found!\n";
      #exit;
   }
}

#################################################################################
sub grabTimestep
{
   my ($object) = @_;
   my $val = -1;
   my $unit = "";
   open OBJ,$object;
   while ($line = <OBJ>)
   {
      @data = split ' ',$line;
      if ($#data > 0)
      {
         if ($data[0] eq 'dt')
         {
            $val = $data[2];
            $unit = $data[3];         
            $unit =~ s/\;//g;  # remove semicolon
         }
      }
   }
   return ($val,$unit);
}

#################################################################################
sub grabECGPrintRate
{
   my ($object) = @_;
   my $val = -1;
   open OBJ,$object;
   $inSensor = 0;
   while ($line = <OBJ>)
   {
      @data = split ' ',$line;
      if ($#data > 1)
      {
         if ($data[0] eq 'method' && $data[2] =~ "VoronoiCoarsening")
         {
            $inSensor = 1;
         }
         if ($inSensor == 1)
         {
            if ($data[0] eq '}') { $inSensor = 0; }
            if ($data[0] eq 'printRate')
            {
               $val = $data[2];
               $val =~ s/\;//g;  # remove semicolon
            }
         }
      }
   }
   return ($val);
}

#################################################################################
sub grabMinStimPeriod
{
   my ($object) = @_;
   my $val = -1;
   open OBJ,$object;
   $inStimulus = 0;
   while ($line = <OBJ>)
   {
      @data = split ' ',$line;
      if ($#data > 0)
      {
         if ($data[1] eq 'STIMULUS')
         {
            $inStimulus = 1;
         }
         if ($inStimulus == 1)
         {
            if ($data[0] eq '}') { $inStimulus = 0; }
            elsif ($data[0] eq 'period')
            {
               $data[2] =~ s/\;//g;  # remove semicolon
               if ($data[2] < $val || $val < 0)
               {
                  $val = $data[2];
               }
            }
         }
      }
   }
   return ($val);
}

#################################################################################
sub plot7PanelGraceBatch
{
   my @inputsets = @_;
   my $difference = 0;

   my $beatnum = $inputsets[0];
   my $graphnbeats = $inputsets[1];

   my $nLeads = 7;
   my $linewidth = 2.5;
   my $controlline = 2.0;
   my $textsize = 1.0;
   my $pagesize = 800;   # graph size in pixels
   my $pagesizey = $pagesize*1.07;	# increase pages size in y direction to accomodate showing title and subtitle
   my $graphspacing = 10; # space between graphs
   my $leftmargin = 60;
   my $rightmargin = 30;
   my $margin = $leftmargin+$rightmargin;
   my $graphsize = ($pagesize-$margin)/3 - $graphspacing;
   my $xtickmajor = $BCL*$graphnbeats/5;	# only want 5 major ticks regardless of how many beats graphing
   my $ytickmajor = int ($ymax - $ymin)/4;
   if ($ytickmajor < 1) { $ytickmajor = 1; }

   my $fulldir = `pwd`;
   chomp $fulldir;
   $fulldir = join '',$fulldir,'/',$ecgdir;
   my @dirsplit = split '/',$fulldir;
   my $thisdir = $dirsplit[$#dirsplit];
   if ($thisdir eq '.')
   {
       $thisdir = $dirsplit[$#dirsplit-1];
   }

   my $filebase = "$ecgdir/ecg-$thisdir-7panelecg";

   my $batchfile = join '',$filebase,'.bat';
   open BATCH,">$batchfile";
   print BATCH "\# xmgrace batch file for plotting ecg data\n\n";
   print BATCH "\n";
   print BATCH "page size $pagesize $pagesizey\n";		# make page size in y direction larger to accomodate title and subtitle
#   print BATCH "view 0.05, 0.05, 0.95, 0.95\n";
#   print BATCH "view 0.0, 0.0, 1.0, 1.0\n";

   # compute differences for lead 1 :  assuming these files are entries 7 and 8 in electrode_list.txt
   {
      my $ecgfile1 = $efiles[6];
      my $ecgfile2 = $efiles[7];
      my $diffFile = join '',$ecgdir,'/electrode_lead1.dat';
      open ECG1,$ecgfile1;
      open ECG2,$ecgfile2;
      open ECGDIFF,">$diffFile";
      while (my $line1 = <ECG1>)
      {
         my $line2 = <ECG2>;
         chomp $line1;
         chomp $line2;
         my $diffval = $line2 - $line1;
         print ECGDIFF "$diffval\n";
      }
      close ECGDIFF;
      close ECG1;
      close ECG2;

      # copy diff file name to efiles[6]
      $efiles[6] = $diffFile;
   }

   # compute differences for lead 1 :  assuming these files are entries 7 and 8 in electrode_list.txt
   if ($controlFound == 1 && $plotControl == 1)
   {
      my $ecgfile1 = $cfiles[6];
      my $ecgfile2 = $cfiles[7];
      my $diffFile = join '',$ecgdir,'/electrode_lead1.control.dat';
      open ECG1,$ecgfile1;
      open ECG2,$ecgfile2;
      open ECGDIFF,">$diffFile";
      while (my $line1 = <ECG1>)
      {
         my $line2 = <ECG2>;
         chomp $line1;
         chomp $line2;
         my $diffval = $line2 - $line1;
         print ECGDIFF "$diffval\n";
      }
      close ECGDIFF;
      close ECG1;
      close ECG2;

      # copy diff file name to efiles[6]
      $cfiles[6] = $diffFile;
   }

   # format each set
   my $xmgrargs = '';
   my $xmax = 0;
   my $xmin = 0;
   for (my $ii=0; $ii<$nLeads; $ii++)
   {
      my $leadname = $electrodeName[$ii];

      if ($beatnum <= 1)
      {
      	$xmin = 0;
	$xmax = $BCL;
      }
      else
      {
      	$xmin = ($startBeat-1)*$BCL;
	if ($pout) {print "\nbeatnum=$beatnum, BCL=$BCL, xmin=$xmin\n";}
      	my $tmpmax = $BCL*($graphnbeats);
      	if ($tmpmax > $xmax) { $xmax = $tmpmax; }
      }
      if ($pout) {print "\nxmax = $xmax and xmin=$xmin\n";}
      print BATCH "\# format data (multiply by dt to rescale data from steps to seconds)\n";
      print BATCH "g$ii.s0.x = g$ii.s0.x + 1\n";
      print BATCH "g$ii.s0.x = g$ii.s0.x * $rescaleFactor\n";
      print BATCH "g$ii.s0.x = g$ii.s0.x - $xmin\n";
      print BATCH "g$ii.s0 legend \"$leadname\"\n";
      print BATCH "g$ii.s0 line linewidth $linewidth\n";
      print BATCH "g$ii.s0 errorbar linewidth $linewidth\n";
      print BATCH "g$ii.s0 errorbar riser linewidth $linewidth\n";
      if ($controlFound == 1 && $plotControl == 1)
      {
         print BATCH "\# format data (multiply by dt to rescale data from steps to seconds)\n";
         print BATCH "g$ii.s1.x = g$ii.s1.x + 1\n";
         print BATCH "g$ii.s1.x = g$ii.s1.x * $rescaleFactor\n";
         print BATCH "g$ii.s1.x = g$ii.s1.x - $xmin\n";
         print BATCH "g$ii.s1 line linewidth $controlline\n";
         print BATCH "g$ii.s1 errorbar linewidth $controlline\n";
         print BATCH "g$ii.s1 errorbar riser linewidth $controlline\n";
      }

      $xmgrargs = join ' ',$xmgrargs,'-graph ',$ii,' -block',$efiles[$ii],'-bxy index:1';
      if ($controlFound == 1 && $plotControl == 1)
      {
         $xmgrargs = join ' ',$xmgrargs,'-block',$cfiles[$ii],'-bxy index:1';
      }
   }
   $xmax *= 0.99;

   # graph layout   X: 0 = left, 1 = middle, 2 = right, Y: 0 = bottom, 1 = middle, 2 = top
   $xcoord[0] = 1;  $ycoord[0] = 2;
   $xcoord[1] = 1;  $ycoord[1] = 1;
   $xcoord[2] = 1;  $ycoord[2] = 0;
   $xcoord[3] = 2;  $ycoord[3] = 2;
   $xcoord[4] = 2;  $ycoord[4] = 1;
   $xcoord[5] = 2;  $ycoord[5] = 0;
   $xcoord[6] = 0;  $ycoord[6] = 2;

   # format graphs
   for (my $ii=0; $ii<$nLeads; $ii++)
   {
      $graph_xmin = ($leftmargin + $xcoord[$ii]*$graphspacing + $xcoord[$ii]*$graphsize)/$pagesize;
      $graph_xmax = $graph_xmin + $graphsize/$pagesize;
      $graph_ymin = ($leftmargin + $ycoord[$ii]*$graphspacing + $ycoord[$ii]*$graphsize)/$pagesize;
      $graph_ymax = $graph_ymin + $graphsize/$pagesize;
      $legendx = $graph_xmin + 0.70*($graph_xmax-$graph_xmin);
      $legendy = $graph_ymin + 0.95*($graph_ymax-$graph_ymin);
      print BATCH "with g$ii\n";
      print BATCH "  view xmin $graph_xmin\n";
      print BATCH "  view xmax $graph_xmax\n";
      print BATCH "  view ymin $graph_ymin\n";
      print BATCH "  view ymax $graph_ymax\n";
      print BATCH "  world xmin 0\n";
      print BATCH "  world xmax $xmax\n";
      print BATCH "  world ymin $ymin\n";
      print BATCH "  world ymax $ymax\n";
      print BATCH "  xaxis  tick major $xtickmajor\n";
      print BATCH "  xaxis  tick minor ticks 4\n";
      print BATCH "  yaxis  tick major $ytickmajor\n";
      print BATCH "  yaxis  tick minor ticks 4\n";
      print BATCH "  legend $legendx, $legendy\n";
      print BATCH "  legend box linestyle 0\n";
      print BATCH "  legend char size 1\n";
      print BATCH "  legend length 0\n";
      if ($ycoord[$ii] == 0 || ($xcoord[$ii] == 0 && $ycoord[$ii] == 2))  # label this x axis
      {
         print BATCH "  xaxis label \"time ( $unit )\"\n";
         print BATCH "  xaxis  label char size $textsize\n";
      }
      else
      {
         print BATCH "  xaxis ticklabel char size 0\n";
      }
      if ( $xcoord[$ii] == 0 || ($xcoord[$ii] == 1 && $ycoord[$ii] < 2))  # label this y axis
      {
         print BATCH "  yaxis label \"ecg ( mV )\"\n";
         print BATCH "  yaxis  label char size $textsize\n";
      }
      else
      {
         print BATCH "  yaxis ticklabel char size 0\n";
      }
      print BATCH "  \# line widths\n";
      print BATCH "  frame linewidth $linewidth\n";
      print BATCH "  xaxis bar linewidth $linewidth\n";
      print BATCH "  xaxis tick major linewidth 1.5\n";
      print BATCH "  xaxis tick minor linewidth 1.5\n";
      print BATCH "  yaxis bar linewidth $linewidth\n";
      print BATCH "  yaxis tick major linewidth 1.5\n";
      print BATCH "  yaxis tick minor linewidth 1.5\n";
      if ($xcoord[$ii] == 1 && $ycoord[$ii] == 2)
      {
         print BATCH "  title \"Showing $graphNBeats heart beat(s), starting at heart beat # $startBeat\"\n";
         print BATCH "  subtitle \"$fulldir\"\n";
      }
      
   }
   print BATCH "\n";
   close BATCH;

   $pngfile = join '',$filebase,'.png';

   if ($interactive == 1)
   {
      system("xmgrace -batch $batchfile $xmgrargs");
   }
   else
   {
      system("xmgrace $xmgrargs -batch $batchfile -nosafe -hardcopy -printfile $pngfile -hdevice PNG");
   }

   unlink($batchfile);

}
