#! /usr/bin/perl
#
# Combine the per-task electrode files written by the Cardioid ECG code when running 
# in parallel, with ngroups > 1
#
# written by Erik Draeger, 02/13/2014

if ($#ARGV < 0) {
  print "syntax:  combine_electrodes.pl [electrode#* file(s)]\n";
  exit;
}

$ecnt = -1;
$lastbase = '';
FILE: foreach $efile (@ARGV[0..$#ARGV]) 
{
   @fsplit = split '\.',$efile;
   if ($#fsplit == 1)
   {
      if ($fsplit[0] ne $lastbase) {
         $ecnt++;
         $filebase[$ecnt] = $fsplit[0];
         $lastbase = $fsplit[0];
         $filenames[$ecnt][0] = $efile;
         $filesuffix[$ecnt][0] = $fsplit[1];
         $nfiles[$ecnt] = 1;
      }
      else
      {
         $n = $nfiles[$ecnt];
         $filenames[$ecnt][$n] = $efile;
         $filesuffix[$ecnt][$n] = $fsplit[1];
         $nfiles[$ecnt]++;
      }
   }
   else
   {
      print "Could not figure out $efile, skipping it.\n";
      next FILE;
   }
}

for ($ee = 0; $ee <= $ecnt; $ee++)
{
   print "file base = $filebase[$ee]\n";

   for ($ff=0; $ff<$nfiles[$ee]; $ff++)
   {
      $tosort[$ff] = $filesuffix[$ee][$ff];
   }
   #@sorted = sort @tosort;
   @sorted = sort {$a <=> $b} @tosort;

   open OUT,">$filebase[$ee]";   
   for ($ff=0; $ff<$nfiles[$ee]; $ff++)
   {
      $readfile = join '',$filebase[$ee],'.',$sorted[$ff];
      open IN, $readfile or die "Could not open $readfile!\n";
      while ($line = <IN>)
      {
         @data = split ' ',$line;
         print OUT "$data[1]\n";
      }
      close IN;
   }
   close OUT;
}

