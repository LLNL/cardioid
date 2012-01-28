#! /usr/bin/perl

# calcCV.pl calculates the conduction velocity by comparing voltage traces and calculating 
# the time difference between pulses
#
# written by Erik Draeger, LLNL, 12/5/2011

if ($#ARGV != 1 && $#ARGV != 2) {
   print "syntax:  calcCV [file1] [file2] [(optional) distance]\n";
   exit;
}

$file1 = $ARGV[0];
$file2 = $ARGV[1];
if (!-e $file1) {
   print "File $file1 not found!\n";
   exit;
}
if (!-e $file2) {
   print "File $file2 not found!\n";
   exit;
}

if ($#ARGV == 2) {
   $distance = $ARGV[2];
}

# loop through each file, calculate the times at which the voltage passed through zero
# (from above and from below)

my @upzeroes1;
my @downzeroes1;
my @upzeroes2;
my @downzeroes2;


open FILE1, "$file1";
$first = 1;
$lastval = 0.0;
$lasttime = 0.0;
LOOP1: while ($line = <FILE1>) {
   $firstchar = substr($line,0,1);
   if ($firstchar eq '#') { next LOOP1; }   # skip comment lines

   @data = split ' ',$line;
   if ($first == 1) {
      $first = 0;
      $lastval = $data[1];
      $lasttime = $data[0];
   }
   elsif ($lastval*$data[1] <= 0.0) {
      if ($lastval < 0.0) {
         $upzeroes1[$#upzeroes1+1] = $data[0];
      }
      elsif ($data[1] == 0.0 && $lastval == 0.0) {
         print "Warning:  both voltages are zero!\n";
         exit;
      }
      elsif ($data[1] <= 0.0) {
         $downzeroes1[$#upzeroes1+1] = $data[0];
      }
   }
   $lastval = $data[1];
}
close FILE1;

open FILE2, "$file2";
$first = 1;
$lastval = 0.0;
$lasttime = 0.0;
LOOP2: while ($line = <FILE2>) {
   $firstchar = substr($line,0,1);
   if ($firstchar eq '#') { next LOOP2; }   # skip comment lines

   @data = split ' ',$line;
   if ($first == 1) {
      $first = 0;
      $lastval = $data[1];
      $lasttime = $data[0];
   }
   elsif ($lastval*$data[1] <= 0.0) {
      if ($lastval < 0.0) {
         $upzeroes2[$#upzeroes2+1] = $data[0];
      }
      elsif ($data[1] == 0.0 && $lastval == 0.0) {
         print "Warning:  both voltages are zero!\n";
         exit;
      }
      elsif ($data[1] <= 0.0) {
         $downzeroes2[$#upzeroes2+1] = $data[0];
      }
   }
   $lastval = $data[1];
}
close FILE2;

print "$#downzeroes1 down zeroes, $#upzeroes1 up zeroes found for $file1\n";
print "$#downzeroes2 down zeroes, $#upzeroes2 up zeroes found for $file2\n";

print "Up pulses: \n";
for ($i=0; $i<$#upzeroes2; $i++) {
   if ($#upzeroes1 >= $i) {
      $deltat = $upzeroes2[$i]-$upzeroes1[$i];
      if ($#ARGV == 2) {
         $velocity = $distance/$deltat;
         printf "pulse $i, deltat = %0.4f, velocity = %0.4f\n",$deltat,$velocity;
      }
      else {
         printf "pulse $i, deltat = %0.4f\n",$deltat;
      }
   }
}


