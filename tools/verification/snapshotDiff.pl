#! /usr/bin/perl

# snapshotDiff.pl reads through anatomy files and computes the differences in voltages and (optionally)
# derivatives to get a quantitative measure of the similarity between two runs
#
# written by Erik Draeger, LLNL, 3/14/2012

if ($#ARGV != 1) {
   print "syntax:  snapshotDiff.pl [snapshot directory 1] [snapshot directory 2]\n";
   exit;
}

# strict_checking = 2:  print warning and exit if any unsigned data doesn't match exactly
#                 = 1:  print warning if any unsigned data doesn't match exactly, don't exit
#                 = 0:  ignore unsigned data
$strict_checking = 2;   

$dir1 = $ARGV[0];
$dir2 = $ARGV[1];

opendir(DIR1,$dir1) or die "Error opening $dir1:  $!";
$ii = 0;
while (my $file = readdir(DIR1)) {
   if ($file =~ /anatomy\#/) { $anatfiles1[$ii++] = $file; }
}
closedir(DIR1);
opendir(DIR2,$dir2) or die "Error opening $dir2:  $!";
$ii = 0;
while (my $file = readdir(DIR2)) {
   if ($file =~ /anatomy\#/) { $anatfiles2[$ii++] = $file; }
}
closedir(DIR2);

$nfiles1 = $#anatfiles1 + 1;
$nfiles2 = $#anatfiles2 + 1;
if ($#anatfiles1 != $#anatfiles2) { 
   print "Directory $dir1 has different number of anatomy files ($nfiles1) than directory $dir2 ($nfiles2)!\n";
   exit;
}
if ($#anatfiles1 < 0) {
   print "No anatomy files found!\n";
   exit;
}

@sortfiles1 = sort @anatfiles1;
@sortfiles2 = sort @anatfiles2;

$nrecords = -1;
$totlines = 0;
for ($ifile = 0; $ifile < $nfiles1; $ifile++) {
   open ANAT1, "$dir1/$sortfiles1[$ifile]";
   open ANAT2, "$dir2/$sortfiles2[$ifile]";
   $linecnt = 0;


   print "Reading $sortfiles1[$ifile]...\n";

   if ($ifile == 0) {   # read file header
      $headerend = 0;
      
      while ($headerend == 0) {
         $line1 = <ANAT1>;
         $line2 = <ANAT2>;
         $linecnt++;
         if ($line1 ne $line2) { 
            print "WARNING: Header mismatch at line $linecnt!\n"; 
            print "    $dir1/$sortfiles1[$ifile]:  $line1";
            print "    $dir2/$sortfiles2[$ifile]:  $line2";
         }
         if ($line1 =~ /field_names/) { 
            chomp $line1; chop $line1;  # remove newline char, semicolon
            @field_names = split ' ',$line1;
            shift @field_names;  # remove field_names
            shift @field_names;  # remove equal sign
            #$nfields = $#field_names+1;
         }
         if ($line1 =~ /field_types/) { 
            chomp $line1; chop $line1;  # remove newline char, semicolon
            @field_types = split ' ',$line1;
            shift @field_types;  # remove field_types
            shift @field_types;  # remove equal sign
         }
         if ($line1 =~ /nrecords/) { 
            chomp $line1; chop $line1;  # remove newline char, semicolon
            @tmpdata = split ' ',$line1;
            $nrecords = $tmpdata[2];
         }
         if ($line1 =~ /}/) { 
            if ($nrecords < 0) { 
               print "WARNING:  nrecords not found!  Averages will be wrong.\n";
            }
            $headerend = 1; 
            $line1 = <ANAT1>; # skip blank line after header
            $line2 = <ANAT2>; # skip blank line after header
         }
      }
   }

   for ($jj=0; $jj<=$#field_types; $jj++) {
      $maxDiff[$jj] = 0;
      $avgDiff[$jj] = 0;
   }

   while ($line1 = <ANAT1>) {
      $line2 = <ANAT2>;
      $linecnt++;
      @data1 = split ' ',$line1;
      @data2 = split ' ',$line2;

      if ($#data1 != $#field_names || $#data1 != $#field_types || $#data2 != $#data1 ) {
         print "Data mismatch at line $linecnt:  ndata1 = $#data1, ndata2 = $#data2, #field_names = $#field_names, #field_types = $#field_types\n";
         exit;
      }

      for ($jj=0; $jj<=$#data1; $jj++) {
         if ($field_types[$jj] eq 'u' && $strict_checking > 0) { 
            if ($data1[$jj] != $data2[$jj]) { 
               print "Unsigned data field $field_names[$jj] does not agree at line $linecnt!  Exiting...\n";
               print "   line1:  $line1";
               print "   line2:  $line2";
               if ($strict_checking == 2) {
                  exit;
               }
            }
         }
         elsif ($field_types[$jj] eq 'f') {
            $diff = abs($data1[$jj] - $data2[$jj]);
            $avgDiff[$jj] += $diff;
            if ($diff > $maxDiff[$jj]) { $maxDiff[$jj] = $diff; }
         }
      }

   }
   close ANAT1;
   close ANAT2;
   $totlines += $linecnt;
}

print "\nFinished!\n";
print "   $totlines total lines read from $nfiles1 anatomy files.\n";
for ($jj=0; $jj<=$#field_types; $jj++) {
   if ($field_types[$jj] eq 'f') {
      $avgDiff[$jj] /= $nrecords;
      printf "   $field_names[$jj]:  maximum difference = %0.12e, average difference = %0.12e\n",$maxDiff[$jj],$avgDiff[$jj];
   }
}

