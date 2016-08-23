#! /usr/bin/perl
#
# snapshotFilter.pl reads Cardioid snapshot files and filters out 
# points within a bounding box of a specific cell
#
# written by Erik Draeger, LLNL, 5/29/2012

if ($#ARGV != 0) {
   print "syntax:  snapshotDiff.pl [snapshot directory]\n";
   exit;
}

# coordinates of specific cell
$xx = 565;
$yy = 444;
$zz = 693;

# filter bounding box size (-bbx:+bbx, -bby:+bby, -bbz:+bbz)
$bbx = 50;
$bby = 50;
$bbz = 50;

$minx = $xx - $bbx;
$maxx = $xx + $bbx;
$miny = $yy - $bby;
$maxy = $yy + $bby;
$minz = $zz - $bbz;
$maxz = $zz + $bbz;

print "Filter all cells around grid point at ($xx, $yy, $zz), using bounding box $minx:$maxx, $miny:$maxy, $minz:$maxz.\n";

$dir = $ARGV[0];
$outdir = "snapshot.filter";
if (! -e $outdir)
{
   system("mkdir $outdir");
}

opendir(DIR,$dir) or die "Error opening $dir:  $!";
$ii = 0;
while (my $file = readdir(DIR)) {
   if ($file =~ /anatomy\#/) { $anatfiles[$ii++] = $file; }
   if ($file =~ /state\#/) { $statefiles[$ii++] = $file; }
}
closedir(DIR);

### READ AND CONVERT ANATOMY FILES FIRST ###

$nfiles = $#anatfiles + 1;
if ($nfiles <= 1) {
   print "No anatomy files found.\n";
}

@sortfiles = sort @anatfiles;

$nrecords = -1;
$reccnt = 0;
for ($ifile = 0; $ifile < $nfiles; $ifile++) {
   open ANAT, "$dir/$sortfiles[$ifile]";
   open OUT, ">$outdir/$sortfiles[$ifile]";
   $linecnt = 0;

   print "Reading $sortfiles[$ifile]...\n";

   if ($ifile == 0) {   # read file header
      $headerend = 0;
      
      while ($headerend == 0) {
         $line = <ANAT>;
         print OUT $line;
         $linecnt++;
         if ($line =~ /field_names/) { 
            chomp $line; chop $line;  # remove newline char, semicolon
            @field_names = split ' ',$line;
            shift @field_names;  # remove field_names
            shift @field_names;  # remove equal sign
            #$nfields = $#field_names+1;
         }
         if ($line =~ /field_types/) { 
            chomp $line; chop $line;  # remove newline char, semicolon
            @field_types = split ' ',$line;
            shift @field_types;  # remove field_types
            shift @field_types;  # remove equal sign
         }
         if ($line =~ /nx/ && $line =~ /ny/ && $line =~ /nz/) { 
            @data = split ' ',$line;
            $nx = $data[2];
            chop $nx;  # remove semicolon
            $ny = $data[5];
            chop $ny;  # remove semicolon
            $nz = $data[8];
            chop $nz;  # remove semicolon
            print "Grid read from header:  $nx x $ny x $nz\n";
         }
         if ($line =~ /nrecords/) { 
            chomp $line; chop $line;  # remove newline char, semicolon
            @tmpdata = split ' ',$line;
            $nrecords = $tmpdata[2];
         }
         if ($line =~ /}/) { 
            $headerend = 1; 
            $line = <ANAT>; # skip blank line after header
            print OUT $line;
         }
      }
   }

   while ($line = <ANAT>) {
      $linecnt++;
      @data = split ' ',$line;

      for ($jj=0; $jj<=$#data; $jj++) {
         if ($field_names[$jj] eq 'gid') {
            $gid = $data[$jj];
            $tmpgid = $gid;
            $thisx = $gid%$nx;
            $tmpgid = int ($gid/$nx);
            $thisy = $tmpgid%$ny;
            $thisz = int ($tmpgid / $ny);

            if ($thisx <= $maxx && $thisx >= $minx && 
                $thisy <= $maxy && $thisy >= $miny &&
                $thisz <= $maxz && $thisz >= $minz)
            {
               print OUT $line;
               $reccnt++;
            }
         }
      }

   }
   close ANAT;
   close OUT;
}

if ($nfiles > 1) {
   print "\nFinished!\n";
   print "   $reccnt coordinates filtered out (edit nrecord in header)\n";
}

### READ AND CONVERT STATE FILES ###

$nfiles = $#statefiles + 1;
if ($nfiles <= 1) {
   print "No state files found.\n";
}

@sortfiles = sort @statefiles;

$nrecords = -1;
$reccnt = 0;
for ($ifile = 0; $ifile < $nfiles; $ifile++) {
   open STATE, "$dir/$sortfiles[$ifile]";
   open OUT, ">$outdir/$sortfiles[$ifile]";
   $linecnt = 0;

   print "Reading $sortfiles[$ifile]...\n";

   if ($ifile == 0) {   # read file header
      $headerend = 0;
      
      while ($headerend == 0) {
         $line = <STATE>;
         print OUT $line;
         $linecnt++;
         if ($line =~ /field_names/) { 
            chomp $line; chop $line;  # remove newline char, semicolon
            @field_names = split ' ',$line;
            shift @field_names;  # remove field_names
            shift @field_names;  # remove equal sign
            #$nfields = $#field_names+1;
         }
         if ($line =~ /field_types/) { 
            chomp $line; chop $line;  # remove newline char, semicolon
            @field_types = split ' ',$line;
            shift @field_types;  # remove field_types
            shift @field_types;  # remove equal sign
         }
         if ($line =~ /nx/ && $line =~ /ny/ && $line =~ /nz/) { 
            @data = split ' ',$line;
            $nx = $data[2];
            chop $nx;  # remove semicolon
            $ny = $data[5];
            chop $ny;  # remove semicolon
            $nz = $data[8];
            chop $nz;  # remove semicolon
            print "Grid read from header:  $nx x $ny x $nz\n";
         }
         if ($line =~ /nrecords/) { 
            chomp $line; chop $line;  # remove newline char, semicolon
            @tmpdata = split ' ',$line;
            $nrecords = $tmpdata[2];
         }
         if ($line =~ /}/) { 
            $headerend = 1; 
            $line = <STATE>; # skip blank line after header
            print OUT $line;
         }
      }
   }

   while ($line = <STATE>) {
      $linecnt++;
      @data = split ' ',$line;

      for ($jj=0; $jj<=$#data; $jj++) {
         if ($field_names[$jj] eq 'gid') {
            $gid = $data[$jj];
            $tmpgid = $gid;
            $thisx = $gid%$nx;
            $tmpgid = int ($gid/$nx);
            $thisy = $tmpgid%$ny;
            $thisz = int ($tmpgid / $ny);

            if ($thisx <= $maxx && $thisx >= $minx && 
                $thisy <= $maxy && $thisy >= $miny &&
                $thisz <= $maxz && $thisz >= $minz)
            {
               print OUT $line;
               $reccnt++;
            }
         }
      }

   }
   close STATE;
   close OUT;
}

if ($nfiles > 1) {
   print "\nFinished!\n";
   print "   $reccnt coordinates filtered out (edit nrecord in header)\n";
}
