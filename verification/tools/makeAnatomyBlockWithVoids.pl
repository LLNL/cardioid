#! /usr/bin/perl
#
# generate Cardioid anatomy input files for block geometry w. voids
#
# written by Erik Draeger, LLNL, 8/17/2012
#
# note:  any anatomies generated before 9/13/2012 had conductivity tensors == 0.0

if ($#ARGV != 4) {
   print "syntax:  makeAnatomyBlockAndVoids.pl [nx] [ny] [nz] [cellType] [nvoids]\n";
   exit;
}

$maxconductivity = 0.02;

# initialize random number generator
srand(88731);

$nfiles = 16;

$nx = $ARGV[0];
$ny = $ARGV[1];
$nz = $ARGV[2];
$cellType = $ARGV[3];
$nVoids = $ARGV[4];

$ncells = $nx*$ny*$nz;
$nloc = $ncells/$nfiles;

print "Calling makeAnatomyBlockAndVoids:  $nx x $ny x $nz, celltype $cellType, nvoids = $nVoids\n";


# print header
$dx = 0.1;
$dy = 0.1;
$dz = 0.1;
$filecnt = 0;
$filename = sprintf("anatomy\#%0.6i",$filecnt);
open ANAT,">$filename";

print ANAT "anatomy FILEHEADER {\n";
print ANAT "  datatype = FIXRECORDASCII\;\n";
print ANAT "  nfiles = $nfiles\;  nrecord = $ncells\;  lrec = 85\;  endian_key = 875770417\;\n";
print ANAT "  nfields = 8\;\n";
print ANAT "  field_names = gid cellType sigma11 sigma12 sigma13 sigma22 sigma23 sigma33\;\n";
print ANAT "  field_types = u u f f f f f f\;\n";
print ANAT "  nx = $nx\; ny = $ny\; nz = $nz\;\n";
print ANAT "  field_units = 1 1 mS/mm mS/mm mS/mm mS/mm mS/mm mS/mm\;\n";
print ANAT "  dx = $dx\; dy = $dy\; dz = $dz\;\n";
print ANAT "  }\n\n";

# create voids
$mindim = $nx;
if ($ny < $mindim) { $mindim = $ny; }
if ($nz < $mindim) { $mindim = $nz; }
$maxradius = 0.2;  # fraction of min dimension
$minradius = 0.05;  # fraction of min dimension

for ($vv=0; $vv<$nVoids; $vv++) 
{
   $voidRadius[$vv] = ($minradius + $maxradius*rand())*$mindim;
   $rx = $nx*rand(); $ry = $ny*rand(); $rz = $nz*rand(); 
   $voidMinX[$vv] = $rx - $voidRadius[$vv];
   $voidMaxX[$vv] = $rx + $voidRadius[$vv];
   $voidMinY[$vv] = $ry - $voidRadius[$vv];
   $voidMaxY[$vv] = $ry + $voidRadius[$vv];
   $voidMinZ[$vv] = $rz - $voidRadius[$vv];
   $voidMaxZ[$vv] = $rz + $voidRadius[$vv];

   printf "Creating void of radius %0.2f, centered at %0.1f %0.1f %0.1f\n",$voidRadius[$vv],$rx,$ry,$rz;

}

for ($k=0; $k<$nz; $k++) {
   for ($j=0; $j<$ny; $j++) {
      for ($i=0; $i<$nx; $i++) {
         $gid = $i + $j*$nx + $k*$nx*$ny;

         if ($gid > ($filecnt+1)*$nloc)
         {
            $filecnt++;
            close ANAT;
            $filename = sprintf("anatomy\#%0.6i",$filecnt);
            open ANAT,">$filename";
         }

         $type = $cellType;
         if ($cellType == "random")
         {
            $rtmp = 100. + 3.*rand();
            $type = int $rtmp;
         }

         # if cell is in a void, set type to 0
         for ($vv=0; $vv<$nVoids; $vv++) 
         {
            if ($i > $voidMinX[$vv] && $i < $voidMaxX[$vv] 
                && $j > $voidMinY[$vv] && $j < $voidMaxY[$vv] 
                && $k > $voidMinZ[$vv] && $k < $voidMaxZ[$vv])
            {
               $type = 0;
            }
         }
         $sigxx = $maxconductivity*rand();
         $sigyy = $maxconductivity*rand();
         $sigzz = $maxconductivity*rand();
         $sigxy = $maxconductivity*rand();
         $sigyz = $maxconductivity*rand();
         $sigxz = $maxconductivity*rand();
         #$line = sprintf("          $gid $type 0.000 0.000 0.000 0.000 0.000 0.000");
         $line = sprintf("          $gid $type %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f",$sigxx,$sigyy,$sigzz,$sigxy,$sigyz,$sigxz);
         printf ANAT "%-84s\n",$line;
      }
   }
}

close ANAT;
print "\n";
