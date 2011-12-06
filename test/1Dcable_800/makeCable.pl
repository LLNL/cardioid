#! /usr/bin/perl

# makeCable generates the initial anatomy file for Cardioid to
# simulate a 1D cable of cells, buffered on all sides by cells with
# type = 0
#
# written by Erik Draeger, LLNL, 11/17/2011

if ($#ARGV != 1) {
  print "syntax:  makeCable [ncells in 1D] [cell type: endo,mid, or epi]\n";
  exit;
}

$ncellz = $ARGV[0];
$celltype = lc($ARGV[1]);
if (!($celltype eq "epi" || $celltype eq "mid" || $celltype eq "endo")) {
  print "Cell type must be epi, endo or mid!\n";
  exit;
}
$cellnum = -1;
if ($celltype eq "epi") { $cellnum = 102; }
elsif ($celltype eq "mid") { $cellnum = 101; }
elsif ($celltype eq "endo") { $cellnum = 100; }

$nx = 3;
$ny = 3;
$nz = $ncellz + 2;
$npts = $nx*$ny*$nz;

# print header
print "anatomy FILEHEADER \{\n";
print "  datatype = FIXRECORDASCII\;\n";
print "  nfiles = 1\;  nrecord = $npts\;  lrec = 35\;  endian_key = 875770417\;\n";
print "  nfields = 4\;\n";
print "  field_names = gid cellType theta phi\;\n";
print "  field_types = u u u u\;\n";
print "  nx = $nx\; ny = $ny\; nz = $nz\;\n";
print "\}\n\n";

# print data
$zero = 0;
$gid = 0;
for ($k=0; $k<$nz; $k++) {
  for ($j=0; $j<$ny; $j++) {
    for ($i=0; $i<$nx; $i++) {
      $num = 0;
      if ($i == 1 && $j == 1 && $k > 0 && $k != ($nz-1)) {
        $num = $cellnum;
      }
      printf " %12u   %4u   %4u   %4u\n",$gid,$num,$zero,$zero;
      $gid++;
    }
  }
}

