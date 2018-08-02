#!/usr/bin/perl -w

use strict;
use warnings;

sub usage 
{
   print <<"_HERE";
$0: outfile.cc Model1.cc [Model2.cc ...]

Scans the models for scan* functions, and includes them in a function
to populate the reactionFactory

_HERE

   die(@_) if @_;
}

usage("Need at least 1 argument") if scalar(@ARGV) < 1;

my $outfilename = shift @ARGV;

my @filenames = @ARGV;


my @reactions;
foreach my $filename (@filenames)
{
   open(my $file, $filename) or die "Could not open file '$filename' for reading: $!";
   while (my $line = <$file>)
   {
      if ($line =~ m,REACTION_FACTORY\((\w+)\)\(,)
      {
         push @reactions, $1;
      }
   }
   close($file);
}

open(my $outfile, ">", $outfilename) or die "Can't open file '$outfilename' for writing: $!";

print $outfile <<'_HERE';
#include "reactionFactory.hh"

_HERE

foreach my $reaction (@reactions)
{
   print $outfile "REACTION_FACTORY($reaction)(OBJECT* obj, const double dt, const int numPoints, const ThreadTeam& group);\n"
}

print $outfile <<'_HERE';

void registerBuiltinReactions()
{
_HERE

foreach my $reaction (@reactions)
{
   print $outfile qq|   registerReactionFactory("$reaction", reactionFactoryFor$reaction);\n|;
}
print $outfile <<'_HERE';
}

_HERE
