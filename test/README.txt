runner.py - thing that does the running, uses argparse, run with no arguments for help
numCompare.py - does numerical comparisons on text files in a smart way, looking for differences in numerical error.
loadTest - This is like virtual environment for a test.  Read the beginning of the file for the arguments.

The basic idea is that you have different build profiles (bluegene, gcc , clang, etc) , and different testing profiles (sequential, mpi, valgrind).  These can be combined combinatorically.  

A test writes a "definition" file that specifies bash functions you can run.  Running a test means that the test writes out a file called "results".  If "results" existis and is empty, that means the test passed.

"save" copies files over from a scratch space to the results file, which can be checked into git to supply things like md5 sums or whatever you'd like.
