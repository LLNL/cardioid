This is a quick and dirty test that runs on 32 nodes and help detect
whether Cardioid is working.

It tests the four main use cases:

* Production mode with rrg
* Production mode with cellML
* Dev mode with fgr simd
* Dev mode with threads

All of the last three should give the same answer within about 2e-11 max
error.  The rrg of course has a different answer.  You can compare to a
reference answer in /usr/gapps/emhm/referenceResults using the snapshot
compare tool that Eirk wrote.

Inputs and batch files are provided for bgq, bgp, and peloton.



Limitations
-----------

There are certain to be many cases these tests don't cover.  In
particular:

1.  The dev code specifies mod=0 so we aren't testing the pade
    approximates.
