ibmAnatomy format definition
----------------------------

The ibm anatomy format is a binary format.  It consists of a 274 byte
header followed by the data.  There is one byte of data for each cell in
the grid.  The indexing is implicit with the x axis changing fastest
(i.e., x is in the inner loop.)

The structure of the header is mostly unknown.  All I can tell from the
code is that the "shared memory value" begins at file position 22 (i.e.,
the 23rd byte).  It fills four bytes.  The next three four byte
quantities (integers) are the x, y, and z sizes of the grid.



Other Notes
-----------



All operations on anatomy files should be coded with 64-bit file
operations.  This means setting the _LARGEFILE64_SOURCE macro.


From the BlueBeats.cpp reader code, it appears that the cell types are
"magic numbers" defined in some other file.  The following meanings can
be guessed:

0:       Cells on outer face of grid
9:       Blood/bath
30-35:   Tissue
75-77:   Tissue
100-102: Tissue
 
30, 31, 75, 100:  TT type 0 (Endo)
76, 101:          TT type 1 (Mid)
77, 102:          TT type 2 (Epi)


Note that zero is assigned to all cells on the boundary regardless of
what is read from the file.
