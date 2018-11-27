# SimUtil

SimUtil is a utility library that makes it easy to quickly create new
High Performance Computing projects for scientific applications. The
library has proven to scale to over a million concurrent processes and
has been used in multiple Gordon Bell submissions and finalist codes.

The library features routines for

* MPI-aware hierarchical input deck processing
* Unit parsing and dimensional analysis
* PIO output file formats that scale to millions of processors
* Other miscellaneous routines commonly needed in HPC scientific codes.
 

## Requirements

SimUtil only depends on C99 and MPI.  The build system requires CMake
and [BLT](https://github.com/LLNL/blt), unless you wish to replace our
build system with your own.

## Getting Started

We recommend placing SimUtil directly in your build tree using
```git subtree```:

    git subtree add --prefix=simutil https://github.com/llnl/simutil.git master

From here, it is a simple matter to integrate SimUtil with your
existing build system

### Using CMake + BLT

The recommended way to include SimUtil in your library is through
[BLT](https://github.com/LLNL/blt). Include BLT as another subtree:

    git subtree add --prefix=blt https://github.com/LLNL/blt.git master

Then, place the following lines in your top level CMakeLists.txt:

    if (NOT BLT_LOADED)
       include(${CMAKE_CURRENT_SOURCE_DIR}/blt/SetupBLT.cmake)
    endif()

    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/simutil/src)

### Using a regular Makefile

Symlink `include/*.h` and `src/*.c` to the directory of your choice
and add these files to your makefile.  You also want add the following
flags:

    CFLAGS += -DDiff_Weight_Type_Single -DWITH_PIO -DWITH_MPI


## Contributing

Please submit any bugfixes or feature improvements as [pull
requests](https://help.github.com/articles/using-pull-requests/).

Authors
----------------

Many thanks go to SimUtil's
[contributors](https://github.com/llnl/simutil/graphs/contributors).

* James P Glosli
* Tomas Oppelstrup
* Shiv Sundram
* Xiaohua Zhang
* David F Richards
* Robert Blake
* Michael P Surh
* Fred Streitz

License
----------------

Simutil is distributed under the terms of both the MIT license. All
new contributions must be made under this license.

See [LICENSE](https://github.com/llnl/simutil/blob/master/LICENSE) and
[NOTICE](https://github.com/llnl/simutil/blob/master/NOTICE) for
details.

`SPDX-License-Identifier: (MIT)`

``LLNL-CODE-762538``
