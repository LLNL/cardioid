# Cardioid

Cardioid is a cardiac multiscale simulation suite spanning from subcellular mechanisms up to simulations of organ-level clinical phenomena. The suite contains tools for simulating cardiac electrophysiology, cardiac mechanics, torso-ECGs, cardiac meshing and fiber generation tools.

## Requirements

Cardioid can be broken down into many separate simulators, each with their own dependencies.

* Cardioid electrics: C99, C++, OpenMP, MPI, and Lapack. It requires a valid perl installation to build. It van also take advantage of CUDA, in which case it also depends on NVTX and NVRTC.

* Cardioid fibers: Depends on C++ and MFEM.

* Cardioid mechanics: Depends on C++ and MFEM.

* Cardioid ecg: Depends on C++ and MFEM.

Some requirements are bundled into the build system, like libkdtree
and simUtil.

## Getting Started

We recommend compiling Cardioid using [Spack](https://github.com/LLNL/spack).  

To build Cardioid with spack, follow the following steps:

* Make a YAML file for your particular cluster telling the system where to find MPI, what compilers to use, where to find lapack, etc. Examples can be found in arch/*.yaml. [Read here](https://spack.readthedocs.io/en/latest/tutorial_environments.html#spack-yaml) for more information on the syntax for this YAML file.

* Clone [Spack](https://spack.io) and set it up:
   ```
   git clone https://github.com/spack/spack.git
   . spack/share/spack/setup-env.sh
   ```
   
* Install your YAML file as a spack environment.
   ```
   spack env create YOURENV arch/YOURENV.yaml
   ```
   
* Activate your environment
   ```
   spack env activate YOURENV
   ```

* Build the cardioid dependencies
   ```
   spack install mfem+hypre+lapack
   ```
   
* Install the dependencies into a directory called deps
   ```
   spack view symlink deps mfem+hypre+lapack
   ```

* Build the rest of Cardioid, using the default settings
   ```
   make build
   ```

## Building without spack


Cardioid is built with CMake, using the BLT make system.  A separate .cmake file is supplied for a variety of architectures.  Please feel free to make your own .cmake architecture for your particular cluster if needed.  Example .cmake architectures can be found in "arch/*.cmake"  Each architecture can be given a separate name, allowing multiple different versions of the code to be built on the same system.

If no architecture file is supplied, BLT will try to pick sane defaults.  To build everything, you will need to install MFEM manually to a directory of your choice and tell Cardioid where that installation lives (through the MFEM_DIR variable)

Builds are performed in `build/<arch>` .  Executables are installed in `build/<arch>/bin`.

## Contributing

Please submit any bugfixes or feature improvements as [pull requests](https://help.github.com/articles/using-pull-requests/).

Authors
----------------

Many thanks go to Cardioid's [contributors](https://github.com/llnl/cardioid/graphs/contributors).

* James P Glosli
* Tomas Oppelstrup
* Xiaohua Zhang
* David F Richards
* Robert Blake
* Erik Draeger
* Jean-Luc Fattebert
* Jamie Bramwell
* Arthur Mirin
* Sebastian Laudenschlager
* Viacheslav Gurev
* Jeremy Rice
* Changhoan Kim
* and many more...

License
----------------

Cardioid is distributed under the terms of both the MIT license. All new contributions must be made under this license.

See [LICENSE](https://github.com/llnl/simutil/blob/master/LICENSE) and [NOTICE](https://github.com/llnl/simutil/blob/master/NOTICE) for details.

`SPDX-License-Identifier: (MIT)`

``LLNL-CODE-764041``

