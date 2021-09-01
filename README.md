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

* Clone [Spack](https://spack.io) and set it up:
   ```
   git clone https://github.com/spack/spack.git
   . spack/share/spack/setup-env.sh
   ```
   
* Make a spack environment in order to build cardioid.
  ```
  mkdir myenv
  cd myenv
  spack env create -d .
  spacktivate .
  ```
* Load any necessary modules in order to build cardioid (cmake, perl, some sort of mpi, optional cuda)
* Add the necessary dependencies in order to build cardioid.
  ```
  spack compiler find
  spack external find cmake perl
  spack external find --not-buildable openmpi mvapich2 mpich spectrum-mpi cuda
  ```
* Add cardioid to the environment
  ```
  spack add cardioid
  # or
  spack add cardioid+cuda
  ```
* Add cardioid as a develop package
  ```
  spack develop cardioid@develop
  ```
  This will check out the cardioid source code into this directory
* Concretize everything
  ```
  spack concretize -f
  ```
* Build!
  ```
  spack install
  ```
  The binary will be in `.spack-env/view/bin/cardioid`
* Once you've done this, you can change directory to `cardioid/spack-build-*` and type `make` to build the source code manually.  All the configure options will be set up for you.

* To remove your changes and build from a fresh source tree again, just use undevelop and reinstall with spack
  ```
  spack undevelop cardioid
  spack install
  ```

## Building without spack


Cardioid is built with CMake, using the BLT make system.  A separate .cmake file is supplied for a variety of architectures.  Please feel free to make your own .cmake architecture for your particular cluster if needed.  Example .cmake architectures can be found in "arch/*.cmake"  Each architecture can be given a separate name, allowing multiple different versions of the code to be built on the same system.

If no architecture file is supplied, BLT will try to pick sane defaults.  To build everything, you will need to install MFEM manually to a directory of your choice and tell Cardioid where that installation lives (through the MFEM_DIR variable)

Builds are performed in `build/<arch>` .  Executables are installed in `build/<arch>/bin`.

## Contributing

Please submit any bugfixes or feature improvements as [pull requests](https://help.github.com/articles/using-pull-requests/).

## Authors

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

## Learn more

Cardioid's history goes back a few years -- it was a finalist for the 2012 Gordon Bell Prize -- but only now is the code available as open source. Initially developed by a team of LLNL and IBM scientists, Cardioid divides the heart into a large number of manageable subdomains. This replicates the electrophysiology of the human heart, accurately simulating the activation of each heart muscle cell and cell-to-cell electric coupling.

- Video: [The Cardioid Project: Simulating the Human Heart](https://computation.llnl.gov/cardioid-project-simulating-human-heart-0)
- *Science & Technology Review* article: [Venturing into the Heart of High-Performance Computing Simulations](https://str.llnl.gov/Sep12/streitz.html)
- *Science & Technology Review* article: [Reaching for New Computational Heights with Sequoia](https://str.llnl.gov/july-2013/mccoy) - Cardioid helped set speed records for the Sequoia supercomputer by clocking in at nearly 12 petaflops while scaling with better than 90% parallel efficiency across all 1,572,864 cores.

## License

Cardioid is distributed under the terms of the MIT license. All new contributions must be made under this license.

See [LICENSE](https://github.com/llnl/cardioid/blob/master/LICENSE) and [NOTICE](https://github.com/llnl/cardioid/blob/master/NOTICE) for details.

`SPDX-License-Identifier: MIT`

``LLNL-CODE-764041``
