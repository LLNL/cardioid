# INSTALLATION ON PERSONAL MAC MACHINE

## INITIAL SETUP
Go to https://access.llnl.gov/cgi-bin/client-download.cgi
Click on OSX
Install the disk image

## MOVE TO DESIRED DIRECTORY
~~~
mkdir cardioid && cd cardioid
echo 'export CARDIOID='$PWD >> ~/.bash_profile && source ~/.bash_profile
~~~

## GET OPEN-MPI LIBRARIES
Install Xcode and Homebrew first
~~~
brew install open-mpi
~~~

## GET HYPRE
~~~
wget https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/download/hypre-2.10.0b.tar.gz
tar -xvf hypre-2.10.0b.tar.gz
cd hypre-2.10.0b/src
./configure
make -j4 && cd -
rm *.tar.gz
~~~

## GET METIS
~~~
wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz
tar -xvf metis-4.0.3.tar.gz
ln -sf metis-4.0.3 metis-4.0
cd metis-4.0
make -j4 && cd -
rm *.tar.gz
~~~

## GET MFEM
~~~
wget --trust-server-names http://goo.gl/Vrpsns
tar -xvf mfem-3.3.tgz
cd mfem-3.3
make parallel -j4
rm ../*.tgz
~~~

## GET FIB-GEN
Open Cisco AnyConnect Secure Mobililty Client using Spotlight search
Connect to vpn.llnl.gov
Group: llnl-vpnc
Username: hafez1
Enter Password and Connect
~~~
LCUSERNAME=hafez1;
git clone https://$LCUSERNAME@lc.llnl.gov/bitbucket/scm/cardioid/fib-gen.git
cd fib-gen
echo 'export TMPDIR=/tmp' >> ~/.bash_profile && source ~/.bash_profile
brew install gnu-sed
tar -xf ddcMD_files_r2265.tgz
gsed -i -e 's/matinv/matinv_g/g' ddcMD_files/src/three_algebra.c
gsed -i -e 's/matinv/matinv_g/g' ddcMD_files/src/three_algebra.h
# make sure the $MFEM_DIR point to the right MFEM installation path
cd ../mfem-3.3
export MFEM_DIR=`pwd`
cd -
make -j4
make surface
make -C ply2vtk
cd $CARDIOID
ln -sf mfem-3.3/fib-gen
~~~

---

## GO TO README FOR USAGE
