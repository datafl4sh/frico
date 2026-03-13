# FRICO - Friendly Radiation Integral COde

FRICO is a code implementing the Method of Moments (known as "Boundary Element Method" outside the Computational Electromagnetics community) as described in the "NASA Contractor Report 189594" titled "MOM3D Method of Moments Code - Theory Manual" by John F. Shaeffer.

This is a research tool for the Computational Electromagnetics and Amateur Radio communities, but it is currently in an early development stage. Currently it can only compute impedance and radiation patterns for antennas, in the future it is planned to add monostatic and bistatic radar cross section capabilities.

FRICO was developed because of a lack of free and open-source 3D MoM solvers and it integrates nicely with modern scientific computing tools, in particular:
 
 * [GMSH](https://gmsh.info/) for geometry definition and mesh generation
 * [VisIt](https://sd.llnl.gov/simulation/computer-codes/visit) and [SILO](https://software.llnl.gov/Silo/) for data visualization
 * [Gnuplot](http://www.gnuplot.info/) for plots and radiation diagrams
 * [Eigen](https://libeigen.gitlab.io/) for linear algebra
 * [HighFive](https://bluebrain.github.io/HighFive/) for HDF5 I/O

## System requirements and installation

FRICO runs on Unix systems, Windows is explicitly NOT supported. It will
probably run fine under WSL, but there is no guarantee about this.

The official way of obtaining FRICO is by downloading the source from this repository and compiling it on your machine. In order to successfully compile FRICO, apart from the packages listed above, you need a compiler supporting
C++23 and a fairly recent version of CMake.

### Installing prerequisites

On a Debian system, the prerequisites for FRICO can be installed with

    apt install libgmsh-dev gmsh libsilo-dev gnuplot libeigen3-dev libhdf5-dev

The visualization software VisIt is not available in the Debian packages and it is downloadable [here](https://sd.llnl.gov/simulation/computer-codes/visit).

### Compiling and installing FRICO

As FRICO is based on CMake, the steps to compile it are the usual ones:

    git clone https://github.com/datafl4sh/frico.git
    cd frico
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make -j `nproc`
    make install   # May need to be run with `sudo`

## Running FRICO

FRICO is a text-based application and the general workflow is the
following:

 * With GMSH, prepare the geometry you want to simulate
 * Run FRICO on the GMSH geometry
 * Use VisIt to visualize the computed fields, or use GNUPLOT (or any other
   plotting tool) to visualize the radiation pattens and frequency sweep
   results.

## Documentation

The authoritative documentation for FRICO is its man page and it is located in the `man` directory in the source tree. If you installed FRICO in a system location, just say `man frico`, otherwise go to the `man` directory in the source tree and say `man ./frico.1`. HTML and PDF versions of the man page is available under `doc`.

In `share/examples` you can find some usage examples of FRICO. Given that FRICO is a code in its early stages, you are encouraged to compare those examples with the results of other codes, like NEC or FEKO.

## References

* NASA Contractor Report 189594, _MOM3D Method of Moments Code - Theory Manual_,  John F. Shaeffer
* _Field computation by Moment Methods_, Roger Harrington
* _The Method of Moments in Electromagnetics_, Walton C. Gibson K4LLA

## Closing remarks

This software is proudly Friulian-made by IV3IWE, so that's the real reason of its name. For more Friulian-flavoured software you can also take a look to [MUSET](https://github.com/rvicedomini/muset).