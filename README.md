# FRICO - Friendly Radiation Integral COde

FRICO is a code implementing the Method of Moments (known as "Boundary Element Method" outside the Computational Electromagnetics community) as described in the "NASA Contractor Report 189594" titled "MOM3D Method of Moments Code - Theory Manual" by John F. Shaeffer.

This is a research tool for the Computational Electromagnetics and Amateur Radio communities, but it is currently in an early development stage.

FRICO was developed because of a lack of free and open-source 3D MoM solvers and it integrates nicely with modern scientific computing tools, in particular:
 
 * [GMSH](https://gmsh.info/) for geometry definition and mesh generation
 * [VisIt](https://sd.llnl.gov/simulation/computer-codes/visit) and [SILO](https://software.llnl.gov/Silo/) for data visualization
 * [Gnuplot](http://www.gnuplot.info/) for plots and radiation diagrams
 * [Eigen](https://libeigen.gitlab.io/) for linear algebra
 * [HighFive](https://bluebrain.github.io/HighFive/) for HDF5 I/O


