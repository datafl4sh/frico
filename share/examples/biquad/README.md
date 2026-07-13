# Biquad antenna simulation with FRICO

The model for the biquad antenna is contained in the file `biquad.geo`, which is a geometry specification suitable for GMSH. You can inspect the model by running

    gmsh biquad.geo

If you look at the `.geo` file, you will find a variable called `MHz` which is used to specify the center frequency of the antenna (in MHz of course). There is also a variable `W` that allows to specify the radiator width.

## Simulating the biquad

Once you are satisfied with your geometry, you can simulate the biquad. For
example, to run a simulation at the frequency of 2.45 GHz, you call FRICO as
follows:

    frico -g biquad.geo -f 2.45e9 -A -s src -t -o biquad.silo -O biquad.h5 -P yz:0.5:0.5:0.005

Let's break down this command line:

* `-g` specifies the geometry file, in this case `biquad.geo`
* `-f` specifies the frequency in Hz, in this case 2.45e9, or 2.45 GHz
* `-A` enables the Approximate mode: integrals are computed as described in the MOM3D manual, without using quadratures.
* `-s` specifies the name of the Physical Curve where to apply the source, in this case `src`
* `-t` (optional) enables text file data output
* `-o` (optional) specifies the name for the SILO output file
* `-O` (optional) specifies the name for the HDF5 output file
* `-P` (optional) evaluates the fields on a rectangle of size 0.5 x 0.5, centered on the origin and aligned with the yz plane. The fields are evaluated on a triangular mesh with size 0.005.


## Visualizing the results of the computation

If everything worked, you will get three files:
* `biquad.silo` (if `-o` was specified)
* `biquad.h5` (if `-O` was specified)
* `polar_0.txt`
* `port_sweep.txt`

### Fields

The computed fields are contained in `biquad.silo`, which is a file that has
to be opened with [VisIt](https://sd.llnl.gov/simulation/computer-codes/visit).
The main variables are `J_real_magnitude` and `J_imag_magnitude` under "Pseudocolor" and `J_real`/`J_imag` under "Vectors". You will find other variables, but this will change in the future.

### Gain patterns

The file `polar_0.txt` contains gain patterns in the XY, YZ and XZ planes. In particular the first column is the angle in degrees whereas the remaining columns are the linear and dB gain in the XY, YZ and XZ planes respectively.
To plot the gain patterns, there is a Gnuplot script named `polar_dB.plot` available in `share/gnuplot` in the FRICO source tree and you can launch it as

    gnuplot --persist -e "filename='polar_0.txt'" ../../gnuplot/polar_dB.plot

Of course, the paths must be adjusted if you launch the script from outside the directory of this example. Alternatively, you can use your preferred programming language or library to plot the contents of the file.

### Impedance and SWR

The file `port_sweep.txt` contains the impedance and SWR values computed at the excitation port. For SWR computation, the system impedance is assumed to be 50 Ohm, this can be changed by using the `-Z` flag on the command line.
The first column contains the sequential number of the frequency point, in this case there will be only one line since we ran the simulation at a single frequency. The following columns are frequency, resistance, reactance and SWR respectively.

## Frequency sweeps

FRICO supports frequency sweeps with the `-R` command line parameter as follows:

    frico -g biquad.geo -R 2.4e9:10e6:2.5e9 -A -s src -o biquad.silo -P yz:0.5:0.5:0.005

The argument of `-R` is a range in the form of `<start>:<step>:<end>`, so the command line above will do a frequency sweep from 2.4 GHz to 2.5 GHz at steps of 10 MHz. In this case, you will get a `biquad.silo` file where all the fields are organized in subfolders and a number of `polar_N.txt` files, one for each simulated frequency. They contents follow the logic described above.

The file `port_sweep.txt` in this case gets more interesting, because in this case it will include multiple lines, one for each frequency. The contents will of course include resistance, reactance and SWR formatted as described above. A Gnuplot script is provided also to plot this information and is launched with

    gnuplot --persist -e "filename='port_sweep.txt'" ../../gnuplot/impedance.plot

Of course, also in this case, you can use your preferred tool to plot this data.