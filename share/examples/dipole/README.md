# Half-wavelength dipole simulation with FRICO

The model for the half-wavelength dipole is contained in the file `dipole.geo`, which is a geometry specification suitable for GMSH. You can inspect the model by running

    gmsh dipole.geo

The dipole is modeled as a rectangular strip, as FRICO does not support thin wires yet.

In order to run a frequency sweep from 250 MHz to 350 MHz at steps of 10 MHz you can launch

    frico -g dipole.geo -R 250e6:10e6:350e6 -A -s src -t -o dipole.silo

This will produce, as explained in the example for the biquad antenna:

* `dipole.silo` file with the fields at the different frequencies
* `polar_N.txt` files with the gain data at the different frequencies
* the file `port_sweep.txt` with the impedance data at the input port

## Comparison with NEC

This example includes also a NEC file for the same dipole, configured for a frequency sweep from 250 MHz to 350 MHz at steps of 10 MHz. The NEC simulation can be run with `xnec2c` as

    xnec2c dipole.nec

On this example, FRICO and NEC should agree about the input parameters of the antenna (resistance and reactance) and the gain. At a frequency of 300 MHz you should observe an impedance of around 73 Ohm and a gain of about 2.15 dB. The dipole in this example is slightly shortened to remove the reactive part of the impedance.

In the `comparison` subdirectory of this example you can find the comparison of the results from NEC and from FRICO. They can be plotted with

    gnuplot --persist compare.plot


## References
* Balanis 3rd Edition section 4.6.