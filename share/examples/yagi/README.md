# 5-element yagi simulation with FRICO

The model for the 5-element yagi is contained in the file `yagi.geo`, which is a geometry specification suitable for GMSH. You can inspect the model by running

    gmsh dipole.geo

The antenna elements are modeled as rectangular strips, as FRICO does not support thin wires yet.

In order to run a frequency sweep from 280 MHz to 320 MHz at steps of 1 MHz you can launch

    frico -g yagi.geo -R 280e6:1e6:320e6 -A -s src

This will produce, as explained in the example for the biquad antenna:

* 41 `default_N.silo` files with the fields at the different frequencies
* 41 `polar_N.txt` files with the gain data at the different frequencies
* the file `port_sweep.txt` with the impedance data at the input port

## Comparison with NEC

This example includes also a NEC file for the same antenna, configured for a frequency sweep from 280 MHz to 320 MHz at steps of 1 MHz. The NEC simulation can be run with `xnec2c` as

    xnec2c dipole.nec

On this example, FRICO and NEC should agree about the input parameters of the antenna (resistance and reactance) and the gain. At a frequency of 300 MHz you should observe an impedance of around 50-j20 Ohm and a gain of about 9.81 dB. The elements in FRICO are quite thicker compared to those of NEC, this is because the different geometry (wires in NEC, rectangular strips in FRICO).

In the `comparison` subdirectory of this example you can find the comparison of the results from NEC and from FRICO. They can be plotted with

    gnuplot --persist compare.plot

There is also a file named `yagi.m` to generate automatically the NEC cards if you want to change frequency.

## References
* Nuova Elettronica, Le Antenne Riceventi e Trasmittenti