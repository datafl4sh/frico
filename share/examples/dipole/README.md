# Half-wavelength dipole simulation with FRICO

The model for the half-wavelength dipole is contained in the file `dipole.geo`, which is a geometry specification suitable for GMSH. You can inspect the model by running

    gmsh dipole.geo

This example includes also a NEC file for the same dipole, configured for a frequency sweep from 250 MHz to 350 MHz at steps of 10 MHz. The NEC simulation can be run with `xnec2c` as

    xnec2c dipole.nec

On this example, FRICO and NEC should agree about the input parameters of the antenna (resistance and reactance) and the gain. At a frequency of 300 MHz you should observe an impedance of around 73 Ohm and a gain of about 2.15 dB. The dipole in this example is slightly shortened to remove the reactive part of the impedance.

## References
* Balanis 3rd Edition section 4.6.