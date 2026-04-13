# Radar cross section computations

## Bistatic RCS
The file `plt-55.geo` can be used to reproduce the results contained in
Section 13, page 13-4 of the MOM3D manual concerning the Bistatic RCS of a
square metallic plate of 5.5 inches of side length. The plate is illuminated
at 3 GHz and the RCS is computed at 0 and 45 degrees.

The bistatic RCS computation is selected via the `-b` flag, which takes three
colon-separated parameters: the distance of the radar, the theta angle (from z to x) and the phi angle (from x to y).

Here is the command line to run the Bistatic RCS simulation at 3.3 GHz with
source at 10 meters from the target, zero degrees of theta and zero degrees of
phi (therefore the source is at the point [0,0,10])

    frico -g plt-55.geo -f 3e9 -b 10:0:0 -A

Once the computation is finished, you should obtain a file named `bistatic_rcs_0.txt`: if you plot the first and second, and the first and third columns on a
polar plot you should get the Figure 13-4 (a) of the MOM3D manual. 

In a similar way it is possible to compute the bistatic RCS with the source at
45 degrees of theta (therefore at [10/sqrt(2), 0, 10/sqrt(2)]) with the command line

    frico -g plt-55.geo -f 3e9 -b 10:45:0 -A

Once the computation is finished, you should obtain a file named `bistatic_rcs_0.txt`: if you plot the first and second, and the first and third columns on a
polar plot you should get the Figure 13-4 (b) of the MOM3D manual. 