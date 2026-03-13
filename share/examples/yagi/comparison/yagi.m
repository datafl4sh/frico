freqMHz = 300;

width = 0.005;

Reflen  = 150/freqMHz;
space1  = 0.20*300/freqMHz;
Diplen  = 142/freqMHz;
space2  = 0.09*300/freqMHz;
Dir1len = 133/freqMHz;
space3  = 0.15*300/freqMHz;
Dir2len = 131/freqMHz;
space4  = 0.2*300/freqMHz;
Dir3len = 129/freqMHz;

zref = -space1;
zdip = 0;
zdir1 = zdip+space2;
zdir2 = zdir1+space3;
zdir3 = zdir2+space4;

M = [ -Reflen/2,   0,     zref,   Reflen/2;
 -Diplen/2,   0,     zdip,   Diplen/2;
-Dir1len/2,   0,    zdir1,  Dir1len/2;
-Dir2len/2,   0,    zdir2,  Dir2len/2;
-Dir3len/2,   0,    zdir3,  Dir3len/2];

for i = 1:5
    row = M(i, :);
    x1 = row(1);
    y1 = 0;
    z1 = row(3);
    x2 = row(4);
    y2 = 0;
    z2 = row(3);
    radius = 0.01;
    fprintf("GW %d 15 %g %g %g %g %g %g %g\n", ...
        i, x1, y1, z1, x2, y2, z2, radius);
end