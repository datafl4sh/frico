% FRICO - Friendly Radiation Integral COde
% 
% Copyright (c) 2025, 2026 Matteo Cicuttin - IV3IWE
% Politecnico di Torino
% Dipartimento di Scienze Matematiche "G. L. Lagrange"
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.
% 
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% This file demonstrates how to obtain gain diagrams and SWR-vs-frequency
% plots from the HDF5 files written by FRICO.

clear;
clf;

hdf5_file = "biquad.h5";

sweep_steps = h5read(hdf5_file, "/frico/antenna_analysis/sweep_steps");

freq = zeros(1, sweep_steps);
swr  = zeros(1, sweep_steps);

for ii = 0:sweep_steps-1
    pathbase = strcat("/frico/antenna_analysis/", num2str(ii));
    freq(ii+1) = h5read(hdf5_file, strcat(pathbase, "/frequency"));
    swr(ii+1) = h5read(hdf5_file, strcat(pathbase, "/port/swr"));
    
    %%%%%
    subplot(2,3,1)
    title("XY plane gain pattern")
    gain_XY = h5read(hdf5_file, strcat(pathbase, "/gain_XY"));
    polarplot(deg2rad(0:359), 10*log10(gain_XY));
    rlim([-30,15]);
    hold on;

    %%%%%
    subplot(2,3,2)
    title("XZ plane gain pattern")
    gain_XZ = h5read(hdf5_file, strcat(pathbase, "/gain_XZ"));
    polarplot(deg2rad(0:359), 10*log10(gain_XZ));
    rlim([-30,15]);
    hold on;
    %%%%%
    subplot(2,3,3)
    title("YZ plane gain pattern")
    gain_YZ = h5read(hdf5_file, strcat(pathbase, "/gain_YZ"));
    polarplot(deg2rad(0:359), 10*log10(gain_YZ));
    rlim([-30,15]);
    hold on;
end

subplot(2,3,[4,5,6])
plot(freq, swr);
title("SWR plot");