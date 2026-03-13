// FRICO - Friendly Radiation Integral COde
// 
// Copyright (c) 2025, 2026 Matteo Cicuttin - IV3IWE
// Politecnico di Torino
// Dipartimento di Scienze Matematiche "G. L. Lagrange"
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
// 
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

SetFactory("OpenCASCADE");

freqMHz = 300;

width = 0.006;

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

Rectangle(1) = {  -Reflen/2,   -width/2,     zref,    Reflen,  width};
Rectangle(2) = {  -Diplen/2,   -width/2,     zdip,  Diplen/2,  width};
Rectangle(3) = {          0,   -width/2,     zdip,  Diplen/2,  width};
Rectangle(4) = { -Dir1len/2,   -width/2,    zdir1,   Dir1len,  width};
Rectangle(5) = { -Dir2len/2,   -width/2,    zdir2,   Dir2len,  width};
Rectangle(6) = { -Dir3len/2,   -width/2,    zdir3,   Dir3len,  width};

Coherence;
MeshSize{ PointsOf{ Surface{:}; } } = 0.01;

Physical Curve("src", 1000) = {26};
