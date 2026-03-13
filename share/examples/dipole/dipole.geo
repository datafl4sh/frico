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

MHz = 300;

c = 299.792458;
lambda = c/MHz;
// For the 0.47 below see "Antenna Theory by Balanis, 3rd edition, section 4.6"
arm_len = 0.47*lambda/2;

width = 0.01;

Rectangle(1) = {     0,    -width/2, 0,  arm_len, width};
Rectangle(2) = { -arm_len, -width/2, 0,  arm_len, width};

Coherence;
MeshSize{ PointsOf{ Surface{:}; } } = width/2;


// Define the Physical curve where the delta-gap is applied
Physical Curve("src", 1000) = {4};
