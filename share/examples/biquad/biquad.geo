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

// This script creates the geometry for a biquad antenna. You can specify
// the center frequency in the MHz variable, of course in MHz. There is also
// the variable W that allows to specify the radiator width. At the end of
// the file you can tune the mesh sizes. Excitation is applied at Physical
// Curve with tag 1000.

///////////////////////////////////////////////////////////////////////////
// Antenna center frequency
MHz = 2450;

///////////////////////////////////////////////////////////////////////////
// From now on there is only code to create the geometry
SetFactory("OpenCASCADE");
c = 299.792458;
lambda = c/MHz;

L = 0.255*lambda;   // Arm length
W = 0.005;          // Radiator width
D = 0.12*lambda;    // Radiator distance from reflector

// Reflector dimensions
wh = lambda*1.2;
wv = lambda*1.1;

// Build the biquad
Louter = L + W/2;
Linner = L - W/2;

trans = Linner/2 + W/4;

Rectangle(1) = {-Louter/2, -Louter/2, 0, Louter, Louter, 0};
Rectangle(2) = {-Linner/2, -Linner/2, 0, Linner, Linner, 0};
Rectangle(3) = {-Louter/2, -Louter/2, 0, Louter, Louter, 0};
Rectangle(4) = {-Linner/2, -Linner/2, 0, Linner, Linner, 0};

Rectangle(5) = {-wh/2, -wv/2, -D, wh, wv, 0};

BooleanDifference(10) = { Surface{1}; Delete; }{ Surface{2}; Delete; };
BooleanDifference(11) = { Surface{3}; Delete; }{ Surface{4}; Delete; };

Rotate {{0, 0, 1}, {0, 0, 0}, Pi/4} { Surface{10}; }
Translate {-trans*Sqrt(2), 0, 0} { Surface{10}; }
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/4} { Surface{11}; }
Translate {trans*Sqrt(2), 0, 0} { Surface{11}; }

BooleanUnion(20) = { Surface{10}; Delete; }{ Surface{11}; Delete; };

Line(100) = {29, 36};
Line{100} In Surface{20};

// Define the Physical curve where the delta-gap is applied
Physical Curve("src", 1000) = {100};

// Set the mesh sizes
MeshSize{ PointsOf{ Surface{5}; } } = lambda/10;
MeshSize{ PointsOf{ Surface{20}; } } = W/2;
MeshSize{ PointsOf{ Curve{100}; } } = W/12;

Mesh 2;


