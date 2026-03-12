# FRICO - Friendly Radiation Integral COde
# 
# Copyright (c) 2025, 2026 Matteo Cicuttin - IV3IWE
# Politecnico di Torino
# Dipartimento di Scienze Matematiche "G. L. Lagrange"
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

set term qt size 1280,500
set multiplot title "Radiation patterns" layout 1,3

unset xtics
unset ytics
set rtics 10
set mrtics 1
set size square
set angle degrees
set grid polar
set polar grid theta [0:360]

set title "XY plane gain pattern"
set rrange [-50:20]
plot filename using 1:3 w l notitle

set title "YZ plane gain pattern"
set rrange [-50:20]
plot filename using 1:5 w l notitle

set title "XZ plane gain pattern"
set rrange [-50:20]
plot filename using 1:7 w l notitle
