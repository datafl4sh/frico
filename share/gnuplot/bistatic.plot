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
set multiplot title "Bistatic RCS" layout 1,2

max(a,b) = (a < b) ? b : a

unset xtics
unset ytics
set rtics 10
set mrtics 1
set size square
set angle degrees
set grid polar
set polar grid theta [0:360]

rmin = -50
rmax = 0
set rrange [rmin:rmax]
set title "Elevation cut (theta)"
plot filename using 1:(max(rmin,10*log10($2))) w l title "Vertical", \
     filename using 1:(max(rmin,10*log10($3))) w l title "Horizontal" 

set title "Azimuth cut (phi)"
plot filename using 1:(max(rmin,10*log10($4))) w l title "Vertical", \
     filename using 1:(max(rmin,10*log10($5))) w l title "Horizontal"
