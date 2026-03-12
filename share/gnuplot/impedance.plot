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

set term qt size 800,800
set multiplot title "Impedance and SWR plots" layout 2,1

set grid

set title "Impedance"
plot filename using (($2)/1e6):3 w l title "Re(Z)", \
     filename using (($2)/1e6):4 w l title "Im(Z)"

set title "SWR"
plot filename using (($2)/1e6):5 w l notitle
