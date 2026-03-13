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
set multiplot title "FRICO vs. NEC comparison - Impedance and SWR plots" layout 2,1
set key bottom right

set title "Impedance"
plot 'yagi_frico.txt' using (($2)/1e6):3 w lp t 'R frico', \
     'yagi_frico.txt' using (($2)/1e6):4 w lp t 'X frico', \
     'yagi_nec.txt' using 1:2 w lp t 'R nec', \
     'yagi_nec.txt' using 1:3 w lp t 'X nec'

set title "SWR"
plot 'yagi_frico.txt' using (($2)/1e6):5 w lp t 'SWR frico', \
     'yagi_nec.txt' using 1:6 w lp t 'SWR nec'

