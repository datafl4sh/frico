/*
 * FRICO - Friendly Radiation Integral COde
 *
 * Copyright (c) 2025,2026 Matteo Cicuttin - IV3IWE
 * Politecnico di Torino
 * Dipartimento di Scienze Matematiche "G. L. Lagrange"
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <QApplication>

#include <highfive/highfive.hpp>

#include "polarwidget.h"

int main(int argc, char *argv[])
{
    HighFive::File file("biquad.h5", HighFive::File::ReadOnly);

    auto dataset = file.getDataSet("/frico/antenna_analysis/0/gain_YZ");
    auto data = dataset.read<std::vector<double>>();
    dataset.read(data);
    
    QApplication app(argc, argv);

    PolarWidget pw;
    pw.setWindowTitle("FRICO");
    pw.resize(500, 500);
    pw.show();
    pw.setData(data);

    return app.exec();
}