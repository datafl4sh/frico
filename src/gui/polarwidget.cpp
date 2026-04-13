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

#include "polarwidget.h"
#include <QtCharts/QValueAxis>
#include <QtMath>

PolarWidget::PolarWidget(QWidget *parent)
    : QChartView(parent),
      m_chart(new QPolarChart()),
      m_series(new QLineSeries())
{
    setRenderHint(QPainter::Antialiasing);

    m_chart->legend()->hide();
    m_chart->addSeries(m_series);
    m_chart->setTitle("Gain (dBi), YZ plane");
    //m_chart->legend()->show();

    setupAxes();

    setChart(m_chart);
}

void PolarWidget::setupAxes()
{
    QValueAxis *angularAxis = new QValueAxis();
    angularAxis->setRange(0, 360);
    angularAxis->setTickCount(13);
    angularAxis->setLabelFormat("%d");

    QValueAxis *radialAxis = new QValueAxis();
    radialAxis->setRange(-40, 20);
    radialAxis->setTickCount(6);

    m_chart->addAxis(angularAxis, QPolarChart::PolarOrientationAngular);
    m_chart->addAxis(radialAxis, QPolarChart::PolarOrientationRadial);

    m_series->attachAxis(angularAxis);
    m_series->attachAxis(radialAxis);
    //m_series->setName("xxx");
}

void PolarWidget::setData(const std::vector<double> &data)
{
    m_series->clear();

    for (int i = 0; i < 360; i++) {
        m_series->append({double(i), 10*log10(data[i])});
    }
}