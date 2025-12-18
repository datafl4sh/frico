/*
 * MoMmy - My experimental Method of Moments code
 *
 * Copyright (c) 2025, Matteo Cicuttin - IV3IWE
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

#pragma once

#include <iostream>

namespace mommy {

class point
{
    double x_;
    double y_;
    double z_;

public:
    point()
        : x_(0.0), y_(0.0), z_(0.0)
    {}

    point(double x, double y, double z)
        : x_(x), y_(y), z_(z)
    {}

    point(const point& other) = default;

    inline double x() const { return x_; }
    inline double y() const { return y_; }
    inline double z() const { return z_; }

    inline point operator-() const {
        return point(-x_, -y_, -z_);
    }

    inline point& operator+=(const point& other) {
        x_ += other.x_;
        y_ += other.y_;
        z_ += other.z_;
        return *this;
    }

    inline point operator+(const point& other) const {
        point ret = *this;
        ret += other;
        return ret;
    }

    point& operator-=(const point& other) {
        (*this) += -other;
        return *this;
    }

    point operator-(const point& other) const {
        return (*this) + (-other);
    }

    point& operator*=(double s) {
        x_ *= s;
        y_ *= s;
        z_ *= s;
        return *this;
    }

    point operator*(double s) const {
        point ret = *this;
        ret *= s;
        return ret;
    }

    point& operator/=(double s) {
        x_ /= s;
        y_ /= s;
        z_ /= s;
        return *this;
    }

    point operator/(double s) const {
        point ret = *this;
        ret /= s;
        return ret;
    }
};

inline point
operator*(double s, const point& p) {
    return p*s;
}

inline std::ostream&
operator<<(std::ostream& os, const point& p)
{
    os << "(" << p.x() << ", " << p.y() << ", " << p.z() << ")";
    return os;
}


}