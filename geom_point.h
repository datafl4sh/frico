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
#include <cmath>

#include <Eigen/Dense>
namespace mommy {

class vec3
{
    double x_;
    double y_;
    double z_;

public:
    vec3()
        : x_(0.0), y_(0.0), z_(0.0)
    {}

    inline vec3(double x, double y, double z)
        : x_(x), y_(y), z_(z)
    {}

    vec3(const vec3& other) = default;

    inline double x() const { return x_; }
    inline double y() const { return y_; }
    inline double z() const { return z_; }

    inline vec3 operator-() const {
        return {-x_, -y_, -z_};
    }

    inline vec3& operator+=(const vec3& other) {
        x_ += other.x_;
        y_ += other.y_;
        z_ += other.z_;
        return *this;
    }

    inline vec3 operator+(const vec3& other) const {
        return { x_ + other.x_, y_ + other.y_, z_ + other.z_ };
    }

    vec3& operator-=(const vec3& other) {
        (*this) += -other;
        return *this;
    }

    vec3 operator-(const vec3& other) const {
        return { x_ - other.x_, y_ - other.y_, z_ - other.z_ };
    }

    vec3& operator*=(double s) {
        x_ *= s;
        y_ *= s;
        z_ *= s;
        return *this;
    }

    vec3 operator*(double s) const {
        return { x_ * s, y_ * s, z_ * s };
    }

    vec3& operator/=(double s) {
        x_ /= s;
        y_ /= s;
        z_ /= s;
        return *this;
    }

    vec3 operator/(double s) const {
        return { x_ / s, y_ / s, z_ / s };
    }

    Eigen::Matrix<double, 3, 1> to_eigen() const {
        return {x_, y_, z_};
    }
};

inline vec3
operator*(double s, const vec3& p) {
    return p*s;
}

inline double
dot(const vec3& a, const vec3& b)
{
    return a.x()*b.x() + a.y()*b.y() + a.z()*b.z();
}

inline vec3
cross(const vec3& a, const vec3& b)
{
    return {
        a.y()*b.z() - a.z()*b.y(),
        a.z()*b.x() - a.x()*b.z(),
        a.x()*b.y() - a.y()*b.x()
    };
}

inline double
norm(const vec3& v)
{
    return std::sqrt( v.x()*v.x() + v.y()*v.y() + v.z()*v.z() );
    //return std::hypot(v.x(), v.y(), v.z());
}

inline std::ostream&
operator<<(std::ostream& os, const vec3& p)
{
    os << "(" << p.x() << ", " << p.y() << ", " << p.z() << ")";
    return os;
}

using point = vec3;

} // namespace mommy