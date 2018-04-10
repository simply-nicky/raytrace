#include <headers/include.hpp>
#include <headers/geometry.hpp>

namespace raytrace{ namespace geometry{

// Point
std::ostream & operator<< (std::ostream & os, const Point & pt) noexcept
{
    os << pt.x_ << "," << pt.y_ << "," << pt.z_;
    return os;
}


double distance (const Point & pt1, const Point & pt2) noexcept
{
    return sqrt(pow(pt1.x_ - pt2.x_, 2) + pow(pt1.y_ - pt2.y_, 2) + pow(pt1.z_ - pt2.z_, 2));
}


// Vector
Vector Vector::operator/ (double k) const
{
    if (k == 0.0)
        throw std::invalid_argument ("Vector::operator/ : division by zerro encountered");
    return Vector(this->vx_ / k, this->vy_ / k, this->vz_ / k);
}

Vector VectorProduct(const Vector & a, const Vector & b) noexcept
{
    return Vector(a.vy_ * b.vz_ - a.vz_ * b.vy_, a.vz_ * b.vx_ - a.vx_ * b.vz_, a.vx_ * b.vy_ - a.vy_ * b.vx_);
}

std::ostream & operator<< (std::ostream & os, const Vector & v) noexcept
{
    os << v.theta() << "," << v.phi();
    return os;
}


// Line
std::ostream & operator<< (std::ostream & os, const Line & line) noexcept
{
    os << line.point() << "," << line.direction();
    return os;
}

}}