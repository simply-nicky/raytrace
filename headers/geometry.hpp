#ifndef GEOMETRY_RAYTRACE_
#define GEOMETRY_RAYTRACE_
#include <headers/include.hpp>

namespace raytrace { namespace geometry {

    class Point
    {
        private:
            static constexpr double TOL = 1e-8;
            double x_, y_, z_;
        public:
            explicit Point (double x = 0.0, double y = 0.0, double z = 0.0) noexcept : x_(x), y_(y), z_(z) {}
            friend std::ostream & operator<< (std::ostream & os, const Point & pt) noexcept;
            friend double distance (const Point & pt1, const Point & pt2) noexcept;
            bool operator== (const Point & pt) const noexcept {return distance(*this, pt) < TOL; }
            bool operator!= (const Point & pt) const noexcept {return distance(*this, pt) > TOL; }
            double x() const noexcept {return x_; }
            double y() const noexcept {return y_; }
            double z() const noexcept {return z_; }
    };

    class Vector
    {
        private:
            double vx_, vy_, vz_;
        public:
            Vector() = default;
            explicit Vector (double theta, double phi = Constants::pi / 2.0) noexcept
            : vx_(std::cos(theta) * std::cos(phi)), vy_(std::cos(theta) * std::sin(phi)), vz_(std::sin(theta)) {}
            Vector (double vx, double vy, double vz) noexcept : vx_(vx), vy_(vy), vz_(vz) {}
            Vector (const Point & pt1, const Point & pt2) noexcept : vx_(pt2.x() - pt1.x()), vy_(pt2.y() - pt1.y()), vz_(pt2.z() - pt1.z()) {} 
            Vector operator+(const Vector & v) const noexcept {return Vector(this->vx_ + v.vx_, this->vy_ + v.vy_, this->vz_ + v.vz_); }
            Vector operator-(const Vector & v) const noexcept {return Vector(this->vx_ - v.vx_, this->vy_ - v.vy_, this->vz_ - v.vz_); }
            Vector operator-() const noexcept {return Vector(-this->vx_, -this->vy_, -this->vz_); }
            double operator*(const Vector & v) const noexcept {return this->vx_ * v.vx_ + this->vy_ * v.vy_ + this->vz_ * v.vz_; }
            Vector operator*(double k) const noexcept {return Vector(k * this->vx_, k * this->vy_, k * this->vz_); }
            Vector operator/(double k) const;
            friend Vector operator*(double k, const Vector & v) noexcept {return Vector(k * v.vx_, k * v.vy_, k * v.vz_); }
            friend std::ostream & operator<< (std::ostream & os, const Vector & v) noexcept;
            friend Vector VectorProduct (const Vector & a, const Vector & b) noexcept;
            double theta() const noexcept {return atan(this->vz_ / sqrt(this->vx_ * this->vx_ + this->vy_ * this->vy_)); }
            double phi() const noexcept {return atan2(this->vy_, this->vx_); }
            double Norm() const noexcept {return sqrt(vx_ * vx_ + vy_ * vy_ + vz_ * vz_); }
            Vector Normalize() {return *this / this->Norm(); }
            double x() const noexcept {return vx_; }
            double y() const noexcept {return vy_; }
            double z() const noexcept {return vz_; }
            friend double VectorAngle (const Vector & v1, const Vector & v2) {return acos(v1 * v2 / v1.Norm() / v2.Norm()); }
    };

    class RotationMatrix
    {
        private:
            Vector m_ [3];
        public:
            RotationMatrix (double theta) noexcept
            : m_{Vector (1.0, 0.0, 0.0), Vector (0.0, cos(theta), -sin(theta)), Vector (0.0, sin(theta), cos(theta))} {}
            RotationMatrix (double theta, Vector v) noexcept
            : m_
            {
                Vector
                {
                    cos(theta) + pow(v.Normalize().x(), 2) * (1 - cos(theta)),
                    v.Normalize().x() * v.Normalize().y() * (1 - cos(theta)) - v.Normalize().z() * sin(theta),
                    v.Normalize().x() * v.Normalize().z() * (1 - cos(theta)) + v.Normalize().y() * sin(theta)
                },
                Vector
                {
                    v.Normalize().y() * v.Normalize().x() * (1 - cos(theta)) + v.Normalize().z() * sin(theta),
                    cos(theta) + pow(v.Normalize().y(), 2) * (1 - cos(theta)),
                    v.Normalize().y() * v.Normalize().z() * (1 - cos(theta)) - v.Normalize().x() * sin(theta)
                },
                Vector
                {
                    v.Normalize().z() * v.Normalize().x() * (1 - cos(theta)) - v.Normalize().y() * sin(theta),
                    v.Normalize().z() * v.Normalize().y() * (1 - cos(theta)) + v.Normalize().x() * sin(theta),
                    cos(theta) + pow(v.Normalize().z(), 2) * (1 - cos(theta))
                }
            } {}       
            Vector operator* (const Vector & v) const noexcept {return Vector (v * m_[0], v * m_[1], v * m_[2]); }
    };

    class Sphere
    {
        public:
            enum part {LOWER, UPPER};
        private:
            Point pt0_;
            double r_, d_;                           //mm
            part part_;
        public:
            explicit Sphere (double x0, double y0, double z0,
            double r, double d, part p) noexcept : pt0_(x0, y0, z0), r_(r), d_(d), part_(p) {}
            Sphere (const Point & pt0, double r, double d, part p) noexcept : pt0_(pt0), r_(r), d_(d), part_(p) {}
            double psi () const {return 2 * asin(d_ / 2.0 / r_); }
            Vector NormVec (const Point & pt) const {return Vector(pt0_, pt).Normalize(); }
            Point center() const noexcept {return pt0_; }
            double radius() const noexcept {return r_; }
            double diameter() const noexcept {return d_; }
            part part() const noexcept {return part_; }
    };

    class SphPoint
    {
        private:
            Point pt_;
            Sphere sph_;
        public:
            SphPoint(double x, double y, double z, const Sphere & sph) : pt_(x, y, z), sph_(sph) {}
            const Point & point() const {return pt_; }
            const Sphere & sphere() const {return sph_; }
            Vector NormVec() const {return sph_.NormVec(pt_); }
    };

    class Line
    {
        private:
            Point pt0_;
            Vector v_;
        protected:
            void SetPoint(const Point & pt) noexcept {pt0_ = pt; }
            void SetVector(const Vector & v) noexcept {v_ = v; }
        public:
            Line() = default;
            virtual ~Line() = default;
            Line (const Point & pt, const Vector & v) noexcept : pt0_(pt), v_(v) {}
            Line (double x, double y, double z, double theta, double phi) noexcept : pt0_(x, y, z), v_(theta, phi) {}
            friend std::ostream & operator<< (std::ostream & os, const Line & line) noexcept;
            Point point() const noexcept {return pt0_; }
            Vector direction() const noexcept {return v_; }
            virtual void message() const noexcept {std::cout << "Line\n"; }
    };
    
}}

#endif
