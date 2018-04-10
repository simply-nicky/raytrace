#ifndef TRACE_RAYTRACE_
#define TRACE_RAYTRACE_
#include <headers/rng.hpp>
#include <headers/geometry.hpp>
#include <boost/iterator/indirect_iterator.hpp>

namespace raytrace {

    using geometry::Point;
    using geometry::Vector;
    using geometry::Sphere;
    using geometry::Line;
    using geometry::SphPoint;

    struct ExpGeometry
    {
        static constexpr float rho = 1000.0;                           //mm
        static constexpr float d = 60.0;                               //mm
        static constexpr float l1 = 195.0;                             //mm
        static constexpr float l2 = 185.0;                             //mm
        static constexpr float L = 100.0;                              //mm
    };

    class Surface
    {
        private:
            setup::Spline delta_;
            setup::Spline gamma_;                                       //Permettivity = 1 - delta + i * gamma
            const double RMSHeight_;                                    //RMS Height
            const double CorrLength_;                                   //Correlation length
            const double alpha_;                                        //alpha parameter in PSD ABC model
        public:
            Surface(double RMSHeight = Constants::RMSHeight,
            double CorrLength = Constants::CorrLength, double alpha = Constants::alpha, const std::string & str = "Al2O3.txt");
            std::complex<double> permettivity (double wl = Constants::WL) const;
            double Rf (double ang, double wl = Constants::WL) const;
            double TIS (double ang, double wl = Constants::WL) const;
            double Indicatrix1D (double th, double th0, double wl = Constants::WL) const;
            double Indicatrix2D (double th, double phi, double th0, double wl = Constants::WL) const;
            double PSD1D (double p) const noexcept;                              //PSD 1D ([micrometer ^ -1]) [micrometer ^ 3]
            double PSD2D (double p1, double p2) const noexcept;                  //PSD 2D ([micrometer ^ -1]) [micrometer ^ 4]
            double CritAng(double wl = Constants::WL) const;            //Critical angle
            double dmax (const Sphere & sph, double wl = Constants::WL) const {return sph.radius() * pow(CritAng(wl), 2) / 2.0; }
    };

    class ExpSetup
    {
        private:
            double RMSHeight_, CorrLength_, wl_, inc_ang_;
            Surface surf_;
            std::vector<Sphere> sphs_;
        public:
            ExpSetup (const std::vector<Sphere> & sphs, double inc_ang, double RMSHeight, double CorrLength, double wl) :
            surf_(Surface(RMSHeight, CorrLength)), sphs_(sphs), inc_ang_(inc_ang), RMSHeight_(RMSHeight), CorrLength_(CorrLength), wl_(wl) {}
            ExpSetup (double inc_ang = 0.0, double RMSHeight = Constants::RMSHeight, double CorrLength = Constants::CorrLength) :
            ExpSetup(std::vector<Sphere>{Sphere(0.0, 0.0, ExpGeometry::rho, ExpGeometry::rho, ExpGeometry::d, Sphere::LOWER)}, inc_ang, RMSHeight, CorrLength, Constants::WL) {}
            const Sphere & substrate() const noexcept {return *sphs_.cbegin(); }
            const std::vector<Sphere> & spheres() const noexcept {return sphs_; }
            const Surface & surface() const noexcept {return surf_; }
            const double IncAngle() const noexcept {return inc_ang_; }
            const double RMSHeight() const noexcept {return RMSHeight_; }
            const double CorrLength() const noexcept {return CorrLength_; }
            const double wavelength() const noexcept {return wl_; }
            void add_sphere(const Sphere & sph) noexcept {sphs_.emplace_back(sph); }
            void add_sphere(double x, double y, double radius, double diameter);
            void add_sphere(const std::vector<Sphere> & sphs) noexcept {sphs_.insert(sphs_.end(), sphs.cbegin(), sphs.cend()); }
    };

    class Beam : public Line
    {
        private:
            static constexpr double intersect_limit = 1e-3;
        public:
            using Line::Line;
            double Delta (const Sphere & sph) const noexcept;
            std::vector<SphPoint> Intersect(const Sphere & sph) const;
            std::vector<SphPoint> Intersect(const std::vector<Sphere> & sphs) const;
            bool is_intersect (const Sphere & sph) const noexcept;
            bool is_intersect (const std::vector<Sphere> & sphs) const noexcept;
            double IncAng (const SphPoint & spt) const;
            Vector SpecVec (const SphPoint & spt) const;
            virtual Vector ScatVec(const SphPoint & pt, const Surface & surf) const = 0;
    };

    class Beam2D : public Beam
    {
        public:
            using Beam::Beam;
            virtual Vector ScatVec(const SphPoint & spt, const Surface& surf) const override;
    };

    class SphBeam2D : public Beam2D
    {
        public:
            using Beam2D::Beam2D;
            SphBeam2D (const Sphere & sph, double inc_ang);
            SphBeam2D (const ExpSetup & setup) : SphBeam2D(setup.substrate(), setup.IncAngle()) {}
    };

    class PlaneBeam2D : public Beam2D
    {
        public:
            using Beam2D::Beam2D;
            PlaneBeam2D (const Sphere & sph, const Surface & surf, double inc_ang);
            PlaneBeam2D (const ExpSetup & setup) : PlaneBeam2D(setup.substrate(), setup.surface(), setup.IncAngle()) {}
    };

    class Beam3D : public Beam
    {
        public:
            using Beam::Beam;
            virtual Vector ScatVec(const SphPoint & spt, const Surface& surf) const override;
    };

    class SphBeam3D : public Beam3D
    {
        public:
            using Beam3D::Beam3D;
            SphBeam3D(const Sphere & sph, double inc_ang);
            SphBeam3D(const ExpSetup & setup) : SphBeam3D(setup.substrate(), setup.IncAngle()) {}
    };

    class PlaneBeam3D : public Beam3D
    {
        public:
            using Beam3D::Beam3D;
            PlaneBeam3D(const Sphere & sph, double inc_ang);
            PlaneBeam3D(const ExpSetup & setup) : PlaneBeam3D(setup.substrate(), setup.IncAngle()) {}
    };

    class Trace
    {
        private:
            std::vector<std::unique_ptr<Beam>> trace_;
            bool is_transmitted_;
        public:
            using iterator = boost::indirect_iterator<std::vector<std::unique_ptr<Beam>>::iterator>;
            using const_iterator = boost::indirect_iterator<std::vector<std::unique_ptr<Beam>>::const_iterator>;
            iterator begin() {return iterator(trace_.begin()); }
            iterator end() {return iterator(trace_.end()); }
            const_iterator begin() const {return const_iterator(trace_.begin()); }
            const_iterator end() const {return const_iterator(trace_.end()); }

            template <typename T>
            Trace (const T & line, const Surface & surf, const std::vector<Sphere> & sphs, double wl)
            {
                trace_.emplace_back(std::make_unique<T>(line));
                is_transmitted_ = trace_.back()->is_intersect(sphs); 
                if (is_transmitted_)
                {
                    double sw = setup::dist(setup::gen);
                    auto pts = trace_.back()->Intersect(sphs);
                    auto next_spt_ = *std::min_element(pts.cbegin(), pts.cend(), [](const SphPoint & a, const SphPoint & b) {return a.point().y() < b.point().y(); });
                    is_transmitted_ = sw < surf.Rf(trace_.back()->IncAng(next_spt_), wl);
                    while (is_transmitted_)
                    {
                        if (sw < surf.TIS(trace_.back()->IncAng(next_spt_), wl))
                            trace_.emplace_back(std::make_unique<T>(next_spt_.point(), trace_.back()->ScatVec(next_spt_, surf)));
                        else
                            trace_.emplace_back(std::make_unique<T>(next_spt_.point(), trace_.back()->SpecVec(next_spt_)));
                        pts = trace_.back()->Intersect(sphs);
                        if (pts.size() == 0)
                            break;
                        next_spt_ = *std::min_element(pts.cbegin(), pts.cend(), [](const SphPoint & a, const SphPoint & b) {return a.point().y() < b.point().y(); });
                        sw = setup::dist(setup::gen);
                        is_transmitted_ = sw < surf.Rf(trace_.back()->IncAng(next_spt_), wl);
                    }  
                }
            }

            template <typename T>
            Trace (const T & line, const ExpSetup & setup) : Trace(line, setup.surface(), setup.spheres(), setup.wavelength()) {}

            template <typename T>
            Trace (const T & line, const Surface & surf, const Sphere & sph, double wl) : Trace(line, surf, std::vector<Sphere> {sph}, wl) {}

            bool is_transmitted() const noexcept {return is_transmitted_; }
            int size() const noexcept {return trace_.size(); }
    };

    template <typename T>
    Trace make_trace(const ExpSetup & setup)
    {
        return Trace (T(setup), setup);
    }

}

#endif