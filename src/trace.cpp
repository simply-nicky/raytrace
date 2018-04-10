#include <headers/trace.hpp>

using namespace raytrace;


// Surface
Surface::Surface (double RMSHeight, double CorrLength, double alpha, const std::string & str) : RMSHeight_(RMSHeight), CorrLength_(CorrLength), alpha_(alpha)
{
    std::ifstream fin;
    fin.open(str);
    if(!fin.good())
    {
        std::cout << "Surface : could not open permettivity file\n";
        exit(EXIT_FAILURE);
    }
    std::vector<double> wl;
    std::vector<double> delta;
    std::vector<double> gamma;
    while (fin.good())
    {
        double temp;
        fin >> temp;
        wl.emplace_back(temp);
        fin >> temp;
        delta.emplace_back(temp);
        fin >> temp;
        gamma.emplace_back(temp);
    }
    if (!fin.good())
    {
        if (fin.eof())
        {
            delta_ = setup::Spline (wl, delta);
            gamma_ = setup::Spline (wl, gamma);
        }
        else if (fin.fail())
        {
            std::cerr << "Surface : input terminated by data mismatch.\n";
            exit(EXIT_FAILURE);
        }
        else
        {
            std::cerr << "Surface : input terminated by unknown reason.\n";
            exit(EXIT_FAILURE);
        }
    }
}

std::complex<double> Surface::permettivity (double wl) const
{
    return std::complex<double> (1.0 - delta_.funcval(wl), - gamma_.funcval(wl));
}

double Surface::Rf (double ang, double wl) const
{
    if (ang < 0.0 || ang > Constants::pi / 2.0)
        throw std::out_of_range ("Surface::Rf : angle is out of range");
    std::complex<double> rf ((sin(ang) - sqrt(permettivity(wl) - cos(ang) * cos(ang))) / (sin(ang) + sqrt(permettivity(wl) - cos(ang) * cos(ang))));
    return abs(rf) * abs(rf);
}

double Surface::TIS (double ang, double wl) const
{
    if (ang < 0.0 || ang > Constants::pi / 2.0)
        throw std::out_of_range ("Surface::TIS : angle is out of renge");
    return Rf(ang, wl) * (1.0 - std::exp(- pow(4 * Constants::pi * RMSHeight_ * sin(ang) / wl, 2)));
}

double Surface::PSD1D (double p) const noexcept
{
    return 2.0 / sqrt(Constants::pi) * pow(10.0, -6) * tgamma(alpha_ + 0.5) / tgamma(alpha_) * pow(RMSHeight_, 2)
    * CorrLength_ / pow(1.0 + p * p * CorrLength_ * CorrLength_, alpha_ + 0.5);
}

double Surface::PSD2D (double p1, double p2) const noexcept
{
    return pow(10.0, -6) * pow(CorrLength_, 2) * pow(RMSHeight_, 2) * alpha_ / Constants::pi / pow(1.0 + (p1 * p1 + p2 * p2) * pow(CorrLength_, 2), 1.0 + alpha_);
}

double Surface::Indicatrix1D (double th, double th0, double wl) const
{
    if (th < 0 || th > Constants::pi / 2)
        throw std::out_of_range ("Surface::Indicatrix1D : angle th is out of range");
    if (th0 < 0 || th0 > Constants::pi / 2)
        throw std::out_of_range ("Surface::Indicatrix1D : angle th0 is out of range");
    double k = 2.0 * Constants::pi / wl;
    std::complex<double> t (2.0 * sin(th) / (sin(th) + sqrt(std::complex<double>(real(permettivity(wl)) - cos(th) * cos(th)))));
    std::complex<double> t0 (2.0 * sin(th0) / (sin(th0) + sqrt(std::complex<double>(real(permettivity(wl)) - cos(th0) * cos(th0)))));
    return pow(k, 3) * pow(abs(1.0 - real(permettivity(wl))), 2) * pow(abs(t * t0), 2) / 16.0 / Constants::pi / sin(th0)
    / sqrt(cos(th0) * cos(th)) * pow(10.0, 9) * PSD1D(1 / wl * pow(10.0, 3) * abs(cos(th0) - cos(th)));
}

double Surface::Indicatrix2D (double th, double phi, double th0, double wl) const
{
    if (th < 0 || th > Constants::pi / 2)
        throw std::out_of_range ("Surface::Indicatrix1D : angle th is out of range");
    if (phi < 0 || phi > Constants::pi / 2)
        throw std::out_of_range ("Surface::Indicatrix1D : angle phi is out of range");
    if (th0 < 0 || th0 > Constants::pi / 2)
        throw std::out_of_range ("Surface::Indicatrix1D : angle th0 is out of range");       
    double k = 2.0 * Constants::pi / wl;
    std::complex<double> t (2.0 * sin(th) / (sin(th) + sqrt(std::complex<double>(real(permettivity(wl)) - cos(th) * cos(th)))));
    std::complex<double> t0 (2.0 * sin(th0) / (sin(th0) + sqrt(std::complex<double>(real(permettivity(wl)) - cos(th0) * cos(th0)))));
    return pow(k, 4) * pow(abs(1.0 - real(permettivity(wl))), 2) * pow(abs(t * t0), 2) / pow(4 * Constants::pi, 2) / sin(th0)
    / sqrt(cos(th0) * cos(th)) * pow(10.0, 12) * PSD2D(1 / wl * pow(10.0, 3) * (cos(th) * cos(phi) - cos(th0)), 1 / wl * pow(10.0, 3) * cos(th) * cos(phi));
}

double Surface::CritAng(double wl) const
{
    return sqrt(1.0 - real(permettivity(wl)));
}


// ExpSetup
void ExpSetup::add_sphere(double x, double y, double radius, double diameter) noexcept
{

}


// Beam
double Beam::Delta (const Sphere & sph) const noexcept
{
    return pow((point().x() - sph.center().x()) * cos(direction().theta()) * cos(direction().phi()) + (point().y() - sph.center().y()) * cos(direction().theta()) * sin(direction().phi()) + (point().z() - sph.center().z()) * sin(direction().theta()), 2)
    - pow(point().x() - sph.center().x(), 2) - pow(point().y() - sph.center().y(), 2) - pow(point().z() - sph.center().z(), 2) + sph.radius() * sph.radius();
}

std::vector<SphPoint> Beam::Intersect (const Sphere & sph) const
{
    std::vector<SphPoint> pts
    {
        SphPoint
        {
            point().x() - cos(direction().theta()) * cos(direction().phi()) * ((point().x() - sph.center().x()) * cos(direction().theta()) * cos(direction().phi()) + (point().y() - sph.center().y()) * cos(direction().theta()) * sin(direction().phi()) + (point().z() - sph.center().z()) * sin(direction().theta()) - sqrt(Delta(sph))),
            point().y() - cos(direction().theta()) * sin(direction().phi()) * ((point().x() - sph.center().x()) * cos(direction().theta()) * cos(direction().phi()) + (point().y() - sph.center().y()) * cos(direction().theta()) * sin(direction().phi()) + (point().z() - sph.center().z()) * sin(direction().theta()) - sqrt(Delta(sph))),
            point().z() - sin(direction().theta()) * ((point().x() - sph.center().x()) * cos(direction().theta()) * cos(direction().phi()) + (point().y() - sph.center().y()) * cos(direction().theta()) * sin(direction().phi()) + (point().z() - sph.center().z()) * sin(direction().theta()) - sqrt(Delta(sph))),
            sph
        },
        SphPoint
        {
            point().x() - cos(direction().theta()) * cos(direction().phi()) * ((point().x() - sph.center().x()) * cos(direction().theta()) * cos(direction().phi()) + (point().y() - sph.center().y()) * cos(direction().theta()) * sin(direction().phi()) + (point().z() - sph.center().z()) * sin(direction().theta()) + sqrt(Delta(sph))),
            point().y() - cos(direction().theta()) * sin(direction().phi()) * ((point().x() - sph.center().x()) * cos(direction().theta()) * cos(direction().phi()) + (point().y() - sph.center().y()) * cos(direction().theta()) * sin(direction().phi()) + (point().z() - sph.center().z()) * sin(direction().theta()) + sqrt(Delta(sph))),
            point().z() - sin(direction().theta()) * ((point().x() - sph.center().x()) * cos(direction().theta()) * cos(direction().phi()) + (point().y() - sph.center().y()) * cos(direction().theta()) * sin(direction().phi()) + (point().z() - sph.center().z()) * sin(direction().theta()) + sqrt(Delta(sph))),
            sph
        }
    };
    pts.erase(
        std::remove_if(
            pts.begin(), 
            pts.end(), 
            [this] (const SphPoint & spt) {return 
                spt.point() == point()
                ||
                (spt.sphere().part() == Sphere::LOWER && Vector(spt.sphere().center(), spt.point()).theta() > (spt.sphere().psi() / 2.0 - Constants::pi / 2.0))
                ||
                (spt.sphere().part() == Sphere::UPPER && Vector(spt.sphere().center(), spt.point()).theta() < (Constants::pi / 2.0 - spt.sphere().psi() / 2.0)); }
            ),
        pts.cend()
    );
    return pts;
}

std::vector<SphPoint> Beam::Intersect(const std::vector<Sphere> & sphs) const
{
    std::vector<SphPoint> spts;
    for(const auto & sph_ : sphs)
    {
        auto pts_ = Intersect(sph_);
        spts.insert(spts.end(), std::make_move_iterator(pts_.cbegin()), std::make_move_iterator(pts_.cend()));
    }
    return spts;
}

bool Beam::is_intersect(const Sphere & sph) const noexcept
{
    return Delta(sph) > intersect_limit;
}

bool Beam::is_intersect(const std::vector<Sphere> & sphs) const noexcept
{
    int is_i_ {};
    for(const auto & sph : sphs)
        is_i_ += is_intersect(sph);
    return is_i_;
}

double Beam::IncAng (const SphPoint & spt) const
{
    return abs(Constants::pi / 2.0 - VectorAngle(spt.NormVec(), direction()));
}

Vector Beam::SpecVec (const SphPoint & spt) const
{
    return direction() - 2 * (spt.NormVec() * direction()) * spt.NormVec();
}


// Beam2D
SphBeam2D::SphBeam2D (const Sphere & sph, double inc_ang)
{
    SetPoint(
        Point {
            sph.center().x(),
            sph.center().y() - ExpGeometry::l1,
            sph.center().z() - sph.diameter() / 2.0 / tan(sph.psi() / 2.0) + tan(sph.psi() / 2.0 + inc_ang) * (ExpGeometry::l1 - sph.diameter() / 2.0)
        });
        double th1 = -sph.psi() / 2.0 - inc_ang;
        double th2 = atan(-tan(sph.psi() / 2.0 + inc_ang) * (ExpGeometry::l1 - sph.diameter() / 2.0) / (ExpGeometry::l1 + sph.diameter() / 2.0));
        SetVector(Vector((th2 - th1) * setup::dist(setup::gen) + th1));
}

PlaneBeam2D::PlaneBeam2D (const Sphere & sph, const Surface & surf, double inc_ang)
{
    SetPoint (Point {
        sph.center().x(),
        sph.center().y() - ExpGeometry::l1,
        sph.center().z() - sph.diameter() / 2.0 / tan(sph.psi() / 2.0) + tan(sph.psi() / 2.0 + inc_ang) * (ExpGeometry::l1 - sph.diameter() / 2.0) + surf.dmax(sph) * setup::dist(setup::gen)
    });
    SetVector(Vector(-sph.psi() / 2.0 - inc_ang));
}

Vector Beam2D::ScatVec (const SphPoint & spt, const Surface & surf) const
{
    setup::RNG ScatAng ([&] (double x) {return surf.Indicatrix1D(x, IncAng(spt)); });
    return geometry::RotationMatrix(ScatAng(setup::gen)) * (direction() - (spt.NormVec() * direction()) * spt.NormVec());
}


// Beam3D
SphBeam3D::SphBeam3D (const Sphere & sph, double inc_ang)
{
    double x0 = ExpGeometry::L * setup::dist(setup::gen) - ExpGeometry::L / 2;
    SetPoint(
    Point{
        sph.center().x() + x0,
        sph.center().y() - ExpGeometry::l1,
        sph.center().z() - sph.diameter() / 2.0 / tan(sph.psi() / 2.0) + tan(sph.psi() / 2.0 + inc_ang) * (ExpGeometry::l1 - sph.diameter() / 2.0)
    });
    double phi = Constants::pi / 2 - (2 * asin(sph.diameter() / 2.0 / sqrt(pow(ExpGeometry::l1, 2) + x0 * x0)) * setup::dist(setup::gen) + atan(-x0 / ExpGeometry::l1) - asin(sph.diameter() / 2.0 / sqrt(pow(ExpGeometry::l1, 2) + x0 * x0)));
    double th1 = atan((ExpGeometry::l1 - sph.diameter() / 2) * tan(sph.psi() / 2 + inc_ang) / (x0 * cos(phi) - ExpGeometry::l1 * sin(phi) + sqrt(pow(sph.diameter() / 2, 2) - pow(x0 * sin(phi) + ExpGeometry::l1 * cos(phi), 2))));
    double th2 = atan((ExpGeometry::l1 - sph.diameter() / 2) * tan(sph.psi() / 2 + inc_ang) / (x0 * cos(phi) - ExpGeometry::l1 * sin(phi) - sqrt(pow(sph.diameter() / 2, 2) - pow(x0 * sin(phi) + ExpGeometry::l1 * cos(phi), 2))));
    SetVector(Vector((th2 - th1) * setup::dist(setup::gen) + th1, phi));
}

PlaneBeam3D::PlaneBeam3D (const Sphere & sph, double inc_ang)
{
    double r = sph.diameter() / 2.0 * setup::dist(setup::gen);
    double phi = 2 * Constants::pi * setup::dist(setup::gen);
    SetPoint(
    Point{
        sph.center().x() - r * cos(phi),
        sph.center().y() - ExpGeometry::l1,
        sph.center().z() - sph.diameter() / 2.0 / tan(sph.psi() / 2.0) + tan(sph.psi() / 2.0 + inc_ang) * (ExpGeometry::l1 - r * sin(phi))
    });
    SetVector(Vector(-sph.psi() / 2.0 - inc_ang));
}

Vector Beam3D::ScatVec (const SphPoint & spt, const Surface & surf) const
{
    setup::RNG Theta ([&] (double x) {return surf.Indicatrix1D(x, IncAng(spt)); });
    double theta = Theta(setup::gen);
    setup::RNG Phi ([&] (double x) {return surf.Indicatrix2D(theta, x, IncAng(spt)); });
    Vector tau = direction() - (spt.NormVec() * direction()) * spt.NormVec();
    return geometry::RotationMatrix(Phi(setup::gen), spt.NormVec()) * (geometry::RotationMatrix(theta, VectorProduct(spt.NormVec(), tau)) * tau);
}