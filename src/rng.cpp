#include <headers/include.hpp>
#include <headers/rng.hpp>

namespace raytrace{ namespace setup{

// SplineSet
std::ostream & operator<< (std::ostream & os, const SplineSet & spl) noexcept
{
    os << "{" << spl.x << ", " << spl.y << ", " << spl.k <<"}";
    return os;
}


// FindRoot
Spline::FindRoot::FindRoot (const Spline & spl, double y)
{
    int i = 0;
    while (spl[i].y < y)
        i++;
    if(spl[i].y == y)
    {
        x = spl[i].x;
        err = 0;
    }
    else
    {
        double x1 = spl[i].x;
        double x0 = spl[i - 1].x;
        x = x1;
        err = x1 - x0;
        i = 0;
        while (err > TOL && i < MAX_IT)
        {
            double x2 = x1 - (spl.funcval(x1) - y) * (x1 - x0) / ((spl.funcval(x1) - y) - (spl.funcval(x0) - y));
            if ((spl.funcval(x2) - y) * (spl.funcval(x0) - y) < 0)
                x1 = x2;
            else
                x0 = x2;
            err = abs(x - x2);
            x = x2;
            i++;
        }
        if (i >= MAX_IT)
            throw std::runtime_error ("FindRoot : no convergence");
    }
}


// Spline
Spline::Spline (const std::vector<double> & x, const std::vector<double> & y)
{
    if (x.size() != y.size())
        throw std::invalid_argument("Spline (std::vector<double> & x, std::vector<double> & y) : array sizes are not equal");
    if(!is_sorted(x.begin(), x.end()))
        throw std::invalid_argument("Spline (std::vector<double> & x, std::vector<double> & y) : x array is not sorted");
    if(!is_sorted(y.begin(), y.end()))
        throw std::invalid_argument("Spline (std::vector<double> & x, std::vector<double> & y) : y array is not sorted");
    int n = x.size() - 1;
    std::vector<double> dx, dy;
    for (int i = 0; i < n; i ++)
    {
        dx.emplace_back(x[i + 1] - x[i]);
        dy.emplace_back(y[i + 1] - y[i]);
    }
    std::vector<double> a (n + 1);
    std::vector<double> b (n + 1);
    std::vector<double> c (n + 1);
    std::vector<double> d (n + 1);
    std::vector<double> alpha (n + 1);
    std::vector<double> betta (n + 1);
    a[0] = 0.0;
    b[0] = 2 / dx[0];
    c[0] = 1 / dx[0];
    d[0] = 3 * dy[0] / pow(dx[0], 2);
    alpha[0] = c[0] / b[0];
    betta[0] = d[0] / b[0];
    a[n] = 1 / dx[n - 1];
    b[n] = 2 / dx[n - 1];
    c[n] = 0.0;
    d[n] = 3 * dy[n - 1] / pow(dx[n - 1], 2);
    for (int i = 1; i < n; i++)
    {
        a[i] = 1 / dx[i - 1];
        b[i] = 2 / dx[i - 1] + 2 / dx[i];
        c[i] = 1 / dx[i];
        d[i] = 3 * (dy[i - 1] / pow(dx[i - 1], 2) + dy[i] / pow(dx[i], 2));
    }
    for (int i = 1; i < n + 1; i++){
        alpha[i] = c[i] / (b[i] - a[i] * alpha[i - 1]);
        betta[i] = (d[i] - a[i] * betta[i - 1]) / (b[i] - a[i] * alpha[i - 1]);
    }
    SplineVec::emplace_back(SplineSet(x[n], y[n], betta[n]));
    for (int j = n - 1; j >= 0; j--)
        SplineVec::emplace_back(SplineSet(x[j], y[j], betta[j] - alpha[j] * SplineVec::operator[](n - 1 - j).k));
    std::reverse(SplineVec::begin(), SplineVec::end());
}

double Spline::funcval (double x) const
{
    if (x < SplineVec::begin()->x || x > SplineVec::back().x)
        throw std::out_of_range ("Spline::funcval : argument is out of range");
    auto pt2 = lower_bound(SplineVec::begin(), SplineVec::end(), SplineSet(x));
    if (x == pt2->x)
        return pt2->y;
    auto pt1 = pt2 - 1;
    double t = (x - pt1->x) / (pt2->x - pt1->x);
    double a = pt1->k * (pt2->x - pt1->x) - (pt2->y - pt1->y);
    double b = - pt2->k * (pt2->x - pt1->x) + (pt2->y - pt1->y);
    return (1.0 - t) * pt1->y + t * pt2->y + t * (1.0 - t) * (a * (1.0 - t) + b * t);
}

double Spline::deriv (double x) const
{
    if (x < SplineVec::begin()->x || x > SplineVec::back().x)
        throw std::out_of_range ("Spline::deriv : argument is out of range");
    auto pt2 = lower_bound(SplineVec::begin(), SplineVec::end(), SplineSet(x));
    if (x == pt2->x)
        return pt2->k;
    auto pt1 = pt2 - 1;
    double t = (x - pt1->x) / (pt2->x - pt1->x);
    double a = pt1->k * (pt2->x - pt1->x) - (pt2->y - pt1->y);
    double b = - pt2->k * (pt2->x - pt1->x) + (pt2->y - pt1->y);
    return (pt2->y - pt1->y) / (pt2->x - pt1->x) + (1.0 - 2.0 * t) * (a * (1.0 - t) + b * t) / (pt2->x - pt1->x) + t * (1.0 - t) * (b - a) / (pt2->x - pt1->x);
}

double Spline::arg (double y) const
{
    if (y < SplineVec::begin()->y || y > SplineVec::back().y)
        throw std::out_of_range ("Spline::arg : argument is out of range");
    return FindRoot(*this, y).x;
}


// RNG
// RNG::RNG (std::function<double (double)> pdf, double low, double high, unsigned int resolution) : low_(low), high_(high), res_(resolution)
// {
//     if (low_ >= high_)
//         throw std::invalid_argument ("RNG (pdf, low, high, resolution) : low value is greater than high value");
//     std::vector<double> val_;
//     std::vector<double> cdf_;
//     std::vector<double> pdf_;
//     double step = (high_ - low_) / res_;
//     for (unsigned int i = 0; i < res_; i++)
//         pdf_.emplace_back(pdf (low_ + i * step));
//     double pdf_max = *std::max_element(pdf_.begin(), pdf_.end());
//     double x = low_;
//     double y = 0.0;
//     while(x < high_)
//     {
//         double dx = step / (1.0 + 5.0 * std::tanh(2 * pdf(x) / pdf_max));
//         y += pdf(x) * dx;
//         val_.emplace_back(x);
//         cdf_.emplace_back(y);
//         x += dx;
//     }
//     cdf_.emplace_back(y + pdf(high_) * (high_ - val_.back()));
//     val_.emplace_back(high_);
//     std::transform(cdf_.begin(), cdf_.end(), cdf_.begin(), [&cdf_] (double x) {return x / cdf_.back(); });
//     CDF = setup::Spline(val_, cdf_);
// }

}}