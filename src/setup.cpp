#include <headers/include.hpp>
#include <headers/setup.hpp>

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
    assert(x.size() == y.size());
    assert(is_sorted(x.begin(), x.end()));
    assert(is_sorted(y.begin(), y.end()));
    int n = x.size() - 1;
    std::vector<double> dx, dy;
    for (int i = 0; i < n; i ++)
    {
        dx.emplace_back(x[i + 1] - x[i]);
        dy.emplace_back(y[i + 1] - y[i]);
    }
    std::vector<double> a (n + 1), b (n + 1), c (n + 1), d (n + 1), alpha (n + 1), betta (n + 1);
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
    assert(x >= SplineVec::begin()->x && x <= SplineVec::back().x);
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
    assert (x >= SplineVec::front().x && x <= SplineVec::back().x);
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
    assert(y >= SplineVec::begin()->y && y <= SplineVec::back().y);
    return FindRoot(*this, y).x;
}

}}