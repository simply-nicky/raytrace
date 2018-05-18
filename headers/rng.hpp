#ifndef RNG_RAYTRACE_
#define RNG_RAYTRACE_
#include <headers/include.hpp>

namespace raytrace { namespace setup {

    struct SplineSet                                                // {x, y = f(x), k = y'(x)}
    {
            double x;
            double y;
            double k;
            SplineSet(double x1 = 0.0, double x2 = 0.0, double x3 = 0.0) noexcept : x(x1), y(x2), k(x3) {}
            friend std::ostream & operator<< (std::ostream &, const SplineSet &) noexcept;
            friend bool operator< (const SplineSet & a, const SplineSet & b) noexcept {return (a.x < b.x); }
            friend bool operator> (const SplineSet & a, const SplineSet & b) noexcept {return (a.x > b.x); }
    };
           
    class Spline : private std::vector<SplineSet>
    {
        private:
            typedef std::vector<SplineSet> SplineVec;
            struct FindRoot
            {
                double x;
                double err;
                static constexpr int MAX_IT = 1000;
                static constexpr double TOL = 1e-14;
                FindRoot (const Spline & spl, double y = 0);
            };
        public:
            using std::vector<SplineSet>::back;
            using std::vector<SplineSet>::begin;
            using std::vector<SplineSet>::end;
            using std::vector<SplineSet>::operator[];
            using std::vector<SplineSet>::size;
            Spline() = default;
            Spline(const std::vector<double> & x, const std::vector<double> & y);
            double funcval (double x) const;
            double deriv (double x) const;
            double arg (double y) const;
    };

    static std::uniform_real_distribution<double> dist;
    static auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::mt19937 gen (seed);

    class RNG
    {
        private:
            static constexpr int RES = 2000;
            double low_, high_;
            unsigned int res_;
            Spline CDF;
        public:
            template <typename fn>
            RNG (fn && pdf, double low, double high, unsigned int resolution = RES) : low_(low), high_(high), res_(resolution)
            {
                if (low_ >= high_)
                    throw std::invalid_argument ("RNG (pdf, low, high, resolution) : low value is greater than high value");
                std::vector<double> val_;
                std::vector<double> cdf_;
                std::vector<double> pdf_;
                double step = (high_ - low_) / res_;
                for (unsigned int i = 0; i < res_; i++)
                    pdf_.emplace_back(pdf(low_ + i * step));
                double pdf_max = *std::max_element(pdf_.begin(), pdf_.end());
                double x = low_;
                double y = 0.0;
                while(x < high_)
                {
                    double dx = step / (1.0 + 5.0 * std::tanh(2 * pdf(x) / pdf_max));
                    y += pdf(x) * dx;
                    val_.emplace_back(x);
                    cdf_.emplace_back(y);
                    x += dx;
                }
                cdf_.emplace_back(y + pdf(high_) * (high_ - val_.back()));
                val_.emplace_back(high_);
                std::transform(cdf_.begin(), cdf_.end(), cdf_.begin(), [&cdf_] (double x) {return x / cdf_.back(); });
                CDF = setup::Spline(val_, cdf_);
            }
            
            template <typename Generator>
            double operator() (Generator && g) {return CDF.arg(dist(std::forward<Generator>(g))); }
    };

    template<class Func, class InputIt, class ... InputIts>
    void transform(Func func, InputIt first, InputIt last, InputIts... firsts)
    {
        while (first != last)
            func(*first++, (*firsts++)...);
    }

}}

#endif