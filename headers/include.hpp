#ifndef INCLUDE_RAYTRACE_
#define INCLUDE_RAYTRACE_

#include <complex>
#include <vector>
#include <iostream>
#include <random>
#include <functional>
#include <chrono>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <ctime>
#include <future>
#include <atomic>
#include <memory>
#include <type_traits>

using std::atan;
using std::atan2;
using std::acos;
using std::asin;
using std::sqrt;
using std::pow;
using std::sin;
using std::cos;
using std::tgamma;
using std::abs;
using std::real;

namespace raytrace {

    struct Constants
    {
        static constexpr double pi = 3.1415926535897932;
        static constexpr double RMSHeight = 0.34815693559743155;    //nm
        static constexpr double CorrLength = 2.6948738995209154;    //micrometer
        static constexpr double alpha = 0.5029559338070592;         //PSD parameter for sapphire
        static constexpr double WL = 0.162;                         //nm
    };
}

#endif