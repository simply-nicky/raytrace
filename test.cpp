#include <headers/make_trace.hpp>
#include <headers/trace.hpp>

using namespace raytrace;

int main()
{
    try
    {
        Sphere sph {0.0, 0.0, 0.0, 1.0, 0.5, Sphere::LOWER};

    }
    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
    }
    return 0;
}