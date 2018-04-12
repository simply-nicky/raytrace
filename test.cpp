#include <headers/make_trace.hpp>
#include <headers/trace.hpp>

using namespace raytrace;

int main()
{
    try
    {
        tracing tr1 (tracing::TRACE);
        auto t1 = std::chrono::high_resolution_clock::now();
        tr1.run<SphBeam3D>(1000000);
        auto t2 = std::chrono::high_resolution_clock::now();
        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms\n";
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
    }
    return 0;
}