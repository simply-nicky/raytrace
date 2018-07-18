#include <headers/trace.hpp>
#include <headers/make_trace.hpp>

using namespace raytrace;

int main()
{
    try
    {       
        fs::path dir = fs::absolute("CorrLength_Series_R1000");
        fs::create_directory(dir);
        tracing tr;
        tr.root_path() = dir;
        for(int i = 1; i < 10; i++)
        {
            tr.surface().CorrLength() = 0.15 * i;
            tr.run<PlaneBeam3D, TracePoint>(100000000, true);
        }
    }
    catch (std::exception & e)
    {
        std::cout << e.what() << std::endl;
    }
    return 0;
}