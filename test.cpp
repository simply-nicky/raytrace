#include <headers/trace.hpp>
#include <headers/make_trace.hpp>

using namespace raytrace;

int main()
{
    try
    {       
        fs::path dir = fs::absolute("CorrLength_Series_R1000_4");
        fs::create_directory(dir);
        tracing tr;
        tr.root_path() = dir;
        tr.surface().RMSHeight() = 2.0;
        for(int i = 1; i < 10; i++)
        {
            tr.surface().CorrLength() = 5.0 * i;
            tr.run<PlaneBeam3D, TracePoint>(200000000, true);
        }
        for(int i = 0; i <= 15; i++)
        {
            tr.surface().CorrLength() = 50.0 + 10.0 * i;
            tr.run<PlaneBeam3D, TracePoint>(200000000, true);
        }
    }
    catch (std::exception & e)
    {
        std::cout << e.what() << std::endl;
    }
    return 0;
}