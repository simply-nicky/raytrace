#include <headers/trace.hpp>
#include <headers/make_trace.hpp>

using namespace raytrace;

int main()
{
    try
    {       
        fs::path dir = fs::absolute("RMSHeight_Series_3D_R2000");
        fs::create_directory(dir);
        tracing tr;
        tr.root_path() = dir;
        tr.substrate() = Sphere(0.0, 0.0, 2000.0, 2000.0, 60, Sphere::LOWER);

        for(int i = 0; i < 5; i++)
        {
            tr.surface().RMSHeight() = 0.3 * i;
            tr.run<PlaneBeam3D, TracePoint>(200000000, true);
        }


        // tr.substrate() = Sphere(0.0, 0.0, 750.0, 750.0, 60.0, Sphere::LOWER);
        // tr.run<PlaneBeam3D, TracePoint>(200000000, true);
        // double d = 0.00310473;
        // double h = 0.25 * d;
        // while(h <= 8 * d)
        // {
        //     tr.add_sphere(0, 0, h, 10 * d);
        //     tr.run<PlaneBeam3D, TracePoint>(400000000, true);
        //     tr.reset_sphere();
        //     if(h < 2 * d)
        //         h += 0.25 * d;
        //     else if(h < 6 * d)
        //         h += 0.5 * d;
        //     else
        //         h += d;
        // }
    }
    catch (std::exception & e)
    {
        std::cout << e.what() << std::endl;
    }
    return 0;
}