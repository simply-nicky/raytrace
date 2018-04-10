#include <headers/trace.hpp>
#include <headers/make_trace.hpp>
#include <headers/include.hpp>

int main (int ac, char * av[])
{
    using namespace raytrace;
    try
    {
        auto check = [](const std::string & str)
        {
            return [&str](const std::vector<double> & vec)
            {
                for(const auto & x : vec)
                    if(x < 0.0)
                        throw std::invalid_argument(str + std::string(" must be a positive number."));
            };
        };
        po::options_description generic("Generic options");
        generic.add_options()
            ("help,h", "This usage infromation.")
            ("help-config-file", "Print config file manual")
            ("rays,N", po::value<int>()->value_name(" ")->notifier([](const int n){if(n < 0) throw std::invalid_argument("The number of rays must be a positive number."); }), "Number of rays to trace.")
            ("verbose,v", po::value<int>()->value_name(" ")->default_value(1), "Verbosity of output.")
            ("config-filename", po::value<std::string>()->default_value(std::string("config.ini")), "Config filename.")
            ;
        po::options_description config("Configuration");
        config.add_options()
            ("BeamMode", po::value<std::string>()->value_name(" ")->default_value("Plane")->notifier([](const std::string & str){if(str != "Spherical" && str != "Plane") throw std::invalid_argument("Beam mode must be Spherical or Plane."); }), "Incident beam mode: Spherical or Plane.")
            ("Dimensions", po::value<std::string>()->value_name(" ")->default_value("3D")->notifier([](const std::string & str){if(str != "2D" && str != "3D") throw std::invalid_argument("Dimensions must be 2D or 3D."); }), "3D or 2D Ray-tracing.")
            ("IncAngle", po::value<float>()->value_name(" ")->default_value(0), "Incident beam angle.")
            ("RMSHeight",po::value<float>()->value_name(" ")->default_value(Constants::RMSHeight, "0.35 nm")->notifier([](const double val){if(val < 0.0) throw std::invalid_argument("RMSHeight must be a positive number."); }),"Root-mean-square height of surface roughness.")
            ("CorrLength",po::value<float>()->value_name(" ")->default_value(Constants::CorrLength, "2.7 um")->notifier([](const double val){if(val < 0.0) throw std::invalid_argument("Correlation length must be a positive number."); }),"Correlation length of surface roughness.")
            ("ImpMode", po::bool_switch()->default_value(false), "Ray tracing with imperfections, config file is required.")
            ;
        po::options_description config_file("Configuration files");
        config_file.add_options()
            ("Sphere.height", po::value<std::vector<double>>()->notifier(check("height")), "Imperfection sphere height.")
            ("Sphere.diameter", po::value<std::vector<double>>()->notifier(check("diameter")), "Imperfection sphere diameter.")
            ("Sphere.location.x", po::value<std::vector<double>>(), "Imperfection sphere location: x coordinate")
            ("Sphere.location.y", po::value<std::vector<double>>(), "Imperfection sphere location: y coordinate")
            ;
        po::options_description all("Ray tracing on spherical surfaces.\nAllowed options");
        all.add(generic).add(config);      
        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, all), vm);
        if(vm.count("ImpMode"))
        {
            std::ifstream config_stream {vm["config-filename"].as<std::string>().c_str()};
            if(config_stream)
                po::store(po::parse_config_file(config_stream, config_file), vm);
        }
        po::notify(vm);

        if(ac == 1 || vm.count("help"))
            std::cout << all;
        else if(vm.count("help-config-file"))
            std::cout << config_file;
        else
        {
            trace_parser tr (vm);
            tr.run();
        }
    }
    catch (const std::exception & e)
    {
        std::cerr << e.what() << std::endl;
    }
}