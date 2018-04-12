#include <headers/make_trace.hpp>

namespace raytrace{


// exml_tree
exml_tree::exml_tree() noexcept
{
    tree_.put("Expression.<xmlattr>.xmlns:mathematica", "http://www.wolfram.com/XML/");
    tree_.put("Expression.<xmlattr>.xmlns", "http://www.wolfram.com/XML/");
}

exml_tree::exml_tree(double x) noexcept : exml_tree()
{
    tree_.put("Expression.Number", x);
    tree_.put("Expression.Number.<xmlattr>.Precision", "");
}

exml_tree::exml_tree(int n) noexcept : exml_tree()
{
    tree_.put("Expression.Number", n);
}

exml_tree::exml_tree (const Line & line) noexcept : exml_tree(
    std::vector<double> 
    {
        line.point().x(), 
        line.point().y(), 
        line.point().z(), 
        line.direction().theta(), 
        line.direction().phi()
    }) {}

exml_tree::exml_tree(const Trace & trace) noexcept : exml_tree()
{
    tree_.put("Expression.Function.Symbol", "List");
    tree_.put("Expression.Function.Function.Symbol", "List");
    for(const auto & pt : trace)
        tree_.add_child("Expression.Function.Function.Function", exml_tree(pt).ptree().get_child("Expression.Function"));
    tree_.add_child("Expression.Function.Number", exml_tree(trace.is_transmitted()).ptree().get_child("Expression.Number"));
}


// xm_info
xml_info::xml_info(const Sphere & sph) noexcept
{
    tree_.put("Center[mm]", sph.center());
    tree_.put("Radius[mm]", sph.radius());
    tree_.put("Diameter[mm]", sph.diameter());
}

xml_info::xml_info(const ExpSetup & setup) noexcept
{
    tree_.put("Material", "");
    tree_.put("Material.Root-mean-square-height[nm]", setup.RMSHeight());
    tree_.put("Material.Correlation-length[um]", setup.CorrLength());
    tree_.put("Wavelength[nm]", setup.wavelength());
    tree_.put("Incident-angle", setup.IncAngle());
    tree_.put("Imperfections-mode", setup.spheres().size() > 1);
    tree_.put_child("Spheres", xml_info(setup.spheres()).ptree());
    tree_.put("Experiment-geometry", "");
    tree_.put("Experiment-geometry.source-to-mirror-distance[mm]", setup.source_distance());
    tree_.put("Experiment-geometry.mirror-to-detector-distance[mm]", setup.detector_distance());
    tree_.put("Experiment-geometry.source-length[mm]", setup.source_length());
}

// base_setup
std::string base_setup::set_date(const time_t t) const noexcept
{
    char buffer [80];
    std::strftime(buffer, sizeof(buffer), "%d %b %Y %I", std::localtime(&t));
    return std::string(buffer);    
}

fs::path base_setup::make_path(const fs::path & path) const
{
    int i = 0;
    auto temp = fs::absolute(path, p_);
    while(fs::exists(temp))
        temp = fs::absolute(path, p_).string() + "_" + std::to_string(i++);
    return temp;
}

fs::path base_setup::make_path(const fs::path & path, const fs::path & ext) const
{
    int i = 0;
    auto temp = fs::absolute(path, p_).replace_extension(ext);
    while(fs::exists(temp))
    {
        temp = fs::absolute(path, p_).string() + "_" + std::to_string(i++);
        temp.replace_extension(ext);
    }
    return temp;
}

base_setup::base_setup(const fs::path & p, const time_t t) : date_(set_date(t)), p_(p)
{
    p_ = make_path(date_);
    fs::create_directory(p_);
}

void base_setup::write_xml(const exml_tree & tree, const std::string & name) const
{
    pt::xml_parser::write_xml(make_path("traces", "xml").string(), tree.ptree(), std::locale(), pt::xml_writer_settings<std::string>(' ', 4));
}

void base_setup::write_xml(const pt::ptree & tree, const std::string & name) const
{
    pt::xml_parser::write_xml(make_path(name, "xml").string(), tree, std::locale(), pt::xml_writer_settings<std::string>(' ', 4));
}


// trace_parser
tracing::mode trace_parser::mode_(const std::string & str) const
{
    if(str == "Trace")
        return tracing::TRACE;
    else if(str == "Detector")
        return tracing::DET;
    else
        throw std::invalid_argument("trace_parser::mode : invalid argument");
}

trace_parser::trace_parser(const po::variables_map & vm, const fs::path & path)
: vm_(vm), tracing(mode_(vm_["write-mode"].as<std::string>()), vm_["IncAngle"].as<float>(), vm_["RMSHeight"].as<float>(), vm_["CorrLength"].as<float>(), path)
{
    tracing::add_sphere(vm_["Sphere.location.x"].as<std::vector<double>>(), vm_["Sphere.location.y"].as<std::vector<double>>(), vm_["Sphere.height"].as<std::vector<double>>(), vm_["Sphere.diameter"].as<std::vector<double>>());
}

void trace_parser::run() const
{
    if(vm_["BeamMode"].as<std::string>() == "Plane" && vm_["Dimensions"].as<std::string>() == "2D")
        tracing::run<PlaneBeam2D>(vm_["rays"].as<int>());
    else if(vm_["BeamMode"].as<std::string>() == "Spherical" && vm_["Dimensions"].as<std::string>() == "2D")
        tracing::run<SphBeam2D>(vm_["rays"].as<int>());       
    else if(vm_["BeamMode"].as<std::string>() == "Plane" && vm_["Dimensions"].as<std::string>() == "3D")
        tracing::run<PlaneBeam3D>(vm_["rays"].as<int>());       
    else if(vm_["BeamMode"].as<std::string>() == "Spherical" && vm_["Dimensions"].as<std::string>() == "3D")
        tracing::run<SphBeam3D>(vm_["rays"].as<int>());  
    else
        throw std::invalid_argument("trace_parser::run : invalid argument.");     
}

}