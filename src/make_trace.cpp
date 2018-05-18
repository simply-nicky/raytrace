#include <headers/make_trace.hpp>

namespace raytrace{


// exml_tree
exml_tree::exml_tree() noexcept
{
    tree_.put("Expression.<xmlattr>.xmlns:mathematica", "http://www.wolfram.com/XML/");
    tree_.put("Expression.<xmlattr>.xmlns", "http://www.wolfram.com/XML/");
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
    tree_.add("Expression.Function.Number", int(trace.is_transmitted()));
    tree_.add("Expression.Function.Number", int(trace.is_scattered()));
}

exml_tree::exml_tree(const TracePoint & tr_pt) noexcept : exml_tree()
{
    tree_.put("Expression.Function.Symbol", "List");
    tree_.put_child("Expression.Function.Function", exml_tree(tr_pt.det_point).ptree().get_child("Expression.Function"));
    tree_.add("Expression.Function.Number", int(tr_pt.is_transmitted));
    tree_.add("Expression.Function.Number", int(tr_pt.is_scattered));
}

exml_tree::exml_tree(const Sphere & sph) noexcept : exml_tree()
{
    tree_.put("Expression.Function.Symbol", "List");
    tree_.put_child("Expression.Function.Function", exml_tree(sph.center()).ptree().get_child("Expression.Function"));
    tree_.add("Expression.Function.Number", sph.radius());
    tree_.add("Expression.Function.Number", sph.diameter());
}


// info_tree
info_tree::info_tree(const Sphere & sph) noexcept
{
    tree_.put("Center[mm]", sph.center());
    tree_.put("Radius[mm]", sph.radius());
    tree_.put("Diameter[mm]", sph.diameter());
}

info_tree::info_tree(const ExpSetup & setup) noexcept
{
    tree_.put("Material", "");
    tree_.put("Material.Root-mean-square-height[nm]", setup.surface().RMSHeight());
    tree_.put("Material.Correlation-length[um]", setup.surface().CorrLength());
    tree_.put("Wavelength[nm]", setup.wavelength());
    tree_.put("Incident-angle", setup.IncAngle());
    tree_.put("Imperfections-mode", setup.spheres().size() > 1);
    tree_.put_child("Spheres", info_tree(setup.spheres()).ptree());
    tree_.put("Experiment-geometry", "");
    tree_.put("Experiment-geometry.source-to-mirror-distance[mm]", setup.source_distance());
    tree_.put("Experiment-geometry.mirror-to-detector-distance[mm]", setup.detector_distance());
    tree_.put("Experiment-geometry.source-length[mm]", setup.source_length());
}


// base_setup
std::string base_setup::date() const noexcept
{
    auto t = std::time(nullptr);
    char buffer [80];
    std::strftime(buffer, sizeof(buffer), "%d %b %Y %H %M", std::localtime(&t));
    return std::string(buffer);    
}

fs::path base_setup::make_path(const fs::path & path) const
{
    int i = 0;
    auto temp = fs::absolute(path, root_path_);
    while(fs::exists(temp))
    {
        temp = fs::absolute(path, root_path_).string() + "_" + std::to_string(i / 10) + std::to_string(i % 10);
        i++;
    }
    return temp;
}

fs::path base_setup::write_path(const fs::path & filename, const fs::path & ext) const
{
    int i = 0;
    auto temp = fs::absolute(filename, current_path_).replace_extension(ext);
    while(fs::exists(temp))
    {
        temp = fs::absolute(filename, current_path_).string() + "_" + std::to_string(i / 10) + std::to_string(i % 10);
        temp.replace_extension(ext);
        i++;
    }
    return temp;
}

void base_setup::make_dir() const
{
    current_path_ = make_path(date());
    fs::create_directory(current_path_);
}

void base_setup::write_xml(const exml_tree & tree, const std::string & name) const
{
    pt::xml_parser::write_xml(write_path(name, "xml").string(), tree.ptree(), std::locale(), pt::xml_writer_settings<std::string>(' ', 4));
}

void base_setup::write_xml(const pt::ptree & tree, const std::string & name) const
{
    pt::xml_parser::write_xml(write_path(name, "xml").string(), tree, std::locale(), pt::xml_writer_settings<std::string>(' ', 4));
}


// tracing
void tracing::log_wolfram() const
{
    base_setup::write_xml(exml_tree(ExpSetup::spheres()), "log_wolfram");
}

void tracing::add_result(std::vector<Trace> & vec, Trace & trace) const
{
    vec.emplace_back(std::move(trace));
}

void tracing::add_result(std::vector<TracePoint> & vec, const Trace & trace) const
{
    vec.emplace_back(trace.det_res(ExpSetup::detector_distance()));
}

void tracing::write_result(std::vector<Trace> & vec) const
{
    base_setup::write_xml(exml_tree(vec), "traces");
}

void tracing::write_result(std::vector<TracePoint> & vec) const
{
    base_setup::write_xml(exml_tree(vec), "det_points");
}


// trace_parser
trace_parser::trace_parser(const po::variables_map & vm, const fs::path & path)
: vm_(vm), tracing(vm["inc-angle"].as<float>(), vm["rms-height"].as<float>(), vm["corr-length"].as<float>(), path)
{
    if(vm_.count("imp-mode"))
        tracing::add_sphere(vm_["Sphere.location.x"].as<std::vector<double>>(), vm_["Sphere.location.y"].as<std::vector<double>>(), vm_["Sphere.height"].as<std::vector<double>>(), vm_["Sphere.radius"].as<std::vector<double>>());
}

void trace_parser::run() const
{
    if(vm_["beam-mode"].as<std::string>() == "Plane" && vm_["dimensions"].as<std::string>() == "2D" && vm_["write-mode"].as<std::string>() == "Trace")
        tracing::run<PlaneBeam2D, Trace>(vm_["rays"].as<int>(), vm_["verbose"].as<int>());
    else if(vm_["beam-mode"].as<std::string>() == "Spherical" && vm_["dimensions"].as<std::string>() == "2D" && vm_["write-mode"].as<std::string>() == "Trace")
        tracing::run<SphBeam2D, Trace>(vm_["rays"].as<int>(), vm_["verbose"].as<int>());       
    else if(vm_["beam-mode"].as<std::string>() == "Plane" && vm_["dimensions"].as<std::string>() == "3D" && vm_["write-mode"].as<std::string>() == "Trace")
        tracing::run<PlaneBeam3D, Trace>(vm_["rays"].as<int>(), vm_["verbose"].as<int>());       
    else if(vm_["beam-mode"].as<std::string>() == "Spherical" && vm_["dimensions"].as<std::string>() == "3D" && vm_["write-mode"].as<std::string>() == "Trace")
        tracing::run<SphBeam3D, Trace>(vm_["rays"].as<int>(), vm_["verbose"].as<int>());
    else if(vm_["beam-mode"].as<std::string>() == "Plane" && vm_["dimensions"].as<std::string>() == "2D" && vm_["write-mode"].as<std::string>() == "Detector")
        tracing::run<PlaneBeam2D, TracePoint>(vm_["rays"].as<int>(), vm_["verbose"].as<int>());
    else if(vm_["beam-mode"].as<std::string>() == "Spherical" && vm_["dimensions"].as<std::string>() == "2D" && vm_["write-mode"].as<std::string>() == "Detector")
        tracing::run<SphBeam2D, TracePoint>(vm_["rays"].as<int>(), vm_["verbose"].as<int>());       
    else if(vm_["beam-mode"].as<std::string>() == "Plane" && vm_["dimensions"].as<std::string>() == "3D" && vm_["write-mode"].as<std::string>() == "Detector")
        tracing::run<PlaneBeam3D, TracePoint>(vm_["rays"].as<int>(), vm_["verbose"].as<int>());       
    else if(vm_["beam-mode"].as<std::string>() == "Spherical" && vm_["dimensions"].as<std::string>() == "3D" && vm_["write-mode"].as<std::string>() == "Detector")
        tracing::run<SphBeam3D, TracePoint>(vm_["rays"].as<int>(), vm_["verbose"].as<int>());  
    else
        throw std::invalid_argument("trace_parser::run : invalid argument.");     
}

}