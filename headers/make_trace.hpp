#ifndef MAKE_TRACE_
#define MAKE_TRACE_
#include <headers/trace.hpp>
#include <boost/core/demangle.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>


namespace raytrace {

    namespace fs = boost::filesystem;
    namespace pt = boost::property_tree;
    namespace po = boost::program_options;

    class exml_tree
    {
        private:
            pt::ptree tree_;
        public:
            exml_tree () noexcept;
            exml_tree (int) noexcept;
            exml_tree (double) noexcept;
            
            template <typename T>
            exml_tree (const std::vector<T> & vec) noexcept(std::is_nothrow_copy_assignable<T>::value) : exml_tree() 
            {
                tree_.put("Expression.Function.Symbol", "List");
                for (const auto & ptr : vec)
                {
                    pt::ptree child = exml_tree(ptr).ptree().get_child("Expression");
                    for(const auto & v : child)
                        if (v.first != "<xmlattr>")
                            tree_.add_child("Expression.Function." + v.first, child.get_child(v.first));
                }
            }

            exml_tree (const Point & pt) : exml_tree(std::vector<double> {pt.x(), pt.y(), pt.z()}) {}
            exml_tree (const Vector & v) : exml_tree(std::vector<double> {v.theta(), v.phi()}) {}
            exml_tree (const Line & line);
            exml_tree (const Trace & trace);
            const pt::ptree& ptree() const noexcept {return tree_; }
    };

    class base_setup
    {
        private:
            std::string date_;
            fs::path p_;
            std::string set_date(const time_t t) const noexcept;
            fs::path make_path(const fs::path & path) const;
            fs::path make_path(const fs::path & path, const fs::path & ext) const;
        public:
            base_setup(const fs::path & p = fs::current_path(), const time_t t = std::time(nullptr));
            void write_xml(const exml_tree & tree, const std::string & name) const;
            void write_xml(const pt::ptree & tree, const std::string & name) const;
            const std::string & date() const noexcept {return date_; }
    };
    
    class tracing : private base_setup, private ExpSetup
    {
        private:
            mutable std::atomic<unsigned long> rays_ {0};
            unsigned int threads_num_ {std::thread::hardware_concurrency()};
            unsigned int write_partition_ {100};
            unsigned long rays_max_;

            template <typename T>
            void do_trace_block() const
            {
                while (rays_ < rays_max_)
                {
                    std::vector<Trace> vec_trace_ {};
                    while(vec_trace_.size() < write_partition_ && ++rays_ < rays_max_)
                    {
                        auto trace_ = make_trace<T>(static_cast<ExpSetup>(*this));
                        if(trace_.is_transmitted())
                            vec_trace_.emplace_back(std::move(trace_));
                    }
                    if(vec_trace_.size() != 0)
                        base_setup::write_xml(exml_tree(vec_trace_), "traces");
                }
            }

            template <typename T>
            void log() const
            {
                pt::ptree log;
                log.put("Title", "Ray tracing on spherical surfaces by Nikolay Ivanov");
                log.put("Date", base_setup::date());
                log.put("Number_of_rays", rays_max_);
                log.put("Incident_angle", ExpSetup::IncAngle());
                log.put("Beam_mode", boost::core::demangle(typeid(T).name()));
                log.put("Root-mean-square_height", ExpSetup::RMSHeight());
                log.put("Correlation_length", ExpSetup::CorrLength());
                base_setup::write_xml(log, "info");
            }
        public:
            tracing(unsigned long rays_max, double inc_ang = 0.0, double RMSHeight = Constants::RMSHeight, double CorrLength = Constants::CorrLength, const fs::path & path = fs::current_path())
            : rays_max_{rays_max}, base_setup(path), ExpSetup(inc_ang, RMSHeight, CorrLength) {}
            void set_partition(unsigned int n) noexcept {write_partition_ = n; }
            void set_threads(unsigned int n) noexcept {threads_num_ = n; }
            void add_sphere(double x, double y, double height, double diameter) noexcept;
            tracing(const tracing & tr) : rays_(tr.rays_.load()), threads_num_(tr.threads_num_), rays_max_(tr.rays_max_), write_partition_(tr.write_partition_) {} 

            template <typename T>
            void run() const
            {
                log<T>();
                std::vector<std::future<void>> futures;
                for(int i = 0; i < threads_num_; i++)
                    futures.emplace_back(std::async(&tracing::do_trace_block<T>, this));
                for(auto & e : futures)
                    e.get();
                std::cout << rays_ << std::endl;
            }
    };

    class trace_parser : private tracing
    {
        private:
            po::variables_map vm_;
        public:
            trace_parser(const po::variables_map & vm, const fs::path & path = fs::current_path())
            : vm_(vm), tracing(vm_["rays"].as<int>(), vm_["IncAngle"].as<float>(), vm_["RMSHeight"].as<float>(), vm_["CorrLength"].as<float>(), path) {}
            void run() const;
    };

}

#endif