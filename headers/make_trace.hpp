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

            exml_tree (const Point & pt) noexcept : exml_tree(std::vector<double> {pt.x(), pt.y(), pt.z()}) {}
            exml_tree (const Vector & v) noexcept : exml_tree(std::vector<double> {v.theta(), v.phi()}) {}
            exml_tree (const Line & line) noexcept;
            exml_tree (const Trace & trace) noexcept;
            exml_tree (const Sphere & sph) noexcept;
            const pt::ptree & ptree() const noexcept {return tree_; }
    };

    class xml_info
    {
        private:
            pt::ptree tree_;
        public:
            template <typename T>
            xml_info(const std::vector<T> & vec) noexcept(std::is_nothrow_copy_assignable<T>::value)
            {
                tree_.put("List", "");
                for(const auto & x : vec)
                {
                    pt::ptree child = xml_info(x).ptree();
                    for(const auto & v : child)
                        tree_.add_child("List." + v.first, child.get_child(v.first));
                }
            }

            xml_info(const Sphere & sph) noexcept;
            xml_info(const ExpSetup & setup) noexcept;
            const pt::ptree & ptree() const noexcept {return tree_; }
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
        public:
            enum mode {TRACE, DET};
        private:
            mode m_;
            std::vector<std::string_view> mode_name_ {"Trace", "Detector"};
            unsigned int threads_num_ {std::thread::hardware_concurrency()};
            unsigned int write_partition_ {500};
            struct run_pars
            {
                run_pars(unsigned long rays) : rays_max_(rays) {}
                mutable std::atomic<unsigned long> rays_ {0};
                unsigned long rays_max_;                
            };

            template <typename T, typename Func>
            void do_trace_block(Func func, const run_pars & pars_) const
            {
                using res_t = typename setup::function_traits<Func>::return_type;
                while (pars_.rays_ < pars_.rays_max_)
                {
                    std::vector<res_t> vec_trace_ {};
                    while(vec_trace_.size() < write_partition_ && ++pars_.rays_ < pars_.rays_max_)
                    {
                        auto trace_ = make_trace<T>(static_cast<ExpSetup>(*this));
                        if(trace_.is_transmitted())
                            vec_trace_.emplace_back(func(trace_));
                    }
                    if(vec_trace_.size() != 0)
                        base_setup::write_xml(exml_tree(vec_trace_), "traces");
                }
            }

            template <typename T>
            void log(const run_pars & pars) const
            {
                pt::ptree log;
                log.put("Title", "Ray tracing on spherical surfaces by Nikolay Ivanov");
                log.put("Date", base_setup::date());
                log.put("Number_of_rays", pars.rays_max_);
                log.put("Write-mode", mode_name_[m_]);
                log.put("Incident_angle", ExpSetup::IncAngle());
                log.put("Beam_mode", boost::core::demangle(typeid(T).name()));
                log.put("Experimental-setup", "");
                log.put_child("Experimental-setup", xml_info(static_cast<ExpSetup>(*this)).ptree());
                base_setup::write_xml(log, "info");
            }
        public:
            using ExpSetup::add_sphere;

            tracing(mode m, double inc_ang = 0.0, double RMSHeight = Constants::RMSHeight, double CorrLength = Constants::CorrLength, const fs::path & path = fs::current_path())
            : m_(m), base_setup(path), ExpSetup(inc_ang, RMSHeight, CorrLength) {}
            void set_partition(unsigned int n) noexcept {write_partition_ = n; }
            void set_threads(unsigned int n) noexcept {threads_num_ = n; }

            template <typename T>
            void run(unsigned long rays) const
            {
                run_pars pars (rays);
                auto trace_ = [](Trace & trace_){return std::move(trace_); };
                auto det_ = [this](const Trace & trace_){return trace_.det_res(static_cast<ExpSetup>(*this)); };
                log<T>(pars);
                std::vector<std::future<void>> futures;
                if(m_ == TRACE)
                    for(int i = 0; i < threads_num_; i++)
                        futures.emplace_back(std::async(&tracing::do_trace_block<T, decltype(trace_)>, this, trace_, std::ref(pars)));
                else
                    for(int i = 0; i < threads_num_; i++)
                        futures.emplace_back(std::async(&tracing::do_trace_block<T, decltype(det_)>, this, det_, std::ref(pars)));
                for(auto & e : futures)
                    e.get();
            }
    };

    class trace_parser : private tracing
    {
        private:
            po::variables_map vm_;
            tracing::mode mode_(const std::string & str) const;
        public:
            trace_parser(const po::variables_map & vm, const fs::path & path = fs::current_path());
            void run() const;
    };

}

#endif