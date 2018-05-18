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
            
            template <
                typename T,
                typename = std::enable_if_t<std::is_arithmetic<std::remove_reference_t<T>>::value>
            >
            exml_tree(T && n) : exml_tree()
            {
                tree_.put("Expression.Number", std::forward<T>(n));
            }

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
            exml_tree (const TracePoint & tr_pt) noexcept;
            exml_tree (const Sphere & sph) noexcept;
            const pt::ptree & ptree() const noexcept {return tree_; }
    };

    class info_tree
    {
        private:
            pt::ptree tree_;
        public:
            template <typename T>
            info_tree(const std::vector<T> & vec) noexcept(std::is_nothrow_copy_assignable<T>::value)
            {
                tree_.put("List", "");
                for(const auto & x : vec)
                {
                    pt::ptree child = info_tree(x).ptree();
                    for(const auto & v : child)
                        tree_.add_child("List." + v.first, child.get_child(v.first));
                }
            }

            info_tree(const Sphere & sph) noexcept;
            info_tree(const ExpSetup & setup) noexcept;
            const pt::ptree & ptree() const noexcept {return tree_; }
    };

    class base_setup
    {
        private:
            fs::path root_path_;
            mutable fs::path current_path_;
            fs::path make_path(const fs::path & path) const;
            fs::path write_path(const fs::path & filename, const fs::path & ext) const;
        public:
            base_setup(const fs::path & p = fs::current_path()) : root_path_(p), current_path_(p) {}
            void make_dir() const;
            void write_xml(const exml_tree & tree, const std::string & name) const;
            void write_xml(const pt::ptree & tree, const std::string & name) const;
            std::string date() const noexcept;
            const fs::path & current_path() const noexcept {return current_path_; }
            const fs::path & root_path() const noexcept {return root_path_; }
            fs::path & root_path() noexcept {return root_path_; }
    };
    
    class tracing : private base_setup, private ExpSetup
    {
        private:
            unsigned int threads_num_ {std::thread::hardware_concurrency()};
            struct run_pars
            {
                run_pars(unsigned long rays) : rays_max_(rays), write_partition_(100 * int(sqrt(rays / 1e6) + 1)) {}
                mutable std::atomic<unsigned long> rays_ {0};
                unsigned long rays_max_;
                unsigned int write_partition_;                
            };

            void add_result(std::vector<Trace> & vec, Trace & trace) const;
            void add_result(std::vector<TracePoint> & vec, const Trace & trace) const;
            void write_result(std::vector<Trace> & vec) const;
            void write_result(std::vector<TracePoint> & vec) const;

            template <typename BeamType, typename ResType>
            void do_trace_block(const run_pars & pars_) const
            {
                while (pars_.rays_ < pars_.rays_max_)
                {
                    std::vector<ResType> vec_ {};
                    while(vec_.size() < pars_.write_partition_ && ++pars_.rays_ < pars_.rays_max_)
                    {
                        auto trace_ = make_trace<BeamType>(static_cast<ExpSetup>(*this));
                        if(trace_.is_transmitted())
                            add_result(vec_, trace_);
                    }
                    if(vec_.size() != 0)
                        write_result(vec_);
                }
            }

            template <typename BeamType>
            void log(const run_pars & pars) const
            {
                pt::ptree log;
                log.put("Title", "Ray tracing on spherical surfaces by Nikolay Ivanov");
                log.put("Date", base_setup::date());
                log.put("Number-of-rays", pars.rays_max_);
                log.put("Beam-mode", boost::core::demangle(typeid(BeamType).name()));
                log.put("Experimental-setup", "");
                log.put_child("Experimental-setup", info_tree(static_cast<ExpSetup>(*this)).ptree());
                base_setup::write_xml(log, "info");
            }
            void log_wolfram() const;
        public:
            using ExpSetup::add_sphere;
            using ExpSetup::spheres;
            using ExpSetup::reset_sphere;
            using ExpSetup::substrate;
            using base_setup::root_path;
            using ExpSetup::IncAngle;
            using ExpSetup::surface;

            tracing(double inc_ang = 0.0, double RMSHeight = Constants::RMSHeight, double CorrLength = Constants::CorrLength, const fs::path & path = fs::current_path())
            : base_setup(path), ExpSetup(inc_ang, RMSHeight, CorrLength) {}
            void set_threads(unsigned int n) noexcept {threads_num_ = n; }

            template <typename BeamType, typename ResType>
            void run(unsigned long rays, bool verbose = false) const
            {
                auto t_ = std::chrono::high_resolution_clock::now();
                run_pars pars (rays);
                base_setup::make_dir();
                if(verbose)
                    std::clog << "Running ray tracing of " << rays << " rays\n"
                            << "Writing results to " << base_setup::current_path() << std::endl;
                log<BeamType>(pars);
                if(std::is_same<Trace, ResType>::value)
                    log_wolfram();
                std::vector<std::future<void>> futures;
                for(int i = 0; i < threads_num_; i++)
                    futures.emplace_back(std::async(&tracing::do_trace_block<BeamType, ResType>, this, std::ref(pars)));
                for(auto & e : futures)
                    e.get();
                if(verbose)
                    std::clog << "Ray tracing has finished\nElapsed time: "
                            << std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::high_resolution_clock::now() - t_).count() << " s\n";
            }
    };

    class trace_parser : private tracing
    {
        private:
            po::variables_map vm_;
        public:
            trace_parser(const po::variables_map & vm, const fs::path & path = fs::current_path());
            void run() const;
    };

}

#endif