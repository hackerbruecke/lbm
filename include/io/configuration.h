#pragma once

#include <string>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <omp.h>

namespace lbm
{
namespace io
{
namespace po = boost::program_options;

class Config
{
    // TODO: Change to private and grant friend access to mpi functions
public:
    std::string _collision_model;
    std::string _input_file;
    std::string _output_dir { "vtk" };
    std::string _output_filename { "output" };
    std::string _scenario_xml;
    uint32_t _timesteps { 0 };
    uint32_t _timesteps_per_plot { 0 };
    uint32_t _omp_threads { 1 };
    uint32_t _iproc { 1 };
    uint32_t _jproc { 1 };
    uint32_t _kproc { 1 };
    double _tau { 1.0 };

public:
    auto iproc() const -> decltype(_iproc)
    {
    	return _iproc;
    }
    auto jproc() const -> decltype(_jproc)
	{
		return _jproc;
	}
    auto kproc() const -> decltype(_kproc)
	{
		return _kproc;
	}
    auto collision_model() const -> decltype(_collision_model)
    {
        return _collision_model;
    }
    auto input_file() const -> decltype(_input_file)
    {
        return _input_file;
    }
    auto output_dir() const -> decltype(_output_dir)
    {
        return _output_dir;
    }
    auto output_filename() const -> decltype(_output_filename)
    {
        return _output_filename;
    }
    void set_output_filename(const std::string& filename)
    {
        _output_filename = filename;
    }
    auto timesteps() const -> decltype(_timesteps)
    {
        return _timesteps;
    }
    auto timesteps_per_plot() const -> decltype(_timesteps_per_plot)
    {
        return _timesteps_per_plot;
    }
    auto tau() const -> decltype(_tau)
    {
        return _tau;
    }
    auto omp_threads() const -> decltype(_omp_threads)
    {
        return _omp_threads;
    }
    auto scenario_xml() const -> decltype(_scenario_xml)
    {
        return _scenario_xml;
    }

    Config() {}
    Config(int argc, char** argv)
    {
        // General options
        po::options_description general("General options");
        general.add_options()
                ("help,h", "Show help message")
                ("input-file,i", po::value<std::string>(&_input_file),
                        "Input file containing configuration values");

        // Configuration options
        po::options_description config("Configuration options");
        config.add_options()
        ("collision-model,c", po::value<std::string>(&_collision_model)->default_value("bgk"),
                "Collision operator model")
        ("tau", po::value<double>(&_tau)->required(),
                "Relaxation factor for BGK collision operator. Must be in (0.5, 2.0)")
        ("timesteps,t", po::value<std::uint32_t>(&_timesteps)->required(),
                "Number of steps to perform")
        ("timesteps-per-plot", po::value<std::uint32_t>(&_timesteps_per_plot)->required(),
                "Number of timesteps after which an output file is written")
        ("scenario-file", po::value<std::string>(&_scenario_xml)->required(),
                "XML file containing scenario to simulate")
        ("omp-threads", po::value<uint32_t>(&_omp_threads),
                "Number of OpenMP threads to use")
        ("output-dir", po::value<std::string>(&_output_dir),
                "Output directory for plots")
        ("iproc", po::value<uint32_t>(&_iproc),
                "Number of processes in x-direction (MPI)")
        ("jproc", po::value<uint32_t>(&_jproc),
                "Number of processes in y-direction (MPI)")
        ("kproc", po::value<uint32_t>(&_kproc),
                "Number of processes in z-direction (MPI)")
        ;// End of options
        // Store command line options of general and config options
        po::options_description args;
        args.add(general).add(config);

        // Positional options: One Filename
        po::positional_options_description positional;
        positional.add("input-file", 1);

        // Store options in map
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(args)
                .allow_unregistered().positional(positional).run(), vm);

        std::cout << "LBM simulation by Krivokapic, Mody, Malcher" << std::endl;
        if (vm.count("help") || argc == 1) {
            std::cout << "Usage: lbm_isotropy [file] additional-options..." << std::endl;
            std::cout << general << std::endl;
            std::cout << config << std::endl;
            std::exit(1);
        }
        if (vm.count("input-file")) {
            std::cout << "Reading configuration file..." << std::endl;
            _input_file = vm["input-file"].as<std::string>();
            po::store(po::parse_config_file<char>(_input_file.c_str(), config), vm);
            po::notify(vm);
        }
        po::notify(vm);

        assert(_timesteps > 0);
        assert(_tau > 0.5 && _tau < 2.0);
        assert(_omp_threads > 0);
        // TODO: Currently only bgk is supported. Remove this assert to enable support for
        // further collision operators (e.g. mrt)
        assert(_collision_model == "bgk");
        omp_set_num_threads(_omp_threads);

        namespace fs = boost::filesystem;
        fs::path dir(_output_dir);
        if (!fs::exists(dir)) {
            // Create directory
            std::cout << "Output directory \"" + _output_dir + "\" does not exist. Creating."
                    << std::endl;
            fs::create_directory(dir);
        } else {
            // Clean existing directory
            fs::directory_iterator end;
            for (fs::directory_iterator file(dir); file != end; ++file) {
                fs::remove(file->path());
            }
        }
    }

//    friend void lbm::mpi::broadcast(lbm::io::Config& cfg, std::array<size_t, 3>& domain_length, MPI_Comm comm);
};

auto operator <<(std::ostream& lhs, const lbm::io::Config& cfg) -> decltype(lhs)
{
    lhs << "Configuration:" << '\n'
        << "> Configuration file:     " << cfg.input_file() << '\n'
        << "> Output directory:       " << cfg.output_dir() << '\n'
        << "> Output filename:        " << cfg.output_filename() << '\n'
        << "> Collision model:        " << cfg.collision_model() << '\n'
        << "> Tau:                    " << cfg.tau() << '\n'
        << "> Timesteps:              " << cfg.timesteps() << '\n'
        << "> Timesteps per plot:     " << cfg.timesteps_per_plot() << '\n'
        << "> Scenario file:          " << cfg.scenario_xml() << '\n'
        << "> === Parallel settings ===\n"
        << "> Number of OMP threads:  " << omp_get_max_threads() << '\n'
        << "> iproc:                  " << cfg.iproc() << '\n'
        << "> jproc:                  " << cfg.jproc() << '\n'
        << "> kproc:                  " << cfg.kproc() << '\n'
        ;
    return lhs;
}

} //namespace io
} //namespace lbm
