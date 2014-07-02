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
    std::string _lattice_model;
    std::string _collision_model;
    std::string _input_file;
    std::string _input_vtk;
    std::string _output_dir { "output" };
    std::string _output_filename { "output" };
    uint64_t _timesteps { 0 };
    uint64_t _timesteps_per_plot { 0 };
    double _tau { 1.0 };
    uint32_t _omp_threads { 1 };

public:
    auto lattice_model() const -> decltype(_lattice_model)
    {
        return _lattice_model;
    }
    auto collision_model() const -> decltype(_collision_model)
    {
        return _collision_model;
    }
    auto input_file() const -> decltype(_input_file)
    {
        return _input_file;
    }
    auto input_vtk() const -> decltype(_input_vtk)
    {
        return _input_vtk;
    }
    auto output_dir() const -> decltype(_output_dir)
    {
        return _output_dir;
    }
    auto output_filename() const -> decltype(_output_filename)
    {
        return _output_filename;
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

    Config(int argc, char** argv)
    {
        // General options
        po::options_description general("General options");
        general.add_options()("help,h", "Show help message")("input-file,i",
                po::value<std::string>(&_input_file), "Input file containing configuration values");

        // Configuration options
        po::options_description config("Configuration options");
        config.add_options()("lattice-model,l",
                po::value<std::string>(&_lattice_model)->default_value("d3q19"),
                "Lattice model, one of d3q15, d3q19, dq327")("collision-model,c",
                po::value<std::string>(&_collision_model)->default_value("bgk"),
                "Collision operator model")("tau", po::value<double>(&_tau),
                "Relaxation factor for BGK collision operator")("timesteps,t",
                po::value<std::uint64_t>(&_timesteps), "Number of steps to perform")(
                "timesteps-per-plotting", po::value<std::uint64_t>(&_timesteps_per_plot),
                "Number of timesteps after which an output file is written")("input-vtk",
                po::value<std::string>(&_input_vtk),
                "Structured points VTK file containing fluid cell mask")("omp-threads",
                po::value<uint32_t>(&_omp_threads), "Number of OpenMP threads to use")("output-dir",
                po::value<std::string>(&_output_dir), "Output directory for plots")(
                "output-filename", po::value<std::string>(&_output_filename),
                "Base name of output file. Will be amended with timestep number and file type");

        // Store command line options of general and config options
        po::options_description args;
        args.add(general).add(config);

        // Positional options: One Filename
        po::positional_options_description positional;
        positional.add("input-file", 1);

        // Store options in map
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(args).positional(positional).run(),
                vm);
        po::notify(vm);

        if (vm.count("help") || argc == 1) {
            std::cout << general << std::endl;
            std::cout << config << std::endl;
            std::exit(1);
        }
        if (vm.count("input-file")) {
            po::store(po::parse_config_file<char>(_input_file.c_str(), config), vm);
        }
        po::notify(vm);

        assert(_tau > 0.5 && _tau < 2.0);
        assert(_omp_threads > 0);
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
};

std::ostream& operator <<(std::ostream& lhs, const lbm::io::Config& cfg)
{
    lhs << "=== LBM configuration ===" << '\n'
        << "Input file:            " << cfg.input_file() << '\n'
        << "Output directory:      " << cfg.output_dir() << '\n'
        << "Number of OMP threads: " << omp_get_max_threads() << '\n'
        << "Lattice model:         " << cfg.lattice_model() << '\n'
        << "Collision model:       " << cfg.collision_model() << '\n'
        << "Tau:                   " << cfg.tau() << '\n'
        << "Timesteps:             " << cfg.timesteps() << '\n'
        << "Timesteps per plot:    " << cfg.timesteps_per_plot() << '\n'
        << "Input VTK file:        " << cfg.input_vtk();
    return lhs;
}

} //namespace io
} //namespace lbm
