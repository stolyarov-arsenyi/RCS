#ifndef  RCS_SOLVER_H
#define  RCS_SOLVER_H


#include "Solver/Logger.h"
#include "Solver/Input.h"
#include "Solver/Solution.h"
#include "Solver/Integrator.h"


template <class Scalar, class Length>

struct Solver
{
    static constexpr Scalar size_to_memory = sizeof (Co <Scalar>) / (1024.0 * 1024.0 * 1024.0);

    std::string conf_name;
    std::string mesh_name;

    Length block_size;

    Mesh <Scalar> mesh;

    Input <Scalar> input;

    BlockSystem <Co <Scalar>, Length> system;

    Logger logger;


    Solver (const std::string & conf_name, const std::string & mesh_name) : conf_name(conf_name), mesh_name(mesh_name)
    {
        logger = Logger(mesh_name + "." + conf_name + "." + "log");
    }

    auto date () -> std::string
    {
        char str [32];

        auto time = std::time(nullptr);

        strftime(str, 32, "%a %b %d %H:%M:%S %Y", localtime(& time));

        return { str };
    }

    void init ()
    {
        logger << date() << " - Parsing input" << std::endl;

        input = Input <Scalar> :: parse(conf_name);

        logger << date() << " - Parsing mesh" << std::endl;

        mesh = Mesh <Scalar> :: load_from_stl(mesh_name);


        auto threads = input["threads"][0];

        omp_set_num_threads((int) threads);


        auto memory = input["block_size"][0];

        block_size = (Length) Sqrt(memory / size_to_memory);


        logger << " - Mesh size: " << mesh.size() << " edges" << std::endl;

        logger << " - Matrix size: " << Pow(mesh.size(), 2.0) * size_to_memory << " Gb" << std::endl;
    }

    void solve ()
    {
        std::stringstream solution_name;

        solution_name << mesh_name << "." << conf_name << ".csv";

        std::ofstream solution_file(solution_name.str());

        solution_file << "wavelength" << ","
                      << "pol_inc"    << ","
                      << "pol_sca"    << ","
                      << "azi_inc"    << ","
                      << "azi_sca"    << ","
                      << "alt_inc"    << ","
                      << "alt_sca"    << ","
                      << "sm"         << ","
                      << "db"         << std::endl;


        for (auto wavelength : input["wavelength"])
        {
            logger << " - Wavelength: " << wavelength << "m" << std::endl;

            auto k = 2.0 * M_PI / wavelength;


            std::vector <Solution <Scalar>> solutions;

            for (auto pol_inc: input["pol_inc"])
            for (auto pol_sca: input["pol_sca"])
            for (auto azi_inc: input["azi_inc"])
            for (auto azi_sca: input["azi_sca"])
            for (auto alt_inc: input["alt_inc"])
            for (auto alt_sca: input["alt_sca"])
            {
                Antenna <Scalar> antenna_inc(alt_inc, azi_inc, pol_inc);
                Antenna <Scalar> antenna_sca(alt_sca, azi_sca, pol_sca);

                solutions.push_back({ antenna_inc, antenna_sca, {} });
            }

            logger << " - Column size: " << mesh.size() * solutions.size() * size_to_memory << " Gb" << std::endl;

            std::stringstream system_name;

            system_name << mesh_name << "." << wavelength;

            system.set_block_size(block_size);

            system.set_matrix_size(mesh.size());

            system.set_column_size(solutions.size());

            system.set_name(system_name.str());

            logger << date() << " - Preparing matrix" << std::endl;

            system.prepare_matrix([&] (Length col, Length row, Co <Scalar> & val)
            {
                auto i_0_0 = Integrator <Scalar> :: efie(mesh.edge(col).face[0], mesh.edge(row).face[0], k);
                auto i_0_1 = Integrator <Scalar> :: efie(mesh.edge(col).face[0], mesh.edge(row).face[1], k);
                auto i_1_0 = Integrator <Scalar> :: efie(mesh.edge(col).face[1], mesh.edge(row).face[0], k);
                auto i_1_1 = Integrator <Scalar> :: efie(mesh.edge(col).face[1], mesh.edge(row).face[1], k);

                val = im(k) * (i_0_0 - i_0_1 - i_1_0 + i_1_1);
            });

            logger << date() << " - Preparing column" << std::endl;

            system.prepare_column([&] (Length col, Length row, Co <Scalar> & val)
            {
                auto i_0 = Integrator <Scalar> :: integral(mesh.edge(row).face[0], solutions[col].inc.v_front, k);
                auto i_1 = Integrator <Scalar> :: integral(mesh.edge(row).face[1], solutions[col].inc.v_front, k);

                val = ((i_0 - i_1), solutions[col].inc.v_polar);
            });

            logger << date() << " - Preparing solution" << std::endl;

            system.prepare_solution();

            logger << date() << " - Processing solution" << std::endl;

            system.process_solution([&] (Length col, Length row, const Co <Scalar> & val)
            {
                auto i_0 = Integrator <Scalar> :: integral(mesh.edge(row).face[0], solutions[col].sca.v_front, k);
                auto i_1 = Integrator <Scalar> :: integral(mesh.edge(row).face[1], solutions[col].sca.v_front, k);

                solutions[col].field.E_sca += val * (i_0 - i_1);
            });


            for (auto & solution : solutions)
            {
                solution.field.rcs_sm = AbsSquared((solution.field.E_sca * im(k), solution.sca.v_polar))
                                      / NormSquared(solution.inc.v_polar) / 4.0 / M_PI;

                solution.field.rcs_db = 10.0 * Log10(solution.field.rcs_sm);
            }

            for (auto & solution : solutions)
            {
                solution_file << wavelength            << ","
                              << solution.inc.pol      << ","
                              << solution.sca.pol      << ","
                              << solution.inc.azi      << ","
                              << solution.sca.azi      << ","
                              << solution.inc.alt      << ","
                              << solution.sca.alt      << ","
                              << solution.field.rcs_sm << ","
                              << solution.field.rcs_db << std::endl;
            }
        }

        logger << date() << " - Finished" << std::endl;
    }

};


#endif //RCS_SOLVER_H
