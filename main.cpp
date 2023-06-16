#include <functional>
#include <stdexcept>
#include <iostream>
#include <numeric>
#include <sstream>
#include <fstream>
#include <vector>
#include <memory>
#include <cmath>
#include <regex>


#include "Complex.h"
#include "Vector.h"
#include "BlockSystem.h"


#define MKL_Complex16 Co <double>


#include <mkl.h>
#include <omp.h>


#include "Specializations.h"
#include "Quadrature.h"
#include "Mesh.h"
#include "Solver.h"


template <class Scalar, class Length>

struct Solver
{
    std::string conf_name;
    std::string mesh_name;

    Length block_size;

    Mesh <Scalar> mesh;

    Input <Scalar> input;

    BlockSystem <Co <Scalar>, Length> system;


    Solver (const std::string & conf_name, const std::string & mesh_name) : conf_name(conf_name), mesh_name(mesh_name) {}

    auto date () -> std::string
    {
        char str [32];

        auto time = std::time(nullptr);

        strftime(str, 32, "%a %b %d %H:%M:%S %Y", localtime(& time));

        return { str };
    }

    void init ()
    {
        std::cout << date() << " - Parsing input" << std::endl;

        input = Input <Scalar> :: parse(conf_name);

        std::cout << date() << " - Parsing mesh" << std::endl;

        mesh = Mesh <Scalar> :: load_from_stl(mesh_name);


        auto threads = input["threads"][0];

        omp_set_num_threads((int) threads);


        auto memory = input["block_size"][0];

        auto size_to_memory = sizeof (Co <Scalar>) / Pow(1024.0, 3);

        block_size = (Length) Sqrt(memory / size_to_memory);


        std::cout << " - Mesh size: " << mesh.size() << " edges" << std::endl;

        std::cout << " - Matrix size: " << Pow(mesh.size(), 2.0) * size_to_memory << " Gb" << std::endl;
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
            std::cout << " - Wavelength: " << wavelength << "m" << std::endl;

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


            std::stringstream system_name;

            system_name << mesh_name << "." << wavelength;

            system.set_block_size(block_size);

            system.set_matrix_size(mesh.size());

            system.set_column_size(solutions.size());

            system.set_name(system_name.str());

            std::cout << date() << " - Preparing matrix" << std::endl;

            system.prepare_matrix([&] (Length col, Length row, Co <Scalar> & val)
            {
                auto i_0_0 = Integrator <Scalar> :: efie(mesh.edge(col).face[0], mesh.edge(row).face[0], k);
                auto i_0_1 = Integrator <Scalar> :: efie(mesh.edge(col).face[0], mesh.edge(row).face[1], k);
                auto i_1_0 = Integrator <Scalar> :: efie(mesh.edge(col).face[1], mesh.edge(row).face[0], k);
                auto i_1_1 = Integrator <Scalar> :: efie(mesh.edge(col).face[1], mesh.edge(row).face[1], k);

                val = im(k) * (i_0_0 - i_0_1 - i_1_0 + i_1_1);
            });

            std::cout << date() << " - Preparing column" << std::endl;

            system.prepare_column([&] (Length col, Length row, Co <Scalar> & val)
            {
                auto i_0 = Integrator <Scalar> :: integral(mesh.edge(row).face[0], solutions[col].inc.v_front, k);
                auto i_1 = Integrator <Scalar> :: integral(mesh.edge(row).face[1], solutions[col].inc.v_front, k);

                val = ((i_0 - i_1), solutions[col].inc.v_polar);
            });

            std::cout << date() << " - Preparing solution" << std::endl;

            system.prepare_solution();

            std::cout << date() << " - Processing solution" << std::endl;

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
                solution_file << wavelength
                              << solution.inc.pol
                              << solution.sca.pol
                              << solution.inc.azi
                              << solution.sca.azi
                              << solution.inc.alt
                              << solution.sca.alt
                              << solution.field.rcs_sm
                              << solution.field.rcs_db << std::endl;
            }
        }

        std::cout << date() << " - Finished" << std::endl;
    }

};


using Scalar = double;
using Length = MKL_INT;


int main (int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cout << "Wrong number of arguments" << std::endl;

        return 0;
    }

    std::string conf_name = argv[1];
    std::string mesh_name = argv[2];

    Solver <Scalar, Length> solver(conf_name, mesh_name);

    solver.init();
    solver.solve();


//    double arr [6][6] =
//    {
//        { 3.966261, 2.133467, 1.216731, 5.262742, 4.057779, 8.091164, },
//        { 2.133467, 5.288848, 6.238135, 2.130408, 6.813768, 8.245826, },
//        { 1.216731, 6.238135, 5.880156, 0.952626, 3.437681, 6.434359, },
//        { 5.262742, 2.130408, 0.952626, 6.579653, 2.719752, 2.975403, },
//        { 4.057779, 6.813768, 3.437681, 2.719752, 4.190280, 5.671886, },
//        { 8.091164, 8.245826, 6.434359, 2.975403, 5.671886, 7.223314, },
//    };
//
//
//    BlockSystem <Scalar, Length> block_system;
//
//    block_system.set_name("system");
//
//    block_system.set_block_size(1);
//
//    block_system.set_matrix_size(6);
//
//    block_system.set_column_size(2);
//
//    block_system.prepare_system([&] (Length col, Length row, Scalar & val)
//    {
//        val = arr[row][col];
//    });
//
//    std::cout << block_system.matrix_to_string() << std::endl;
//
//    block_system.prepare_column([&] (Length col, Length row, Scalar & val)
//   {
//       val = 1.0 + (Scalar) row;
//   });
//
//    std::cout << block_system.column_to_string() << std::endl;
//
//    block_system.prepare_solution();
//
//    std::cout << block_system.column_to_string() << std::endl;

    return 0;
}
