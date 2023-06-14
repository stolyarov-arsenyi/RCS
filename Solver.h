#ifndef  RCS_SOLVER_H
#define  RCS_SOLVER_H


#include "Solver/Input.h"
#include "Solver/Solution.h"
#include "Solver/Integrator.h"


template <class Sca, class Len>

struct Solver
{
    static void impedance_matrix_decomposition (const Mesh <Sca> & mesh, const Re <Sca> & wavelength, const std::string & name)
    {
        auto k = 2.0 * M_PI / wavelength;

        BlockMatrix <Co <Sca>, Len> matrix((Len) mesh.size(), (Len) mesh.size());

        for (Len col = 0; col < matrix.grid.cols.size(); col ++)
        for (Len row = 0; row < matrix.grid.rows.size(); col ++)
        {

        }

        matrix.for_each_entry([&] (Len col, Len row, Co <Sca> & val)
        {
            auto i_0_0 = Integrator <Sca> :: efie(mesh.edge(col).face[0], mesh.edge(row).face[0], k);
            auto i_0_1 = Integrator <Sca> :: efie(mesh.edge(col).face[0], mesh.edge(row).face[1], k);
            auto i_1_0 = Integrator <Sca> :: efie(mesh.edge(col).face[1], mesh.edge(row).face[0], k);
            auto i_1_1 = Integrator <Sca> :: efie(mesh.edge(col).face[1], mesh.edge(row).face[1], k);

            val = im(k) * (i_0_0 - i_0_1 - i_1_0 + i_1_1);
        });
    }

    static auto solutions (const Mesh <Sca> & mesh, const Re <Sca> & wavelength, const std::vector <Source <Sca>> & sources, const BlockMatrix <Co <Sca>, Len> & inverse) -> std::vector <Solution <Sca>>
    {
        auto k = 2.0 * M_PI / wavelength;

        std::vector <Solution <Sca>> solutions;

        for (const auto & source : sources)
        {
            auto radians = M_PI / 180.0;

            Field <double> field;

            field.front = { 0.0, 0.0, 1.0 };
            field.E_inc = { 0.0, 1.0, 0.0 };

            field.E_inc = Rotate(field.E_inc, { 0.0,                  0.0, source.pol * radians });
            field.E_inc = Rotate(field.E_inc, { 0.0, source.alt * radians, source.azi * radians });
            field.front = Rotate(field.front, { 0.0, source.alt * radians, source.azi * radians });

            solutions.push_back({ source, field });
        }


        BlockMatrix <Co <Sca>, Len> column((Len) solutions.size(), (Len) mesh.size());

        column.for_each_entry([&] (Len col, Len row, Co <Sca> & val)
        {
            auto i_0 = Integrator <Sca> :: integral(mesh.edge(row).face[0], solutions[col].field.front, k);
            auto i_1 = Integrator <Sca> :: integral(mesh.edge(row).face[1], solutions[col].field.front, k);

            val = ((i_0 - i_1), solutions[col].field.E_inc);
        });

        BlockMatrix <Co <Sca>, Len> result = inverse * column;

        result.for_each_entry([&] (Len col, Len row, Co <Sca> & val)
        {
            auto i_0 = Integrator <Sca> :: integral(mesh.edge(row).face[0], solutions[col].field.front, k);
            auto i_1 = Integrator <Sca> :: integral(mesh.edge(row).face[1], solutions[col].field.front, k);

            solutions[col].field.E_sca += val * (i_0 - i_1);
        });

        for (auto & solution : solutions)
        {
            solution.field.rcs_sm = AbsSquared((solution.field.E_sca * im(k), solution.field.E_inc)) / (NormSquared(solution.field.E_inc) * 4.0 * M_PI);

            solution.field.rcs_db = 10.0 * Log10(solution.field.rcs_sm);
        }

        return solutions;
    }
};


#endif //RCS_SOLVER_H
