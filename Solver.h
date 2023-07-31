#include "Solver/Logger.h"
#include "Solver/Quadrature.h"
#include "Solver/STL.h"
#include "Solver/Face.h"
#include "Solver/Edge.h"
#include "Solver/Mesh.h"
#include "Solver/Input.h"
#include "Solver/Integrator.h"


template <class Scalar, class Length>

class Solver
{
    struct Params
    {
        Scalar wavelength;

        Length omp;

        bool monostatic;

        double block;

        struct
        {
            struct
            {
                Scalar min;
                Scalar max;
                Scalar inc;
            }
            phi,
            the;
        }
        inc,
        sca;

        struct
        {
            Scalar cos;
            Scalar sin;
        }
        tau,
        psi;
    };

    struct Angles
    {
        struct Angle
        {
            Scalar phi;
            Scalar the;
        }
        inc,
        sca;
    };

    struct Fields
    {
        struct Field
        {
            Vector <Re <Scalar>> e;
            Vector <Re <Scalar>> v;
            Vector <Re <Scalar>> h;

            Co <Scalar> p_x;
            Co <Scalar> p_y;

            Vector <Co <Scalar>> p;

            Field (const Vector <Re <Scalar>> & dir, const class Angles::Angle & angle, const Params & params) : e(dir)
            {
                auto radians = M_PI / 180.0;

                auto cos_the = Cos(angle.the * radians);
                auto sin_the = Sin(angle.the * radians);
                auto cos_phi = Cos(angle.phi * radians);
                auto sin_phi = Sin(angle.phi * radians);

                Vector <Re <Scalar>> rot [3] =
                {
                    { + cos_the * cos_phi, - sin_phi, + sin_the * cos_phi },
                    { + cos_the * sin_phi, + cos_phi, + sin_the * sin_phi },
                    { - sin_the          , +     0.0, + cos_the           }
                };

                h = { 0.0, 1.0, 0.0 };

                v = h % e;

                e = { (rot[0], e), (rot[1], e), (rot[2], e) };
                v = { (rot[0], v), (rot[1], v), (rot[2], v) };
                h = { (rot[0], h), (rot[1], h), (rot[2], h) };

                p_x = co(+ params.psi.cos * params.tau.cos, - params.tau.sin * params.psi.sin);
                p_y = co(+ params.psi.sin * params.tau.cos, + params.psi.cos * params.tau.sin);

                p = p_x * v + p_y * h;
            }
        }
        inc,
        sca;

        Fields (const Angles & angles, const Params & params) : inc({ + 0.0, + 0.0, - 1.0 }, angles.inc, params),
                                                                sca({ + 0.0, + 0.0, + 1.0 }, angles.sca, params) {}
    };

    struct Output
    {
        Vector <Co <Scalar>> e_field;

        Re <Scalar> vv = 0.0;
        Re <Scalar> vh = 0.0;
        Re <Scalar> hv = 0.0;
        Re <Scalar> hh = 0.0;
    };


    std::string conf_name;
    std::string mesh_name;

    Mesh <Scalar> mesh;

    BlockSystem <Co <Scalar>, Length> system;

    Params params;

    Logger logger;

    std::vector <Angles> angles;
    std::vector <Fields> fields;
    std::vector <Output> output;

public:

    Solver (std::string conf_name, std::string mesh_name) : conf_name(std::move(conf_name)), mesh_name(std::move(mesh_name)) {}

    void init ()
    {
        Input input = Input::parse(conf_name);

        std::stringstream(input["block"]) >> params.block;

        std::stringstream(input["nopenmp"]) >> params.omp;
        std::stringstream(input["lymda"  ]) >> params.wavelength;

        std::stringstream(input["ptay0"  ]) >> params.tau.cos >> params.tau.sin;
        std::stringstream(input["ppsi0"  ]) >> params.psi.cos >> params.psi.sin;

        std::stringstream(input["phi0"  ]) >> params.inc.phi.min >> params.inc.phi.max;
        std::stringstream(input["phi1"  ]) >> params.sca.phi.min >> params.sca.phi.max;
        std::stringstream(input["tet0"  ]) >> params.inc.the.min >> params.inc.the.max;
        std::stringstream(input["tet1"  ]) >> params.sca.the.min >> params.sca.the.max;

        std::stringstream(input["dphi0" ]) >> params.inc.phi.inc;
        std::stringstream(input["dphi1" ]) >> params.sca.phi.inc;
        std::stringstream(input["dtet0" ]) >> params.inc.the.inc;
        std::stringstream(input["dtet1" ]) >> params.sca.the.inc;

        params.monostatic = input[ "phi1"].empty() || input["dphi1"].empty()
                         || input[ "tet1"].empty() || input["dtet1"].empty();

        if (params.monostatic) params.sca = params.inc;

        std::filesystem::create_directories(logger_name());
        std::filesystem::create_directories(matrix_name());
        std::filesystem::create_directories(column_name());
        std::filesystem::create_directories(result_name());

        logger = Logger(logger_name() + "log.txt");

        logger.date() << " - Preparing mesh" << std::endl;

        logger << " - Try to load parsed mesh" << std::endl;

        if (! mesh.load_binary(mesh_name + ".bin"))
        {
            logger << " - Parsed mesh doesn't exist" << std::endl;

            logger << " - Parsing mesh" << std::endl;

            mesh = Mesh <Scalar> :: load_from_stl(mesh_name);

            logger << " - Saving parsed mesh" << std::endl;

            mesh.save_binary(mesh_name + ".bin");
        }

        logger << " - Mesh size: " << mesh.size() << " edges" << std::endl;
    }

    void solve ()
    {
        auto k = 2.0 * M_PI / params.wavelength;

        omp_set_num_threads(params.omp);


        logger.date() << " - Preparing fields" << std::endl;

        if (params.monostatic)
        {
            for (auto inc_phi : range(params.inc.phi.min, params.inc.phi.max, params.inc.phi.inc))
            for (auto inc_the : range(params.inc.the.min, params.inc.the.max, params.inc.the.inc))

                angles.push_back({{ inc_phi, inc_the }, { inc_phi, inc_the }});
        }
        else
        {
            for (auto inc_phi : range(params.inc.phi.min, params.inc.phi.max, params.inc.phi.inc))
            for (auto inc_the : range(params.inc.the.min, params.inc.the.max, params.inc.the.inc))
            for (auto sca_phi : range(params.sca.phi.min, params.sca.phi.max, params.sca.phi.inc))
            for (auto sca_the : range(params.sca.the.min, params.sca.the.max, params.sca.the.inc))

                angles.push_back({{ inc_phi, inc_the }, { sca_phi, sca_the }});
        }

        for (auto angle : angles)

            fields.emplace_back(angle, params);

        output.resize(fields.size());


        system.set_matrix_name(matrix_name());

        system.set_column_name(column_name());

        if (params.block == 0.0)

            params.block = 4.0;

        system.set_block_memory_usage(params.block);

        system.set_matrix_size(mesh.size());

        system.set_column_size(angles.size());

        logger.date() << " - Preparing matrix" << std::endl << " - Matrix size: " << system.matrix_memory_usage() << " Gb" << std::endl;

        system.prepare_matrix([&] (Length col, Length row, Co <Scalar> & val)
        {
            auto i_0_0 = Integrator <Scalar> :: efie(mesh.edge(col).face[0], mesh.edge(row).face[0], k);
            auto i_0_1 = Integrator <Scalar> :: efie(mesh.edge(col).face[0], mesh.edge(row).face[1], k);
            auto i_1_0 = Integrator <Scalar> :: efie(mesh.edge(col).face[1], mesh.edge(row).face[0], k);
            auto i_1_1 = Integrator <Scalar> :: efie(mesh.edge(col).face[1], mesh.edge(row).face[1], k);

            val = im(k) * (i_0_0 - i_0_1 - i_1_0 + i_1_1);
        });

        logger.date() << " - Preparing column" << std::endl << " - Column size: " << system.column_memory_usage() << " Gb" << std::endl;

        system.prepare_column([&] (Length col, Length row, Co <Scalar> & val)
        {
            auto i_0 = Integrator <Scalar> :: integral(mesh.edge(row).face[0], - fields[col].inc.e, k);
            auto i_1 = Integrator <Scalar> :: integral(mesh.edge(row).face[1], - fields[col].inc.e, k);

            val = ((i_0 - i_1), fields[col].inc.p);
        });

        logger.date() << " - Processing column" << std::endl;

        system.process_column([&] (Length col, Length row, Co <Scalar> & val)
        {
            auto i_0 = Integrator <Scalar> :: integral(mesh.edge(row).face[0], + fields[col].sca.e, k);
            auto i_1 = Integrator <Scalar> :: integral(mesh.edge(row).face[1], + fields[col].sca.e, k);

            output[col].e_field += val * (i_0 - i_1);
        });

        logger.date() << " - Processing result" << std::endl;

        for (std::size_t o = 0; o < output.size(); o ++)
        {
            output[o].e_field *= im(k);

            output[o].hh = AbsSquare((output[o].e_field, fields[o].sca.h)) / AbsSquare(fields[o].inc.p_y) / 4.0 / M_PI;
            output[o].hv = AbsSquare((output[o].e_field, fields[o].sca.h)) / AbsSquare(fields[o].inc.p_x) / 4.0 / M_PI;
            output[o].vh = AbsSquare((output[o].e_field, fields[o].sca.v)) / AbsSquare(fields[o].inc.p_y) / 4.0 / M_PI;
            output[o].vv = AbsSquare((output[o].e_field, fields[o].sca.v)) / AbsSquare(fields[o].inc.p_x) / 4.0 / M_PI;
        }

        logger.date() << " - Saving result" << std::endl;

        save_result_matrix_sm();

        save_result_matrix_db();

        save_result_binary_ef();

        logger.date() << " - Finished" << std::endl;
    }

private:

    void save_result_matrix_sm () const
    {
        std::ofstream file(result_name() + "rcs.sm.csv");

        file << "phi_inc, the_inc, phi_sca, the_sca, hh, hv, vh, vv" << std::endl;

        file << std::fixed << std::showpos;

        for (std::size_t o = 0; o < output.size(); o ++)
        {
            file << angles[o].inc.phi << ","
                 << angles[o].inc.the << ","
                 << angles[o].sca.phi << ","
                 << angles[o].sca.the << ","

                 << output[o].hh << ","
                 << output[o].hv << ","
                 << output[o].vh << ","
                 << output[o].vv << std::endl;
        }
    }

    void save_result_matrix_db () const
    {
        std::ofstream file(result_name() + "rcs.db.csv");

        file << "phi_inc, the_inc, phi_sca, the_sca, hh, hv, vh, vv" << std::endl;

        file << std::fixed << std::showpos;

        for (std::size_t o = 0; o < output.size(); o ++)
        {
            file << angles[o].inc.phi << ","
                 << angles[o].inc.the << ","
                 << angles[o].sca.phi << ","
                 << angles[o].sca.the << ","

                 << 10.0 * Log10(output[o].hh) << ","
                 << 10.0 * Log10(output[o].hv) << ","
                 << 10.0 * Log10(output[o].vh) << ","
                 << 10.0 * Log10(output[o].vv) << std::endl;
        }
    }

    void save_result_binary_ef () const
    {
        std::ofstream file(result_name() + "rcs.ef.bin", std::ios::binary);

        std::uint64_t size = output.size();

        file.write((char *) & size, 8);

        for (std::size_t o = 0; o < output.size(); o ++)
        {
            file.write((char *) & angles[o].inc.phi, sizeof(Scalar));
            file.write((char *) & angles[o].inc.the, sizeof(Scalar));
            file.write((char *) & angles[o].sca.phi, sizeof(Scalar));
            file.write((char *) & angles[o].sca.the, sizeof(Scalar));

            file.write((char *) & output[o].e_field.x.r, sizeof(Scalar));
            file.write((char *) & output[o].e_field.x.i, sizeof(Scalar));
            file.write((char *) & output[o].e_field.y.r, sizeof(Scalar));
            file.write((char *) & output[o].e_field.y.i, sizeof(Scalar));
            file.write((char *) & output[o].e_field.z.r, sizeof(Scalar));
            file.write((char *) & output[o].e_field.z.i, sizeof(Scalar));
        }
    }

    auto range (Scalar min, Scalar max, Scalar inc) const -> std::vector <Scalar>
    {
        std::vector <Scalar> range;

        while (min < max)
        {
            range.push_back(min);

            min += inc;
        }

        range.push_back(max);

        return range;
    }


    auto matrix_name () const -> std::string
    {
        std::stringstream matrix_name;

        matrix_name << "matrix" << "/" << mesh_name << "/" << params.wavelength << "/";

        return matrix_name.str();
    }

    auto column_name () const -> std::string
    {
        std::stringstream column_name;

        column_name << "column" << "/" << mesh_name << "/"

        << "wavelength " << params.wavelength << "/"

        << "phi inc " << params.inc.phi.min << " " << params.inc.phi.max << " " << params.inc.phi.inc << "/"
        << "the inc " << params.inc.the.min << " " << params.inc.the.max << " " << params.inc.the.inc << "/"
        << "phi sca " << params.sca.phi.min << " " << params.sca.phi.max << " " << params.sca.phi.inc << "/"
        << "the sca " << params.sca.the.min << " " << params.sca.the.max << " " << params.sca.the.inc << "/"

        << "tau " << params.tau.cos << " " << params.tau.sin << "/"
        << "psi " << params.psi.cos << " " << params.psi.sin << "/";

        return column_name.str();
    }

    auto result_name () const -> std::string
    {
        std::stringstream result_name;

        result_name << "result" << "/" << mesh_name << "/"

        << "wavelength " << params.wavelength << "/"

        << "phi inc " << params.inc.phi.min << " " << params.inc.phi.max << " " << params.inc.phi.inc << "/"
        << "the inc " << params.inc.the.min << " " << params.inc.the.max << " " << params.inc.the.inc << "/"
        << "phi sca " << params.sca.phi.min << " " << params.sca.phi.max << " " << params.sca.phi.inc << "/"
        << "the sca " << params.sca.the.min << " " << params.sca.the.max << " " << params.sca.the.inc << "/"

        << "tau " << params.tau.cos << " " << params.tau.sin << "/"
        << "psi " << params.psi.cos << " " << params.psi.sin << "/";

        return result_name.str();
    }

    auto logger_name () const -> std::string
    {
        std::stringstream result_name;

        result_name << "logger" << "/" << mesh_name << "/"

        << "wavelength " << params.wavelength << "/"

        << "phi inc " << params.inc.phi.min << " " << params.inc.phi.max << " " << params.inc.phi.inc << "/"
        << "the inc " << params.inc.the.min << " " << params.inc.the.max << " " << params.inc.the.inc << "/"
        << "phi sca " << params.sca.phi.min << " " << params.sca.phi.max << " " << params.sca.phi.inc << "/"
        << "the sca " << params.sca.the.min << " " << params.sca.the.max << " " << params.sca.the.inc << "/"

        << "tau " << params.tau.cos << " " << params.tau.sin << "/"
        << "psi " << params.psi.cos << " " << params.psi.sin << "/";

        return result_name.str();
    }
};