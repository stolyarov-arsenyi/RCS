#ifndef  RCS__BLOCKMATRIX_H
#define  RCS__BLOCKMATRIX_H


#include "Matrix/Lapack.h"
#include "Matrix/Pivots.h"
#include "Matrix/Matrix.h"


template <class Scalar, class Length>

struct BlockSystem
{
    static Length block_size;

    std::string name;

    struct
    {
        Length size;

        struct
        {
            std::vector <Length> size;
        }
        grid;
    }
    matrix;

    struct
    {
        Length cols;
        Length rows;

        struct
        {
            std::vector <Length> cols;
            std::vector <Length> rows;
        }
        grid;
    }
    column;


//    struct Matrix : :: Matrix <Scalar, Length>
//    {
//        std::string name;
//
//        Matrix (Length cols, Length rows, std::string name) : name(std::move(name))
//        {
//            Matrix::cols = cols;
//            Matrix::rows = rows;
//        }
//
//        bool load ()
//        {
//            Matrix::matrix.data.resize(Matrix::cols() * Matrix::rows());
//
//            auto data = (char *) Matrix::matrix.data.data();
//
//            auto size = (std::streamsize) (Matrix::cols() * Matrix::rows() * sizeof(Scalar));
//
//            return std::ifstream(name, std::ios::binary).read(data, size).good();
//        }
//
//        bool save () const
//        {
//            auto data = (const char *) Matrix::matrix.data.data();
//
//            auto size = (std::streamsize) (Matrix::cols() * Matrix::rows() * sizeof(Scalar));
//
//            return std::ofstream(name, std::ios::binary).write(data, size).good();
//        }
//
//        bool exist () const
//        {
//            return std::ifstream(name).is_open();
//        }
//    };
//
//    struct Pivots : :: Pivots <Length>
//    {
//        std::string name;
//
//        Pivots (Length size, std::string name) : name(std::move(name))
//        {
//            Pivots::size = size;
//        }
//
//        bool load ()
//        {
//            Pivots::pivots.data.resize(Pivots::size());
//
//            auto data = (char *) Pivots::pivots.data.data();
//
//            auto size = (std::streamsize) (Pivots::size() * sizeof(Length));
//
//            return std::ifstream(name, std::ios::binary).read(data, size).good();
//        }
//
//        bool save () const
//        {
//            auto data = (const char *) Pivots::pivots.data.data();
//
//            auto size = (std::streamsize) (Pivots::size() * sizeof(Scalar));
//
//            return std::ofstream(name, std::ios::binary).write(data, size).good();
//        }
//
//        bool exist () const
//        {
//            return std::ifstream(name).is_open();
//        }
//    };


    template <class ... Num>

    bool matrix_exist (Length c, Length r, Num ... n) const
    {
        auto name = matrix_name(c, r, n ...);

        return std::ifstream(name).is_open();
    }


    template <class ... Num>

    auto matrix_name (Length c, Length r, Num ... n) const -> std::string
    {
        std::stringstream string;

        string << name << "." << c << "." << r;

        for (auto num : { n ... })

            string << "." << num;

        string << ".matrix";

        return string.str();
    }

    auto pivots_name (Length n) const -> std::string
    {
        std::stringstream string;

        string << name << "." << n << ".pivots";

        return string.str();
    }


    template <class ... Num>

    auto matrix_load (Length c, Length r, Num ... n) const -> Matrix <Scalar, Length>
    {
        Matrix <Scalar, Length> mat(matrix.grid.size[c], matrix.grid.size[r]);

        auto data = (char *) mat[0];

        auto size = (std::streamsize) (mat.cols() * mat.rows() * sizeof(Scalar));

        auto name = matrix_name(c, r, n ...);

        std::ifstream(name, std::ios::binary).read(data, size).good();

        return mat;
    }

    auto pivots_load (Length n) const -> Pivots <Length>
    {
        Pivots <Length> piv(matrix.grid.size[n]);

        auto data = (char *) & piv[0];

        auto size = (std::streamsize) (piv.size() * sizeof(Scalar));

        auto name = pivots_name(n);

        std::ifstream(name, std::ios::binary).read(data, size).good();

        return piv;
    }


    template <class ... Num>

    void matrix_save (const Matrix <Scalar, Length> & mat, Length c, Length r, Num ... n) const
    {
        auto data = (const char *) mat[0];

        auto size = (std::streamsize) (mat.cols() * mat.rows() * sizeof(Scalar));

        auto name = matrix_name(c, r, n ...);

        std::ofstream(name, std::ios::binary).write(data, size).good();
    }

    void pivots_save (const Pivots <Length> & piv, Length n) const
    {
        auto data = (const char *) piv[0];

        auto size = (std::streamsize) (piv.cols() * piv.rows() * sizeof(Scalar));

        auto name = pivots_name(n);

        std::ofstream(name, std::ios::binary).write(data, size).good();
    }


    void prepare_system (const std::function <void (Length, Length, Scalar &)> & func)
    {
        Length N = matrix.grid.size.size();

        for (Length n = 0; n < N; n ++)
        {
            Length step = n;

            while (step > 0)
            {
                if (matrix_exist(n, n, step))

                    break;

                else step --;
            }

            auto A00 = matrix_load(n, n, step);


            if (step == 0 && ! matrix_exist(n, n, step))
            {
                /** Fill matrix **/

                Length offset = matrix.grid.size[n] * block_size;

                for (Length col = 0; col < A00.cols(); col ++)
                for (Length row = 0; row < A00.rows(); row ++)

                    func(offset + col, offset + row, A00[col][row]);

                matrix_save(A00, n, n, step);
            }

            /** A00(n + 1) = A00(n) - L10(n) * U01(n) **/

            while (step < n)
            {
                step ++;

                auto U01 = matrix_load(n, step - 1, step);
                auto L10 = matrix_load(step - 1, n, step);

                A00 -= L10 * U01;

                matrix_save(A00, n, n, step);
            }

            step ++;

            /** A00 = P00 * L00 * U00 **/

            if (! matrix_exist(n, n, step))
            {
                auto P00 = A00.lu_factorization();

                matrix_save(A00, n, n, step);

                pivots_save(P00, n);
            }

            /** Solve L00 * U01 = P00 * A01 **/

            for (Length col = n + 1; col < N; col ++)
            {
                if (matrix_exist(n, n, step))
                {
                    /** Fill matrix **/

                    Length offset = matrix.grid.size[n] * block_size;

                    for (Length col = 0; col < A00.cols(); col ++)
                    for (Length row = 0; row < A00.rows(); row ++)

                        func(offset + col, offset + row, A00[col][row]);

                    matrix_save(A00, n, n, step);
                }
            }
        }
    }

    void prepare_column (const std::function <void (Length, Length, Scalar &)> & func)
    {

    }

    void prepare_solution ()
    {

    }

    void process_solution (const std::function <void (Length, Length, Scalar &)> & func)
    {

    }
};


template <class Scalar, class Length>

Length BlockSystem <Scalar, Length> :: block_size = 0;


#endif //RCS__BLOCKMATRIX_H