#ifndef  RCS__BLOCKSYSTEM_H
#define  RCS__BLOCKSYSTEM_H


#include "Matrix/Lapack.h"
#include "Matrix/Pivots.h"
#include "Matrix/Matrix.h"


template <class Scalar, class Length>

struct BlockSystem
{
    Length block_size;

    std::string name;

    struct
    {
        Length size;

        std::vector <Length> grid;
    }
    matrix,
    column;

    void set_name (const std::string & name)
    {
        this -> name = name;
    }

    void set_block_size (Length size)
    {
        block_size = size;
    }

    void set_matrix_size (Length size)
    {
        matrix.size = size;

        matrix.grid.resize(size / block_size, block_size);

        if (size %= block_size) matrix.grid.push_back(size);
    }

    void set_column_size (Length size)
    {
        column.size = size;

        column.grid.resize(size / block_size, block_size);

        if (size %= block_size) column.grid.push_back(size);
    }


    bool matrix_exist (Length c, Length r, Length n) const
    {
        auto name = matrix_name(c, r, n);

        return std::ifstream(name).is_open();
    }

    bool column_exist (Length c, Length r, Length n) const
    {
        auto name = column_name(c, r, n);

        return std::ifstream(name).is_open();
    }

    bool pivots_exist (Length n) const
    {
        auto name = pivots_name(n);

        return std::ifstream(name).is_open();
    }


    auto matrix_name (Length c, Length r, Length n) const -> std::string
    {
        std::stringstream string;

        string << name << "." << block_size << "." << c << "." << r << "." << n << ".matrix";

        return string.str();
    }

    auto column_name (Length c, Length r, Length n) const -> std::string
    {
        std::stringstream string;

        string << name << "." << block_size << "." << c << "." << r << "." << n << ".column";

        return string.str();
    }

    auto pivots_name (Length n) const -> std::string
    {
        std::stringstream string;

        string << name << "." << block_size << "." << n << ".pivots";

        return string.str();
    }


    auto matrix_load (Length c, Length r, Length n) const -> Matrix <Scalar, Length>
    {
        Matrix <Scalar, Length> mat(matrix.grid[c], matrix.grid[r]);

        auto data = (char *) mat[0];

        auto size = (std::streamsize) (mat.cols() * mat.rows() * sizeof(Scalar));

        auto name = matrix_name(c, r, n);

        std::ifstream(name, std::ios::binary).read(data, size).good();

        return mat;
    }

    auto column_load (Length c, Length r, Length n) const -> Matrix <Scalar, Length>
    {
        Matrix <Scalar, Length> col(column.grid[c], matrix.grid[r]);

        auto data = (char *) col[0];

        auto size = (std::streamsize) (col.cols() * col.rows() * sizeof(Scalar));

        auto name = column_name(c, r, n);

        std::ifstream(name, std::ios::binary).read(data, size).good();

        return col;
    }

    auto pivots_load (Length n) const -> Pivots <Length>
    {
        Pivots <Length> piv(matrix.grid[n]);

        auto data = (char *) & piv[0];

        auto size = (std::streamsize) (piv.size() * sizeof(Length));

        auto name = pivots_name(n);

        std::ifstream(name, std::ios::binary).read(data, size).good();

        return piv;
    }


    void matrix_save (const Matrix <Scalar, Length> & mat, Length c, Length r, Length n) const
    {
        auto data = (const char *) mat[0];

        auto size = (std::streamsize) (mat.cols() * mat.rows() * sizeof(Scalar));

        auto name = matrix_name(c, r, n);

        std::ofstream(name, std::ios::binary).write(data, size).good();
    }

    void column_save (const Matrix <Scalar, Length> & col, Length c, Length r, Length n) const
    {
        auto data = (const char *) col[0];

        auto size = (std::streamsize) (col.cols() * col.rows() * sizeof(Scalar));

        auto name = column_name(c, r, n);

        std::ofstream(name, std::ios::binary).write(data, size).good();
    }

    void pivots_save (const Pivots <Length> & piv, Length n) const
    {
        auto data = (const char *) & piv[0];

        auto size = (std::streamsize) (piv.size() * sizeof(Length));

        auto name = pivots_name(n);

        std::ofstream(name, std::ios::binary).write(data, size).good();
    }


    bool matrix_remove (Length c, Length r, Length n)
    {
        auto name = matrix_name(c, r, n);

        return std::remove(name.c_str());
    }

    bool column_remove (Length c, Length r, Length n)
    {
        auto name = column_name(c, r, n);

        return std::remove(name.c_str());
    }


    void prepare_A00 (Length n, const std::function <void (Length, Length, Scalar &)> & func)
    {
        Length step = n + 1;

        while (step > 0)
        {
            if (matrix_exist(n, n, step))

                break;

            else step --;
        }

        /** If A00 exist **/

        if (step == n + 1)
        {
            if (pivots_exist(n))

                return;

            /** If P00 doesn't exist - compute A00 from scratch**/

            matrix_remove(n, n, step);

            step = 0;
        }

        /** If A00(0) doesn't exist **/

        if (step == 0 && ! matrix_exist(n, n, step))
        {
            /** Fill matrix **/

            matrix_update(n, n, step, func);
        }

        auto A00 = matrix_load(n, n, step);

        /** A00(n + 1) = A00(n) - L10(n) * U01(n) **/

        while (++ step < n + 1)
        {
            auto U01 = matrix_load(n, step - 1, step);
            auto L10 = matrix_load(step - 1, n, step);

            A00 -= L10 * U01;

            matrix_remove(n, n, step - 1);

            matrix_save(A00, n, n, step);
        }

        /** A00 = P00 * L00 * U00 **/

        auto P00 = A00.lu_factorization();

        matrix_remove(n, n, step - 1);

        matrix_save(A00, n, n, step);

        pivots_save(P00, n);
    };

    void prepare_U01 (Length c, Length n, const std::function <void (Length, Length, Scalar &)> & func)
    {
        Length step = n + 1;

        while (step > 0)
        {
            if (matrix_exist(c, n, step))

                break;

            else step --;
        }

        if (step == n + 1)

            return;

        /** If A01(0) doesn't exist **/

        if (step == 0 && ! matrix_exist(c, n, step))
        {
            /** Fill matrix **/

            matrix_update(c, n, step, func);
        }

        auto A01 = matrix_load(c, n, step);

        /** A01(n + 1) = A01(n) - L10(n) * U01(n) **/

        while (++ step < n + 1)
        {
            auto U01 = matrix_load(c, step - 1, step);
            auto L10 = matrix_load(step - 1, n, step);

            A01 -= L10 * U01;

            matrix_remove(c, n, step - 1);

            matrix_save(A01, c, n, step);
        }

        /** Solve P00 * A01 = L00 * U01 **/

        auto A00 = matrix_load(n, n, step);

        auto P00 = pivots_load(n);

        A01 = P00 * A01;

        A01 = SolveTriangular(A00, UPLO::Lower, DIAG::Unit, SIDE::Left, TRAN::None);

        matrix_remove(c, n, step - 1);

        matrix_save(A01, c, n, step);
    };

    void prepare_L10 (Length n, Length r, const std::function <void (Length, Length, Scalar &)> & func)
    {
        Length step = n + 1;

        while (step > 0)
        {
            if (matrix_exist(n, r, step))

                break;

            else step --;
        }

        if (step == n + 1)

            return;

        /** If A10(0) doesn't exist **/

        if (step == 0 && ! matrix_exist(n, r, step))
        {
            /** Fill matrix **/

            matrix_update(n, r, step, func);
        }

        auto A10 = matrix_load(n, r, step);

        /** A10(n + 1) = A10(n) - L10(n) * U01(n) **/

        while (++ step < n + 1)
        {
            auto U01 = matrix_load(n, step - 1, step);
            auto L10 = matrix_load(step - 1, r, step);

            A10 -= L10 * U01;

            matrix_remove(n, r, step - 1);

            matrix_save(A10, n, r, step);
        }

        /** Solve A10 = L10 * U00 **/

        auto A00 = matrix_load(n, n, step);

        A10 = SolveTriangular(A00, UPLO::Upper, DIAG::NonUnit, SIDE::Right, TRAN::None);

        matrix_remove(n, r, step - 1);

        matrix_save(A10, n, r, step);
    };


    void prepare_matrix (const std::function <void (Length, Length, Scalar &)> & func)
    {
        Length N = matrix.grid.size();

        for (Length n = 0; n < N; n ++)
        {
            prepare_A00(n, func);

            for (Length c = n + 1; c < N; c ++)

                prepare_U01(c, n, func);

            for (Length r = n + 1; r < N; r ++)

                prepare_L10(n, r, func);
        }
    }

    void prepare_column (const std::function <void (Length, Length, Scalar &)> & func)
    {
        for (Length c = 0; c < column.grid.size(); c ++)
        for (Length r = 0; r < matrix.grid.size(); r ++)
        {
            column_update(c, r, 0, func);
        }
    }


    void prepare_solution ()
    {
        Length N = matrix.grid.size();
        Length C = column.grid.size();

        /** Forward substitution **/

        for (Length n = 0; n < N; n ++)
        {
            auto L = matrix_load(n, n, n + 1);

            auto P = pivots_load(n);

            for (Length c = 0; c < C; c ++)
            {
                auto Y = column_load(c, n, 0);

                Y = P * Y;

                Y = SolveTriangular(L, UPLO::Lower, DIAG::Unit, SIDE::Left, TRAN::None);

                column_save(Y, c, n, 0);

                for (Length r = n + 1; r < N; r ++)
                {
                    auto B = column_load(c, r, 0);

                    auto A = matrix_load(n, r, n + 1);

                    B -= A * Y;

                    column_save(B, c, r, 0);
                }
            }
        }

        /** Backward substitution **/

        for (Length n = N - 1; n >= 0; n --)
        {
            auto U = matrix_load(n, n, n + 1);

            for (Length c = 0; c < C; c ++)
            {
                auto X = column_load(c, n, 0);

                X = SolveTriangular(U, UPLO::Upper, DIAG::NonUnit, SIDE::Left, TRAN::None);

                column_save(X, c, n, 0);

                for (Length r = n - 1; r >= 0; r --)
                {
                    auto Y = column_load(c, r, 0);

                    auto A = matrix_load(n, r, r + 1);

                    Y -= A * X;

                    column_save(Y, c, r, 0);
                }
            }
        }
    }

    void process_solution (const std::function <void (Length, Length, const Scalar &)> & func)
    {
        for (Length c = 0; c < column.grid.size(); c ++)
        for (Length r = 0; r < matrix.grid.size(); r ++)
        {
            column_traversal(c, r, 0, func);
        }
    }


    void matrix_traversal (Length c, Length r, Length n, const std::function <void (Length, Length, const Scalar &)> & func) const
    {
        auto matrix = matrix_load(c, r, n);

        Length cols = c * block_size;
        Length rows = r * block_size;

        #pragma omp parallel for default(none) shared(cols, rows, matrix, func)

        for (Length col = 0; col < matrix.cols(); col ++)
        for (Length row = 0; row < matrix.rows(); row ++)

            func(cols + col, rows + row, matrix[col][row]);
    }

    void column_traversal (Length c, Length r, Length n, const std::function <void (Length, Length, const Scalar &)> & func) const
    {
        auto column = column_load(c, r, n);

        Length cols = c * block_size;
        Length rows = r * block_size;

        #pragma omp parallel for default(none) shared(cols, rows, column, func)

        for (Length col = 0; col < column.cols(); col ++)
        for (Length row = 0; row < column.rows(); row ++)

            func(cols + col, rows + row, column[col][row]);
    }


    void matrix_update (Length c, Length r, Length n, const std::function <void (Length, Length, Scalar &)> & func) const
    {
        auto matrix = matrix_load(c, r, n);

        Length cols = c * block_size;
        Length rows = r * block_size;

        #pragma omp parallel for default(none) shared(cols, rows, matrix, func)

        for (Length col = 0; col < matrix.cols(); col ++)
        for (Length row = 0; row < matrix.rows(); row ++)

            func(cols + col, rows + row, matrix[col][row]);

        matrix_save(matrix, c, r, n);
    }

    void column_update (Length c, Length r, Length n, const std::function <void (Length, Length, Scalar &)> & func) const
    {
        auto column = column_load(c, r, n);

        Length cols = c * block_size;
        Length rows = r * block_size;

        #pragma omp parallel for default(none) shared(cols, rows, column, func)

        for (Length col = 0; col < column.cols(); col ++)
        for (Length row = 0; row < column.rows(); row ++)

            func(cols + col, rows + row, column[col][row]);

        column_save(column, c, r, n);
    }


    auto matrix_to_string () const -> std::string
    {
        auto system_size = matrix.size;

        std::vector <Scalar> array(system_size * system_size);

        Length N = matrix.grid.size();

        for (Length n = 0; n < N; n ++)
        {
            matrix_traversal(n, n, n + 1, [&] (Length col, Length row, Scalar & val)
            {
                array[col * system_size + row] = val;
            });

            for (Length c = n + 1; c < N; c ++)
            {
                matrix_traversal(c, n, n + 1, [&] (Length col, Length row, Scalar & val)
                {
                    array[col * system_size + row] = val;
                });
            }

            for (Length r = n + 1; r < N; r ++)
            {
                matrix_traversal(n, r, n + 1, [&] (Length col, Length row, Scalar & val)
                {
                    array[col * system_size + row] = val;
                });
            }
        }

        std::stringstream string;

        string << std::showpos << std::fixed;

        for (Length row = 0; row < system_size; row ++)
        for (Length col = 0; col < system_size; col ++)

            string << array[col * system_size + row] << (col == system_size - 1 ? "\n" : ",");

        return string.str();
    }

    auto column_to_string () const -> std::string
    {
        std::vector <Scalar> array(column.size * matrix.size);

        for (Length c = 0; c < column.grid.size(); c ++)
        for (Length r = 0; r < matrix.grid.size(); r ++)
        {
            column_traversal(c, r, 0, [&] (Length col, Length row, Scalar & val)
            {
                array[col * matrix.size + row] = val;
            });
        }

        std::stringstream string;

        string << std::showpos << std::fixed;

        for (Length row = 0; row < matrix.size; row ++)
        for (Length col = 0; col < column.size; col ++)

            string << array[col * matrix.size + row] << (col == column.size - 1 ? "\n" : ",");

        return string.str();
    }
};


#endif //RCS__BLOCKSYSTEM_H