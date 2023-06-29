#ifndef RCS__BLOCKSYSTEM_H
#define RCS__BLOCKSYSTEM_H


#include "Matrix/Lapack.h"
#include "Matrix/Pivots.h"
#include "Matrix/Matrix.h"


template <class Scalar, class Length>

class BlockSystem
{
    class MatrixDescriptor;

    using Matrix = Matrix <Scalar, Length>;

    using Pivots = Pivots <Length>;


    struct
    {
        std::vector <Length> grid;

        std::string name;

        Length size = 0;
    }
    matrix,
    column;

    struct
    {
        std::vector <std::weak_ptr <MatrixDescriptor>> P;
        std::vector <std::weak_ptr <MatrixDescriptor>> A;
        std::vector <std::weak_ptr <MatrixDescriptor>> B;

        Length size = 0;
    }
    blocks;

public:

    auto matrix_name () const -> std::string
    {
        return matrix.name;
    }

    auto column_name () const -> std::string
    {
        return column.name;
    }


    void set_matrix_name (const std::string & name)
    {
        matrix.name = name;
    }

    void set_column_name (const std::string & name)
    {
        column.name = name;
    }


    auto blocks_size () const -> Length
    {
        return blocks.size;
    }

    auto matrix_size () const -> Length
    {
        return matrix.size;
    }

    auto column_size () const -> Length
    {
        return column.size;
    }


    void set_blocks_size (Length size)
    {
        blocks.size = size;

        set_matrix_size(matrix_size());
        set_column_size(column_size());
    }

    void set_matrix_size (Length size)
    {
        matrix.size = size;

        matrix.grid.resize(size / blocks.size, blocks.size);

        if (size %= blocks.size) matrix.grid.push_back(size);


        blocks.P.resize(matrix.grid.size());

        blocks.A.resize(matrix.grid.size() * matrix.grid.size());
    }

    void set_column_size (Length size)
    {
        column.size = size;

        column.grid.resize(size / blocks.size, blocks.size);

        if (size %= blocks.size) column.grid.push_back(size);


        blocks.B.resize(column.grid.size() * matrix.grid.size());
    }


    void set_block_memory_usage (double gb)
    {
        auto size_to_gb = sizeof (Scalar) / (1024.0 * 1024.0 * 1024.0);

        set_blocks_size(std::sqrt(gb / size_to_gb));
    }

    auto matrix_memory_usage () const -> double
    {
        auto size_to_gb = sizeof (Scalar) / (1024.0 * 1024.0 * 1024.0);

        return matrix_size() * matrix_size() * size_to_gb;
    }

    auto column_memory_usage () const -> double
    {
        auto size_to_gb = sizeof (Scalar) / (1024.0 * 1024.0 * 1024.0);

        return column_size() * matrix_size() * size_to_gb;
    }


    void prepare_matrix (const std::function <void (Length, Length, Scalar &)> & func)
    {
        Length N = matrix.grid.size();

        for (Length n = 0; n < N; n ++)
        {
            prepare_matrix_block(n, n, func);

            for (Length c = n + 1; c < N; c ++)

                prepare_matrix_block(c, n, func);

            for (Length r = n + 1; r < N; r ++)

                prepare_matrix_block(n, r, func);
        }
    }

    void prepare_column (const std::function <void (Length, Length, Scalar &)> & func)
    {
        Length C = column.grid.size();
        Length R = matrix.grid.size();

        for (Length c = 0; c < C; c ++)
            for (Length r = R; r > 0; r --)

                prepare_column_block_x(c, r - 1, func);

        for (Length c = 0; c < C; c ++)
            for (Length r = 0; r < R; r ++)
            {
                Length step = 0;

                while (step < R + 1)
                {
                    if (column_block(c, r).exist(step))

                        column_block(c, r).clear(step);

                    step ++;
                }
            }
    }

    void process_column (const std::function <void (Length, Length, Scalar &)> & func)
    {
        Length C = column.grid.size();
        Length R = matrix.grid.size();

        for (Length c = 0; c < C; c ++)
            for (Length r = 0; r < R; r ++)
            {
                auto C = column_block(c, r);

                C.init(R + 1);

                C.iterate(func);
            }
    }


    auto matrix_to_string () -> std::string
    {
        std::vector <Scalar> array(matrix.size * matrix.size);

        Length N = matrix.grid.size();

        for (Length n = 0; n < N; n ++)
        {
            for (Length c = n + 1; c < N; c ++)
            {
                auto U = matrix_block(c, n);

                U.init(n + 1);

                U.iterate([&] (Length col, Length row, Scalar & val)
                          {
                              array[col * matrix.size + row] = val;
                          });
            }

            for (Length r = n + 1; r < N; r ++)
            {
                auto L = matrix_block(n, r);

                L.init(n + 1);

                L.iterate([&] (Length col, Length row, Scalar & val)
                          {
                              array[col * matrix.size + row] = val;
                          });
            }

            auto A = matrix_block(n, n);

            A.init(n + 1);

            A.iterate([&] (Length col, Length row, Scalar & val)
                      {
                          array[col * matrix.size + row] = val;
                      });
        }

        std::stringstream string;

        string << std::showpos << std::fixed;

        for (Length row = 0; row < matrix.size; row ++)
            for (Length col = 0; col < matrix.size; col ++)

                string << array[col * matrix.size + row] << (col == matrix.size - 1 ? "\n" : ",");

        return string.str();
    }

    auto column_to_string () -> std::string
    {
        std::vector <Scalar> array(column.size * matrix.size);

        Length C = column.grid.size();
        Length R = matrix.grid.size();

        for (Length c = 0; c < C; c ++)
            for (Length r = 0; r < R; r ++)
            {
                auto C = column_block(c, r);

                C.init(R + 1);

                C.iterate([&] (Length col, Length row, Scalar & val)
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

private:

    void prepare_matrix_block (Length c, Length r, const std::function <void (Length, Length, Scalar &)> & func)
    {
        auto A = matrix_block(c, r);

        /** Determine diagonal index **/
        auto n = std::min(c, r);

        auto step = n + 1;

        /** Determine missing steps **/
        while (step > 0)
        {
            if (A.exist(step))

                break;

            step --;
        }

        /** If last step exist **/
        if (step == n + 1)
        {
            /** If it's 'L' or 'U' block - exit **/
            if (c != r)

                return;

            /** If it's diagonal block - check if it's pivots exist **/
            if (pivots_block(n).exist())

                return;

            /** If not - compute from diagonal block from scratch **/
            A.clear(step);

            step = 0;
        }

        A.init(step);

        /** If block doesn't exist fill it with 'func' **/
        if (step == 0 && ! A.exist(step))
        {
            A.iterate(func);

            A.update(step);
        }

        /** A(c,r) = L(n,r) * U(c, n) **/
        while (++ step < n + 1)
        {
            auto L = matrix_block(step - 1, r);
            auto U = matrix_block(c, step - 1);

            L.init(step);
            U.init(step);

            A -= L * U;

            A.update(step);
        }

        /** Diagonal block **/
        if (c == r)
        {
            auto P = pivots_block(n);

            P = A.lu_factorization();

            P.save();
        }
        else

            /** 'L' column **/
        if (c > r)
        {
            auto L = matrix_block(n, n);
            auto P = pivots_block(n);

            L.init(step);
            P.init();

            A = P * A;

            A = SolveTriangular(L, UPLO::Lower, DIAG::Unit, SIDE::Left, TRAN::None);
        }
        else

            /** 'U' row **/
        if (c < r)
        {
            auto U = matrix_block(n, n);

            U.init(step);

            A = SolveTriangular(U, UPLO::Upper, DIAG::NonUnit, SIDE::Right, TRAN::None);
        }

        A.update(step);
    }

    void prepare_column_block_x (Length c, Length r, const std::function <void (Length, Length, Scalar &)> & func)
    {
        auto B = column_block(c, r);

        Length N = matrix.grid.size();

        Length step = N + 1;

        while (step > r + 1)
        {
            if (B.exist(step))

                break;

            step --;
        }

        if (step == N + 1)

            return;

        if (step == r + 1 && ! B.exist(step))
        {
            for (Length row = 0; row < r + 1; row ++)

                prepare_column_block_y(c, row, func);
        }

        B.init(step);

        while (++ step < N + 1)
        {
            auto X = column_block(c, step - 1);

            auto U = matrix_block(step - 1, r);

            X.init(N + 1);

            U.init(r + 1);

            B -= U * X;

            B.update(step);
        }

        auto LU = matrix_block(r, r);

        LU.init(r + 1);

        B = SolveTriangular(LU, UPLO::Upper, DIAG::NonUnit, SIDE::Left, TRAN::None);

        B.update(step);
    }

    void prepare_column_block_y (Length c, Length r, const std::function <void (Length, Length, Scalar &)> & func)
    {
        auto B = column_block(c, r);

        Length N = matrix.grid.size();

        Length step = r + 1;

        while (step > 0)
        {
            if (B.exist(step))

                break;

            step --;
        }

        /** If forward substitution completed - exit **/

        if (step == r + 1)

            return;

        B.init(step);

        if (step == 0 && ! B.exist(step))
        {
            B.iterate(func);

            B.update(step);
        }

        while (++ step < r + 1)
        {
            auto Y = column_block(c, step - 1);

            auto L = matrix_block(step - 1, r);

            Y.init(step);

            L.init(step);

            B -= L * Y;

            B.update(step);
        }

        auto L = matrix_block(r, r);

        auto P = pivots_block(r);

        L.init(step);

        P.init();

        B = P * B;

        B = SolveTriangular(L, UPLO::Lower, DIAG::Unit, SIDE::Left, TRAN::None);

        B.update(step);
    }


    struct MatrixDescriptor : Matrix
    {
        std::string prefix;
        std::string suffix;

        Length col = 0;
        Length row = 0;

        struct
        {
            Length cols = 0;
            Length rows = 0;
        }
        size,
        grid,
        bias;


        using Matrix::operator  =;

        using Matrix::operator -=;


        auto name (Length num) const -> std::string
        {
            auto w = std::to_string(std::min(grid.cols, grid.rows)).length();

            std::stringstream string;

            string << std::setfill('0');

            string << prefix << ".";

            string << "c" << std::setw(w) << col << "."
                   << "r" << std::setw(w) << row << "."
                   << "n" << std::setw(w) << num << ".";

            string << suffix;

            return string.str();
        }


        void init (Length n)
        {
            Matrix::matrix.cols = size.cols;

            Matrix::matrix.rows = size.rows;

            Matrix::matrix.data.resize(Matrix::matrix.cols * Matrix::matrix.rows);

            load(n);
        }

        void update (Length n)
        {
            save(n);

            clear(n - 1);
        }


        bool load (Length n)
        {
            auto data = (char *) Matrix::matrix.data.data();

            auto size = (std::streamsize) (Matrix::matrix.data.size() * sizeof(Scalar));

            return std::ifstream(name(n), std::ios::binary).read(data, size).good();
        }

        bool save (Length n) const
        {
            auto data = (const char *) Matrix::matrix.data.data();

            auto size = (std::streamsize) (Matrix::matrix.data.size() * sizeof(Scalar));

            return std::ofstream(name(n), std::ios::binary).write(data, size).good();
        }

        bool clear (Length n) const
        {
            return std::remove(name(n).c_str());
        }

        bool exist (Length n) const
        {
            auto cols = size.cols;

            auto rows = size.rows;

            auto size = (std::ifstream::pos_type) (cols * rows * sizeof(Scalar));

            std::ifstream file(name(n), std::ios::ate | std::ios::binary);

            return file.is_open() && size == file.tellg();
        }


        void iterate (const std::function <void (Length, Length, Scalar &)> & func)
        {
#pragma omp parallel for default(none) shared(size, bias, func)

            for (Length col = 0; col < size.cols; col ++)
                for (Length row = 0; row < size.rows; row ++)

                    func(bias.cols + col, bias.rows + row, (* this)[col][row]);
        }
    };

    struct PivotsDescriptor : Pivots
    {
        std::string prefix;
        std::string suffix;

        Length grid = 0;
        Length diag = 0;
        Length size = 0;


        using Pivots::operator =;


        auto name () const -> std::string
        {
            auto w = std::to_string(grid).length();

            std::stringstream string;

            string << std::setfill('0');

            string << prefix << ".";

            string << "d" << std::setw(w) << diag << "." << "pivots";

            return string.str();
        }


        void init ()
        {
            Pivots::pivots.size = size;

            Pivots::pivots.data.resize(Pivots::pivots.size);

            load();
        }


        bool load ()
        {
            auto data = (char *) Pivots::pivots.data.data();

            auto size = (std::streamsize) (Pivots::pivots.data.size() * sizeof(Length));

            return std::ifstream(name(), std::ios::binary).read(data, size).good();
        }

        bool save () const
        {
            auto data = (const char *) Pivots::pivots.data.data();

            auto size = (std::streamsize) (Pivots::pivots.data.size() * sizeof(Length));

            return std::ofstream(name(), std::ios::binary).write(data, size).good();
        }

        bool clear () const
        {
            return std::remove(name().c_str());
        }

        bool exist () const
        {
            std::ifstream file(name(), std::ios::ate | std::ios::binary);

            return file.is_open() && size * sizeof(Length) == file.tellg();
        }
    };


    auto matrix_block (Length c, Length r) -> MatrixDescriptor
    {
        MatrixDescriptor descriptor;

        descriptor.prefix = matrix_name() + "b" + std::to_string(blocks_size());
        descriptor.suffix = "matrix";

        descriptor.col = c;
        descriptor.row = r;

        descriptor.size.cols = matrix.grid[c];
        descriptor.size.rows = matrix.grid[r];

        descriptor.grid.cols = (Length) matrix.grid.size();
        descriptor.grid.rows = (Length) matrix.grid.size();

        descriptor.bias.cols = c * blocks_size();
        descriptor.bias.rows = r * blocks_size();

        return descriptor;
    }

    auto column_block (Length c, Length r) -> MatrixDescriptor
    {
        MatrixDescriptor descriptor;

        descriptor.prefix = column_name() + "b" + std::to_string(blocks_size());
        descriptor.suffix = "column";

        descriptor.col = c;
        descriptor.row = r;

        descriptor.size.cols = column.grid[c];
        descriptor.size.rows = matrix.grid[r];

        descriptor.grid.cols = (Length) column.grid.size();
        descriptor.grid.rows = (Length) matrix.grid.size();

        descriptor.bias.cols = c * blocks_size();
        descriptor.bias.rows = r * blocks_size();

        return descriptor;
    }

    auto pivots_block (Length d) -> PivotsDescriptor
    {
        PivotsDescriptor descriptor;

        descriptor.prefix = matrix_name() + "b" + std::to_string(blocks_size());
        descriptor.suffix = "pivots";

        descriptor.diag = d;
        descriptor.size = matrix.grid[d];
        descriptor.grid = matrix.grid.size();

        return descriptor;
    };
};


#endif //RCS__BLOCKSYSTEM_H
