#ifndef  RCS__BLOCKMATRIX_H
#define  RCS__BLOCKMATRIX_H


#include "Matrix/Lapack.h"
#include "Matrix/Pivots.h"
#include "Matrix/Matrix.h"


template <class Scalar, class Length>

struct BlockMatrix
{
    struct Block
    {
        static Length size;

        static bool rename (const std::string & name_old, const std::string & name_new)
        {
            return std::rename(name_old.c_str(), name_new.c_str());
        }

        static bool remove (const std::string & name)
        {
            return std::remove(name.c_str());
        }

        static bool exist (const std::string & name)
        {
            return std::ifstream(name).is_open();
        }

        static auto suffix_matrix (Length col, Length row) -> std::string
        {
            return "." + std::to_string(col) + "." + std::to_string(row) + ".matrix";
        }

        static auto suffix_pivots (Length col) -> std::string
        {
            return "." + std::to_string(col) + ".pivots";
        }
    };


    std::uintmax_t addr = (std::uintmax_t) this;

    std::string name;

    const Length cols; // Total matrix columns
    const Length rows; // Total matrix rows


    struct
    {
        std::vector <Length> cols; // Block columns
        std::vector <Length> rows; // Block rows
    }
    grid;


    BlockMatrix (Length cols, Length rows, std::string name = "") : cols(cols), rows(rows), name(std::move(name))
    {
        if (Block::size <= 0)
        {
            std::stringstream string;

            string << "Invalid block size " << Block::size;

            throw std::runtime_error(string.str());
        }

        /* Divide grid into blocks of size Block::size */

        grid.cols.resize(cols / Block::size, Block::size);
        grid.rows.resize(rows / Block::size, Block::size);

        /* Add column and row with odd size */

        if (cols %= Block::size) grid.cols.push_back(cols);
        if (rows %= Block::size) grid.rows.push_back(rows);
    }

  ~ BlockMatrix ()
    {
        /* If name doesn't set, remove from disk */

        if (name.empty())
        {
            for (Length col = 0; col < grid.cols.size(); col ++)
            for (Length row = 0; row < grid.rows.size(); row ++)

                Block::remove(prefix() + Block::suffix_matrix(col, row));

            for (Length col = 0; col < grid.cols.size(); col ++)

                Block::remove(prefix() + Block::suffix_pivots(col));
        }
    }


    void for_each_entry (const std::function <void (Length, Length, const Scalar &)> & func) const
    {
        for (Length grid_col = 0; grid_col < grid.cols.size(); grid_col ++)
        for (Length grid_row = 0; grid_row < grid.rows.size(); grid_row ++)
        {
            auto block = load_matrix(grid_col, grid_row);

            auto col = grid_col * Block::size;
            auto row = grid_row * Block::size;

            #pragma omp parallel for default(none) shared(col, row, block, func)

            for (Length c = 0; c < block.cols(); c ++)
            for (Length r = 0; r < block.rows(); r ++)

                func(col + c, row + r, block[c][r]);
        }
    }

    void for_each_entry (const std::function <void (Length, Length, Scalar &)> & func)
    {
        for (Length grid_col = 0; grid_col < grid.cols.size(); grid_col ++)
        for (Length grid_row = 0; grid_row < grid.rows.size(); grid_row ++)
        {
            auto block = load_matrix(grid_col, grid_row);

            auto col = grid_col * Block::size;
            auto row = grid_row * Block::size;

            #pragma omp parallel for default(none) shared(col, row, block, func)

            for (Length c = 0; c < block.cols(); c ++)
            for (Length r = 0; r < block.rows(); r ++)

                func(col + c, row + r, block[c][r]);

            save_matrix(block, grid_col, grid_row);
        }
    }


    void lu_factorization ()
    {
        /** |A00 A01|   |L00  0 |   |U00 U01| **/
        /** |A10 A11| = |L10 L11| x | 0  U11| **/

        /** 1) A00  = L00 * U00                    **/
        /** 2) A10  = L10 * U00                    **/
        /** 3) A01  = L00 * U01                    **/
        /** 4) A11' = L11 * U11 = A11 - L10 * U01  **/
        /** 5) Repeat steps 1-4 recursive for A11' **/

        Length N = std::min(grid.cols.size(), grid.rows.size());

        for (Length n = 0; n < N; n ++)
        {
            /** A00 = P00 * L00 * U00 **/

            auto A00 = load_matrix(n, n);

            auto P00 = A00.lu_factorization();

            save_matrix(A00, n, n);
            save_pivots(P00, n);

            /** Solve P00 * L00 * U01 = A01 **/

            for (Length c = n + 1; c < N; c ++)
            {
                auto U01 = load_matrix(c, n);

                U01 = P00 * U01;

                U01 = SolveTriangular(A00, Lower, Unit, Left, None);

                save_matrix(U01, c, n);
            }

            /** Solve L10 * U00 = A10 **/

            for (Length r = n + 1; r < N; r ++)
            {
                auto L10 = load_matrix(n, r);

                L10 = SolveTriangular(A00, Upper, NonUnit, Right, None);

                save_matrix(L10, n, r);
            }

            /** Free unused memory **/

            A00 = {};
            P00 = {};

            /** A11' = A11 - L10 * U01 **/

            for (Length c = n + 1; c < N; c ++)
            for (Length r = n + 1; r < N; r ++)
            {
                auto A11 = load_matrix(c, r);
                auto L10 = load_matrix(n, r);
                auto U01 = load_matrix(c, n);

                A11 -= L10 * U01;

                save_matrix(A11, c, r);
            }
        }
    }

    void solve_for (BlockMatrix & column)
    {
        /** Example of solving L * U * x = b, where L * U is factorization of 3x3 block matrix **/

        /** Step 1: Solve L * y = b **/

        /** Pxx - permutation matrix for Lxx coming from PLU decomposition, such that P * L * U = A **/

        /** |L11        |   |Y1|   |P11 * B1| **/
        /** |L21 L22    | x |Y2| = |P22 * B2| **/
        /** |L31 L32 L33|   |Y3|   |P33 * B3| **/

        /** Solve L11 * Y1 = B1                       **/
        /** Solve L22 * Y2 = B2 - L21 * Y1            **/
        /** Solve L33 * Y3 = B3 - L31 * Y1 - L32 * Y2 **/

        /** Step 2: Solve U * x = y **/

        /** |U11 U12 U13|   |X1|   |Y1| **/
        /** |    U22 U23| x |X2| = |Y2| **/
        /** |        U33|   |X3|   |Y3| **/

        /** Solve U33 * Y3 = B3                       **/
        /** Solve U22 * Y2 = B2 - L23 * Y3            **/
        /** Solve U11 * Y1 = B1 - L13 * Y3 - L12 * Y2 **/


        for (Length col = 0; col < column.grid.cols.size(); col ++)
        {
            /** Solve L * y = b **/

            for (Length row = 0; row < column.grid.rows.size(); row ++)
            {
                auto B = column.load_matrix(col, row);

                for (Length k = 0; k < row; k ++)
                {
                    auto L =        load_matrix(k, row);
                    auto Y = column.load_matrix(col, k);

                    B -= L * Y;
                }

                auto L = load_matrix(row, row);
                auto P = load_pivots(row);

                B = P * B;

                B = SolveTriangular(L, Lower, Unit, Left, None);

                column.save_matrix(B, col, row);
            }

            /** Solve U * x = y **/

            for (Length row = column.grid.rows.size() - 1; row >= 0; row --)
            {
                auto Y = column.load_matrix(col, row);

                for (Length k = row + 1; k < grid.rows.size(); k ++)
                {
                    auto U =        load_matrix(k, row);
                    auto X = column.load_matrix(col, k);

                    Y -= U * X;
                }

                auto U = load_matrix(row, row);
                auto P = load_pivots(row);

                Y = SolveTriangular(U, Upper, NonUnit, Left, None);

                column.save_matrix(Y, col, row);
            }
        }
    }


    auto load_matrix (Length col, Length row) const -> Matrix <Scalar, Length>
    {
        Matrix <Scalar, Length> matrix(grid.cols[col], grid.rows[row]);

        auto data = (char *) matrix[0];

        auto size = (std::streamsize) (matrix.cols() * matrix.rows() * sizeof(Scalar));

        auto name = prefix() + Block::suffix_matrix(col, row);

        if (Block::exist(name))

            std::ifstream(name, std::ios::binary).read(data, size);

        return matrix;
    }

    auto load_pivots (Length col) -> Pivots <Length>
    {
        Pivots <Length> pivots(grid.cols[col]);

        auto data = (char *) & pivots[0];

        auto size = (std::streamsize) (pivots.size() * sizeof(Length));

        auto name = prefix() + Block::suffix_pivots(col);

        if (Block::exist(name))

            std::ifstream(name, std::ios::binary).read(data, size);

        return pivots;
    }


    void save_matrix (const Matrix <Scalar, Length> & matrix, Length col, Length row)
    {
        auto data = (const char *) matrix[0];

        auto size = (std::streamsize) (matrix.cols() * (std::streamsize) matrix.rows() * sizeof(Scalar));

        auto name = prefix() + Block::suffix_matrix(col, row);

        std::ofstream(name, std::ios::binary).write(data, size);
    }

    void save_pivots (const Pivots <Length> & pivots, Length col)
    {
        auto data = (const char *) & pivots[0];

        auto size = (std::streamsize) (pivots.size() * sizeof(Length));

        auto name = prefix() + Block::suffix_pivots(col);

        std::ofstream(name, std::ios::binary).write(data, size);
    }


    auto prefix () const -> std::string
    {
        std::stringstream stream;

        if (name.empty())

            stream << std::hex << addr << std::dec;

        else

            stream << name;

        stream << "." << Block::size;

        return stream.str();
    }


    friend std::ostream & operator << (std::ostream & stream, const BlockMatrix & block_matrix)
    {
        std::vector <Scalar> array(block_matrix.cols * block_matrix.rows);

        block_matrix.for_each_entry([&] (Length col, Length row, const Scalar & val)
        {
            array[col * block_matrix.rows + row] = val;
        });

        std::stringstream string;

        string << std::showpos << std::fixed;

        for (Length row = 0; row < block_matrix.rows; row ++)
        for (Length col = 0; col < block_matrix.cols; col ++)

            string << array[col * block_matrix.rows + row] << (col == block_matrix.cols - 1 ? "\n" : ",");

        return stream << string.str();
    }
};


template <class Scalar, class Length>

Length BlockMatrix <Scalar, Length> :: Block :: size;


#endif //RCS__BLOCKMATRIX_H