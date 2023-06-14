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





//    static void rename(const std::string & name_old, const std::string & name_new)
//    {
//        std::rename(name_old.c_str(), name_new.c_str());
//    }
//
//    static bool save (const Matrix <Scalar, Length> & matrix, const std::string & block_name)
//    {
//        auto data = (const char *) matrix[0];
//
//        auto size = (std::streamsize) matrix.cols() * (std::streamsize) matrix.rows() * (std::streamsize) sizeof(Scalar);
//
//        return std::ofstream(block_name).write(data, size).good();
//    }
//
//    static bool save (const Pivots <Length> & pivots, const std::string & block_name)
//    {
//        auto data = (const char *) & pivots[0];
//
//        auto size = (std::streamsize) pivots.size() * (std::streamsize) sizeof(Length);
//
//        return std::ofstream(block_name).write(data, size).good();
//    }
//
//    static bool load (Matrix <Scalar, Length> & matrix, const std::string & block_name)
//    {
//        auto data = (char *) matrix[0];
//
//        auto size = (std::streamsize) matrix.cols() * (std::streamsize) matrix.rows() * (std::streamsize) sizeof(Scalar);
//
//        return std::ifstream(block_name).read(data, size).good();
//    }
//
//    static bool load (Pivots <Length> & pivots, const std::string & block_name)
//    {
//        auto data = (char *) & pivots[0];
//
//        auto size = (std::streamsize) pivots.size() * (std::streamsize) sizeof(Length);
//
//        return std::ifstream(block_name).read(data, size).good();
//    }
//
//
//    auto load_matrix ()
//
//
//    template <class ... Int>
//
//    auto block_name(const std::string & suffix, const std::string & prefix, Int ... ints)-> std::string
//    {
//        std::stringstream string;
//
//        string << name << "." << prefix << ".";
//
//        for (auto i : { ints ... })
//
//            string << i << ".";
//
//        string << suffix;
//
//        return string.str();
//    }
//
//    bool exist(const std::string & block_name)
//    {
//        return std::ifstream(block_name).is_open();
//    }


//    template <class  ... Args> void make(Args ...);


    void prepare_system (const std::function <void (Length, Length, Scalar &)> & func)
    {
        Length N = matrix.grid.size.size();

        std::function <void (Length, Length, Length)> A;
        std::function <void (Length, Length        )> L;
        std::function <void (Length, Length        )> U;
        std::function <void (Length                )> P;

        A = [&] (Length c, Length r, Length n)
        {
            if (n == 0)
            {
                Matrix <Scalar, Length> A_c_r(matrix.grid.size[c], matrix.grid.size[r]);

                Length cols = std::accumulate(& matrix.grid.size[0], & matrix.grid.size[c], 0);
                Length rows = std::accumulate(& matrix.grid.size[0], & matrix.grid.size[r], 0);

                for (Length col = 0; col < A_c_r.cols(); col ++)
                for (Length row = 0; row < A_c_r.rows(); row ++)

                    func(cols + col, rows + row, A_c_r[col][row]);

                save(A_c_r, block_name("matrix", "A", c, r, n));
            }
            else
            {
                if (! exist("matrix", "U", c, r - 1))

                    U(c, r - 1);

                if (! exist("matrix", "L", c - 1, r))

                    L(c - 1, r);

                if (! exist("matrix", "A", c, r, n))

                    A(c, r, n);


                //A -= L * U;
            }
        };
        L = [&] (Length c, Length r)
        {
            if (! exist("pivots", "P", c))

                P(c);

            if (! exist("matrix", "A", c, r))

                A(c, r);

            make("L", c, r);
        };
        U = [&] (Length c, Length r)
        {
            if (! exist("pivots", "P", r))

                P(r);

            if (! exist("matrix", "A", c, r))

                A(c, r);

            make("U", c, r);
        };
        P = [&] (Length n)
        {
            if (! exist("matrix", "A", n, n, n))

                A(n, n, n);

            make("LU", n, n, n);
        };

        A(0, 0, 0);
        L(0, 0);
        U(0, 0);

        for (Length n = N - 1; n > 0; n ++)
        {
            if (! exist("pivots", "P", n))

                P(n);
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