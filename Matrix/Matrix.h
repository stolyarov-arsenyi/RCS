#ifndef  RCS__MATRIX_H
#define  RCS__MATRIX_H


enum SIDE : char
{
    Left  = 'L',
    Right = 'R',
};

enum UPLO : char
{
    Upper = 'U',
    Lower = 'L',
};

enum TRAN : char
{
    None = 'N',
    Tran = 'T',
    Conj = 'C',
};

enum DIAG : char
{
    Unit    = 'U',
    NonUnit = 'N',
};


template <class Scalar, class Length>

class Matrix
{
    struct
    {
        std::vector <Scalar> data;

        Length cols;
        Length rows;
    }
    matrix;

public:

    Matrix () = default;

    Matrix (Length cols, Length rows) : matrix { {}, cols, rows }
    {
        matrix.data.resize(cols * rows);
    }


    auto cols () const -> Length
    {
        return matrix.cols;
    }

    auto rows () const -> Length
    {
        return matrix.rows;
    }


    auto operator [] (Length col) const -> const Scalar *
    {
        return matrix.data.data() + col * rows();
    }

    auto operator [] (Length col) -> Scalar *
    {
        return matrix.data.data() + col * rows();
    }

private:

    struct MulMat
    {
        const Matrix & mat_a;
        const Matrix & mat_b;
    };

    struct PivMat
    {
        const Pivots <Length> & piv;
        const Matrix          & mat;
    };

    struct SolTri
    {
        const Matrix & mat;

        UPLO uplo;
        DIAG diag;
        SIDE side;
        TRAN tran;
    };

public:

    Matrix & operator -= (const MulMat & mul_mat)
    {
        const Matrix & mat_a = mul_mat.mat_a;
        const Matrix & mat_b = mul_mat.mat_b;

        char trans_a = 'N';
        char trans_b = 'N';

        Length m = mat_a.rows();
        Length n = mat_b.cols();
        Length k = mat_a.cols();

        Length lda = m;
        Length ldb = k;
        Length ldc = m;

        Scalar alpha = - 1.0;
        Scalar  beta = + 1.0;

        const Scalar * a = mat_a[0];
        const Scalar * b = mat_b[0];

        Scalar * c = (* this)[0];

        Lapack <Scalar, Length> :: gemm(& trans_a, & trans_b, & m, & n, & k, & alpha, a, & lda, b, & ldb, & beta, c, & ldc);

        return * this;
    }

    Matrix & operator  = (const PivMat & piv_mat)
    {
        if (this != & piv_mat.mat)

            * this = piv_mat.mat;

        Length n = cols();

        Length lda = rows();

        Length k1 = 1;

        Length k2 = rows();

        Length incx = 1;

        Scalar * a = (* this)[0];

        const Length * ipiv = & piv_mat.piv[0];

        Lapack <Scalar, Length> :: laswp(& n, a, & lda, & k1, & k2, ipiv, & incx);

        return * this;
    }

    Matrix & operator  = (const SolTri & sol_tri)
    {
        char side   = (char) sol_tri.side;
        char uplo   = (char) sol_tri.uplo;
        char transa = (char) sol_tri.tran;
        char diag   = (char) sol_tri.diag;

        Length m = rows();
        Length n = cols();

        Length lda = sol_tri.side == Left ? m : n;
        Length ldb = m;

        Scalar alpha = 1.0;

        const Scalar * a = sol_tri.mat[0];

        Scalar * b = (* this)[0];

        Lapack <Scalar, Length> :: trsm(& side, & uplo, & transa, & diag, & m, & n, & alpha, a, & lda, b, & ldb);

        return * this;
    }


    auto lu_factorization () -> Pivots <Length>
    {
        Length info;

        Length m = rows();
        Length n = cols();

        Length lda = m;

        Scalar * a = (* this)[0];

        Pivots <Length> pivots(std::min(m, n));

        Lapack <Scalar, Length> :: getrf(& m, & n, a, & lda, & pivots[0], & info);

        if (info)
        {
            std::stringstream string;

            string << "Invalid matrix factorization. INFO: " << info;

            throw std::invalid_argument(string.str());
        }

        return pivots;
    }


    friend auto SolveTriangular (const Matrix & matrix, UPLO uplo, DIAG diag, SIDE side, TRAN tran) -> SolTri
    {
        return { matrix, uplo, diag, side, tran };
    }


    friend auto operator * (const Pivots <Length> & pivots, const Matrix & matrix) -> PivMat
    {
        return { pivots, matrix };
    }

    friend auto operator * (const Matrix & matrix_a, const Matrix & matrix_b) -> MulMat
    {
        return { matrix_a, matrix_b };
    }

    friend std::ostream & operator << (std::ostream & stream, const Matrix & matrix)
    {
        std::stringstream string;

        string << std::fixed << std::showpos;

        for (Length row = 0; row < matrix.rows(); row ++)
            for (Length col = 0; col < matrix.cols(); col ++)

                string << matrix[col][row] << (col == matrix.cols() - 1 ? "\n" : ",");

        return stream << string.str();
    }
};


#endif //RCS__MATRIX_H
