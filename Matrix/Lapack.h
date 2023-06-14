#ifndef  RCS__LAPACK_H
#define  RCS__LAPACK_H


template <class Scalar, class Length>

struct Lapack
{
    template <class ... Args>

    static void gemm (Args ... args);


    template <class ... Args>

    static void laswp (Args ... args);


    template <class ... Args>

    static void getrf (Args ... args);


    template <class ... Args>

    static void trsm (Args ... args);
};


#endif //RCS__LAPACK_H
