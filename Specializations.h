#ifndef BLOCKLU__SPECIALIZATIONS_H
#define BLOCKLU__SPECIALIZATIONS_H


template <>

template <class ... Args>

void Lapack <double, MKL_INT> :: gemm (Args ... args)
{
    dgemm(args ...);
}


template <>

template <class ... Args>

void Lapack <double, MKL_INT> :: laswp (Args ... args)
{
    dlaswp(args ...);
}


template <>

template <class ... Args>

void Lapack <double, MKL_INT> :: getrf (Args ... args)
{
    dgetrf(args ...);
}


template <>

template <class ... Args>

void Lapack <double, MKL_INT> :: trsm (Args ... args)
{
    dtrsm(args ...);
}


template <>

template <class ... Args>

void Lapack <MKL_Complex16, MKL_INT> :: gemm (Args ... args)
{
    zgemm(args ...);
}


template <>

template <class ... Args>

void Lapack <MKL_Complex16, MKL_INT> :: laswp (Args ... args)
{
    zlaswp(args ...);
}


template <>

template <class ... Args>

void Lapack <MKL_Complex16, MKL_INT> :: getrf (Args ... args)
{
    zgetrf(args ...);
}


template <>

template <class ... Args>

void Lapack <MKL_Complex16, MKL_INT> :: trsm (Args ... args)
{
    ztrsm(args ...);
}


#endif //BLOCKLU__SPECIALIZATIONS_H
