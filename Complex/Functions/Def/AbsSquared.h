template <class X>

auto AbsSquare (const Re <X> & x) -> Re <X>
{
    return x * x;
}


template <class X>

auto AbsSquare (const Co <X> & x) -> Re <X>
{
    return x.r * x.r + x.i * x.i;
}


template <class X>

auto AbsSquare (const X & x) -> decltype (re(x))
{
    return AbsSquare(re(x));
}