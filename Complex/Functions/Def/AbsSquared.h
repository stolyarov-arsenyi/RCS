template <class X>

auto AbsSquared (const Re <X> & x) -> Re <X>
{
    return x * x;
}


template <class X>

auto AbsSquared (const Co <X> & x) -> Re <X>
{
    return x.r * x.r + x.i * x.i;
}


template <class X>

auto AbsSquared (const X & x) -> decltype (re(x))
{
    return AbsSquared(re(x));
}