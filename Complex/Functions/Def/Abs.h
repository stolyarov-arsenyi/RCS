template <class X>

auto Abs (const Re <X> & x) -> Re <X>
{
    return (X) x < 0.0 ? - x : + x;
}


template <class X>

auto Abs (const Co <X> & x) -> Re <X>
{
    return Sqrt(x.r * x.r + x.i * x.i);
}


template <class X>

auto Abs (const X & x) -> decltype (re(x))
{
    return Abs(re(x));
}
