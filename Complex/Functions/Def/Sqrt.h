template <class X>

auto Sqrt (const Re <X> & x) -> Re <X>
{
    return std::sqrt(x);
}


template <class X>

auto Sqrt (const Co <X> & x) -> Co <X>
{
    Re <X> h = Hypot(x.r, x.i);

    return { Sqrt((h + x.r) / 2.0), Sgn(x.i) * Sqrt((h - x.r) / 2.0) };
}


template <class X>

auto Sqrt (const Im <X> & x) -> Im <X>
{
    return Sgn(x) * Sqrt(Abs(x) / 2.0);
}


template <class X>

auto Sqrt (const X & x) -> decltype (re(x))
{
    return Sqrt(re(x));
}