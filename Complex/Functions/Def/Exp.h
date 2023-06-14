template <class X>

auto Exp (const Re <X> & x) -> Re <X>
{
    return std::exp(x);
}


template <class X>

auto Exp (const Im <X> & x) -> Co <X>
{
    return { Cos(re(x)), Sin(re(x)) };
}


template <class X>

auto Exp (const Co <X> & x) -> Co <X>
{
    return Exp(re(x.r)) * Exp(im(x.i));
}


template <class X>

auto Exp (const X & x) -> decltype (re(x))
{
    return Exp(re(x));
}