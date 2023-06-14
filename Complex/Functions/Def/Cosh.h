template <class X>

auto Cosh (const Re <X> & x) -> Re <X>
{
    return std::cosh(x);
}


template <class X>

auto Cosh (const X & x) -> decltype (re(x))
{
    return Cosh(re(x));
}