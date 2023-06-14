template <class X>

auto Cbrt (const Re <X> & x) -> Re <X>
{
    return std::cbrt(x);
}


template <class X>

auto Cbrt (const X & x) -> decltype (re(x))
{
    return Cbrt(re(x));
}