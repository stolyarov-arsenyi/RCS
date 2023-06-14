template <class X>

auto Expm1 (const Re <X> & x) -> Re <X>
{
    return std::expm1(x);
}


template <class X>

auto Expm1 (const X & x) -> decltype (re(x))
{
    return Expm1(re(x));
}