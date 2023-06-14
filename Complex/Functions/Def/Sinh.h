template <class X>

auto Sinh (const Re <X> & x) -> Re <X>
{
    return std::sinh(x);
}


template <class X>

auto Sinh (const X & x) -> decltype (re(x))
{
    return Sinh(re(x));
}