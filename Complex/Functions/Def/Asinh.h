template <class X>

auto Asinh (const Re <X> & x) -> Re <X>
{
    return std::asinh(x);
}


template <class X>

auto Asinh (const X & x) -> decltype (re(x))
{
    return Asinh(re(x));
}