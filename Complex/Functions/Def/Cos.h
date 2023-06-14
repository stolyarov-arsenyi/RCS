template <class X>

auto Cos (const Re <X> & x) -> Re <X>
{
    return std::cos(x);
}


template <class X>

auto Cos (const X & x) -> decltype (re(x))
{
    return Cos(re(x));
}