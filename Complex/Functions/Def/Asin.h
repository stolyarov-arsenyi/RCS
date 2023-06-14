template <class X>

auto Asin (const Re <X> & x) -> Re <X>
{
    return std::asin(x);
}


template <class X>

auto Asin (const X & x) -> decltype (re(x))
{
    return Asin(re(x));
}