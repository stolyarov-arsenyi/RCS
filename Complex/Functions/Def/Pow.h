template <class X, class Y>

auto Pow (const Re <X> & x, const Y & y) -> Re <X>
{
    return std::pow((X) x, y);
}


template <class X, class Y>

auto Pow (const X & x, const Y & y) -> decltype (std::pow(x, y))
{
    return std::pow(x, y);
}