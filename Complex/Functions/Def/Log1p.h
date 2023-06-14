template <class X>

auto Log1p (const Re <X> & x) -> Re <X>
{
    return std::log1p(x);
}


template <class X>

auto Log1p (const X & x) -> decltype (re(x))
{
    return Log1p(re(x));
}