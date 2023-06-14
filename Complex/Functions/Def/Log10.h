template <class X>

auto Log10 (const Re <X> & x) -> Re <X>
{
    return std::log10(x);
}


template <class X>

auto Log10 (const X & x) -> decltype (re(x))
{
    return Log10(re(x));
}