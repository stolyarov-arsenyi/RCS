template <class X>

auto Log2 (const Re <X> & x) -> Re <X>
{
    return std::log2(x);
}


template <class X>

auto Log2 (const X & x) -> decltype (re(x))
{
    return Log2(re(x));
}