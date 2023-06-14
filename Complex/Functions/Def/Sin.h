template <class X>

auto Sin (const Re <X> & x) -> Re <X>
{
    return std::sin(x);
}


template <class X>

auto Sin (const X & x) -> decltype (re(x))
{
    return Sin(re(x));
}