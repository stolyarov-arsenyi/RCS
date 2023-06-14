template <class X>

auto Atan (const Re <X> & x) -> Re <X>
{
    return std::atan(x);
}


template <class X>

auto Atan (const X & x) -> decltype (re(x))
{
    return Atan(re(x));
}