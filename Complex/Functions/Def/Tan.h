template <class X>

auto Tan (const Re <X> & x) -> Re <X>
{
    return std::tan(x);
}


template <class X>

auto Tan (const X & x) -> decltype (re(x))
{
    return Tan(re(x));
}