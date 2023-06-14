template <class X>

auto Acos (const Re <X> & x) -> Re <X>
{
    return std::acos(x);
}


template <class X>

auto Acos (const X & x) -> decltype (re(x))
{
    return Acos(re(x));
}