template <class X>

auto Acosh (const Re <X> & x) -> Re <X>
{
    return std::acosh(x);
}


template <class X>

auto Acosh (const X & x) -> decltype (re(x))
{
    return Acosh(re(x));
}