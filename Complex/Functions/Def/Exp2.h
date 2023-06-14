template <class X>

auto Exp2 (const Re <X> & x) -> Re <X>
{
    return std::exp2(x);
}


template <class X>

auto Exp2 (const X & x) -> decltype (re(x))
{
    return Exp2(re(x));
}