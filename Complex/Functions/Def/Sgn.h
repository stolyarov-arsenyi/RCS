template <class X>

auto Sgn (const X & x) -> decltype (re(x))
{
    return re(x) < 0.0 ? - 1.0 : + 1.0;
}