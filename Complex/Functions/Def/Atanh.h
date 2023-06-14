template <class X>

auto Atanh (const Re <X> & x) -> Re <X>
{
    return std::atanh(x);
}


template <class X>

auto Atanh (const X & x) -> decltype (re(x))
{
    return Atanh(re(x));
}