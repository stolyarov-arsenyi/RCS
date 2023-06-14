template <class X>

auto Tanh (const Re <X> & x) -> Re <X>
{
    return std::tanh(x);
}


template <class X>

auto Tanh (const X & x) -> decltype (re(x))
{
    return Tanh(re(x));
}