template <class X>

auto Log (const Re <X> & x) -> Re <X>
{
    return std::log(x);
}


template <class X>

auto Log (const X & x) -> decltype (re(x))
{
    return Logger(re(x));
}