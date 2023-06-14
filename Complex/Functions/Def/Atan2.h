template <class X>

auto Atan2 (const Re <X> & y, const Re <X> & x) -> Re <X>
{
    return std::atan2(y, x);
}


template <class X, class Y, class R = decltype (re(X() + Y()))>

auto Atan2 (const X & z, const Y & y) -> R
{
    return Atan2(R(z), R(y));
}