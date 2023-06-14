template <class X>

auto Hypot (const Re <X> & x, const Re <X> & y) -> Re <X>
{
    return std::hypot(x, y);
}


template <class X, class Y, class R = decltype (re(X()) + re(Y()))>

auto Hypot (const X & x, const Y & y) -> R
{
    return Hypot(R(x), R(y));
}