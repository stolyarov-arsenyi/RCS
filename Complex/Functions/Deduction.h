template <class X>

auto re (const X & x) -> Re <X>
{
    return x;
}

template <class X>

auto re (const Re <X> & x) -> Re <X>
{
    return x;
}

template <class X>

auto re (const Im <X> & x) -> Re <X>
{
    return x;
}


template <class X>

auto im (const X & x) -> Im <X>
{
    return x;
}

template <class X>

auto im (const Re <X> & x) -> Im <X>
{
    return x;
}

template <class X>

auto im (const Im <X> & x) -> Im <X>
{
    return x;
}



template <class R, class I = R>

auto co (const R & r, const I & i = {}) -> decltype (re(r) + im(i))
{
    return re(r) + im(i);
}