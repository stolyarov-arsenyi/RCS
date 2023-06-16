template <class R>

Vector <R> operator + (const Vector <R> & r)
{
    return { + r.x, + r.y, + r.z };
}

template <class R>

Vector <R> operator - (const Vector <R> & r)
{
    return { - r.x, - r.y, - r.z };
}


template <class L, class R>

auto operator == (const Vector <L> & l, const Vector <R> & r) -> decltype (l.x == r.x)
{
    return l.x == r.x && l.y == r.y && l.z == r.z;
}

template <class L, class R>

auto operator + (const Vector <L> & l, const Vector <R> & r) -> Vector <decltype (l.x + r.x)>
{
    return { l.x + r.x, l.y + r.y, l.z + r.z };
}

template <class L, class R>

auto operator - (const Vector <L> & l, const Vector <R> & r) -> Vector <decltype (l.x - r.x)>
{
    return { l.x - r.x, l.y - r.y, l.z - r.z };
}

template <class L, class R>

auto operator , (const Vector <L> & l, const Vector <R> & r) -> decltype (l.x * ~ r.x)
{
    return { l.x * ~ r.x + l.y * ~ r.y + l.z * ~ r.z };
}

template <class L, class R>

auto operator % (const Vector <L> & l, const Vector <R> & r) -> Vector <decltype (l.x * r.x)>
{
    return { l.y * r.z - l.z * r.y, l.z * r.x - l.x * r.z, l.x * r.y - l.y * r.x };
}

template <class L, class R>

auto operator * (const L & l, const Vector <R> & r) -> Vector <decltype (r.x * l)>
{
    return { r.x * l, r.y * l, r.z * l };
}

template <class L, class R>

auto operator * (const Vector <L> & l, const R & r) -> Vector <decltype (l.x * r)>
{
    return { l.x * r, l.y * r, l.z * r };
}

template <class L, class R>

auto operator / (const Vector <L> & l, const R & r) -> Vector <decltype (l.x / r)>
{
    return { l.x / r, l.y / r, l.z / r };
}