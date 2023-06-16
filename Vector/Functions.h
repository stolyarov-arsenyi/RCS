template <class X, class Y, class Z>

auto vector (const X & x, const Y & y, const Z & z) -> Vector <decltype (x + y + z)>
{
    return { x, y, z };
}


template <class X>

auto Norm (const Vector <X> & v) -> decltype (AbsSquared(v.x))
{
    return Sqrt(AbsSquared(v.x) + AbsSquared(v.y) + AbsSquared(v.z));
}

template <class X>

auto NormSquared (const Vector <X> & v) -> decltype (AbsSquared(v.x))
{
    return AbsSquared(v.x) + AbsSquared(v.y) + AbsSquared(v.z);
}


template <class X>

auto Unit (const Vector <X> & v) -> Vector <X>
{
    return v / Norm(v);
}


template <class X>

auto Rotate (Vector <Re <X>> vec, const Vector <Re <X>> & rot) -> Vector <Re <X>>
{
    Re <X> cos_x = Cos(rot.x);
    Re <X> sin_x = Sin(rot.x);

    Re <X> cos_y = Cos(rot.y);
    Re <X> sin_y = Sin(rot.y);

    Re <X> cos_z = Cos(rot.z);
    Re <X> sin_z = Sin(rot.z);

    vec = { vec.x, vec.y * cos_x - vec.z * sin_x, vec.y * sin_x + vec.z * cos_x };
    vec = { vec.x * cos_y + vec.z * sin_y, vec.y, vec.z * cos_y - vec.x * sin_y };
    vec = { vec.x * cos_z - vec.y * sin_z, vec.x * sin_z + vec.y * cos_z, vec.z };

    return vec;
}