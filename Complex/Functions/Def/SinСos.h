template <class X>

auto SinCos (const Re <X> & x, Re <X> & sin, Re <X> & cos) -> void
{
    std::sincos(x, sin, cos);
}


template <class A, class X>

auto SinCos (const A & arg, Re <X> & sin, Re <X> & cos) -> void
{
    SinCos <X> (arg, sin, cos);
}