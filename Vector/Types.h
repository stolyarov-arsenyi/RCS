template <class X>

struct Vector
{
    X x, y, z;

    Vector () = default;

    Vector (const X & x, const X & y, const X & z) : x(x), y(y), z(z) {}

    template <class ... C>

    friend auto operator << (std::basic_ostream <C ...> & stream, const Vector & v) -> std::basic_ostream <C ...> &
    {
        std::stringstream string;

        string << std::showpos << std::fixed;

        string << "[" << v.x << ", " << v.y << ", " << v.z << "]";

        return stream << string.str();
    }
};