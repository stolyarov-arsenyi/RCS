template <class X>

struct Co
{
    Re <X> r;
    Re <X> i;

    template <class R = X, class I = X>

    Co (const R & r = {}, const I & i = {}) : r(r), i(i) {}

    Co (const Im <X> & i) : Co({}, i) {}

    template <class ... C>

    friend auto operator << (std::basic_ostream <C ...> & stream, const Co & c) -> std::basic_ostream <C ...> &
    {
        std::stringstream string;

        string << std::showpos << std::fixed;

        string << "(" << c.r << ", " << c.i << ")";

        return stream << string.str();
    }
};