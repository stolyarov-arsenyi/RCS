template <class X>

struct Re
{
    X r;

    template <class R = X>

    Re (const R & r = {}) : r(r) {}

    operator X () const
    {
        return r;
    }

    template <class ... C>

    friend auto operator << (std::basic_ostream <C ...> & stream, const Re & r) -> std::basic_ostream <C ...> &
    {
        std::stringstream string;

        string << std::showpos << std::fixed;

        string << r.r;

        return stream << string.str();
    }
};
