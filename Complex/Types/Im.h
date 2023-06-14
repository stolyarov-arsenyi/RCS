template <class X>

struct Im
{
    X i;

    template <class I = X>

    Im (const I & i = {}) : i(i) {}

    operator X () const
    {
        return i;
    }

    template <class ... C>

    friend auto operator << (std::basic_ostream <C ...> & stream, const Im & i) -> std::basic_ostream <C ...> &
    {
        std::stringstream string;

        string << std::showpos << std::fixed;

        string << co(0.0, i);

        return stream << string.str();
    }
};
