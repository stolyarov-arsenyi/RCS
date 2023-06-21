#ifndef RCS_SOLVER_LOGGER_H
#define RCS_SOLVER_LOGGER_H


struct Logger : std::ofstream
{
    using std::ofstream::ofstream;

    template <class X>

    Logger & operator << (X x)
    {
        ((std::ofstream &) * this)  << x;

        std::cout << x;

        return * this;
    }

    Logger & operator << (std::ostream & (* func) (std::ostream &))
    {
        func((std::ofstream &) * this);

        func(std::cout);

        return * this;
    }

    Logger & date ()
    {
        char str [32];

        auto time = std::time(nullptr);

        strftime(str, 32, "%a %b %d %H:%M:%S %Y", localtime(& time));

        return (* this) << str;
    }
};


#endif //RCS_SOLVER_LOGGER_H
