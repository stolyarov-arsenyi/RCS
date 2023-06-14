#ifndef  RCS__PIVOTS_H
#define  RCS__PIVOTS_H


template <class Length>

class Pivots
{
    struct
    {
        std::vector <Length> data;

        Length size;
    }
    pivots;

public:

    Pivots () = default;

    Pivots (Length size) : pivots { {}, size }
    {
        pivots.data.resize(size);
    }


    auto size () const -> Length
    {
        return pivots.size;
    }

    auto operator [] (Length idx) const -> const Length &
    {
        return pivots.data[idx];
    }

    auto operator [] (Length idx) -> Length &
    {
        return pivots.data[idx];
    }


    friend std::ostream & operator << (std::ostream & stream, const Pivots & pivots)
    {
        std::stringstream string;

        string << std::fixed << std::showpos;

        for (Length idx = 0; idx < pivots.size(); idx ++)

            string << pivots[idx] << (idx == pivots.size() - 1 ? "\n" : ",");

        return stream << string.str();
    }
};


#endif //RCS__PIVOTS_H
