template <class X>

struct Edge
{
    const Face <X> face [2];

    Edge (const Face <X> & f_0, const Face <X> & f_1) : face { f_0, f_1 } {}
};