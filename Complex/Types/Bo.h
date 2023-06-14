template <class X>

struct Bo
{
    bool b;

    Bo (const bool & b) : b(b) {}

    operator bool () const
    {
        return b;
    }
};