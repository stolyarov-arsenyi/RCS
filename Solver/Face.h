template <class X>

struct Face
{
    std::shared_ptr <Vector <Re <X>>> vertices [3];

    Vector <Re <X>> r [16];
    Vector <Re <X>> n;

    Re <X> l;
    Re <X> s;

    Face (std::shared_ptr <Vector <Re <X>>> v_0,
          std::shared_ptr <Vector <Re <X>>> v_1,
          std::shared_ptr <Vector <Re <X>>> v_2) : vertices { v_0, v_1, v_2 }
    {
        n = Unit((v(1) - v(0)) % (v(2) - v(0)));
        s = Norm((v(1) - v(0)) % (v(2) - v(0)) / 2.0);
        l = Norm((v(1) - v(0)));

        auto v_0_2 = v(0) - v(2);
        auto v_1_2 = v(1) - v(2);

        for (int q = 0; q < 16; q ++)

            r[q] = Quadrature <X> :: x[q] * v_0_2 + Quadrature <X> :: y[q] * v_1_2;
    }

    auto v (int i) const -> const Vector <Re <X>> &
    {
        return * vertices[(3 + i % 3) % 3];
    }
};

