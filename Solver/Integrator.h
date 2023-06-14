#ifndef RCS_SOLVER_INTEGRATOR_H
#define RCS_SOLVER_INTEGRATOR_H


template <class X>


class Integrator
{
    static bool singular (const Face <X> & m, const Face <X> & n)
    {
        return & m.v(0) == & n.v(0) || & m.v(0) == & n.v(1) || & m.v(0) == & n.v(2)
            || & m.v(1) == & n.v(0) || & m.v(1) == & n.v(1) || & m.v(1) == & n.v(2)
            || & m.v(2) == & n.v(0) || & m.v(2) == & n.v(1) || & m.v(2) == & n.v(2);
    }

    template <class K>

    static auto efie_norm (const Face <X> & m, const Face <X> & n, const K & k) -> Co <X>
    {
        Co <X> integral;

        auto k_i = 4.0 / (k * k);

        auto m_2_n_2 = m.v(2) - n.v(2);

        for (int o = 0; o < 16 ; o ++)
        {
            for (int i = 0; i < 16; i ++)
            {
                auto R = Norm(m.r[o] - n.r[i] + m_2_n_2);

                auto e = Exp(im(- k * R)) / R;

                integral += e * (Quadrature <X> :: z[o] * Quadrature <X> :: z[i] * ((m.r[o], n.r[i]) - k_i));
            }
        }

        return integral * m.l * n.l / (4.0 * M_PI);
    }

    template <class K>

    static auto efie_sing (const Face <X> & m, const Face <X> & n, const K & k) -> Co <X>
    {
        Co <X> integral;

        auto k_i = 4.0 / (k * k);

        auto m_2_n_2 = m.v(2) - n.v(2);

        for (int o = 0; o < 16; o ++)
        {
            for (int i = 0; i < 16; i ++)
            {
                auto R = Norm(m.r[o] - n.r[i] + m_2_n_2);

                auto e = R < 1e-6 ? co(0.0, - k) : (Exp(im(- k * R)) - 1.0) / R;

                integral += e * (Quadrature <X> :: z[o] * Quadrature <X> :: z[i] * ((m.r[o], n.r[i]) - k_i));
            }

            auto m_r = m.r[o] + m.v(2);

            auto p_m = m_r - (m_r, n.n) * n.n;

            for (int i = 0; i < 3; i ++)
            {
                auto p_0 = n.v(i + 1) - (n.v(i + 1), n.n) * n.n;
                auto p_1 = n.v(i - 1) - (n.v(i - 1), n.n) * n.n;
                auto p_2 = n.v(    2) - (n.v(    2), n.n) * n.n;

                auto I = Unit(p_1 - p_0);

                auto u = I % n.n;

                auto d = Abs((m_r - n.v(i - 1), n.n));

                auto l_0 = ((p_0 - p_m), I);
                auto l_1 = ((p_1 - p_m), I);

                auto P_m = Abs(((p_1 - p_m), u));

                auto P_0 = NormSquared(p_0 - p_m);
                auto P_1 = NormSquared(p_1 - p_m);

                auto R_2 = (P_m * P_m + d * d);

                auto R_0 = Sqrt(P_0 + d * d);
                auto R_1 = Sqrt(P_1 + d * d);

                auto ln = Log((R_1 + l_1) / (R_0 + l_0));

                auto atan_0 = Atan((P_m * l_0) / (R_2 + d * R_0));
                auto atan_1 = Atan((P_m * l_1) / (R_2 + d * R_1));

                auto P = ((p_1 - p_m - l_1 * I), u) / P_m;

                auto E = P_m < 1e-6 ? re(0.0) : P * (P_m * ln - d * (atan_1 - atan_0));

                auto M = (p_m - p_2) * E + u * (((R_2 < 1e-12 ? re(0.0) : R_2 * ln) - R_0 * l_0 + R_1 * l_1) / 2.0);

                integral += (Quadrature <X> :: z[o] * ((m.r[o], M) - E * k_i)) / (2.0 * n.s);
            }
        }

        return integral * m.l * n.l / (4.0 * M_PI);
    }

public:

    static auto integral (const Face <X> & face, const Vector <Re <X>> & e_s, const Re <X> & k) -> Vector <Co <X>>
    {
        auto e_0 = (e_s, face.v(0));
        auto e_1 = (e_s, face.v(1));
        auto e_2 = (e_s, face.v(2));

        auto e_0_1 = e_0 - e_1;
        auto e_1_2 = e_1 - e_2;
        auto e_0_2 = e_0 - e_2;

        if (Abs(e_0_2) < 1e-6)
        {
            if (Abs(e_0_1) < 1e-6)

                return Exp(im(k * e_0)) * face.l / 6.0 * (face.v(0) + face.v(1) - 2.0 * face.v(2));

            auto exp_1 = Exp(im(k * e_1));
            auto exp_2 = Exp(im(k * e_2));

            auto a_k   =   k * e_1_2;
            auto a_k_2 = a_k * a_k;
            auto a_k_3 = a_k * a_k_2;

            auto I_x = (im(         2.0) * exp_1 + co(2.0 * a_k, a_k_2 - 2.0) * exp_2) / (2.0 * a_k_3);
            auto I_y = (co(- a_k, - 2.0) * exp_1 + co(    - a_k,         2.0) * exp_2) / (      a_k_3);

            return face.l * ((face.v(0) - face.v(2)) * I_x + (face.v(1) - face.v(2)) * I_y);
        }

        else

        if (Abs(e_0_1) < 1e-6)
        {
            auto exp_1 = Exp(im(k * e_1));
            auto exp_2 = Exp(im(k * e_2));

            auto a_k   =   k * e_1_2;
            auto a_k_2 = a_k * a_k;
            auto a_k_3 = a_k * a_k_2;

            auto I = (im(- 2.0) * exp_2 + co(2.0 * a_k, 2.0 - a_k_2) * exp_1) / (2.0 * a_k_3);

            return I * face.l * (face.v(0) + face.v(1) - 2.0 * face.v(2));
        }

        else

        if (Abs(e_1_2) < 1e-6)
        {
            auto exp_0 = Exp(im(k * e_0));
            auto exp_2 = Exp(im(k * e_2));

            auto a_k   =   k * e_0_2;
            auto a_k_2 = a_k * a_k;
            auto a_k_3 = a_k * a_k_2;

            auto I_x = (co(- a_k, - 2.0) * exp_0 + co(    - a_k,         2.0) * exp_2) / (      a_k_3);
            auto I_y = (im(         2.0) * exp_0 + co(2.0 * a_k, a_k_2 - 2.0) * exp_2) / (2.0 * a_k_3);

            return face.l * ((face.v(0) - face.v(2)) * I_x + (face.v(1) - face.v(2)) * I_y);
        }

        auto exp_0 = Exp(im(k * e_0));
        auto exp_1 = Exp(im(k * e_1));
        auto exp_2 = Exp(im(k * e_2));

        auto k_3 = k * k * k;

        auto I_x = ( (im(1.0) * exp_1 - co(e_0_1 * k, + 1.0) * exp_0) / (e_0_1 * e_0_1) -
                     (im(1.0) * exp_2 - co(e_0_2 * k, + 1.0) * exp_0) / (e_0_2 * e_0_2) ) / (e_1_2 * k_3);

        auto I_y = ( (im(1.0) * exp_0 + co(e_0_1 * k, - 1.0) * exp_1) / (e_0_1 * e_0_1) -
                     (im(1.0) * exp_2 - co(e_1_2 * k, + 1.0) * exp_1) / (e_1_2 * e_1_2) ) / (e_0_2 * k_3);

        return face.l * ((face.v(0) - face.v(2)) * I_x + (face.v(1) - face.v(2)) * I_y);
    }

    template <class K>

    static auto efie (const Face <X> & m, const Face <X> & n, const K & k) -> Co <X>
    {
        if (singular(m, n))

            return efie_sing(m, n, k);

        return efie_norm(m, n, k);
    }
};


#endif //RCS_SOLVER_INTEGRATOR_H
