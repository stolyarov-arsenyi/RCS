#ifndef RCS_SOLVER_SOLUTION_H
#define RCS_SOLVER_SOLUTION_H


template <class X>

struct Antenna
{
    Re <X> alt;
    Re <X> azi;
    Re <X> pol;

    Vector <Re <X>> v_polar;
    Vector <Re <X>> v_front;

    Antenna (const Re <X> & alt, const Re <X> & azi, const Re <X> & pol) : alt(alt), azi(azi), pol(pol)
    {
        auto radians = M_PI / 180.0;

        v_front = { 0.0, 0.0, 1.0 };
        v_polar = { 0.0, 1.0, 0.0 };

        v_polar = Rotate(v_polar, { 0.0,           0.0, pol * radians });
        v_polar = Rotate(v_polar, { 0.0, alt * radians, azi * radians });
        v_front = Rotate(v_front, { 0.0, alt * radians, azi * radians });
    }
};


template <class X>

struct Field
{
    Vector <Co <X>> E_sca;

    Re <X> rcs_db;
    Re <X> rcs_sm;
};


template <class X>

struct Solution
{
    Antenna <X> inc;
    Antenna <X> sca;

    Field <X> field;
};


#endif //RCS_SOLVER_SOLUTION_H
