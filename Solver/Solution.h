#ifndef RCS_SOLVER_SOLUTION_H
#define RCS_SOLVER_SOLUTION_H


template <class X>

struct Source
{
    Re <X> alt;
    Re <X> azi;
    Re <X> pol;
};


template <class X>

struct Field
{
    Vector <Re <X>> front;
    Vector <Re <X>> E_inc;
    Vector <Co <X>> E_sca;

    Re <X> rcs_db;
    Re <X> rcs_sm;
};


template <class X>

struct Solution
{
    Source <X> source;
    Field  <X> field;
};


#endif //RCS_SOLVER_SOLUTION_H
