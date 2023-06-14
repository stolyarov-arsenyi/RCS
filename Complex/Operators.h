#include "Operators/Def/Pos.h"
#include "Operators/Def/Neg.h"
#include "Operators/Def/Con.h"
#include "Operators/Def/Add.h"
#include "Operators/Def/Sub.h"
#include "Operators/Def/Mul.h"
#include "Operators/Def/Div.h"
#include "Operators/Def/Cmp.h"
#include "Operators/Def/Log.h"


template <typename L, typename R> L & operator += (L & l, const R & r) {return l = l + r; }
template <typename L, typename R> L & operator -= (L & l, const R & r) {return l = l - r; }
template <typename L, typename R> L & operator *= (L & l, const R & r) {return l = l * r; }
template <typename L, typename R> L & operator /= (L & l, const R & r) {return l = l / r; }