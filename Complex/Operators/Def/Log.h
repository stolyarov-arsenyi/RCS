template <class X> Bo <X> operator ! (const Bo <X> & r) { return ! X(r); }

template <class X> Bo <X> operator && (const Bo <X> & l, const Bo <X> & r) { return X(l) && X(r); }
template <class X> Bo <X> operator && (const Bo <X> & l, const     X  & r) { return X(l) && X(r); }
template <class X> Bo <X> operator && (const     X  & l, const Bo <X> & r) { return X(l) && X(r); }

template <class X> Bo <X> operator || (const Bo <X> & l, const Bo <X> & r) { return X(l) || X(r); }
template <class X> Bo <X> operator || (const Bo <X> & l, const     X  & r) { return X(l) || X(r); }
template <class X> Bo <X> operator || (const     X  & l, const Bo <X> & r) { return X(l) || X(r); }