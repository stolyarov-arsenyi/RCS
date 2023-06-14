template <class X> Bo <X> operator == (const Re <X> & l, const Re <X> & r) { return X(l) == X(r);  }
template <class X> Bo <X> operator == (const Re <X> & l, const     X  & r) { return X(l) == X(r);  }
template <class X> Bo <X> operator == (const     X  & l, const Re <X> & r) { return X(l) == X(r);  }

template <class X> Bo <X> operator != (const Re <X> & l, const Re <X> & r) { return X(l) != X(r);  }
template <class X> Bo <X> operator != (const Re <X> & l, const     X  & r) { return X(l) != X(r);  }
template <class X> Bo <X> operator != (const     X  & l, const Re <X> & r) { return X(l) != X(r);  }

template <class X> Bo <X> operator >  (const Re <X> & l, const Re <X> & r) { return X(l) >  X(r);  }
template <class X> Bo <X> operator >  (const Re <X> & l, const     X  & r) { return X(l) >  X(r);  }
template <class X> Bo <X> operator >  (const     X  & l, const Re <X> & r) { return X(l) >  X(r);  }

template <class X> Bo <X> operator >= (const Re <X> & l, const Re <X> & r) { return X(l) >= X(r);  }
template <class X> Bo <X> operator >= (const Re <X> & l, const     X  & r) { return X(l) >= X(r);  }
template <class X> Bo <X> operator >= (const     X  & l, const Re <X> & r) { return X(l) >= X(r);  }

template <class X> Bo <X> operator <  (const Re <X> & l, const Re <X> & r) { return X(l) <  X(r);  }
template <class X> Bo <X> operator <  (const Re <X> & l, const     X  & r) { return X(l) <  X(r);  }
template <class X> Bo <X> operator <  (const     X  & l, const Re <X> & r) { return X(l) <  X(r);  }

template <class X> Bo <X> operator <= (const Re <X> & l, const Re <X> & r) { return X(l) <= X(r);  }
template <class X> Bo <X> operator <= (const Re <X> & l, const     X  & r) { return X(l) <= X(r);  }
template <class X> Bo <X> operator <= (const     X  & l, const Re <X> & r) { return X(l) <= X(r);  }