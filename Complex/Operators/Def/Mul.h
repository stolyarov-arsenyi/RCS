template <class X> Co <X> operator * (const Co <X> & l, const Co <X> & r) { return { X(l.r) * X(r.r) - X(l.i) * X(r.i), X(l.r) * X(r.i) + X(l.i) * X(r.r) }; }
template <class X> Co <X> operator * (const Co <X> & l, const Im <X> & r) { return {                 - X(l.i) * X(r  ), X(l.r) * X(r  )                   }; }
template <class X> Co <X> operator * (const Co <X> & l, const Re <X> & r) { return { X(l.r) * X(r  )                  ,                 + X(l.i) * X(r  ) }; }
template <class X> Co <X> operator * (const Co <X> & l, const     X  & r) { return { X(l.r) * X(r  )                  ,                 + X(l.i) * X(r  ) }; }
template <class X> Co <X> operator * (const Im <X> & l, const Co <X> & r) { return {                 - X(l  ) * X(r.i),                 + X(l  ) * X(r.r) }; }
template <class X> Re <X> operator * (const Im <X> & l, const Im <X> & r) { return {                 - X(l  ) * X(r  )                                    }; }
template <class X> Im <X> operator * (const Im <X> & l, const Re <X> & r) { return {                                                    + X(l  ) * X(r  ) }; }
template <class X> Im <X> operator * (const Im <X> & l, const     X  & r) { return {                                                    + X(l  ) * X(r  ) }; }
template <class X> Co <X> operator * (const Re <X> & l, const Co <X> & r) { return { X(l  ) * X(r.r)                  , X(l  ) * X(r.i)                   }; }
template <class X> Im <X> operator * (const Re <X> & l, const Im <X> & r) { return {                                    X(l  ) * X(r  )                   }; }
template <class X> Re <X> operator * (const Re <X> & l, const Re <X> & r) { return { X(l  ) * X(r  )                                                      }; }
template <class X> Re <X> operator * (const Re <X> & l, const     X  & r) { return { X(l  ) * X(r  )                                                      }; }
template <class X> Co <X> operator * (const     X  & l, const Co <X> & r) { return { X(l  ) * X(r.r)                  , X(l  ) * X(r.i)                   }; }
template <class X> Im <X> operator * (const     X  & l, const Im <X> & r) { return {                                    X(l  ) * X(r  )                   }; }
template <class X> Re <X> operator * (const     X  & l, const Re <X> & r) { return { X(l  ) * X(r  )                                                      }; }