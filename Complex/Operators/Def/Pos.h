template <class X> Co <X> operator + (const Co <X> & r) { return { + X(r.r), + X(r.i) }; }
template <class X> Im <X> operator + (const Im <X> & r) { return {           + X(r  ) }; }
template <class X> Re <X> operator + (const Re <X> & r) { return { + X(r  )           }; }