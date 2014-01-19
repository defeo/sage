#include "compositum.h"
#include "nmod_poly_extra.h"
#include <flint/nmod_vec.h>

/*
  res must not alias any x[i]
*/
void _compositum_iso_from_mono(mp_ptr res, mp_srcptr* x,
			       const nmod_poly_t P, 
			       const nmod_poly_t Q, mp_srcptr Qnewton) {
  slong
    k,
    m = nmod_poly_degree(P),
    n = nmod_poly_degree(Q),
    M = m*n;

  mp_ptr y = _nmod_vec_init(M+n-1);
  nmod_poly_trem(y, Qnewton, Q, M+n-1);
  mp_ptr tmp = _nmod_vec_init(M);
  _nmod_vec_zero(res, M);
  for (n--; n >= 0; n--) {
    if (x[n]) {
      nmod_poly_trem(tmp, x[n], P, M);
      k = M;
      __COEFF_PROD(tmp, tmp, y+n, P->mod, k);
      _nmod_vec_add(res, res, tmp, M, P->mod);
    }
  }
  _nmod_vec_clear(y);
  _nmod_vec_clear(tmp);
}

/*
  x must not alias any res[i]
*/
void _compositum_iso_to_dual(nmod_poly_struct** res, const nmod_poly_t P,
			     const nmod_poly_t x, const nmod_poly_t PQ,
			     const nmod_poly_t Q, mp_srcptr Qnewton) {
  slong
    k,
    n = nmod_poly_degree(Q),
    M = nmod_poly_degree(PQ);
  mp_ptr y = _nmod_vec_init(M+n-1);
  nmod_poly_trem(y, Qnewton, Q, M+n-1);
  for (n--; n >= 0; n--) {
    nmod_poly_fit_length(res[n], M);
    k = FLINT_MIN(x->length, M);
    res[n]->length = k;
    __COEFF_PROD(res[n]->coeffs, x->coeffs, y+n, P->mod, k);
    _nmod_poly_normalise(res[n]);
    nmod_poly_rem(res[n], res[n], P);
  }
  _nmod_vec_clear(y);
}


/*
  x must not alias any res[i]
*/
void _compositum_iso_to_mono(nmod_poly_struct** res, const nmod_poly_t P,
			     const nmod_poly_t x, const nmod_poly_t PQ,
			     const nmod_poly_t Q) {
  slong
    k,
    n = nmod_poly_degree(Q),
    M = nmod_poly_degree(PQ);
  mp_ptr y = _nmod_vec_init(M);
  _nmod_vec_zero(y, n);
  for (n--; n >= 0; n--) {
    y[n] = WORD(1);
    nmod_poly_trem(y, y, Q, M);
    nmod_poly_fit_length(res[n], M);
    k = FLINT_MIN(x->length, M);
    res[n]->length = k;
    __COEFF_PROD(res[n]->coeffs, x->coeffs, y, P->mod, k);
    _nmod_poly_normalise(res[n]);
    nmod_poly_rem(res[n], res[n], P);
    y[n] = WORD(0);
  }
  _nmod_vec_clear(y);
}
