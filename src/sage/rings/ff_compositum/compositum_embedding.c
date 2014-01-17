#include "compositum.h"
#include "nmod_poly_extra.h"
#include <flint/nmod_vec.h>


/************** EMBEDDING *****************/

/*
  Given x in A=k[X]/P and y in B=k[X]/Q, computes the image of xy in
  the compositum AB. x, y and the result are in dual form.
  
  res must have enough memory to hold deg P x deg Q limbs.

  res must not be equal to P->coeffs. You would be crazy to do so.
 */
void _compositum_embed(mp_ptr res,
		       mp_srcptr x, const nmod_poly_t P, 
		       mp_srcptr y, const nmod_poly_t Q,
		       slong M) {
  mp_ptr tmp = _nmod_vec_init(M);
  nmod_poly_trem(tmp, y, Q, M);
  nmod_poly_trem(res, x, P, M);
  __COEFF_PROD(res, res, tmp, P->mod, M);
  _nmod_vec_clear(tmp);
}



/************** PROJECTION *****************/

/*
  Project x onto the space (k[X]/P)y*, where y* is any element dual to
  y. x is in monomial form, y is in dual form. Result is in monomial
  form.
 */
void _compositum_project(nmod_poly_t res, const nmod_poly_t P,
			 const nmod_poly_t x,
			 mp_srcptr y, const nmod_poly_t Q) {
  slong
    m = nmod_poly_degree(P),
    n = nmod_poly_degree(Q),
    M = FLINT_MIN(x->length, m*n);
  
  mp_ptr tmp = _nmod_vec_init(M);
  nmod_poly_trem(tmp, y, Q, M);
  nmod_poly_fit_length(res, M);
  res->length = M;
  __COEFF_PROD(res->coeffs, x->coeffs, tmp, P->mod, M);
  _nmod_poly_normalise(res);
  _nmod_vec_clear(tmp);

  nmod_poly_rem(res, res, P);
}
