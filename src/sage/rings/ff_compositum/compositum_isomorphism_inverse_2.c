#include "compositum.h"
#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"

/*------------------------------------------------------------------------*/
/* F is an array of size m*n (m=deg(P), n=deg(Q))                         */
/* of the form F0..F_{m-1}, with Fi = coeff(F,X^i) \in Fp[y]/Q(y)         */
/*------------------------------------------------------------------------*/

void _compositum_isomorphism_inverse_2(mp_ptr F, const nmod_poly_t G, 
				       const nmod_poly_t P, mp_srcptr Pt,
				       const nmod_poly_t Q, mp_srcptr Qt,
				       const nmod_poly_t R, mp_srcptr Rt,
				       const nmod_poly_t iR, const nmod_poly_t SR){

  long m = nmod_poly_degree(P);
  long n = nmod_poly_degree(Q);

  mp_ptr ell = _nmod_vec_init(m*n);
  nmod_poly_tmulmod(ell, Rt, G, R, SR);

  _compositum_tisomorphism_2(F, ell, P, Pt, Q, Qt, R, iR, SR);

  _nmod_vec_clear(ell);
}
