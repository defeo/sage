#include "compositum.h"
#include "nmod_poly_extra.h"
#include <flint/nmod_vec.h>

/*
  res must not alias any z[i]
*/
void _compositum_iso_matrix_from_mono(mp_ptr res, mp_srcptr* z,
				      const nmod_poly_t P, mp_srcptr Pnewton,
				      const nmod_poly_t Q, mp_srcptr Qnewton) {
  slong
    i, j, k,
    m = nmod_poly_degree(P),
    n = nmod_poly_degree(Q),
    M = m*n;

  mp_ptr x = _nmod_vec_init(M+m-1);
  nmod_poly_trem(x, Pnewton, P, M+m-1);
  mp_ptr y = _nmod_vec_init(M+n-1);
  nmod_poly_trem(y, Qnewton, Q, M+n-1);
  _nmod_vec_zero(res, M);
  for (k = 0; k < M; k++) {
    for (j = 0; j < n; j++) {
      if (z[j]) {
	for (i = 0; i < m; i++) {
	  NMOD_ADDMUL(res[k], z[j][i], nmod_mul(x[i+k], y[j+k], P->mod), P->mod);
	}
      }
    }
  }
  _nmod_vec_clear(x);
  _nmod_vec_clear(y);
}

/*
  z must not alias any res[i]
*/
void _compositum_iso_matrix_to_dual(nmod_poly_struct** res, mp_srcptr z,
				    const nmod_poly_t P, mp_srcptr Pnewton,
				    const nmod_poly_t Q, mp_srcptr Qnewton) {
  slong
    i, j, k,
    m = nmod_poly_degree(P),
    n = nmod_poly_degree(Q),
    M = m*n;

  mp_ptr x = _nmod_vec_init(M+m-1);
  nmod_poly_trem(x, Pnewton, P, M+m-1);
  mp_ptr y = _nmod_vec_init(M+n-1);
  nmod_poly_trem(y, Qnewton, Q, M+n-1);
  for (j = 0; j < n; j++) {
    _nmod_vec_zero(res[j]->coeffs, m);
    for (k = 0; k < M; k++) {
      for (i = 0; i < m; i++) {
	NMOD_ADDMUL(res[j]->coeffs[i], z[k], nmod_mul(x[i+k], y[j+k], P->mod), P->mod);
      }
    }
  }
  //_nmod_vec_clear(x);
  //_nmod_vec_clear(y);
}
