#include "compositum.h"
#include <math.h>
#include <flint/nmod_poly_mat.h>
#include <flint/nmod_vec.h>
#include "nmod_poly_extra.h"
#include "nmod_poly_mat_extra.h"

/*------------------------------------------------------------------------*/
/* input G is a linear form defined modulo R                              */
/* output F is an array of size m*n (m=deg(P), n=deg(Q))                  */
/* of the form F0..F_{m-1}, with Fi = G(S^i T^j), j=0..n-1y]/Q(y)         */
/* with S=Phi(X), T=Phi(Y), Phi = isomorphism PxQ -> R                    */
/*------------------------------------------------------------------------*/

void _compositum_tisomorphism_2(mp_ptr F, mp_srcptr G,
				const nmod_poly_t P, mp_srcptr Pt,
				const nmod_poly_t Q, mp_srcptr Qt,
				const nmod_poly_t R, const nmod_poly_t iR, const nmod_poly_t SR){
  

  nmod_t mod = P->mod;
  long m = nmod_poly_degree(P);
  long n = nmod_poly_degree(Q);
  long M = m*n;

  long np = n + m - 1;
  long p = ceil(sqrt(np));
  long q = ceil((1.0*np)/p);
  long i, j, k;

  nmod_poly_t iT, iTm, T;
  nmod_poly_struct *TT, *VV;

  nmod_poly_init(iTm, mod.n);
  nmod_poly_init(iT, mod.n);
  nmod_poly_init(T, mod.n);

  mp_ptr 
    X = _nmod_vec_init(M+1),
    uP = _nmod_vec_init(M);

  TT = flint_malloc(sizeof(nmod_poly_struct) * (q+1));
  for (i = 0; i < q+1; i++)
    nmod_poly_init(TT+i, mod.n);

  VV = flint_malloc(sizeof(nmod_poly_struct) * p);
  for (i = 0; i < p; i++)
    nmod_poly_init(VV+i, mod.n);


  //nmod_poly_set_coeff_ui(X, 1, 1);
  //embeddings_embed(T, X, FQ, FP, FR);
  nmod_poly_trem(uP, Pt, P, M);
  nmod_poly_trem(X, Qt, Q, M+1);
  __COEFF_PROD(X+1, X+1, uP, Q->mod, M);
  nmod_poly_convert_from_trace(T, X+1, R, iR);

  nmod_poly_invmod(iT, T, R);
  nmod_poly_powmod_ui_binexp(iTm, iT, m-1, R);

  nmod_poly_mat_t MT, MH, MV;
  nmod_poly_mat_init(MT, q, n, mod.n);
  nmod_poly_mat_init(MV, p, n, mod.n);
  nmod_poly_mat_init(MH, p, q, mod.n);


  nmod_poly_zero(TT);
  nmod_poly_set_coeff_ui(TT, 0, 1);

  for (i = 1; i < q+1; i++) {
    nmod_poly_mulmod(TT+i, TT+(i-1), T, R);
  }

  for (i = 0; i < q; i++)
    for (j = 0; j < n; j++){
      long jm = j*m;
      for (k = 0; k < m; k++)
  	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(MT, i, j), k, nmod_poly_get_coeff_ui(TT+i, k+jm));
    }
  
  mp_ptr TG = _nmod_vec_init(m*n), VTG = _nmod_vec_init(m*n+m-1);
  nmod_poly_tmulmod(TG, G, iTm, R, SR);
  for (i = 0; i < p; i++){
    nmod_poly_trem(VTG, TG, R, m*n+m-1);
    for (j = 0; j < n; j++){
      long jm = j*m;
      for (k = jm; k < jm+2*m-1; k++)
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(MV, i, j), k-jm, VTG[k]);
    }
    nmod_poly_tmulmod(TG, TG, TT+q, R, SR);
  }

  nmod_poly_mat_right_tmul(MH, MV, MT, m-1, m);

  for (i = 0; i < m; i++)
    for (j = 0; j < np; j++)
      if ((i+j-(m-1) >= 0) && (i+j-(m-1) < n))
	F[i*n+i+j-(m-1)] = nmod_poly_get_coeff_ui(MH->entries+j, i);

  
  _nmod_vec_clear(TG);
  _nmod_vec_clear(VTG);
 
  nmod_poly_mat_clear(MT);
  nmod_poly_mat_clear(MH);
  nmod_poly_mat_clear(MV);

  for (i = 0; i < q+1; i++)
    nmod_poly_clear(TT+i);
  flint_free(TT);

  for (i = 0; i < p; i++)
    nmod_poly_clear(VV+i);
  flint_free(VV);

  nmod_poly_clear(iTm);
  nmod_poly_clear(iT);
  nmod_poly_clear(T);

  _nmod_vec_clear(X);
  _nmod_vec_clear(uP);

}
