#ifndef COMPOSITUM__H
#define COMPOSITUM__H

#include <flint/nmod_poly.h>

#define __COEFF_PROD(res, x, y, mod, M) for ((M)--; (M) >= WORD(0); (M)--) (res)[M] = nmod_mul((x)[M], (y)[M], (mod))


/*********** EMBEDDING & PROJECTION ***************/
void _compositum_embed(mp_ptr res,
		       mp_srcptr x, const nmod_poly_t P, 
		       mp_srcptr y, const nmod_poly_t Q,
		       slong M);
void _compositum_project(nmod_poly_t res, const nmod_poly_t P,
			 const nmod_poly_t x,
			 mp_srcptr y, const nmod_poly_t Q);
/************** ISOMORPHISM *****************/
void _compositum_iso_from_mono(mp_ptr res, mp_srcptr* x,
			       const nmod_poly_t P, 
			       const nmod_poly_t Q, mp_srcptr Qnewton);
void _compositum_iso_to_dual(nmod_poly_struct** res, const nmod_poly_t P,
			     const nmod_poly_t x, const nmod_poly_t PQ,
			     const nmod_poly_t Q, mp_srcptr Qnewton);
void _compositum_iso_to_mono(nmod_poly_struct** res, const nmod_poly_t P,
			     const nmod_poly_t x, const nmod_poly_t PQ,
			     const nmod_poly_t Q);
void _compositum_isomorphism_2(nmod_poly_t G, mp_srcptr F, 
			       const nmod_poly_t P, mp_srcptr Pt, 
			       const nmod_poly_t Q, mp_srcptr Qt,
			       const nmod_poly_t R, const nmod_poly_t iR);
void _compositum_tisomorphism_2(mp_ptr F, mp_srcptr G,
				const nmod_poly_t P, mp_srcptr Pt,
				const nmod_poly_t Q, mp_srcptr Qt,
				const nmod_poly_t R, const nmod_poly_t iR, const nmod_poly_t SR);
void _compositum_isomorphism_inverse_2(mp_ptr F, const nmod_poly_t G, 
				       const nmod_poly_t P, mp_srcptr Pt,
				       const nmod_poly_t Q, mp_srcptr Qt,
				       const nmod_poly_t R, mp_srcptr Rt,
				       const nmod_poly_t iR, const nmod_poly_t SR);
void _compositum_iso_matrix_from_mono(mp_ptr res, mp_srcptr* z,
				      const nmod_poly_t P, mp_srcptr Pnewton,
				      const nmod_poly_t Q, mp_srcptr Qnewton);
void _compositum_iso_matrix_to_dual(nmod_poly_struct** res, mp_srcptr z,
				    const nmod_poly_t P, mp_srcptr Pnewton,
				    const nmod_poly_t Q, mp_srcptr Qnewton);
/************** UTILS *****************/
void _compositum_iP(nmod_poly_t res, const nmod_poly_t P);


#endif
