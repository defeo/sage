from sage.libs.flint.nmod_poly cimport *

cdef extern from "flint/nmod_poly.h":
    ctypedef long slong
    void nmod_poly_invmod(nmod_poly_t res, nmod_poly_t x, nmod_poly_t m)
    void nmod_poly_rem(nmod_poly_t res, nmod_poly_t x, nmod_poly_t m)
    void _nmod_poly_normalise(nmod_poly_t x)
    void _nmod_poly_add(mp_ptr res, mp_ptr poly1, slong len1, mp_ptr poly2, slong len2, nmod_t mod)
    void _nmod_poly_sub(mp_ptr res, mp_ptr poly1, slong len1, mp_ptr poly2, slong len2, nmod_t mod)

cdef extern from "nmod_poly_extra.h":
    void nmod_poly_inverse_reverse(nmod_poly_t res, const nmod_poly_t P, const long k)
    void nmod_poly_inverse_reverse_main(nmod_poly_t res, const nmod_poly_t P)
    void nmod_poly_ratrecon(nmod_poly_t num, nmod_poly_t den, const nmod_poly_t s, long d)
    void nmod_poly_minimal_polynomial_sequence(nmod_poly_t res, mp_srcptr val, long d)
    void nmod_poly_trem(mp_ptr res, mp_srcptr ell, const nmod_poly_t P, const long k)
    void nmod_poly_tmulmod(mp_ptr res, mp_srcptr ell, const nmod_poly_t B, const nmod_poly_t P, const nmod_poly_t S)
    void nmod_poly_to_newton_sums(mp_ptr newton, const nmod_poly_t P, long k)
    void nmod_poly_from_newton_sums(nmod_poly_t P, mp_srcptr newton, long d)
    void nmod_poly_convert_from_trace(nmod_poly_t C, mp_srcptr t, const nmod_poly_t M, const nmod_poly_t iM)
    void nmod_poly_convert_from_trace_bi(mp_ptr C, mp_srcptr t, 
                                         const nmod_poly_t M, const nmod_poly_t iM, const nmod_poly_t N, const nmod_poly_t iN)
    
cdef extern from "compositum.h":
    void _compositum_embed(mp_ptr res,
                           mp_srcptr x, const nmod_poly_t P, 
                           mp_srcptr y, const nmod_poly_t Q,
                           slong M)
    void _compositum_project(nmod_poly_t res, const nmod_poly_t P,
                             const nmod_poly_t x,
                             mp_srcptr y, const nmod_poly_t Q)
    void _compositum_iso_from_mono(mp_ptr res, mp_srcptr* x,
                                   const nmod_poly_t P, 
                                   const nmod_poly_t Q, mp_srcptr Qnewton)
    void _compositum_iso_to_dual(nmod_poly_struct** res, const nmod_poly_t P,
                                 const nmod_poly_t x, const nmod_poly_t PQ,
                                 const nmod_poly_t Q, mp_srcptr Qnewton)
    void _compositum_iso_to_mono(nmod_poly_struct** res, const nmod_poly_t P,
                                 const nmod_poly_t x, const nmod_poly_t PQ,
                                 const nmod_poly_t Q)
    void _compositum_isomorphism_2(nmod_poly_t G, mp_srcptr F, 
                                   const nmod_poly_t P, mp_srcptr Pt, 
                                   const nmod_poly_t Q, mp_srcptr Qt,
			       const nmod_poly_t R, const nmod_poly_t iR)
    void _compositum_tisomorphism_2(mp_ptr F, mp_srcptr G,
                                    const nmod_poly_t P, mp_srcptr Pt,
                                    const nmod_poly_t Q, mp_srcptr Qt,
                                    const nmod_poly_t R, const nmod_poly_t iR, const nmod_poly_t SR)
    void _compositum_isomorphism_inverse_2(mp_ptr F, const nmod_poly_t G, 
                                           const nmod_poly_t P, mp_srcptr Pt,
                                           const nmod_poly_t Q, mp_srcptr Qt,
                                           const nmod_poly_t R, mp_srcptr Rt,
                                           const nmod_poly_t iR, const nmod_poly_t SR)
    void _compositum_iso_matrix_from_mono(mp_ptr res, mp_srcptr* z,
                                          const nmod_poly_t P, mp_srcptr Pnewton,
                                          const nmod_poly_t Q, mp_srcptr Qnewton)
    void _compositum_iso_matrix_to_dual(nmod_poly_struct** res, mp_srcptr z,
                                        const nmod_poly_t P, mp_srcptr Pnewton,
                                        const nmod_poly_t Q, mp_srcptr Qnewton)
    void _compositum_iP(nmod_poly_t res, const nmod_poly_t P)
