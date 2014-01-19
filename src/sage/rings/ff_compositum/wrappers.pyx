from sage.rings.ff_compositum.ff_compositum_lib cimport *
from sage.libs.flint.nmod_poly cimport *
from sage.structure.sage_object cimport SageObject

import cython
#from libc.stdio cimport printf
from libc.stdlib cimport malloc, free

cdef class FFDesc(SageObject):
    '''
    Finite field descriptor.
    
    This object wraps four precomputed polynomials defining a finite
    field. Not meant to be exposed to the user.
    '''
    cdef nmod_poly_t M
    cdef nmod_poly_t im
    cdef nmod_poly_t s
    cdef nmod_poly_t mt
    cdef int degree
    
    def __cinit__(self, int n, *args, **kwds):
        nmod_poly_init(self.M, n)
        nmod_poly_init(self.im, n)
        nmod_poly_init(self.s, n)
        nmod_poly_init(self.mt, n)
        self.degree = 0
        
    def __init__(self, int n, modulus=None):
        if modulus:
            self.set_modulus(modulus)

    def __dealloc__(self):
        nmod_poly_clear(self.M)
        nmod_poly_clear(self.im)
        nmod_poly_clear(self.s)
        nmod_poly_clear(self.mt)

    cpdef set_modulus(self, coeffs):
        cdef int i, k
        set_coeffs(self.M, coeffs)
        self.degree = nmod_poly_degree(self.M)

    cdef nmod_poly_struct* iM(self):
        if not self.im.length and self.degree:
            _compositum_iP(self.im, self.M)
        return self.im
    
    cdef nmod_poly_struct* S(self):
        if not self.s.length and self.degree:
            nmod_poly_inverse_reverse_main(self.s, self.M)
        return self.s

    cdef nmod_poly_struct* Mt(self):
        if not self.mt.length and self.degree:
            nmod_poly_fit_length(self.mt, self.degree)
            nmod_poly_to_newton_sums(self.mt.coeffs, self.M, self.degree)
            self.mt.length = self.degree
        return self.mt

    cpdef FFDesc compositum(self, FFDesc other):
        cdef int p = self.M.mod.n
        cdef int d = self.degree * other.degree
        cdef int k = d + 1 + (p < d) * d
        cdef FFDesc comp = FFDesc(p)
        comp.degree = d

        nmod_poly_fit_length(comp.mt, k)
        _compositum_embed(comp.mt.coeffs, 
                          self.Mt().coeffs, self.M,
                          other.Mt().coeffs, other.M,
                          k)
        comp.mt.length = k

        if (k == d + 1):
            nmod_poly_from_newton_sums(comp.M, comp.mt.coeffs, d);
        else:
            nmod_poly_minimal_polynomial_sequence(comp.M, comp.mt.coeffs, d);

        return comp

    cpdef FFElt zero(self):
        cdef FFElt res = FFElt(self)
        return res

    cpdef FFElt one(self):
        return self.mono_basis_element(0)
    
    cpdef FFElt gen(self):
        return self.mono_basis_element(1)

    cpdef FFElt mono_basis_element(self, int i):
        cdef FFElt res = FFElt(self)
        nmod_poly_set_coeff_ui(res.mono, i, 1)
        return res    

    cpdef FFElt dual_basis_element(self, int i):
        cdef FFElt res = FFElt(self)
        set_coeffs(res.dual, [], self.degree)
        nmod_poly_set_coeff_ui(res.dual, i, 1)
        return res

    cpdef FFElt seqelt_mono(self, elts, FFDesc Q):
        cdef FFElt res = FFElt(self)
        cdef FFDesc P = (<FFElt>elts[0]).ff
        cdef mp_srcptr* coeffs = []
        for i in range(Q.degree):
            coeffs[i] = (<FFElt>elts[i]).get_dual().coeffs
        nmod_poly_fit_length(res.dual, self.degree)
        _compositum_iso_from_mono(res.dual.coeffs, coeffs, P.M, Q.M, Q.Mt().coeffs)
        res.dual.length = self.degree
        return res
                                  
    cpdef FFElt seqelt_mono_python(self, elts, FFDesc Q):
        cdef FFElt res = FFElt(self)
        for i in range(Q.degree):
            res += elts[i].embed(Q.mono_basis_element(i), self)
        return res

    cpdef FFElt seqelt_mono_BSGS(self, elts, FFDesc Q):
        cdef FFElt res = FFElt(self)
        cdef FFDesc P = (<FFElt>elts[0]).ff
        cdef int m = P.degree, n = Q.degree
        cdef mp_limb_t* coeffs = <mp_limb_t*>malloc(m*n*sizeof(mp_limb_t))
        cdef mp_limb_t* tmp
        for j in range(n):
            if elts[j]:
                tmp = (<FFElt>elts[j]).get_dual().coeffs
            else:
                tmp = NULL
            for i in range(m):
                coeffs[i*n + j] = tmp[i] if tmp else 0
        _compositum_isomorphism_2(res.mono, coeffs, 
                                  P.M, P.Mt().coeffs,
                                  Q.M, Q.Mt().coeffs,
                                  self.M, self.iM())
        free(coeffs)
        return res
    

    cpdef FFElt seqelt_dual_python(self, elts, FFDesc Q):
        cdef FFElt res = FFElt(self)
        for i in range(Q.degree):
            res += elts[i].embed(Q.dual_basis_element(i), self)
        return res

    def __repr__(self):
        return "FFDesc(%d, (%s))\n  %s\n  %s \n  %s" % (self.M.mod.n, print_coeffs(self.M),
                                                        print_coeffs(self.S()), print_coeffs(self.iM()), print_coeffs(self.Mt()))


cdef class FFElt(SageObject):
    '''
    Elements of FFDesc with lazy double basis representation.
    '''
    cdef FFDesc ff
    cdef nmod_poly_t mono
    cdef nmod_poly_t dual

    def __cinit__(self, FFDesc ff, *args, **kwds):
        nmod_poly_init(self.mono, ff.M.mod.n)
        nmod_poly_init(self.dual, ff.M.mod.n)
    
    def __init__(self, FFDesc ff, coeffs=None, dual_coeffs=None):
        self.ff = ff
        if coeffs is not None:
            self.set_mono(coeffs)
        elif dual_coeffs:
            self.set_dual(dual_coeffs)

    def __dealloc__(self):
        nmod_poly_clear(self.mono)
        nmod_poly_clear(self.dual)

    cpdef set_mono(self, coeffs):
        set_coeffs(self.mono, coeffs)
        nmod_poly_rem(self.mono, self.mono, self.ff.M)
        nmod_poly_zero(self.dual)

    cpdef set_dual(self, coeffs):
        set_coeffs(self.dual, coeffs, self.ff.degree)
        nmod_poly_zero(self.mono)

    cdef nmod_poly_struct* get_mono(self):
        if not self.mono.length and self.dual.length:
            nmod_poly_convert_from_trace(self.mono, self.dual.coeffs, 
                                         self.ff.M, self.ff.iM())
        return self.mono
    
    cdef nmod_poly_struct* get_dual(self):
        if not self.dual.length and self.mono.length:
            nmod_poly_fit_length(self.dual, self.ff.degree)
            nmod_poly_tmulmod(self.dual.coeffs, self.ff.Mt().coeffs,
                              self.mono, self.ff.M, self.ff.S())
            self.dual.length = self.ff.degree
        return self.dual

    def __nonzero__(self):
        return self.mono.length or self.dual.length

    def __richcmp__(FFElt a, FFElt b, int op):
        eq = op in (1,2,5)
        if not a or not b:
            return (a or b) != eq
        elif a.dual.length and b.dual.length:
            return nmod_poly_equal(a.dual, b.dual) == eq
        else:
            return nmod_poly_equal(a.get_mono(), b.get_mono()) == eq

    def __add__(FFElt a, FFElt b):
        if not a:
            return b
        if not b:
            return a
        cdef FFElt res = FFElt(a.ff)
        if a.mono.length and b.mono.length:
            nmod_poly_add(res.mono, a.mono, b.mono)
        if a.dual.length and b.dual.length:
            nmod_poly_fit_length(res.dual, a.ff.degree)
            _nmod_poly_add(res.dual.coeffs, a.dual.coeffs, a.ff.degree,
                           b.dual.coeffs, a.ff.degree, a.ff.M.mod)
            res.dual.length = a.ff.degree
        return res

    def __sub__(FFElt a, FFElt b):
        if not a:
            return b
        if not b:
            return a
        cdef FFElt res = FFElt(a.ff)
        if a.mono.length and b.mono.length:
            nmod_poly_sub(res.mono, a.mono, b.mono)
        if a.dual.length and b.dual.length:
            nmod_poly_fit_length(res.dual, a.ff.degree)
            _nmod_poly_sub(res.dual.coeffs, a.dual.coeffs, a.ff.degree,
                           b.dual.coeffs, a.ff.degree, a.ff.M.mod)
            res.dual.length = a.ff.degree
        return res

    def __mul__(FFElt a, FFElt b):
        cdef FFElt res = FFElt(a.ff)
        if not a or not b:
            return res
        if a.mono.length and not b.mono.length:
            a, b = b, a
        if not a.mono.length and b.mono.length:
            nmod_poly_fit_length(res.dual, a.ff.degree)
            nmod_poly_tmulmod(res.dual.coeffs, a.dual.coeffs, b.mono, a.ff.M, a.ff.S())
            res.dual.length = a.ff.degree
            nmod_poly_zero(res.mono)
        else:
            nmod_poly_mulmod(res.mono, a.get_mono(), b.get_mono(), a.ff.M)
            nmod_poly_zero(res.dual)
        return res

    def __inv__(self):
        cdef FFElt res = FFElt(self.ff)
        if not self:
            raise ZeroDivisionError()
        nmod_poly_invmod(res.mono, self.get_mono(), self.ff.M)
        
    cpdef FFElt embed(FFElt a, b_or_Q, FFDesc R):
        cdef FFElt b, res = FFElt(R)
        if a and b_or_Q:
            if isinstance(b_or_Q, FFDesc):
                b = b_or_Q.one()
            elif isinstance(b_or_Q, FFElt):
                b = b_or_Q
            else:
                raise ValueError

            nmod_poly_fit_length(res.dual, R.degree)
            _compositum_embed(res.dual.coeffs, 
                              a.get_dual().coeffs, a.ff.M,
                              b.get_dual().coeffs, b.ff.M,
                              R.degree)
            res.dual.length = R.degree
        return res

    cpdef FFElt project(FFElt a, FFDesc P, b_or_Q):
        cdef FFElt b, res = FFElt(P)
        if a:
            if isinstance(b_or_Q, FFDesc):
                b = b_or_Q.dual_basis_element(0)
            elif isinstance(b_or_Q, FFElt) and b_or_Q:
                b = b_or_Q
            else:
                raise ValueError
            _compositum_project(res.mono, P.M, a.get_mono(), 
                                b.get_dual().coeffs, b.ff.M)
        return res

    cpdef FFElt trace(self, FFDesc P, FFDesc Q):
        'Relative trace to P'
        cdef FFElt res = FFElt(P)
        if self:
            _compositum_project(res.mono, P.M, self.get_mono(), 
                                Q.Mt().coeffs, Q.M)
        return res

    cpdef eltseq_mono(self, FFDesc P, FFDesc Q):
        '''
        Returns the list of coefficients of self as an element on the
        monomial basis of Q
        '''
        cdef nmod_poly_struct** coeffs = []
        res = [FFElt(P) for i in range(Q.degree)]
        if self:
            for i in range(Q.degree):
                coeffs[i] = <nmod_poly_struct*>(<FFElt>res[i]).mono
            _compositum_iso_to_mono(coeffs, P.M, self.get_mono(), self.ff.M, Q.M)
        return res

    cpdef eltseq_mono_python(self, FFDesc P, FFDesc Q):
        cdef FFElt b
        cdef int i
        res = []
        for i in range(Q.degree):
            b = Q.dual_basis_element(i)
            res.append(self.project(P, b))
        return res

    cpdef eltseq_mono_BSGS(self, FFDesc P, FFDesc Q):
        cdef int m = P.degree, n = Q.degree
        cdef mp_limb_t* coeffs = <mp_limb_t*>malloc(m*n*sizeof(mp_limb_t))
        _compositum_isomorphism_inverse_2(coeffs, self.get_mono(),
                                          P.M, P.Mt().coeffs,
                                          Q.M, Q.Mt().coeffs,
                                          self.ff.M, self.ff.Mt().coeffs,
                                          self.ff.iM(), self.ff.S())
        res = [FFElt(P) for i in range(n)]
        cdef nmod_poly_struct* tmp
        for j in range(n):
            tmp = <nmod_poly_struct*>(<FFElt>res[i]).mono
            nmod_poly_fit_length(tmp, m)
            for i in range(m):
                tmp.coeffs[i] = coeffs[i*n + j]
            tmp.length = m
            _nmod_poly_normalise(tmp)
        free(coeffs)
        return res

    cpdef eltseq_dual(self, FFDesc P, FFDesc Q):
        '''
        Returns the list of coefficients of self as an element on the
        monomial basis of Q
        '''
        cdef nmod_poly_struct** coeffs = []
        res = [FFElt(P) for i in range(Q.degree)]
        if self:
            for i in range(Q.degree):
                coeffs[i] = <nmod_poly_struct*>(<FFElt>res[i]).mono
            _compositum_iso_to_dual(coeffs, P.M, self.get_mono(), self.ff.M, Q.M, Q.Mt().coeffs)
        return res

    cpdef eltseq_dual_python(self, FFDesc P, FFDesc Q):
        cdef FFElt b
        cdef int i
        res = []
        for i in range(Q.degree):
            b = Q.mono_basis_element(i)
            res.append(self.project(P, b))
        return res
    

    def __repr__(self):
        return 'FFElt\n  mono: %s\n  dual: %s' % (print_coeffs(self.mono),
                                                  print_coeffs(self.dual))


cdef void set_coeffs(nmod_poly_t res, coeffs, int stop=0):
    'Helper function to set a polynomial from a list'
    cdef int i, k
    stop = max(len(coeffs), stop)
    nmod_poly_fit_length(res, stop)
    res.length = stop
    for i in range(len(coeffs)):
        res.coeffs[i] = coeffs[i]
    for i in range(len(coeffs), stop):
        res.coeffs[i] = 0

cdef print_coeffs(nmod_poly_t x):
    cdef int i
    return " ".join([str(nmod_poly_get_coeff_ui(x, i)) for i in range(x.length)])
