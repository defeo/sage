from sage.rings.ff_compositum.ff_compositum_lib cimport *
from sage.libs.flint.nmod_poly cimport *
from sage.structure.sage_object cimport SageObject

import cython

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
        cdef FFElt res = FFElt(self)
        res.set_mono([1])
        nmod_poly_set(res.dual, self.Mt())
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
                b = FFElt(b_or_Q)
                b.set_dual([1] + [0] * (b.ff.degree - 2))
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

    def __repr__(self):
        return 'FFElt\n  mono: %s\n  dual: %s' % (print_coeffs(self.mono),
                                                  print_coeffs(self.dual))


cdef void set_coeffs(nmod_poly_t res, coeffs, int stop=0):
    'Helper function to set a polynomial from a list'
    cdef int i, k
    stop = max(len(coeffs), stop)
    nmod_poly_fit_length(res, stop)
    for i in range(len(coeffs)):
        nmod_poly_set_coeff_ui(res, i, coeffs[i])
    for i in range(len(coeffs), stop):
        nmod_poly_set_coeff_ui(res, i, 0)

cdef print_coeffs(nmod_poly_t x):
    cdef int i
    return " ".join([str(nmod_poly_get_coeff_ui(x, i)) for i in range(x.length)])
