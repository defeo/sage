r"""
Implementation of the algebraic closures of finite fields using fast
composita.
"""

import sage.rings.algebraic_closure_finite_field as fpbar
from sage.structure.element cimport Element
from sage.structure.parent cimport Parent
from sage.rings.finite_rings.hom_finite_field cimport FiniteFieldHomomorphism_generic, SectionFiniteFieldHomomorphism_generic
from sage.structure.sage_object import SageObject

class AlgebraicClosureFiniteField_primary_decomp(fpbar.AlgebraicClosureFiniteField_generic):
    """
    Algebraic closure of a finite field, constructed using
    pseudo-Conway polynomials.

    EXAMPLE::

        sage: F = GF(5).algebraic_closure(implementation='primary_decomposition')
        sage: type(F)
        <class 'sage.rings.algebraic_closure_finite_field.AlgebraicClosureFiniteField_primary_decomp_with_category'>

    """
    def __init__(self, base_ring, name, category=None):
        """
        TEST::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField_pseudo_conway
            sage: F = AlgebraicClosureFiniteField_pseudo_conway(GF(3), 'z')
            sage: F
            Algebraic closure of Finite Field of size 3

        """
        from sage.rings.finite_rings.finite_field_base import is_FiniteField
        if not (is_FiniteField(base_ring) and base_ring.is_prime_field()):
            raise NotImplementedError
        self._PR = base_ring.polynomial_ring()
        self._lattice = { }
        fpbar.AlgebraicClosureFiniteField_generic.__init__(self, base_ring, name, category)

    def __getstate__(self):
        """
        Used for pickling.

        TEST::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField_pseudo_conway
            sage: F = AlgebraicClosureFiniteField_pseudo_conway(GF(5), 'z')
            sage: F.__getstate__()
            <class 'sage.rings.finite_rings.conway_polynomials.PseudoConwayLattice'>

        """
        return self._lattice

    def __setstate__(self, state):
        """
        Used for pickling.

        TEST::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField_pseudo_conway
            sage: F = AlgebraicClosureFiniteField_pseudo_conway(GF(5), 'z')
            sage: F.__setstate__(F.__getstate__())

        """
        self._lattice = state

    def __cmp__(self, other):
        """
        Compare ``self`` with ``other``.

        TEST::

            sage: F3 = GF(3).algebraic_closure()
            sage: F3 == F3
            True
            sage: F5 = GF(5).algebraic_closure()
            sage: F3 == F5
            False

        """
        c = fpbar.AlgebraicClosureFiniteField_generic.__cmp__(self, other)
        if c != 0:
            return c
        return cmp(self._lattice, other._lattice)            

    def _get_polynomial(self, n):
        """
        Return the defining polynomial of the unique subfield of
        degree `n` of ``self``.

        EXAMPLE::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField_pseudo_conway
            sage: F = AlgebraicClosureFiniteField_pseudo_conway(GF(5), 'z')
            sage: F._get_polynomial(1)
            x + 3

        """
        try:
            return self._lattice[n].M
        except KeyError:
            f = n.factor()
            if len(f) == 1:
                self._lattice[n] = _field_descriptor(self._PR.irreducible_element(n))
            else:
                from sage.matrix.berlekamp_massey import berlekamp_massey
                from sage.misc.misc_c import prod

                # we solve a knapsack to split n
                # into two halves of roughly the same size
                from sage.numerical.knapsack import knapsack
                _, sq = knapsack([(p.log(prec=10)*e,)*2 + (p,e) for (p,e) in f],
                                 max=n.log(prec=10)/2)
                m = prod(p**e for (_,_,p,e) in sq)
                # we (recursively) compute two subfields
                self._get_polynomial(m)
                self._get_polynomial(n//m)
                P = self._lattice[m]
                Q = self._lattice[n // m]
                b = 2 * P.degree() * Q.degree() + 4
                t = P.trace(b)
                u = Q.trace(b)
                M = berlekamp_massey(_dotprod(t, u).padded_list(b))
                self._lattice[n] = _field_descriptor(M)
            return self._lattice[n].M
            

    def _get_im_gen(self, m, n):
        raise NotImplementedError()

    def inclusion(self, m, n):
        """
        Return the canonical inclusion map from the subfield
        of degree `m` to the subfield of degree `n`.

        EXAMPLE::

            sage: F = GF(3).algebraic_closure()
            sage: F.inclusion(1, 2)
            Ring Coercion morphism:
              From: Finite Field of size 3
              To:   Finite Field in z2 of size 3^2
            sage: F.inclusion(2, 4)
            Ring morphism:
              From: Finite Field in z2 of size 3^2
              To:   Finite Field in z4 of size 3^4
              Defn: z2 |--> 2*z4^3 + 2*z4^2 + 1

        """
        if m == 1:
            return self.base_ring().hom(self._subfield(n))
        elif m.divides(n):
            from sage.categories.homset import Hom
            H = Hom(self._subfield(m), self._subfield(n))
            # this is to make sure that the other field has been computed
            self._subfield(n // m)
            if H.is_aut():
                return H.identity()
            else:
                return FiniteFieldHomomorphism_trace_formulas(H, self._lattice[m], 
                                                              self._lattice[n // m],
                                                              self._lattice[n])
        else:
            raise ValueError("subfield of degree %s not contained in subfield of degree %s" % (m, n))

cdef class FiniteFieldHomomorphism_trace_formulas(FiniteFieldHomomorphism_generic):
    'Trace formula based morphism'

    def __init__(self, parent, P, xP, PxP):
        cdef Element im_gens, ell, ellstar, v, im

        self._P, self._xP, self._PxP = P, xP, PxP

        m = P.degree()
        n = xP.degree()

        u = xP.trace(m*n)
        ell = P.trace(m+1).shift(-1)
        ellstar = P.trem(ell, m*n)

        v = _dotprod(ellstar, u)
        im = self._PxP.convert_from_trace(v)

        im_gens = parent.codomain()(im)
        FiniteFieldHomomorphism_generic.__init__(self, parent, (im_gens,), check=False,
                                                 section_class=SectionFiniteFieldHomomorphism_trace_formulas)

    cpdef Element _call_(self, x):
        cdef Element t, u, ell, ellstar, v, im
        cdef Parent U
    
        m = self._P.degree()
        n = self._xP.degree()
        U = self._P.M.parent()

        t = self._P.trace()
        u = self._xP.trace(m*n)
        ell = self._P.tmulmod(t, U(x.polynomial()))
        ellstar = self._P.trem(ell, m*n)

        v = _dotprod(ellstar, u)
        im = self._PxP.convert_from_trace(v)
        return self.codomain()(im)


cdef class SectionFiniteFieldHomomorphism_trace_formulas(SectionFiniteFieldHomomorphism_generic):
    'Sections of FiniteFieldHomomorphism_trace_formulas'

    cdef Element _tembed(self, ell):
        'Transpose of self._inverse.__call__'
        cdef Element t, u, v, ellS

        h = <FiniteFieldHomomorphism_trace_formulas>self._inverse
        m = h._P.degree()
        n = h._xP.degree()

        t = h._P.trace()
        u = h._xP.trace(m*n)
        v = h._PxP.convert_from_trace(ell)
        ellS = _dotprod(v, u) % h._P.M
        return h._P.tmulmod(t, ellS)

    cpdef Element _call_(self, x):
        cdef Element t, ell, v, im
        cdef Parent U

        h =  <FiniteFieldHomomorphism_trace_formulas>self._inverse
        U = h._P.M.parent()
        m = h._P.degree()
        n = h._xP.degree()

        t = h._PxP.trace()
        ell = h._PxP.tmulmod(t, U(x.polynomial()))
        v = self._tembed(ell)
        im = h._P.convert_from_trace(v) 
        return self.codomain()(im) / n


class _field_descriptor(SageObject):
    from sage.misc.cachefunc import cached_method
    from sage.misc.lazy_attribute import lazy_attribute

    '''
    Data describing a finite field in the lattice
    '''
    def __init__(self, M):
        from sage.rings.power_series_ring import PowerSeriesRing
        R = M.parent()
        self.M = M
        self._PS = PowerSeriesRing(R.base_ring(), name=R.gen().variable_name())
        
    def __getstate__(self):
        return {'M': self.M, '_PS': self._PS}

    def _repr_(self):
        return self.__class__.__name__ + '(' + repr(self.M) + ')'

    def degree(self):
        return self.M.degree()
    
    @cached_method
    def _S(self, m=None):
        'return 1/rev(M) mod x^m'
        if m is None:
            m = self.degree() - 1
        return ( 1 / self._PS(self.M.reverse()).O(m) ).polynomial()

    @lazy_attribute
    def S(self):
        'return 1/rev(P) mod x^m'
        return self._S()

    @lazy_attribute
    def I(self):
        "return 1/M' mod M"
        return self.M.derivative().inverse_mod(self.M)

    @cached_method
    def trace(self, k=None):
        'The k first Newton sums of P (as a polynomial)'
        m = self.M.degree()
        if k is None:
            k = m
        tmp = self.M.derivative().reverse(m-1) * self._S(k) ## m-1 = degree
        return tmp.truncate(k)

    def convert_from_trace(self, t):
        'given traces tr(A B^i) where B is a root of M, recovers the expression A = C(B)'
        m = self.degree()
        N = (self.M.reverse(m) * t).truncate(m)
        Nstar = N.reverse(m-1)
        return (Nstar * self.I) % self.M

    def trem(self, ell, k):
        'Transposed remainder. ell is reduced mod P, outputs k terms'
        m = self.degree()
        G = self._S(k-m)
        A = _tmul(ell, self.M, m, k-m)
        C = (G*A).truncate(k-m)
        return ell - C.shift(m)

    def tmulmod(self, ell, B):
        'Transposed modular multiplication, B.ell mod P'
        m = self.degree()
        A = _tmul(ell, self.M, m, m-1)
        C = (self.S * A).truncate(m-1)
        D = ell - C.shift(m)
        return _tmul(D, B, m-1, m)

cpdef Element _tmul(Element A, Element B, int m, int k):
    'univariate transposed product K^{m+k} -> K^k, with deg(B) <= m'
    return ((A*(B.reverse(m))).truncate(k+m)).shift(-m) ## m = degree

# todo: needs to be sped up
cpdef Element _dotprod(Element A, Element B):
    from itertools import izip_longest
    U = A.parent()
    return U([a*b for a,b in izip_longest(A, B, fillvalue=U.base_ring().zero())])
