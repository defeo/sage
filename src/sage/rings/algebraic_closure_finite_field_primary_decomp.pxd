from sage.structure.element cimport Element
from sage.structure.parent cimport Parent
from sage.rings.finite_rings.hom_finite_field cimport FiniteFieldHomomorphism_generic, SectionFiniteFieldHomomorphism_generic

cdef class FiniteFieldHomomorphism_trace_formulas(FiniteFieldHomomorphism_generic):
    cdef _P, _xP, _PxP
    cpdef Element _call_(self, x)

cdef class SectionFiniteFieldHomomorphism_trace_formulas(SectionFiniteFieldHomomorphism_generic):
    cdef Element _tembed(self, ell)
    cpdef Element _call_(self, x)

cpdef Element _tmul(Element A, Element B, int m, int k)
cpdef Element _dotprod(Element A, Element B)
