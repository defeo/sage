r"""
Algebraic closures of finite fields

Let `\Bold{F}` be a finite field, and let `\overline{\Bold{F}}` be an
algebraic closure of `\Bold{F}`; this is unique up to (non-canonical)
isomorphism.  For every `n\ge 1`, there is a unique subfield
`\Bold{F}_n` of `\overline{\Bold{F}}` such that
`\Bold{F}\subset\Bold{F}_n` and '[\Bold{F}_n:\Bold{F}]=1`.

In Sage, algebraic closures of finite fields are implemented using
compatible systems of finite fields.  The resulting Sage object keeps
track of a finite lattice of the subfields `\Bold{F}_n` and the
embeddings between them.  This lattice is extended as necessary.

The Sage class corresponding to `\overline{\Bold{F}}` can be
constructed from the finite field `\Bold{F}` by using the
:meth:`~sage.rings.finite_rings.finite_field_base.FiniteField.algebraic_closure`
method.  This invokes the ``AlgebraicClosureFiniteField`` factory
object to get unique representation.

The Sage class for elements of `\overline{\Bold{F}}` is
:class:`AlgebraicClosureFiniteFieldElement`.  Such an element is
represented as an element of one of the `\Bold{F}_n`.  This means that
each element `x\in\Bold{F}` has infinitely many different
representations, one for each `n` such that `x` is in `\Bold{F}_n`.

.. NOTE::

    Only prime finite fields are currently accepted as base fields for
    algebraic closures.  To obtain an algebraic closure of a non-prime
    finite field `\Bold{F}`, take an algebraic closure of the prime
    field of `\Bold{F}` and embed `\Bold{F}` into this.

    Algebraic closures of finite fields are currently implemented
    using (pseudo-)Conway polynomials; see
    :class:`AlgebraicClosureFiniteField_pseudo_conway` and the module
    :mod:`~sage.rings.finite_rings.conway_polynomials`.  Other
    implementations may be added by creating appropriate subclasses of
    :class:`AlgebraicClosureFiniteField_generic`.

TEST::

    sage: F = GF(5).algebraic_closure('z')
    sage: TestSuite(F).run()

AUTHORS:

- Peter Bruin (August 2013): initial version

"""

from sage.rings.finite_rings.element_base import is_FiniteFieldElement
from sage.rings.finite_rings.finite_field_base import is_FiniteField
from sage.rings.ring import Field
from sage.structure.element import FieldElement
from sage.structure.factory import UniqueFactory

class AlgebraicClosureFiniteFieldElement(FieldElement):
    """
    Element of an algebraic closure of a finite field.

    EXAMPLE::

        sage: F = GF(3).algebraic_closure('z')
        sage: F.gen(2)
        z2
        sage: type(F.gen(2))
        <class 'sage.rings.algebraic_closure_finite_field.AlgebraicClosureFiniteField_pseudo_conway_with_category.element_class'>

    """
    def __init__(self, parent, value):
        """
        TEST::

            sage: F = GF(3).algebraic_closure('z')
            sage: TestSuite(F.gen(2)).run()

        """
        if is_FiniteFieldElement(value):
            n = value.parent().degree()
        else:
            from sage.rings.integer import Integer
            n = Integer(1)
        self._value = parent._subfield(n).coerce(value)
        self._level = n
        FieldElement.__init__(self, parent)

    def __reduce__(self):
        """
        Used for pickling.

        TESTS::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
            sage: F = AlgebraicClosureFiniteField(GF(2), 'z')
            sage: loads(dumps(F.gen(5))) == F.gen(5)
            True

        """
        return unpickle_AlgebraicClosureFiniteFieldElement, (self.parent(), self._level, self._repr_())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLE::

            sage: F = GF(3).algebraic_closure('z')
            sage: F._repr_()
            'Algebraic closure of Finite Field of size 3'

        """
        return self._value._repr_()

    def __cmp__(self, right):
        """
        Compare ``self`` with ``right``.

        EXAMPLE::

            sage: F = GF(3).algebraic_closure('z')
            sage: F.gen(2) == F.gen(3)
            False

        """
        x, y = self.parent()._coerce_2(self, right)
        return cmp(x, y)

    def _add_(self, right):
        """
        Return ``self`` + ``right``.

        EXAMPLE::

            sage: F = GF(3).algebraic_closure('z')
            sage: F.gen(2) + F.gen(3)
            z6^5 + 2*z6^4 + 2*z6^3 + z6^2 + 2*z6 + 1

        """
        F = self.parent()
        x, y = F._coerce_2(self, right)
        return self.__class__(F, x + y)

    def _sub_(self, right):
        """
        Return ``self`` - ``right``.

        EXAMPLE::

            sage: F = GF(3).algebraic_closure('z')
            sage: F.gen(2) - F.gen(3)
            z6^4 + 2*z6^3 + z6^2 + 2*z6

        """
        F = self.parent()
        x, y = F._coerce_2(self, right)
        return self.__class__(F, x - y)

    def _mul_(self, right):
        """
        Return ``self`` * ``right``.

        EXAMPLE::

            sage: F = GF(3).algebraic_closure('z')
            sage: F.gen(2) * F.gen(3)
            z6^5 + 2*z6^4 + z6^2 + 2

        """
        F = self.parent()
        x, y = F._coerce_2(self, right)
        return self.__class__(F, x * y)

    def _div_(self, right):
        """
        Return ``self`` / ``right``.

        EXAMPLE::

            sage: F = GF(3).algebraic_closure('z')
            sage: F.gen(2) / F.gen(3)
            z6^5 + 2*z6^4 + z6^3 + 1

        """
        F = self.parent()
        x, y = F._coerce_2(self, right)
        return self.__class__(F, x / y)

    def _change_level(self, n):
        """
        Return a representation of ``self`` as an element of the
        subfield of degree ``n`` of the parent, if possible.

        EXAMPLE::

            sage: F = GF(3).algebraic_closure('z')
            sage: z = F.gen(4)
            sage: (z^10)._change_level(6)
            2*z6^5 + 2*z6^3 + z6^2 + 2*z6 + 2

        """
        F = self.parent()
        l = self._level
        m = l.gcd(n)
        xl = self._value
        xm = F.inclusion(m, l).section()(xl)
        xn = F.inclusion(m, n)(xm)
        return self.__class__(F, xn)

    def is_square(self):
        """
        Return ``True`` if ``self`` is a square.

        This always returns ``True``.

        EXAMPLE::

            sage: F = GF(3).algebraic_closure('z')
            sage: F.gen(2).is_square()
            True

        """
        return True

    def sqrt(self):
        """
        Return a square root of ``self``.

        EXAMPLE::

            sage: F = GF(3).algebraic_closure('z')
            sage: F.gen(2).sqrt()
            z4^3 + z4 + 1

        """
        F = self.parent()
        x = self._value
        if x.is_square():
            return self.__class__(F, x.sqrt(extend=False))
        else:
            l = self._level
            x = F.inclusion(l, 2*l)(x)
            return self.__class__(F, x.sqrt(extend=False))

    def nth_root(self, n):
        r"""
        Return the ``n``-th root of ``self``.

        EXAMPLE::

            sage: F = GF(5).algebraic_closure('z')
            sage: t = F.gen(2) + 1
            sage: s = t.nth_root(15); s
            4*z6^5 + 3*z6^4 + 2*z6^3 + 2*z6^2 + 4
            sage: s**15 == t
            True
        """
        from sage.rings.integer import Integer
        F = self.parent()
        x = self._value
        n = Integer(n)
        l = self._level
        # in order to be smart we look for the smallest subfield that actually
        # contains the root.
        #TODO: it certainly can be faster
        for d in n.divisors():
            xx = F.inclusion(l, d*l)(x)
            try:
                y = xx.nth_root(n, extend=False)
            except ValueError:
                continue
            return self.__class__(F, y)

        raise ValueError("this should not happen!")

    def multiplicative_order(self):
        r"""
        Return the multiplicative order of ``self``.

        EXAMPLES::

            sage: K = GF(7).algebraic_closure('t')
            sage: K.gen(5).multiplicative_order()
            16806
            sage: (K.gen(1) + K.gen(2) + K.gen(3)).multiplicative_order()
            7353
        """
        return self._value.multiplicative_order()

    def pth_power(self, k=1):
        r"""
        Return the ``p^k``-th power of ``self``.

        EXAMPLES::

            sage: K = GF(13).algebraic_closure('t')
            sage: t3 = K.gen(3)
            sage: s = 1 + t3 + t3**2
            sage: s.pth_power()
            10*t3^2 + 6*t3
            sage: s.pth_power(2)
            2*t3^2 + 6*t3 + 11
            sage: s.pth_power(3)
            t3^2 + t3 + 1
        """
        return self._value.pth_power(k)

    def pth_root(self, k=1):
        r"""
        Return the ``p^k``-th root of ``self``.

        EXAMPLES::

            sage: K = GF(13).algebraic_closure('t')
            sage: t3 = K.gen(3)
            sage: s = 1 + t3 + t3**2
            sage: s.pth_root()
            2*t3^2 + 6*t3 + 11
            sage: s.pth_root(2)
            10*t3^2 + 6*t3
            sage: s.pth_root(3)
            t3^2 + t3 + 1
        """
        return self._value.pth_root(k)

    def as_finite_field_element(self, minimal=False):
        r"""
        Return ``self`` as a finite field element. Is ``minimal`` is set to
        ``True`` then returns the smallest subfield in which ``self`` is
        contained.

        EXAMPLES::

            sage: F = GF(3).algebraic_closure('t')
            sage: t = F.gen(5)
            sage: t.as_finite_field_element()
            (Finite Field in t5 of size 3^5,
             t5,
             Ring morphism:
              From: Finite Field in t5 of size 3^5
              To:   Algebraic closure of Finite Field of size 3
              Defn: t5 |--> t5)

        By default, the finite field returned is not minimal, but this behavior
        can be modified through the ``minimal`` option::

            sage: s = t+1-t
            sage: s.as_finite_field_element()[0]
            Finite Field in t5 of size 3^5
            sage: s.as_finite_field_element(minimal=True)[0]
            Finite Field of size 3
        """
        if not minimal:
            l = self._level
        else:
            l = self._value.minpoly().degree()

        F,phi = self.parent().subfield(l)
        return (F, F(self._value), phi)

def unpickle_AlgebraicClosureFiniteFieldElement(parent, level, x):
    """
    Unpickle an element `x` of an algebraic closure of a finite field.

    TEST::

        sage: F = GF(7).algebraic_closure('z')
        sage: loads(dumps(F.gen(2))) == F.gen(2)  # indirect doctest
        True

    """
    return parent.coerce(parent._subfield(level)(x))


class AlgebraicClosureFiniteField_generic(Field):
    """
    Algebraic closure of a finite field.
    """
    def __init__(self, base_ring, name, category=None):
        """
        TEST::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField_generic
            sage: F = AlgebraicClosureFiniteField_generic(GF(5), 'z')
            sage: F
            Algebraic closure of Finite Field of size 5

        """
        Field.__init__(self, base_ring=base_ring, names=name,
                       normalize=False, category=category)

    def cardinality(self):
        r"""
        Return infinity.
        """
        from sage.rings.infinity import Infinity
        return Infinity

    def __getstate__(self):
        """
        Used for pickling.

        TEST::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField_generic
            sage: F = AlgebraicClosureFiniteField_generic(GF(5), 'z')
            sage: F.__getstate__() is None
            True

        """
        pass

    def __setstate__(self, state):
        """
        Used for pickling.

        TEST::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField_generic
            sage: F = AlgebraicClosureFiniteField_generic(GF(5), 'z')
            sage: F.__setstate__(None)

        """
        pass

    def __reduce__(self):
        """
        Used for pickling.

        TEST::

            sage: F = GF(5).algebraic_closure('z')
            sage: loads(dumps(F)) == F
            True

        """
        return (AlgebraicClosureFiniteField,
                (self.base_ring(), self.variable_name(), self.category()),
                self.__getstate__())

    def __cmp__(self, other):
        """
        Compare ``self`` with ``other``.

        TEST::

            sage: F3 = GF(3).algebraic_closure('z')
            sage: F3 == F3
            True
            sage: F5 = GF(5).algebraic_closure('z')
            sage: F3 == F5
            False

        """
        if self is other:
            return 0
        c = cmp(type(self), type(other))
        if c != 0:
            return c
        return cmp((self.base_ring(), self.variable_name(), self.category()),
                   (other.base_ring(), other.variable_name(), other.category()))

    def characteristic(self):
        """
        Return the characteristic of ``self``.

        EXAMPLE::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
            sage: p = next_prime(1000)
            sage: F = AlgebraicClosureFiniteField(GF(p), 'z')
            sage: F.characteristic() == p
            True

        """
        return self.base_ring().characteristic()

    Element = AlgebraicClosureFiniteFieldElement

    def _element_constructor_(self, x):
        """
        Construct an element of ``self``.

        TEST::

            sage: F = GF(5).algebraic_closure('z')
            sage: type(F(3))
            <class 'sage.rings.algebraic_closure_finite_field.AlgebraicClosureFiniteField_pseudo_conway_with_category.element_class'>

        """
        if isinstance(x, self.element_class) and x.parent() is self:
            return x
        else:
            return self.element_class(self, x)

    def _coerce_map_from_(self, other):
        """
        Return ``True`` if elements of ``other`` can be coerced into
        ``self``.

        EXAMPLE::

            sage: F = GF(7).algebraic_closure('z')
            sage: F.has_coerce_map_from(Integers())
            True

        """
        if other is self:
            return True
        elif is_FiniteField(other) and self._subfield(other.degree()) is other:
            return True
        elif self._subfield(1).has_coerce_map_from(other):
            return True

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLE::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
            sage: F = AlgebraicClosureFiniteField(GF(5), 'z')
            sage: F._repr_()
            'Algebraic closure of Finite Field of size 5'

        """
        return 'Algebraic closure of %s' % self.base_ring()

    def _coerce_2(self, x, y):
        """
        Coerce `x` and `y` to a common subfield of ``self``.

        TEST::

            sage: F = GF(3).algebraic_closure('z')
            sage: x, y = F._coerce_2(F.gen(2), F.gen(3))
            sage: x.parent()
            Finite Field in z6 of size 3^6
            sage: y.parent()
            Finite Field in z6 of size 3^6

        """
        n = x._level.lcm(y._level)
        mx = self.inclusion(x._level, n)
        my = self.inclusion(y._level, n)
        return mx(x._value), my(y._value)

    def _get_polynomial(self, n):
        """
        Return the polynomial defining the unique subfield of degree
        `n` of ``self``.

        This must be implemented by subclasses.

        EXAMPLE::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField_generic
            sage: F = AlgebraicClosureFiniteField_generic(GF(5), 'z')
            sage: F._get_polynomial(1)
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        raise NotImplementedError

    def _get_im_gen(self, m, n):
        """
        Return the image of ``self.gen(m)`` under the canonical
        inclusion into ``self.subfield(n)``.

        This must be implemented by subclasses.

        EXAMPLE::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField_generic
            sage: F = AlgebraicClosureFiniteField_generic(GF(5), 'z')
            sage: F._get_im_gen(2, 4)
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        raise NotImplementedError

    def _subfield(self, n):
        """
        Return the unique subfield of degree `n` of ``self``.

        EXAMPLE::

            sage: F = GF(3).algebraic_closure('z')
            sage: F._subfield(4)
            Finite Field in z4 of size 3^4

        """
        if n == 1:
            return self.base_ring()
        else:
            from sage.rings.finite_rings.constructor import FiniteField
            return FiniteField(self.base_ring().cardinality() ** n,
                               name=self.variable_name() + str(n),
                               modulus=self._get_polynomial(n))

    def subfield(self, n):
        """
        Return the unique subfield of degree `n` of ``self``
        together with its canonical embedding into ``self``.

        EXAMPLE::

            sage: F = GF(3).algebraic_closure('z')
            sage: F.subfield(1)
            (Finite Field of size 3,
             Ring morphism:
               From: Finite Field of size 3
               To:   Algebraic closure of Finite Field of size 3
               Defn: 1 |--> 1)
            sage: F.subfield(4)
            (Finite Field in z4 of size 3^4,
             Ring morphism:
               From: Finite Field in z4 of size 3^4
               To:   Algebraic closure of Finite Field of size 3
               Defn: z4 |--> z4)

        """
        Fn = self._subfield(n)
        return Fn, Fn.hom((self.gen(n),))

    def inclusion(self, m, n):
        """
        Return the canonical inclusion map from the subfield
        of degree `m` to the subfield of degree `n`.

        EXAMPLE::

            sage: F = GF(3).algebraic_closure('z')
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
            return self._subfield(m).hom((self._get_im_gen(m, n),))
        else:
            raise ValueError("subfield of degree %s not contained in subfield of degree %s" % (m, n))

    def ngens(self):
        """
        Return the number of generators of ``self``, which is
        infinity.

        EXAMPLE::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
            sage: AlgebraicClosureFiniteField(GF(5), 'z').ngens()
            +Infinity

        """
        from sage.rings.infinity import Infinity
        return Infinity

    def gen(self, n):
        """
        Return the `n`-th generator of ``self``.

        EXAMPLE::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
            sage: F = AlgebraicClosureFiniteField(GF(5), 'z')
            sage: F.gen(2)
            z2

        """
        F = self._subfield(n)
        return self(F.gen())

    def _first_ngens(self, n):
        """
        Return the first `n` generators of ``self``.

        EXAMPLE::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
            sage: F = AlgebraicClosureFiniteField(GF(5), 'z')
            sage: F._first_ngens(3)
            (1, z2, z3)

        """
        return tuple([self.gen(i + 1) for i in xrange(n)])

    def algebraic_closure(self):
        """
        Return an algebraic closure of ``self``.

        This always returns ``self``.

        EXAMPLE::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
            sage: F = AlgebraicClosureFiniteField(GF(5), 'z')
            sage: F.algebraic_closure() is F
            True

        """
        return self


class AlgebraicClosureFiniteField_pseudo_conway(AlgebraicClosureFiniteField_generic):
    """
    Algebraic closure of a finite field, constructed using
    pseudo-Conway polynomials.

    EXAMPLE::

        sage: F = GF(5).algebraic_closure('z')
        sage: type(F)
        <class 'sage.rings.algebraic_closure_finite_field.AlgebraicClosureFiniteField_pseudo_conway_with_category'>

    """
    def __init__(self, base_ring, name, category=None):
        """
        TEST::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField_pseudo_conway
            sage: F = AlgebraicClosureFiniteField_pseudo_conway(GF(3), 'z')
            sage: F
            Algebraic closure of Finite Field of size 3

        """
        if not (is_FiniteField(base_ring) and base_ring.is_prime_field()):
            raise NotImplementedError
        from sage.rings.finite_rings.conway_polynomials import PseudoConwayLattice
        self._pseudo_conway_lattice = PseudoConwayLattice(base_ring.characteristic())
        AlgebraicClosureFiniteField_generic.__init__(self, base_ring, name, category)

    def __getstate__(self):
        """
        Used for pickling.

        TEST::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField_pseudo_conway
            sage: F = AlgebraicClosureFiniteField_pseudo_conway(GF(5), 'z')
            sage: F.__getstate__()
            <class 'sage.rings.finite_rings.conway_polynomials.PseudoConwayLattice'>

        """
        return self._pseudo_conway_lattice

    def __setstate__(self, state):
        """
        Used for pickling.

        TEST::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField_pseudo_conway
            sage: F = AlgebraicClosureFiniteField_pseudo_conway(GF(5), 'z')
            sage: F.__setstate__(F.__getstate__())

        """
        self._pseudo_conway_lattice = state

    def __cmp__(self, other):
        """
        Compare ``self`` with ``other``.

        TEST::

            sage: F3 = GF(3).algebraic_closure('z')
            sage: F3 == F3
            True
            sage: F5 = GF(5).algebraic_closure('z')
            sage: F3 == F5
            False

        """
        c = AlgebraicClosureFiniteField_generic.__cmp__(self, other)
        if c != 0:
            return c
        return cmp(self._pseudo_conway_lattice, other._pseudo_conway_lattice)

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
        return self._pseudo_conway_lattice.polynomial(n)

    def _get_im_gen(self, m, n):
        """
        Return the image of ``self.gen(m)`` under the canonical
        inclusion into ``self.subfield(n)``.

        EXAMPLE::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField_pseudo_conway
            sage: F = AlgebraicClosureFiniteField_pseudo_conway(GF(5), 'z')
            sage: F._get_im_gen(2, 4)
            z4^3 + z4^2 + z4 + 3

        """
        p = self.characteristic()
        return self._subfield(n).gen() ** ((p**n - 1)//(p**m - 1))


class AlgebraicClosureFiniteFieldFactory(UniqueFactory):
    """
    Factory for constructing algebraic closures of finite fields.

    EXAMPLE::

        sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
        sage: F = GF(2).algebraic_closure('z')
        sage: F1 = AlgebraicClosureFiniteField(GF(2), 'z')
        sage: F1 is F
        True
        sage: loads(dumps(F)) is F
        True

    """
    def create_key(self, base_ring, name, category=None, implementation=None, **kwds):
        """
        TEST::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
            sage: AlgebraicClosureFiniteField.create_key(GF(3), 'z')
            (Finite Field of size 3, 'z', Category of fields, 'pseudo_conway')

        """
        if category is None:
            from sage.categories.fields import Fields
            category = Fields()
        if implementation is None:
            implementation = 'pseudo_conway'
        return (base_ring, name, category, implementation)

    def create_object(self, version, key, **kwds):
        """
        TEST::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
            sage: key = (GF(3), 'z', Fields(), 'pseudo_conway')
            sage: AlgebraicClosureFiniteField.create_object(0, key)
            Algebraic closure of Finite Field of size 3

        """
        base_ring, name, category, implementation = key
        if implementation == 'pseudo_conway':
            return AlgebraicClosureFiniteField_pseudo_conway(base_ring, name, category, **kwds)
        else:
            raise ValueError('unknown implementation for algebraic closure of finite field: %s'
                             % implementation)


AlgebraicClosureFiniteField = AlgebraicClosureFiniteFieldFactory('sage.rings.algebraic_closure_finite_field.AlgebraicClosureFiniteField')
