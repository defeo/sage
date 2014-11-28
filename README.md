# Fast composita of finite fields in Sage

This is a fork of the
[Sage mathematical software](http://github.com/sagemath/sage/)
containing experimental code developed for the research project
[*"Fast arithmetic for the algebraic closure of finite fields"*](http://github.com/defeo/ff_compositum/).

All the added code is contained in the directory
`src/sage/rings/ff_compositum`, its contents are released in the
public domain. For the copying conditions of Sage, see `COPYING.md`.

## How to use

This code is a proof of concent. It is very user unfriendly and
buggy. If you are interested in using it for your research project,
please [drop me a line](http://defeo.lu), and I will try and help you.

If nevertheless you want to try it:

1. Clone this branch.
2. Compile Sage
   [as usual](http://www.sagemath.org/doc/installation/source.html)
   (this usually means just `make`).
3. Launch sage and `from sage.rings.ff_compositum.all import FFDesc,
   FFElt`.

See the examples in <http://github.com/defeo/ff_compositum/>, have
fun.
