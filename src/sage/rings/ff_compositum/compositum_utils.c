#include "compositum.h"
#include "nmod_poly_extra.h"

void _compositum_iP(nmod_poly_t res, const nmod_poly_t P) {
  nmod_poly_t dP;
  nmod_poly_init(dP, P->mod.n);
  nmod_poly_derivative(dP, P);
  if (!nmod_poly_invmod(res, dP, P)) {
    printf("Exception (_compositum_iP). Ramified extension.\n");
    abort();
  }
  nmod_poly_clear(dP);
}
