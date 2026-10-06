#ifndef SEMSOLVE2_INTEGRALS_H
#define SEMSOLVE2_INTEGRALS_H

#include <basis.h>
#include <potential.h>

// double	 potential_energy(  );
double	 inner_prod( basis_t *, long int, long int );				///< The inner product of two basis functions on the simplex, (1,2)
double	 potential_inner_prod( basis_t *, potential_t *, long int, long int );	///< The potential inner product
double	 kinetic_inner_prod( basis_t *, long int, long int, double * );		///< The inner product a(1,2) = (grad 1, grad 2)

#endif
