#ifndef SEMSOLVE2_DLINALG_ASM_H
#define SEMSOLVE2_DLINALG_ASM_H

#include <dlinalg.h>

#ifdef __cplusplus
extern "C" {
#endif
	double __asm_dot_prod_dim2_amd64( dvector_t *, dvector_t * );
	double __asm_dot_prod_dim2_x86( dvector_t *, dvector_t * );
	double __asm_dot_prod_amd64( dvector_t *, dvector_t * );
	double __asm_dot_prod_x86( dvector_t *, dvector_t * );
	void __asm_matrix_mult_2x2_amd64( dmatrix_t *, dmatrix_t * );		///< Use this to optimize LU decomposition
	void __asm_matrix_mult_2x2_x86( dmatrix_t *, dmatrix_t *, dmatrix_t * );	///< Same, except for x86 machines
#ifdef __cplusplus
}
#endif

#endif
