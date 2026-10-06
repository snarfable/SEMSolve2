#ifndef SEMSOLVE2_LAPACK_H
#define SEMSOLVE2_LAPACK_H

extern "C" void	dgesv_( int *, int *, double *, int *, int *, double *, int *, int * );
extern "C" void	dsysv_( char *, int *, int *, double *, int *, int *, double *, int *, double *, int *, int * );
extern "C" void	dgemm_( char *, char *, int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int * );
extern "C" void	dsymm_( char *, char *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int * );

#endif
