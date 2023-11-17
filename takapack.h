#pragma once

/* ----------------------------------------------------
takapack.h is written by Takashi Ijiri @ RIKEN
This libraly support computing 
 + Dense  linear system (by CG method)
 + Sparse linear system (by CG method)
 + Sparse linear system (by LU factorization)

 This codes is distributed with NYSL (Version 0.9982) license. 

 免責事項   : 本ソースコードによって起きたいかなる障害についても、著者は一切の責任を負いません。
 disclaimer : The author (Takashi Ijiri) is not to be held responsible for any problems caused by this software.
 ----------------------------------------------------*/


void takapack_test();

//CG method (for sparse matrix)
bool takapack_CG_sparse_solve( const int N, 
								const int    *Ap, const int* Ai, const double* Ax, 
								const double *b , 
								double* result, double  threshold);
//CG method (for dense matrix)
void takapack_CG_dense_solve( const int N, double** A, const double* b, double* result);

//LU method
void takapack_traceMat( int N, const int *Ap , const int *Ai , const double *Ax );
void takapack_LU_factorization( const int N, const int *Ap , const int *Ai , const double *Ax, int* &LUp, int* &LUi, double* &LUx, int* &LU_rowFlip);
void takapack_LU_solve        ( const int N, const int *LUp, const int *LUi, const double *LUx, const int *LU_rowFlip, const double *b, double *res);
void takapack_LU_free         (              int *LUp, int *LUi, double *LUx, int *LU_rowFlip );
