
#include "stdio.h"
#include "takapack.h"

int main() {


/*----------------------------
    2  3  0  0  0        8
    3  0  4  0  6       45
A =   0 -1 -3  2  0    b =-3
    0  0  1  0  0        3
    0  4  2  0  1       19
----------------------------*/

    //(umfpackは complessed columnだけど、実装の都合上 Compressed Rowを採用 ())
    // umfpackと併用するときは umfpack_solve の第一引数に UMFPACK_At を食わせて転置すればOK
    //compressed row form of A 
    const int N = 5;
    int    Ap [ ] = { 0,        2,          5,            8,   9,      12} ;
    int    Ai [ ] = { 0,  1,    0, 2, 4,   1,  2, 3,    2,   1, 2, 4   } ;
    double Ax [ ] = { 2., 3.,   3, 4, 6,  -1, -3, 2,    1,   4, 2, 1   } ;

    double b  [ ] = { 8, 45, -3, 3, 19} ;

    double x1 [5] = {0,0,0,0,0};
    
    int    *LUp, *LUi, *LU_rowflip;
    double *LUx;
    
    takapack_LU_factorization( N, Ap, Ai, Ax, LUp, LUi, LUx, LU_rowflip);     
    takapack_LU_solve        ( N, LUp, LUi, LUx, LU_rowflip, b, x1 ); 
    fprintf( stderr, "%f %f %f %f %f\n", x1[0], x1[1], x1[2], x1[3], x1[4]);//適当にprintfしてください...
    takapack_LU_free( LUp, LUi, LUx, LU_rowflip);
}