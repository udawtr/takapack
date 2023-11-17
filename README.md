# takapack

http://takashiijiri.com/study/miscs/takapack.html より、

　
- 疎な連立方程式を解く
  -  与える行列は, umfpackと同様の comressed row form (umfpackではcompressed column formだった)
  - CG法による　　連立方程式 Ax = b の 求解をサポート
  - LU分解による 連立方程式 Ax = b の 求解をサポート
　　
- NYSL（Version 0.9982）ライセンスなので、好きに使ってください（大学のレポートとかの参考になれば良いなと思っています.）．

## 使い方
　
### 1) 解きたい連立方程式 Ax = b の行列Aをcompressed row formで作る。
　　　　配列bもつくる。

```
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
```

### 2) CG法で解くなら CG法の関数を呼ぶ

```
    double x1 [5] = {0,0,0,0,0};
    takapack_CG_sparse_solve(N, Ap, Ai, Ax, b, x2, 0.00000001);
    fprintf( stderr, "%f %f %f %f %f\n", x2[0], x2[1], x2[2], x2[3], x2[4]);}
```

### 3)LU分解で解くなら i)LU分解関数 ii)解く関数 iii)解放する関数 を順に呼ぶ

```
    double x1 [5] = {0,0,0,0,0};
    
    int    *LUp, *LUi, *LU_rowflip;
    double *LUx;
    
    takapack_LU_factorization( N, Ap, Ai, Ax, LUp, LUi, LUx, LU_rowflip);     
    takapack_LU_solve        ( N, LUp, LUi, LUx, LU_rowflip, b, x1 ); 
    fprintf( stderr, "%f %f %f %f %f\n", x1[0], x1[1], x1[2], x1[3], x1[4]);//適当にprintfしてください...
    takapack_LU_free( LUp, LUi, LUx, LU_rowflip);
``` 

## 雑な動作確認方法

```
$ gcc main.cpp takapack.cpp
$ ./a.out
1.000000 2.000000 3.000000 4.000000 5.000000
```