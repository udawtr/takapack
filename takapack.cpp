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




//#include "StdAfx.h"
#include "stdio.h"
#include "string.h"
#include "takapack.h"

#include "math.h"
#include <list>
#include <algorithm>
using namespace std;

#ifdef _DEBUG
#define new DEBUG_NEW
#endif





void takapack_test()
{

/*----------------------------
	2  3  0  0  0        8
	3  0  4  0  6       45
A =	0 -1 -3  2  0    b =-3
	0  0  1  0  0        3
	0  4  2  0  1       19
----------------------------*/

	//(umfpackは complessed columnだけど、実装の都合上 Compressed Rowを採用 ())
	// umfpackと併用するときは umfpack_solve の第一引数に UMFPACK_At を食わせて転置すればOK
	//compressed row form of A 
	const int N = 5;
	int    Ap [ ] = { 0,        2,         5,           8,   9,      12} ;
	int    Ai [ ] = { 0,  1,    0, 2, 4,   1,  2, 3,    2,   1, 2, 4   } ;
	double Ax [ ] = { 2., 3.,   3, 4, 6,  -1, -3, 2,    1,   4, 2, 1   } ;

	takapack_traceMat(N,Ap,Ai,Ax);

	double b  [ ] = { 8, 45, -3, 3, 19} ;
	double x1 [5] = {0,0,0,0,0};
	double x2 [5] = {0,0,0,0,0};;
	{
		int    *LUp, *LUi, *LU_rowflip;
		double *LUx;
	
		takapack_LU_factorization( N, Ap, Ai, Ax, LUp, LUi, LUx, LU_rowflip); 
	
		takapack_traceMat(N,LUp,LUi,LUx);
	
		takapack_LU_solve        ( N, LUp, LUi, LUx, LU_rowflip, b, x1 ); 
		fprintf( stderr, "%f %f %f %f %f\n", x1[0], x1[1], x1[2], x1[3], x1[4]);

		takapack_LU_free( LUp, LUi, LUx, LU_rowflip);
	}
	{
		takapack_CG_sparse_solve(N, Ap, Ai, Ax, b, x2, 0.00000001);
		fprintf( stderr, "%f %f %f %f %f\n", x2[0], x2[1], x2[2], x2[3], x2[4]);

	}
}





void takapack_CG_dense_solve( const int N, double** A, const double* b, double* result)
{
	size_t double_N = sizeof( double ) * N;

	double *cg_r = new double[ N ];
	double *cg_d = new double[ N ];
	double *cg_q = new double[ N ];

	memset( result, 0, double_N );
	int iteration = 0;

	//r = b - Ax (  A is symmetric!!  )
	memcpy( cg_r, b, double_N );
	for( int i = 0; i < N; ++i) 
	for( int j = 0; j < N; ++j) cg_r[ i ] -= A[i][j] * result[ j ];
		
	//d = r
	memcpy( cg_d, cg_r, double_N );

	//deltaNew = r_t * r
	double deltaNew = 0;
	for( int i = 0; i < N; ++i) deltaNew += cg_r[i] * cg_r[i];

	while( iteration < N && deltaNew > 0.000000001)
	{
		//q = Ad
		memset( cg_q, 0, double_N );
		for( int i = 0; i < N; ++i) 
		for( int j = 0; j < N; ++j) cg_q[i] += A[i][j] * cg_d[j];

		//alpha = deltaNew / (d_t * q )
		double alpha = 0;
		for( int i = 0; i < N ; ++i ) alpha += cg_d[i] * cg_q[i];
		alpha = deltaNew / alpha;

		// x = x + alpha * d
		for(int i = 0; i < N; ++i) result[i] += alpha * cg_d[i];

		if( iteration % 30 == 0 ){//r = b - Ax
			memcpy( cg_r, b, double_N );
			for( int i = 0; i < N; ++i) 
			for( int j = 0; j < N; ++j) cg_r[i] -= A[i][j] * result[j];

		}else{//r = r - alpha * q
			for( int i = 0; i < N; ++i) cg_r[i] -= alpha * cg_q[i];		
		}

		//deltaOld = deltaNew
		double deltaOld = deltaNew;

		//deltaNew = r*r
		deltaNew = 0;
		for( int i = 0; i < N; ++i) deltaNew += cg_r[i] * cg_r[i];

		//beta = deltaNew / deltaOld
		double beta = deltaNew / deltaOld;
	
		//d = r + beta + d
		for( int i = 0; i < N; ++i) cg_d[i] = cg_r[i] + beta * cg_d[i];
		++iteration;
	}

	delete[] cg_r;
	delete[] cg_d;
	delete[] cg_q;
}




bool takapack_CG_sparse_solve( const int N, 
							   const int    *Ap, const int* Ai, const double* Ax, 
							   const double *b , 
							   double* result, double  threshold)
{
	size_t double_N = sizeof( double ) * N;	
	double *m_cg_r = new double[ N ];
	double *m_cg_d = new double[ N ];
	double *m_cg_q = new double[ N ];
	const int    *m_Ap_c = Ap;
	const int    *m_Ai_c = Ai;
	const double *m_Ax_c = Ax;

	int iteration = 0;

	//r = b - Ax (  A is symmetric!!  )
	memcpy( m_cg_r, b, double_N );
	memset( result, 0, double_N );

	for( int i = 0        ; i < N          ; ++i)
	for( int j = m_Ap_c[i]; j < m_Ap_c[i+1]; ++j)
			m_cg_r[ i ] -= m_Ax_c[ j ] * result[ m_Ai_c[j] ];
		
	//d = r
	memcpy( m_cg_d, m_cg_r, double_N );

	//deltaNew = r_t * r
	double deltaNew = 0;
	for( int i = 0; i < N; ++i) deltaNew += m_cg_r[i] * m_cg_r[i];


	//deltaZero = deltaNew
	double deltaZero = deltaNew;

	while( ( iteration < N && deltaNew > threshold) )
	{
		//q = Ad
		memset( m_cg_q, 0, double_N );
		for( int i = 0        ; i < N          ; ++i)
		for( int j = m_Ap_c[i]; j < m_Ap_c[i+1]; ++j)
			m_cg_q[ i ] += m_Ax_c[ j ] * m_cg_d[ m_Ai_c[j] ];

	
		//alpha = deltaNew / (d_t * q )
		double alpha = 0;
		for( int i = 0; i < N ; ++i ) alpha += m_cg_d[i] * m_cg_q[i];
		alpha = deltaNew / alpha;


		// x = x + alpha * d
		for(int i = 0; i < N; ++i) result[i] += alpha * m_cg_d[i];

		if( iteration % 20 == 0 ){
			//r = b - Ax
			memcpy( m_cg_r, b, double_N );
			for( int i = 0        ; i < N          ; ++i)
			for( int j = m_Ap_c[i]; j < m_Ap_c[i+1]; ++j)
					m_cg_r[ i ] -= m_Ax_c[ j ] * result[ m_Ai_c[j] ];

		}else{	
			//r = r - alpha * q
			for( int i = 0; i < N; ++i) m_cg_r[i] -= alpha * m_cg_q[i];		
		}

		//deltaOld = deltaNew
		double deltaOld = deltaNew;

		//deltaNew = r*r
		deltaNew = 0;
		for( int i = 0; i < N; ++i) deltaNew += m_cg_r[i] * m_cg_r[i];

		//beta = deltaNew / deltaOld
		double beta = deltaNew / deltaOld;
	
		//d = r + beta + d
		for( int i = 0; i < N; ++i) m_cg_d[i] = m_cg_r[i] + beta * m_cg_d[i];
		
		++iteration;
	}

	delete[] m_cg_r ;
	delete[] m_cg_d ;
	delete[] m_cg_q ;
	return true;
}











void takapack_traceMat( int N, const int *Ap , const int *Ai , const double *Ax )
{
	fprintf( stderr, "\ntakapack_tracemat\n");

	int idx = 0;
	for( int y = 0; y<N; ++y)
	{
		for( int x = 0; x<N; ++x)
		{
			if( x == Ai[idx] ) { fprintf( stderr, "%.4f  ", Ax[idx]); ++idx; }
			else               { fprintf( stderr, "*.****  " );                  }
		}
		fprintf( stderr, "\n");
	}
}







void takapack_LU_factorization( const int N, const int *Ap , const int *Ai , const double *Ax, int* &LUp, int* &LUi, double* &LUx, int* &LU_rowFlip)
{
	typedef list<pair<int,double>>::iterator RowItr;

	//compressed row form用の LU分解
	LU_rowFlip = new int[ N ];
	for( int i=0;i<N; ++i) LU_rowFlip[i] = i;

	int nunZeroEntryNum  = 0;
	list<pair<int,double>> *Row = new list<pair<int,double>> [ N ];
	double   *Myi_tmp = new double[N]; //(これに変えて Axを縦方向につぶしたcompressed column行列を置いてもいいかも)
	double   *Myi     = new double[N]; //(これに変えて Axを縦方向につぶしたcompressed column行列を置いてもいいかも)

	for( int I = 0; I < N; ++I)// for each column(各列(縦)について)
	{
		//1)-----------I列 (縦)を全て集める : Myi -------------------------------------------------------------
#pragma omp parallel for
		for( int y=0; y < N; ++y){ 
			Myi_tmp[y] = 0;
			for( int kk=Ap[y]; kk < Ap[y+1] ; ++kk) if( Ai[kk] == I ){ Myi_tmp[y] = Ax[kk]; break;}
		}
		//row flipを適用
		for( int y=0; y<N; ++y) Myi[y] = Myi_tmp[ LU_rowFlip[y] ];

		//2)-----------Myiの値を計算 (この時Aiにアクセスしてはダメ (計算結果は入っていないから))---------------
		for( int y=0; y <=I; ++y) {
			double v = Myi[y]; 
			for( RowItr it = Row[y].begin(); it != Row[y].end() && it->first < y; ++it) v -= it->second * Myi[ it->first ];
			Myi[y] = v;
		}

		double big = fabs( Myi[I] ), tmp;
		int    piv = I;

		for( int y=I+1; y < N; ++y) {
			double v = Myi[y];
			for( RowItr it = Row[y].begin(); it != Row[y].end() && it->first < I; ++it) v -= it->second * Myi[ it->first ];
			Myi[y] = v;
			if( (tmp = fabs( Myi[y] )) > big ){ big = tmp;  piv = y;  }
		}

		//3)------------係数掛ける-----------------------
		double pivV = Myi[piv],  coef = 1.0 / pivV;
#pragma omp parallel for
		for( int y = I; y < N; ++y) Myi[y] *= coef;
		Myi[piv] = pivV;

		//4)------------軸選択(I<-->piv)----------------
		if( piv != I ) {
			std::swap( Myi       [I], Myi       [piv] );
			std::swap( Row       [I], Row       [piv] ); //listのflipは遅いかも
			std::swap( LU_rowFlip[I], LU_rowFlip[piv] );
		}

		//5) -----------値を代入 LU <-- Myi-------------
		for( int y = 0; y < N; ++y) if( Myi[y] != 0.0 ) {
			++nunZeroEntryNum;
			Row[y].push_back( make_pair( I, Myi[y] ) );
		}
	}

	//最後に Row を Compressed row formにまとめる
	LUp = new int   [ N+1 ];
	LUi = new int   [ nunZeroEntryNum ];
	LUx = new double[ nunZeroEntryNum ];
	
	int index = 0;
	for( int i = 0; i < N; ++i)
	{
		LUp[i] = index;
		for( list<pair<int,double>>::iterator it = Row[i].begin(); it != Row[i].end(); ++it) 
		{
			LUi[index] = it->first;
			LUx[index] = it->second;
			++index;
		}
	}
	LUp[ N ] = index;

	delete[] Myi ;
	delete[] Myi_tmp;
	delete[] Row ;
}




//solve Ax = b --> LUx = b  --> L(Ux=y)=b
void takapack_LU_solve( const int N, const int *LUp, const int *LUi, const double *LUx, const int *LU_rowFlip, const double *b, double *res)
{
	double *f_B = new double[N];//fliped B
	double *f_Y = new double[N];//fliped Y
	for( int i=0; i<N;++i) f_B[i] = b[ LU_rowFlip[i] ];

	//前進代入  L * f_Y = flipB --> f_Y (Lの対角成分は1であることに注意)
	for( int y = 0; y < N; ++y){
		double val = f_B[y];
		for( int kk=LUp[y]; kk < LUp[y+1] && LUi[kk] < y; ++kk) val -= LUx[ kk ] * f_Y[ LUi[kk] ];
		f_Y[y] = val;
	}

	//後退代入 solve U x = flipY (L-1 * flipB) 
	for( int y = N-1; y>=0; --y){
		int centerIdx = 0;
		for( int kk = LUp[y] ; kk < LUp[y+1]; ++kk) if( LUi[kk] == y) { centerIdx = kk; break;} //対角部分に当たるindexを取得
		
		double val = f_Y[y];
		for( int kk = centerIdx + 1 ; kk < LUp[y+1]; ++kk) val -= LUx[kk] * res[ LUi[kk] ];
		res[y] = val / LUx[centerIdx];
	}
	delete[] f_B;
	delete[] f_Y;
}

void takapack_LU_free ( int *LUp, int *LUi, double *LUx, int *LU_rowFlip )
{
	delete[] LUp;
	delete[] LUi;
	delete[] LUx;
	delete[] LU_rowFlip;
}