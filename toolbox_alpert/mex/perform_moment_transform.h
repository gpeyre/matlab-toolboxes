// header file for perform_moment_transform


#include <math.h>
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gram_schmidt.h"
#include "matrix.h"


extern double* v;
extern double* w;
extern double* pos;
extern double* monomials;
extern double** part;
extern int* G;
extern int* si;

#define posx_(i) pos[ 0+(i)*2 ]
#define posy_(i) pos[ 1+(i)*2 ]
#define monomialsx_(i) ((int)monomials[ 0+(i)*2 ])
#define monomialsy_(i) ((int)monomials[ 1+(i)*2 ])

inline
int fncompare(const void * elem1, const void * elem2 )
{
	int n1 = *((int*) elem1)-1;
	int n2 = *((int*) elem2)-1;
	if( posx_(n1)<posx_(n2) )
		return -1;
	else if ( posx_(n1)==posx_(n2) )
		return 0;
	else
		return 1;
}

inline 
int construct_part(const double* pos, double** &part, int* &si, int* &G, int n, int ptmax)
{
	int j = 0;
	double m = (double) n;
	while( m>=ptmax )
	{
		m /= 2; j++;
	}
	j--;
	GW_ASSERT(j>=0);
	int nj = 1<<j;		// 2^j : number of groups
	int R = n % nj;		// mod(n,2^j);
	int L = (n-R)/nj;

	part = new double*[nj];
	si = new int[nj+1];
	G = new int[nj];

	// the array 0:nj-1 to reorder
	int* I = new int[n];
	for( int i=0; i<n; ++i )
		I[i] = i+1;

	// [x,I] = sort(pos(1,:));
	qsort(I, n, sizeof(int), fncompare);

	si[0] = 0;
	for( int k=0; k<nj; ++k )
	{
		if( k<nj-R )
			G[k] = L;
		else
			G[k] = L+1;
		part[k] = new double[(int) G[k]];
		// part{k} = I(si[k]:si[k]+G(k)-1);
		for( int i=0; i<G[k]; ++i )
			part[k][i] = I[si[k]+i];
		si[k+1] = si[k] + G[k]; 
	}

	GW_DELETEARRAY(I);

	return nj;
}

inline 
double monomial(double x, double y, int nx, int ny)
{
	double r = 1;
	for( int i=0; i<nx; ++i )
		r = r*x;
	for( int i=0; i<ny; ++i )
		r = r*y;
	return r;
}

inline
void perform_moment_transform(int n, int P, int k, int dir)
{
	int k2 = k >> 1;
	int destroy_data = 0;
	if( part==NULL )
	{
		// construct part
		P = construct_part(pos, part, si, G, n, k2);
		destroy_data = 1;
	}

	int J = log2(P) + 1;

	matrix_list* Uj = new matrix_list[J];

	// initialization of the moments matrix
	matrix_list Mi(P);
	matrix_list& Ui = Uj[0]; Ui.SetSize(P);
	for( int i=0; i<P; ++i ) 
	{
		double* seli = part[i];
		// current matrix is of size G[i] x 2*k2
		matrix& M = *(new matrix(G[i],k));
		for( int ii=0; ii<G[i]; ++ii )
		{
			int num = (int) seli[ii]-1;
			for( int jj=0; jj<2*k2; ++jj )
			{
				M(ii,jj) = monomial( posx_(num),posy_(num),monomialsx_(jj),monomialsy_(jj) );
			}
		}
		Mi(i) = &M;		// assign to list
		// orthogonalize
		matrix& Q = *(new matrix(G[i],G[i]));
		M.PerformGramSchmidt(Q);
		Ui(i) = &Q;		// assign to list
	}

	// computation of Uj[j] for j>0
	for( int j=1; j<J; ++j )
	{
		// at this scale, we have nj = P/2^(j) groups
		int nj = P >> j;	// P/( 2^j )
		int mj = nj*k;    // total length of the blocks

		matrix_list& Ui = Uj[j-1];
		matrix_list& Uii = Uj[j];	Uii.SetSize(nj);
		// update each sub matrix
		for( int i=0; i<nj; ++i )
		{
			// MM is a (k x k) matrix
			matrix& MM = *(new matrix(k,k));
			// the two orthogonal matrix from coarser level
			matrix& Ui1 = *Ui(2*i);
			matrix& Ui2 = *Ui(2*i+1);
			// the two moment matrix from coarser level
			matrix& Mi1 = *Mi(2*i);
			matrix& Mi2 = *Mi(2*i+1);
			// compute M via : 
			//		MM(0:k2-1,:) = Ui1'(0:k2-1,:)*Mi1
			//		MM(k2:k,:)   = Ui2'(0:k2-1,:)*Mi2
			// recall that Ui1->G[2*i]xG[2*i], Ui2->G[2*i+1]xG[2*i+1],
			//			   Mi1->G[2*i]x k    , Mi2->G[2*i+1]x k
			for( int s=0; s<k2; ++s )
			{
				// MM(s,a) = sum_b{ Ui1(b,s)*Mi1(b,a) }
				// MM(s+k2,a) = sum_b{ Ui2(b,s)*Mi2(b,a) }
				for( int a=0; a<k; ++a )
				{
					MM(s,a) = MM(s+k2,a) = 0;
					for( int b=0; b<Ui1.GetNbrRow(); ++b )
						MM(s,a) += Ui1(b,s)*Mi1(b,a);
					for( int b=0; b<Ui2.GetNbrRow(); ++b )
						MM(s+k2,a) += Ui2(b,s)*Mi2(b,a);
				}
			}
			// orthogonalize MM and store in Ui
			matrix& Q = *(new matrix(k,k));
			MM.PerformGramSchmidt(Q);
			Uii(i) = &Q;

			// store M in Mi[i]
			GW_DELETE( Mi(i) );
			Mi(i) = &MM;
		}
	}

	// initialize the result with v
	memcpy( w, v, n*sizeof(double) );
	double* ww = new double[n];			// to store temporary multiplication

	// Sparse Matrix Multiplication
	for( int jj=0; jj<J; ++jj )
	{
		int j = jj;
		if( dir==-1 )
			j = J-1-j;
		matrix_list& Ui = Uj[j];
		// remember previous value of w
		memcpy( ww, w, n*sizeof(double) );
		if( j==0 )
		{
			/****************************************************/
			// Special treatment for 1st scale
			for( int i=0; i<P; ++i )
			{
				double* selj = part[i];

				///////////////////////////////////////////////////////////
				// to keep : upper part is of size n-P*k^2
				int offsr = si[i]-i*k2;			// offset on row
				int longr = G[i]-k2;			// length on row
				// we keep the G(i)-k^2 last
				matrix& U = *Ui(i);
				// here we have 
				if( dir==1 )
				{
					// forward : w(seli) = U(:,k2:G[i])' * ww(selj)
					for( int s=0; s<longr; ++s )		// s is the number of the vector (column of U)
					{
						int ii = offsr + s;		// seli(s)
						w[ii] = 0;
						for( int t=0; t<G[i]; ++t )	// t is the number of the point (row of U)
						{
							int jj = (int) selj[t]-1;
							w[ii] += U( t, s+k2 )*ww[jj];
						}
					}
				}
				else
				{
					// backward : w(selj) = U(:,k2:G[i]) * ww(seli)
					for( int t=0; t<G[i]; ++t )		// t is the number of the point (row of U)
					{
						int jj = (int) selj[t]-1;
						w[jj] = 0;
						for( int s=0; s<longr; ++s )	// s is the number of the vector (column of U)
						{
							int ii = offsr + s;
							w[jj] += U( t, s+k2 )*ww[ii];
						}
					}
				}


				//////////////////////////////////////////////////////
				// to retransform : lower part is of size P*k2
				offsr = n - P*k2 + i*k2;
				GW_ASSERT( offsr>=0 );
				// here we have : seli = offs+(0:k2-1);
				if( dir==1 )
				{
					// forward : w(seli) = U(:,k2:k)' * ww(selj)
					for( int s=0; s<k2; ++s )		// s is the number of the vector (column of U)
					{
						int ii = offsr + s;		// seli(s)
						w[ii] = 0;
						for( int t=0; t<G[i]; ++t )	// t is the number of the point (row of U)
						{
							int jj = (int) selj[t]-1;
							w[ii] += U(t,s)*ww[jj];
						}
					}
				}
				else
				{
					// backward : w(selj) = w(selj) + U(:,1:k2)' * ww(seli)
					for( int t=0; t<G[i]; ++t )	// t is the number of the point (row of U)
					{
						int jj = (int) selj[t]-1;
						for( int s=0; s<k2; ++s )		// s is the number of the vector (column of U)
						{
							int ii = offsr + s;		// seli(s)
							w[jj] += U(t,s)*ww[ii];
						}
					}
				}
			}
		}
		else
		{
			/****************************************************/
			// Other following scales

			// at this scale, we have nj = P/2^(j) groups
			int nj = P >> j;	// P/( 2^j )
			int mj = nj*k;    // total length of the blocks

			int offsr = n-mj;

			for( int i=0; i<nj; ++i )
			{
				// selj = offs + k*i+(0:k-1);
				// seli = offs + k2*i+(0:k2-1);
				matrix& U = *Ui(i);
				if( dir==1 )
				{
					// w(seli) = U(:,k2:k)' * ww(selj);
					for( int s=0; s<k2; ++s )
					{
						int ii = offsr + k2*i+s;	// seli(s)
						w[ii] = 0;
						for( int t=0; t<k; ++t )
						{
							int jj = offsr + k*i+t;	// selj(t)
							w[ii] += U(t,s+k2)*ww[jj];
						}
					}

					// seli = offs + mj/2+k2*i+(0:k2-1)
					// forward : w(seli) = U(:,0:k2-1)' * ww(selj)
					for( int s=0; s<k2; ++s )		// s is the number of the vector (column of U)
					{
						int ii = offsr + mj/2 + k2*i+s;		// seli(s)
						w[ii] = 0;
						for( int t=0; t<k; ++t )	// t is the number of the point (row of U)
						{
							int jj = offsr + k*i+t;	// selj(t)
							w[ii] += U(t,s)*ww[jj];
						}
					}
				}
				else
				{
					// w(selj) = w(selj) + U(:,k2:k) * ww(seli); 
					for( int t=0; t<k; ++t )
					{
						int jj = offsr + k*i+t;	// selj(t)
						w[jj] = 0;
						for( int s=0; s<k2; ++s )
						{
							int ii = offsr + k2*i+s;	// seli(s)
							w[jj] += U(t,s+k2)*ww[ii];
						}
					}  


					// seli = offs + mj/2+k2*i+(0:k2-1)
					// backward : w(selj) = w(selj) + U(:,0:k2-1) * ww(seli)
					for( int t=0; t<k; ++t )	// t is the number of the point (row of U)
					{
						int jj = offsr + k*i+t;	// selj(t)
						for( int s=0; s<k2; ++s )		// s is the number of the vector (column of U)
						{
							int ii = offsr + mj/2 + k2*i+s;		// seli(s)
							w[jj] += U(t,s)*ww[ii];
						}
					}
				}
			}

		}

	}

	GW_DELETEARRAY(ww);
	if( destroy_data )
	{
		// destroy data we are responsible for
		GW_DELETEARRAY(G);
		GW_DELETEARRAY(si);
		for( int i=0; i<P; ++i )
			GW_DELETEARRAY(part[i]);
		GW_DELETEARRAY(part);
	}
	// delete matrix list array
	GW_DELETEARRAY(Uj);
}