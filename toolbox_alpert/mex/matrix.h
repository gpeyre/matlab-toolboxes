#include <math.h>
#include "config.h"
#ifdef GW_DEBUG
	#include <crtdbg.h>
#endif
// #include "mex.h"


class matrix
{
private:
	double* M_;
	int nb_row_, nb_col_;
public:
	matrix(int n, int m)
	{
		nb_row_ = n;
		nb_col_ = m;
		M_ = new double[n*m];
	}
	matrix( const matrix& A )
	{
		nb_col_ = A.nb_col_;
		nb_row_ = A.nb_row_;
		M_ = new double[nb_row_*nb_col_];
		memcpy( M_, A.M_, nb_row_*nb_col_*sizeof(double) );
	}
	~matrix()
	{
		delete[] M_;
	}
	double operator () (int a, int b) const
	{
		GW_ASSERT( (a<nb_row_) & (a>=0) & (b<nb_col_) & (b>=0) );
		return M_[a+(b)*nb_row_];
	}
	double& operator () (int a, int b)
	{
		GW_ASSERT( (a<nb_row_) & (a>=0) & (b<nb_col_) & (b>=0) );
		return M_[a+(b)*nb_row_];
	}
	int GetNbrRow()
	{
		return nb_row_; 	
	}
	int GetNbrCol()
	{
		return nb_col_; 	
	}
	void PerformGramSchmidt(matrix& Q)
	{
		if( Q.nb_row_==Q.nb_col_ && Q.nb_row_==1 )
		{
			Q(0,0) = 1;
			return;
		}
		if( Q.nb_row_==Q.nb_col_ && Q.nb_row_==2 )
		{
			const double a = 1.0/sqrt(2);
			Q(0,0) = Q(1,0) = Q(0,1) = a;
			Q(1,1) = -a;
			return;
		}
#if 0	// using Matlab QR function
		mxArray* lhs[2];
		mxArray* rhs;
		rhs    = mxCreateDoubleMatrix(nb_row_, nb_col_, mxREAL);
		double* m = mxGetPr(rhs);
		memcpy( m, M_, nb_row_*nb_col_*sizeof(double) );
		int res = mexCallMATLAB(2, lhs, 1, &rhs,  "qr");
		double* q = mxGetPr(lhs[0]);
		memcpy( Q.M_, q, nb_row_*nb_row_*sizeof(double) );
		mxDestroyArray(rhs);
		mxDestroyArray(lhs[0]);
		mxDestroyArray(lhs[1]);
		// gram_schmidt(this->M_, Q.M_, this->GetNbrRow(), this->GetNbrCol());
#else
		// using our own QR function
		PerformQR(Q);
#endif
	}

	double hypot(const double &a, const double &b)
	{
		if( a==0 )
			return GW_ABS(b);
		else
		{
			double c = b/a;
			return a * sqrt(1 + c*c);
		}
	}

	void PerformQR(matrix& Q)
	{
		matrix QR_(*this);
		// double* Rdiag = new double[nb_col_];

		int i=0, j=0, k=0;

		// Main loop.
		for (k = 0; k < Q.nb_col_; k++) 
		{
			// Compute 2-norm of k-th column without under/overflow.
			double nrm = 0;
			for (i = k; i < Q.nb_row_; i++) 
			{
				nrm = hypot(nrm,QR_(i,k));
			}

			if (nrm != 0.0) 
			{
				// Form k-th Householder vector.
				if (QR_(k,k) < 0) 
				{
					nrm = -nrm;
				}
				for (i = k; i < Q.nb_row_; i++) 
				{
					QR_(i,k) /= nrm;
				}
				QR_(k,k) += 1.0;

				// Apply transformation to remaining columns.
				for (j = k+1; j < Q.nb_col_; j++) 
				{
					double s = 0.0; 
					for (i = k; i < Q.nb_row_; i++) 
					{
						s += QR_(i,k)*QR_(i,j);
					}
					s = -s/QR_(k,k);
					for (i = k; i < Q.nb_row_; i++) 
					{
						QR_(i,j) += s*QR_(i,k);
					}
				}
			}
			// Rdiag[k] = -nrm;
		}

		// compute the Q matrix
		i=0; j=0; k=0;

		for( k=Q.nb_col_-1; k >= 0; k--) 
		{
			for( i=0; i < Q.nb_row_; i++ ) 
			{
				Q(i,k) = 0.0;
			}
			Q(k,k) = 1.0;
			for( j=k; j < Q.nb_col_; j++ ) 
			{
				if( QR_(k,k)!=0 ) 
				{
					double s = 0.0;
					for( i=k; i < Q.nb_row_; i++ ) 
					{
						s += QR_(i,k)*Q(i,j);
					}
					s = -s/QR_(k,k);
					for( i=k; i<Q.nb_row_; i++ ) 
					{
						Q(i,j) += s*QR_(i,k);
					}
				}
			}
		}

		// GW_DELETEARRAY(Rdiag);
	}

};


// a list of matrix
class matrix_list
{
private:
	matrix** list_;
	int n_;			// size of the list

public:
	matrix_list(int n = 0)
	{
		if( n>0 )
			SetSize(n);
		else
		{
			list_ = NULL; 
			n_ = -1;
		}
	}
	void SetSize(int n)
	{
		list_ = new matrix*[n];
		memset( list_, 0, n*sizeof(matrix*) );
		n_ = n;
	}
	~matrix_list()
	{
		for( int i=0; i<n_; ++i )
			GW_DELETE(list_[i]);
		GW_DELETE(list_);
	}
	matrix* operator () (int a) const
	{
		GW_ASSERT( (a<n_) & (a>=0) );
		return list_[a];
	}
	matrix*& operator () (int a)
	{
		GW_ASSERT( (a<n_) & (a>=0) );
		return list_[a];
	}
};
