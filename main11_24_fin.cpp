#include <iostream>
 #include <stdio.h>
//#include "stdafx.h"
//#include "targetver.h"


#include <math.h>

#define NR_END 1
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#include <SDKDDKVer.h>
//#include <condefs>
#include <iostream>
#include <iomanip>

#include <time.h>
#include <cmath>
//#include "nrutil.h"

#include <complex>
#include "Eigen/Eigenvalues"

/* #include <stddef.h> */
#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*
using namespace std;


 	double Mu=	0.1,
	  U=1,
	   z=4,
	   t = 0.01,
	   T =0.1;// 1/500;
	   
	    double U3A = 1;
		double UA= U, UB = U, 
		W=0.1,
		 MuA = Mu, MuB = Mu,
			tA= t, tB= t; 
//			mA = 1.5, mB = 1.8;

	double gamma = 1, ksi = 1, eta = 1;
	double * d;
    double * e;
	double ** A;


int n = 5, m = 5;	    
int N = n*m;	    
/* run this program using the console pauser or add your own getch, system("pause") or input loop */
void nrerror(int Number )
// Numerical Recipes standard error handler
{
	fprintf(stderr, "Numerical Recipes run-time error...\n");
	fprintf(stderr, "%s\n", "error_text/");
	fprintf(stderr, "...now exiting to system...\n");
	exit(1);
} 
double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	double **m;

	/* allocate pointers to rows */
	m = (double **)malloc((size_t)((nrow + NR_END) * sizeof(double*)));
	if (!m) nrerror(1);
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl] = (double *)malloc((size_t)((nrow*ncol + NR_END) * sizeof(double)));
	if (!m[nrl]) nrerror(2);
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;
//	zeroM(m);
	/* return pointer to array of pointers to rows */
	return m;
}
double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v = (double *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(double)));
	if (!v) nrerror(3);
//	zeroV(v);
	return v - nl + NR_END;
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}
void MdV(double ** A, double * b, double * out)
{
	int i, j; double value = 0;
	for (j = 0; j<N; j++) {
		for (i = 0, value = 0; i<N; i++)
			value += A[j][i] * b[i];
		out[j] = value;
	}
}
void MdM(double ** C, double ** B, double ** out)
{
	int i, j, k; 
	double value = 0;
	for (i = 1, value = 0; i<N+1; i++)
	{
	
	for (j = 1; j<N+1; j++)
	 {
		for (k = 1; k<N+1; k++)
		{
			value += C[i][k] * B[k][j];
/*			if (B[k][j] != 0)
			{
						cout << C[i][k] << "   " << B[k][j];
				
			}*/
			
		}
					
		out[i][j] = value;
//		if (value != 0) cout << value;
		value = 0;
	}
}
}

void VTdM(double * bT, double ** A, double * out)
{
	int i, j; double value = 0;
	for (j = 0; j<N; j++) {
		for (i = 0, value = 0; i<N; i++)
			value += A[i][j] * bT[i];
		out[j] = value;
	}
}
void VTdV(double *bT, double * c, double & out)
{
	out = 0;
	for (int i = 0; i<N; i++)
		out += bT[i] * c[i];
}
void VmV(double *b, double & s, double *c)
{
	for (int i = 0; i<N; i++)
		b[i] -= s * c[i];
}
void VdS(double *b, double s)
{
	for (int i = 0; i<N; i++)
		b[i] *= s;
}
// Various outputs
void printV(double * b)
{
	cout << " [";
	for (int i = 0; i<N; i++)
		cout << setw(8) << b[i] << "    "; cout << "] ";
}
void printV(double * b, int offset)
{
	cout << " [";
	for (int i = 0 + offset; i<N + offset; i++)
		cout << setw(8) << b[i]; cout << "] ";
}
void zeroM( double ** out)
{
	int i, j; 

	for (j = 0; j<N+1; j++) {
		for (i = 0; i<N+1; i++)
			out[i][j] = 0;
	}
}
void zeroV(double * out)
{
	int j; 
	for (j = 0; j<N+1; j++) {
		out[j] = 0;
	}
}
void TextOut(string Text, int Number)
{

	// #include <stdio.h>

	FILE	 	*stream;
	char str[10];


	char Stri[25];

	//   itoa(Number, Stri, 12);
	//  string strs(Stri);

	//   strcpy(str,("D/Data"+strs+".txt").c_str() ) ;

//	sprintf(Stri, "Data.txt");
	stream =
		fopen("DataEnh.txt",  "a");
fopen(Stri,  "a");

	//	stream = fopen(Stri, , "rw");

	//  "Data.txt"
	char * p;
	p = (char *)Text.c_str();
	if (stream != NULL)
	{
		fprintf(stream, p, "DataEnh.txt");
//		fprintf(stream, p, Stri);
		fclose(stream);
	}
	return;
}

void PrintMatrix(int Dim, double ** H, bool ToFile)
{

	char buf[256];

	//	memset(buf[256], ' ');

	//     TextOut(ST, 0);
	int i, j;
	double Val;
	for (i = 0; i < Dim; i++)
	{
		for (j = 0; j < Dim; j++)
		{ 
			Val = H[i][j];
			sprintf(buf, "%4.4f    ", Val);
			string ST(buf);
			if(ToFile == 0)
			{
				cout << ST;
			}
			else TextOut(ST, 0);
		}
		sprintf(buf, "\n");
		string ST(buf);
		if(ToFile == 0)
		{
				cout << ST;
		}
		else TextOut(ST, 0);
	}
}
void PrintVector(int Dim, double * H, bool ToFile)
{

	char buf[256];

	//	memset(buf[256], ' ');

	//     TextOut(ST, 0);
	int i;
	double Val;
	for (i = 0; i < Dim; i++)
	{

			Val = H[i];
			sprintf(buf, "%4.4f    ", Val);
			string ST(buf);
			if(ToFile == 0)
			{
				cout << ST;
			}
			else TextOut(ST, 0);

	}
}


void tred2(double **a, int n, double d[], double e[])
/*******************************************************************************
Householder reduction of a real, symmetric matrix a[1..n][1..n]. 
On output, a is replaced by the orthogonal matrix Q effecting the
transformation. d[1..n] returns the diagonal elements of the tridiagonal matrix,
and e[1..n] the off-diagonal elements, with e[1]=0. Several statements, as noted
in comments, can be omitted if only eigenvalues are to be found, in which case a
contains no useful information on output. Otherwise they are to be included.
*******************************************************************************/
{
	int l,k,j,i;
	double scale,hh,h,g,f;

	for (i=n;i>=2;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 1) {
			for (k=1;k<=l;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0) /* Skip transformation. */
				e[i]=a[i][l];
			else {
				for (k=1;k<=l;k++) {
					a[i][k] /= scale; /* Use scaled a's for transformation. */
					h += a[i][k]*a[i][k]; /* Form sigma in h. */
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g; /* Now h is equation (11.2.4). */
				a[i][l]=f-g; /* Store u in the ith row of a. */
				f=0.0;
				for (j=1;j<=l;j++) {
					/* Next statement can be omitted if eigenvectors not wanted */
					a[j][i]=a[i][j]/h; /* Store u/H in ith column of a. */
					g=0.0; /* Form an element of A.u in g. */
					for (k=1;k<=j;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<=l;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h; /* Form element of p in temporarily unused element of e. */
					f += e[j]*a[i][j];
				}
				hh=f/(h+h); /* Form K, equation (11.2.11). */
				for (j=1;j<=l;j++) { /* Form q and store in e overwriting p. */
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=1;k<=j;k++) /* Reduce a, equation (11.2.13). */
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
		}
		/* Next statement can be omitted if eigenvectors not wanted */
		d[1]=0.0;
		e[1]=0.0;
		/* Contents of this loop can be omitted if eigenvectors not
		   wanted except for statement d[i]=a[i][i]; */
		for (i=1;i<=n;i++) { /* Begin accumulation of transformation matrices. */
			l=i-1;
		if (d[i]) { /* This block skipped when i=1. */
			for (j=1;j<=l;j++) {
				g=0.0;
				for (k=1;k<=l;k++) /* Use u and u/H stored in a to form P.Q. */
					g += a[i][k]*a[k][j];
				for (k=1;k<=l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i]; /* This statement remains. */
		a[i][i]=1.0; /* Reset row and column of a to identity matrix for next iteration. */
		for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
	}
}

/******************************************************************************/
void tqli(double d[], double e[], int n, double **z)
/*******************************************************************************
QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors
of a real, symmetric, tridiagonal matrix, or of a real, symmetric matrix
previously reduced by tred2 sec. 11.2. On input, d[1..n] contains the diagonal
elements of the tridiagonal matrix. On output, it returns the eigenvalues. The
vector e[1..n] inputs the subdiagonal elements of the tridiagonal matrix, with
e[1] arbitrary. On output e is destroyed. When finding only the eigenvalues,
several lines may be omitted, as noted in the comments. If the eigenvectors of
a tridiagonal matrix are desired, the matrix z[1..n][1..n] is input as the
identity matrix. If the eigenvectors of a matrix that has been reduced by tred2
are required, then z is input as the matrix output by tred2. In either case,
the kth column of z returns the normalized eigenvector corresponding to d[k].
*******************************************************************************/
{
	double pythag(double a, double b);
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;

	for (i=2;i<=n;i++) e[i-1]=e[i]; /* Convenient to renumber the elements of e. */
	e[n]=0.0;
	for (l=1;l<=n;l++) {
		iter=0;
		do {
			for (m=l;m<=n-1;m++) { /* Look for a single small subdiagonal element to split the matrix. */
				dd=fabs(d[m])+fabs(d[m+1]);
				if ((double)(fabs(e[m])+dd) == dd) break;
			}
			if (m != l) {
				if (iter++ == 150)
				{
				 cout << "Too many iterations in tqli";
				 return;
				 } 
				g=(d[l+1]-d[l])/(2.0*e[l]); /* Form shift. */
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g)); /* This is dm - ks. */
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) { /* A plane rotation as in the original QL, followed by Givens */
					f=s*e[i];          /* rotations to restore tridiagonal form.                     */
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) { /* Recover from underflow. */
						d[i+1] -= p;
						e[m]=0.0;
						cout << "!!!";
						break;
						
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					/* Next loop can be omitted if eigenvectors not wanted */
					for (k=1;k<=n;k++) { /* Form eigenvectors. */
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}
				if (r == 0.0 && i >= l)
				{
				
				   cout << "!!!";
				   continue;
				}
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
}

/******************************************************************************/
double pythag(double a, double b)
/*******************************************************************************
Computes (a2 + b2)1/2 without destructive underflow or overflow.
*******************************************************************************/
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+(absb/absa)*(absb/absa));
	else
	{
//	 if	(absb == 0.0) cout << "p1";
		 return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
	 }
}
void SortVector(int Dim, double * H)
{
	//array declaration

	int i,j;
	double temp;
	

	//sorting - ASCENDING ORDER
	for(i=0;i<Dim;i++)
	{		
		for(j=i+1;j<Dim;j++)
		{
			if(H[i]>H[j])
			{
				temp  =H[i];
				H[i]=H[j];
				H[j]=temp;
			}
	}
   }

	
	return ;			
						
}

long double FreeNrg(double mA, double mB, double t_)
{
	tA = t_;

	 MuA = Mu;
	 UA = U;

/*	 tB= t_; 
     MuB = Mu;	 
	 UB = U;*/
	 
	 tB= gamma*t_; 
     MuB = ksi*Mu;	 
	 UB = eta*U;


	

    
    int i, j;
	zeroV(d);
	zeroV(e);
	zeroM(A);
	

	int i2, j2, nn, mm;
	for( i=0; i <n; i++)
	{
	for( i2=0; i2 <n; i2++)
	{	
	 for( j=0; j <m; j++) 	
 	 {
 	 for( j2=0; j2 <m; j2++) 	
 	 {	
		if(i==i2 && j == j2)
		{

				nn =(i*m+j+1);
				mm =(i2*m+j2+1);
					
					A[nn][mm] += -MuA*i - MuB*j+ UA*(i)*(i-1)/2+ UB*(j)*(j-1)/2+W*i*j+U3A/6*(i)*(i-1)*(i-2);
					
		}
		if( j == j2 && i-1==i2)
		{
				nn =(i*m+j+1);
				mm =(i2*m+j2+1);		
					
					A[nn][mm] += -z*tA*mA*(sqrt(i));
					
		}
		if( j == j2 && i+1==i2)
		{
  			nn =(i*m+j+1);
			mm =(i2*m+j2+1);	
					
					A[nn][mm] += -z*tA*mA*(sqrt(i+1));
					
		}
		
		
		
		if( j-1 == j2 && i==i2)
		{
  			nn =(i*m+j+1);
			mm =(i2*m+j2+1);		
					
					A[nn][mm] += -z*tB*mB*(sqrt(j));
					
		}
		if( j+1 == j2 && i==i2)
		{
  			nn =(i*m+j+1);
			mm =(i2*m+j2+1);	
					
					A[nn][mm] += -z*tB*mB*(sqrt(j+1));
					
		}
		
		
     }
	 }	
	}
	}	
	

//diagonalization

///!!!fast alg 
   tred2(A, N, d, e);    
   tqli (d, e, N, A);
 
 
 //  slow alg
/*    
  Eigen::MatrixXd B(N, N);
	
	for( i=0; i <N; i++)
	{
	for( j=0; j <N; j++)
	{
		B(i,j) = A[i+1][j+1] ;
	}
	}	
	
  Eigen::EigenSolver<Eigen::MatrixXd> ces;	
  ces.compute(B);
 // std::cout << "The eigenvalues of a are:" << std::endl << ces.eigenvalues() << std::endl;	
  
  
	for( i=0; i <N; i++)
	{
		d[i]=std::real(ces.eigenvalues()[i] );
	}   
    
    */
    
    
    
    
    

	SortVector(N+1,d); 

	long double Z_=0, fN=0;
	

	fN = tA*z/2*mA*mA+tB*z/2*mB*mB+d[0];
	
	
	return (fN);
}
long double FreeNrgEnhanced(double m, double t_, double Mu_)
{
// not
	double * d = dvector(0, N);


    double * e = dvector(0, N);
    int i;
	double ** A = dmatrix(0, N, 0, N);
	
 	for( i=2; i <N+1; i++)
 	{
 		d[i]= -(i-2)*Mu_ + (i-2)*(i-3)*U/2; // here +1 if compare to Matematica because 0-th rows do not take into account
 		e[i+1]= -z/2*t_*m*sqrt(i-1);       // here +1 if compare to Matematica because 0-th rows do not take into account
	}

    
    tqli (d, e, N, A);


	SortVector(N+1,d);
 


	long double lnZ_, lnZ_add=0;
   for( i=1; i <N; i++)
 	{
 		lnZ_add+=exp(-(d[i]-d[0])/T);
	}

    lnZ_= d[0]+lnZ_add ;

return lnZ_;
    
}


long double lnZ(double m, double t_)
{
	double * d = dvector(0, N);
    double * e = dvector(0, N);
    int i;
	double ** A = dmatrix(0, N, 0, N);
	
 	for( i=2; i <N+1; i++)
 	{
 		d[i]= -(i-2)*Mu + (i-2)*(i-3)*U/2; // here +1 if compare to Matematica because 0-th rows do not take into account
 		e[i+1]= -z/2*t_*m*sqrt(i-1);       // here +1 if compare to Matematica because 0-th rows do not take into account
	}

    
    tqli (d, e, N, A);
    


	SortVector(N+1,d); 

//	T=1;
	long double lnZ_, lnZ_add=0;
    for( i=1; i <N; i++)
 	{
 		lnZ_add+=exp(-(d[i]-d[0])/T);
	}

    lnZ_= d[0]-lnZ_add ;

	
	return (lnZ_);
}



double dF( double t_)
{
	double  fN_1,  fN_2;
	fN_1 = (FreeNrg(-0.1,0.0, t_) - 2*FreeNrg(0, 0.0, t_) + FreeNrg(0.1, 0.0, t_	))*100.0;
//	n = 8; m = 1;
	fN_2 = (FreeNrg(0.1, -0.1, t_) - 2*FreeNrg(0.1, 0, t_) + FreeNrg(0.1, 0.1, t_	))*100.0;
	
	
	return fN_1;
	
}
double dF( double t_, double m)
{
	double  fN_1,  fN_2;
	fN_1 = (FreeNrg(-0.1,m, t_) - 2*FreeNrg(0, m, t_) + FreeNrg(0.1, m, t_	))*100.0;
	
	
	return fN_1;
	
}
double dF2( double t_, double m)
{
	double  fN_1;
	fN_1 = (FreeNrg(m,-0.1, t_) - 2*FreeNrg(m, 0, t_) + FreeNrg(m, 0.1, t_	))*100.0;
	
	
	return fN_1;
	
}





void Eigenproblem()
{
  const int n = 2;
  Eigen::MatrixXd a(n, n);
  //typedef std::complex<double> C;
  a << 1.2, 123,
  31, 124;
    
  Eigen::EigenSolver<Eigen::MatrixXd> ces;
  ces.compute(a);
  std::cout << "The eigenvalues of a are:" << std::endl << ces.eigenvalues() << std::endl;
  std::cout << "The eigenvalues of a are:" <<std::endl << std::real(ces.eigenvalues()[1] )<< std::endl;
  
  
}


double CritLine(double t_)
{

    int n = 2;
    double t2 = 0.2;
    double t1 =   t_;
	double tt;

    while( ( fabs(dF(t1)) > 0.0001 ) && ( n <= 100 ) )
    {
        tt = t1 - (dF(t1) * (t1 - t2)) / (dF(t1) - dF(t2));
        t2 = t1;
        t1 = tt;

        n++;
    }
	return 	tt ;


}
double CritLine(double t_, double m)
{

    int n = 2;
    double t2 = 0.2;
    double t1 =   t_;
	double tt;

    while( ( fabs(dF(t1, m)) > 0.000000001 ) && ( n <= 100 ) )
    {
        tt = t1 - (dF(t1, m) * (t1 - t2)) / (dF(t1, m) - dF(t2, m));
        t2 = t1;
        t1 = tt;

        n++;
    }
	return 	tt ;


}
double CritLineEnh(double t_, bool cross)
{

	double m = 0;
	double fN, tempfN = 10000, tempmB = 0;
    int n = 2;
	int i, j, min_j;
    double t2 = 0.2;
    double t1 =   t_;
	double tt, ttt = 0.04, tttt = 1, Delta=1000, min_Delta = 999;


	for( i=0; i <10; i++)
	{
    while( ( fabs(dF(t1, m)) > 0.000000001 ) && ( n <= 100 ) )
    {
        tt = t1 - (dF(t1, m) * (t1 - t2)) / (dF(t1, m) - dF(t2, m));
        t2 = t1;
        t1 = tt;

        n++;
    }

  	
	  m=0;	  	
	for( j=0; j <40; j++)
	{

		m += 0.05;
		
		fN=FreeNrg(0,m, tt) ;
		if (fN<tempfN)
		{
             if(cross== true && j!=min_j && j!=0
		//  && (min_j+1!=j&& j+1!=min_j && j!=min_j && j!=0)
		 )	
			   {
			   	// uslowie pri kotorom proish faz perehod 1 roda
			   	// v etom sluchae FreeNRG ravny primerno na linii faz perehoda
			   	// naidem t pri kotorom eto vipolnaetsja
			   	
		//	   	  printf( "\n %4.4f %d %d    \n", fN , j,min_j );
			   	  
			   	    while( ttt>=0 )
   				    {
					//  Delta = fabs(FreeNrg(0,0, ttt) - FreeNrg(0.1,0.1,ttt));
					  
					  Delta = fabs(FreeNrg(1,1, ttt) - FreeNrg(0,m,ttt)); 
					  double FN1, FN2;
					  FN1= FreeNrg(1,1, ttt);
					  FN2= FreeNrg(0,m,ttt);
					  if(Delta<min_Delta)
					  {	 
					  
	//				    printf( "\n %9.9f   %4.19f    %4.19f     %4.4f     %4.4f    %4.4f    \n", Delta , FN1, FN2,ttt,m,tt );
					    // ne pracuje%4.4f 
					    // pochemu to chto to = 0
					    
					    
					  	tttt = ttt;
					  	min_Delta = Delta;
					  	
					  }
				
					 
					ttt-=0.0001;  					   
			        }
			     return tttt;
			     
				}	//	 */  
			 

				tempmB = m;	
				tempfN = fN;	
				min_j = j;
		}
	}
	
    m = tempmB;
	}
	
	return 	tt ;
}
double CritLineEnh2(double t_, bool cross)
{

	double m = 0;
	double fN, tempfN = 10000, tempmB = 0;
    int n = 2;
	int i, j, min_j;
    double t2 = 0.2;
    double t1 =   t_;
	double tt;
	double ttt = 0.04, tttt = 1, Delta=1000, min_Delta = 999;

	for( i=0; i <10; i++)
	{
    while( ( fabs(dF2(t1, m)) > 0.000000001 ) && ( n <= 100 ) )
    {
        tt = t1 - (dF2(t1, m) * (t1 - t2)) / (dF2(t1, m) - dF2(t2, m));
        t2 = t1;
        t1 = tt;

        n++;
    }

  	
	  m=0;	
	for( j=0; j <=40; j++)
	{

		m += 0.05;
		
		fN=FreeNrg(m,0, tt) ;
		if (fN<tempfN)
		{
         if(cross== true && j!=min_j && j!=0
		//  && (min_j+1!=j&& j+1!=min_j && j!=min_j && j!=0)
		 )	
			   {
			   	// uslowie pri kotorom proish faz perehod 1 roda
			   	// v etom sluchae FreeNRG ravny primerno na linii faz perehoda
			   	// naidem t pri kotorom eto vipolnaetsja
			   	
		//	   	  printf( "\n %4.4f %d %d    \n", fN , j,min_j );
			   	  
			   	    while( ttt>=0 )
   				    {
					//  Delta = fabs(FreeNrg(0,0, ttt) - FreeNrg(0.1,0.1,ttt));
					  
					  Delta = fabs(FreeNrg(1,1, ttt) - FreeNrg(m,0,ttt)); 
					  double FN1, FN2;
					  FN1= FreeNrg(1,1, ttt);
					  FN2= FreeNrg(m,0,ttt);
					  if(Delta<min_Delta)
					  {	 
					  
	//				    printf( "\n %9.9f   %4.19f    %4.19f     %4.4f     %4.4f    %4.4f    \n", Delta , FN1, FN2,ttt,m,tt );
					    // ne pracuje%4.4f 
					    // pochemu to chto to = 0
					    
					    
					  	tttt = ttt;
					  	min_Delta = Delta;
					  	
					  }
				
					 
					ttt-=0.0001;  					   
			        }
			     return tttt;
			     
				}	//	 */ 
				tempmB = m;	
				tempfN = fN;
				min_j = j;
		}
	}
	
    m = tempmB;
	}
	
	return 	tt ;
}


double OccupationEnh()
{
	double res;


    double  fN1, fN2, dFn_dMu, deltaMu;
	deltaMu = 0.01;

	
	
	fN1 = FreeNrgEnhanced(0.1, t, Mu+deltaMu) ;
	fN2 = FreeNrgEnhanced(0.1, t, Mu) ;

/*	Mu += deltaMu;	
	fN1 = FreeNrg(3, t) ;
	Mu -= deltaMu;		
	fN2 = FreeNrg(3, t) ;	*/
	
	
	
	dFn_dMu = (fN1-fN2)/deltaMu;
	res = -dFn_dMu;
	
//	res = FreeNrg(1, t)/Mu;
	
	return res;
}
int mapping(double t_)
{
	double mA = 0, mB = 0, fN, tempfN = 10000, tempmA = 0, tempmB = 0;
	int i, j;
	
	for( i=0; i <=6; i++)
	{
//		mA = i*0.0125;		//mA +=0.0125
		mB=0;
	for( j=0; j <=6; j++)
	{

//		mB = j*0.0125;
//		printf( "%4.4f     \n", mB  );
		
		fN=FreeNrg(mA,mB, t_) ;
		if (fN<tempfN)  /// bylo <     11/16
		{
				tempmA = mA;
				tempmB = mB;	
				tempfN = fN;
		}
		mB += 0.2;
		
	}
	mA +=  0.2;
	}	
	if (tempmA == 0 && tempmB == 0)	
	{
		return 3;
	}
	if (tempmA != 0 && tempmB == 0)	
	{
		return 1;
	}
	if (tempmA == 0 && tempmB != 0)	
	{
		return 2;
	}
	if (tempmA != 0 && tempmB != 0)	
	{
		return 0;
	}
	else return 4;
	
	
	
}

int main() 
{


	d = dvector(0, N);
    e = dvector(0, N);
	A = dmatrix(0, N, 0, N);


	int i = 0, j = 0;
	int  k = 0;
	long double fN, 	d2fN_dm;
	
	char buf[256];


// Psi Average
/*



	for (i=10; i<120; i++)
	{
		Mu = 0.01*i;

		cout   << Mu << "   " << 	PsiAverage(0.1,0.1) << "\n" ;
				

		
	}// */
	


//	mapping 
/*
    int  map=3, mapt;
//	Mu =  0.75;
	W=0.60;	
	U3A = 0.0;
	
		double t_;

	
	ksi = 0.9;
	gamma = 1.1;
	eta = 0.8;
	
	


	T=0.001;//1.0/500;
	
	for (j=0; j<700; j++)
	{

	Mu =  0.0000+0.005*j;
	for (i=0; i<800; i+=5)
	{
		t_ = 0.000001+ 0.0001*i;	
		
		map= mapping(t_) ;
		if (map == 0) 
		{
			break;	
		}
		cout   << Mu << "   " << t_ << "   "<< map   << 	
		 "\n" ;	
		 
		sprintf(buf, "%4.4f    %4.4f    %d   \n", Mu, t_, map  );
		string ST(buf);			
		TextOut(ST, 0);
	}
    }
	// */

///critical line enhanced
   cout   <<"result crit line" << "\n";


	double t1, t2;
	bool crossing = false;
	
	for (i=0; i<1000; i++)
	{
	Mu =  0.005*i;
//	for (j=0; j<20; j++)
	{
//	W=0.05*j;
	W=0.4;	

//	for (k=0; k<30; k++)
	{
//		U3A = 0.02*k;
		
		U3A = 0.5;
//		ksi = 0.5+0.025*k;
	ksi = 0.8;
 	eta = 0.9;// 0.0+0.025*k;
 	gamma = 1.1;

//	Mu = 10.5;
	double tt1, tt2;

		
//		t1=CritLineEnh(0.1, false);
//		t2=CritLineEnh2(0.1, false);
//     	t1=CritLineEnh(0.1, crossing);
//		t2=CritLineEnh2(0.1, crossing);
        t1=CritLineEnh(0.1, true);
		t2=CritLineEnh2(0.1, true);
		tt1=CritLineEnh(0.1, false);
		tt2=CritLineEnh2(0.1, false);

		
		if (fabs(t1-t2)<0.002)
		{
			 printf( "\n  crossing %9.9f   %4.19f      \n",t1,t2 );
			 crossing = true;

		}
		else
		{
			crossing = false;
		}
//		t1=CritLineEnh(0.1, crossing);
//		t2=CritLineEnh2(0.1, crossing);
	
	// gladkaya shiwka reshenij	
	// wybiraetsja to cho s bolshei t
		if (tt1>t1)
		{
			t1=tt1;
		}
		if (tt2>t2)
		{
			t2=tt2;
		}	
		
		
		cout   << t1 << "   "<< t2 << "   " << Mu << "   "<< W  << "   "<< U3A << "   " << 	dF(t) <<
		"   " << gamma << "   "<< ksi  << "   "<< eta << "   " <<		
		 "\n" ;
				
		sprintf(buf, "%4.4f    %4.4f    %4.4f   %4.4f    %4.4f    %4.4f   %4.4f  %4.4f     %4.4f  \n", Mu, U3A,  t1, t2 , W, dF(t),  gamma , ksi  , eta  );
		string ST(buf);	
		TextOut(ST, 0);
}
		
	}
}
	sprintf(buf, "t     Mu        W       U3A    	dF(t)    gamma      ksi       eta \n" );
		string ST(buf);	
	TextOut(ST, 0);
	
//	*/
// occupation
/*	cout   << " T =  " << T << "   "  << "\n" ;
double n;

	for (i=0; i<50; i++)
	{
		Mu = 0.1*i+0.00;
///	Mu = 3.05;

		n=OccupationEnh();
		cout   << n << "   " << Mu << "   "  << "\n" ;
		
		
		sprintf(buf, "%4.4f    %4.4f    \n", Mu, n  );
		string ST(buf);			
		TextOut(ST, 0);

		
	}
		cout   << " T =  " << T << "   "  << "\n" ;
	
//	TextOut(string Text, int Number);;
	
	
	// vyvesti text v fail i razobratsya s tempetaturoi
	
	*/
    free_dvector(d, 0, N);
    free_dvector(e, 0, N);
    free_dmatrix(A, 0, N, 0, N);
    
	getchar();
	return 0;
	
}
