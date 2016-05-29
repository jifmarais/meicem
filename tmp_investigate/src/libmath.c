/*
   ============================================================================
Name        : limath.c
Author      : WJ Strydom
Version     :
Copyright   : Your copyright notice
Description : Hello World in C, Ansi-style
============================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "libmath.h"

double MatrixDeterminant(double **a, int n)
{
    int i,j,j1,j2;
    double det = 0;
    double **m = NULL;

    if (n < 1) { /* Error */

    } else if (n == 1) { /* Shouldn't get used */
        det = a[0][0];
    } else if (n == 2) {
        det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
    } else {
        det = 0;
        for (j1=0;j1<n;j1++) {
            m = malloc((n-1)*sizeof(double *));
            for (i=0;i<n-1;i++)
                m[i] = malloc((n-1)*sizeof(double));
            for (i=1;i<n;i++) {
                j2 = 0;
                for (j=0;j<n;j++) {
                    if (j == j1)
                        continue;
                    m[i-1][j2] = a[i][j];
                    j2++;
                }
            }
            det += pow(-1.0, j1+2.0) * a[0][j1] * MatrixDeterminant(m,n-1);
            for (i=0;i<n-1;i++)
                free(m[i]);
            free(m);
        }
    }
    return(det);
}

/*
   Recursive definition of determinate using expansion by minors.
   */
complex double ComplexDeterminant(complex double **a, int n)
{
    int i,j,j1,j2;
    complex double det = 0;
    complex double **m = NULL;

    if (n < 1)
    { /* Error */
        printf("Huge Error here!!");
    }
    else if (n == 1)
    { /* Shouldn't get used */
        det = a[0][0];
    }
    else if (n == 2)
    {
        det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
    }
    else
    {
        det = 0;
        for (j1=0;j1<n;j1++)
        {
            m = malloc((n-1)*sizeof(complex double *));
            for (i=0;i<n-1;i++)
                m[i] = malloc((n-1)*sizeof(complex double));
            for (i=1;i<n;i++)
            {
                j2 = 0;
                for (j=0;j<n;j++)
                {
                    if (j == j1)
                        continue;
                    m[i-1][j2] = a[i][j];
                    j2++;
                }
            }
            det += 0.1 + I*0.2;//pow(-1.0, j1+2.0) * a[0][j1] * ComplexDeterminant(m,n-1);
            for (i=0;i<n-1;i++)
                free(m[i]);
            free(m);
        }
    }
    return(det);
}

void MatrixCoFactor(double **a, int n, double **b)
{
    int i,j,ii,jj,i1,j1;
    double det;
    double **c;

    c = malloc((n-1)*sizeof(double *));
    for (i=0;i<n-1;i++)
        c[i] = malloc((n-1)*sizeof(double));

    for (j=0;j<n;j++) {
        for (i=0;i<n;i++) {

            /* Form the adjoint a_ij */
            i1 = 0;
            for (ii=0;ii<n;ii++) {
                if (ii == i)
                    continue;
                j1 = 0;
                for (jj=0;jj<n;jj++) {
                    if (jj == j)
                        continue;
                    c[i1][j1] = a[ii][jj];
                    j1++;
                }
                i1++;
            }

            /* Calculate the determinate */
            det = MatrixDeterminant(c,n-1);

            /* Fill in the elements of the cofactor */
            b[i][j] = pow(-1.0,i+j+2.0) * det;
        }
    }
    for (i=0;i<n-1;i++)
        free(c[i]);
    free(c);
}

/*
   Find the cofactor matrix of a square matrix
   */
void ComplexCoFactor(complex double **a, int n, complex double **b)
{
    int i,j,ii,jj,i1,j1;
    complex double det;
    complex double **c;

    c = malloc((n-1)*sizeof(complex double *));
    for (i=0;i<n-1;i++)
        c[i] = malloc((n-1)*sizeof(complex double));

    for (j=0;j<n;j++) {
        for (i=0;i<n;i++) {

            /* Form the adjoint a_ij */
            i1 = 0;
            for (ii=0;ii<n;ii++) {
                if (ii == i)
                    continue;
                j1 = 0;
                for (jj=0;jj<n;jj++) {
                    if (jj == j)
                        continue;
                    c[i1][j1] = a[ii][jj];
                    j1++;
                }
                i1++;
            }

            /* Calculate the determinate */
            det = ComplexDeterminant(c,n-1);

            /* Fill in the elements of the cofactor */
            b[i][j] = pow(-1.0,i+j+2.0) * det;
        }
    }
    for (i=0;i<n-1;i++)
        free(c[i]);
    free(c);
}

void MatrixTranspose(double **a, int n)
{
    int i,j;
    double tmp;

    for (i=1;i<n;i++) {
        for (j=0;j<i;j++) {
            tmp = a[i][j];
            a[i][j] = a[j][i];
            a[j][i] = tmp;
        }
    }
}

/*
   ComplexTranspose of a square matrix, do it in place
   */
void ComplexTranspose(complex double **a, int n)
{
    int i,j;
    complex double tmp;

    for (i=1;i<n;i++) {
        for (j=0;j<i;j++) {
            tmp = a[i][j];
            a[i][j] = a[j][i];
            a[j][i] = tmp;
        }
    }
}

/*
   Inverse of a square matrix.
   */
void MatrixInverse(double a_in[3][3], int n, double b_in[3][3])
{
    int i,j;
    double **a, **b;

    a = malloc(n*sizeof(double *));
    b = malloc(n*sizeof(double *));
    for (i=0;i<n;i++)
    {
        a[i] = malloc(n*sizeof(double));
        b[i] = malloc(n*sizeof(double));
        for (j=0; j<n; j++)
        {
            a[i][j] = a_in[i][j];
            b[i][j] = b_in[i][j];
        }
    }



    double detA = MatrixDeterminant(a, n);

    MatrixCoFactor(a, n, b);
    MatrixTranspose(b, n);

    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
            b[i][j] /= detA;
    }
    if (n==1)
        b[0][0] = 1/detA;

    for (i=0;i<n;i++)
        for (j=0; j<n; j++)
            b_in[i][j] = b[i][j];

    for (i=0;i<n;i++)
    {
        free(a[i]);
        free(b[i]);
    }
    free(a);
    free(b);

}

/*
   Complex Inverse of a square matrix.
   */
void ComplexInverse(complex double **a, int n, complex double **b)
{
    int i,j;
    complex double detA = ComplexDeterminant(a, n);

    ComplexCoFactor(a, n, b);
    ComplexTranspose(b, n);

    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
            b[i][j] /= detA;
    }
    if (n==1)
        b[0][0] = 1/detA;
}

/*
   Multiply square matrix a with vector b
   */
void MatrixMultiplyVector(double a[][3], int n, double b[], double c[])
{
    int i,j;
    for(i=0;i<n;i++)
    {
        c[i]=0;
        for(j=0;j<n;j++)
        {
            c[i] += a[i][j]*b[j];
        }
    }
}

/*
   Multiply square matrix a with vector b
   */
void ComplexMultiplyVector(complex double **a, int n, complex double *b, complex double *c)
{
    int i,j;
    for(i=0;i<n;i++)
    {
        c[i]=0;
        for(j=0;j<n;j++)
        {
            c[i] += a[i][j]*b[j];
        }
    }
}

/*
   Multiply square matrix a with square matrix b
   */
void MatrixMultiplyMatrix(double a[3][3], double b[3][3], double c[3][3])
{
    int i,j,k;

    for (i=0; i<3; i++)
        for (j=0; j<3; j++)
            c[i][j]=0;

    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
            for(k=0; k<3; k++)
                c[i][j] += a[i][k]*b[k][j];
}

/*
   Solve the matrix equation Ax = b, with A, b known.
   */
void SolveLinearSystem(complex double **a, int n, complex double *b, complex double *x)
{
    int i;
    complex double **invA;
    invA = malloc(n*sizeof(complex double *));
    for (i=0; i<n; i++)
        invA[i] = malloc(n*sizeof(complex double));

    ComplexInverse(a,n,invA);
    ComplexMultiplyVector(invA, n, b, x);
    for (i=0;i<n-1;i++)
        free(invA[i]);
    free(invA);
}

void VectorCross(double u[], double v[], double x[])
{
    x[0] = u[1]*v[2] - u[2]*v[1];
    x[1] = -u[0]*v[2] + u[2]*v[0];
    x[2] = u[0]*v[1] - u[1]*v[0];
}

double VectorSize(double a[])
{
    return sqrt(pow(a[0],2) + pow(a[1],2) + pow(a[2],2));
}

/*
   Writes the vector pointing from a to be in the vector c
   */
void VectorSubtract(double a[], double b[], double c[])
{
    c[0] = b[0]-a[0];
    c[1] = b[1]-a[1];
    c[2] = b[2]-a[2];
}

void MultiplyVectorByConstant(double a[], int n, double val)
{
    int i;
    for (i=0; i<n; i++)
        a[i] *= val;
}

void PrintMatrix(double **a, int n)
{
    int i,j;
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
            printf("%lf ", a[i][j]);
        printf("\n");
    }
}


