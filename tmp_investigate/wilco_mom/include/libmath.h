/*
   ============================================================================
Name        : mathlib.c
Author      : WJ Strydom
Version     :
Copyright   : Your copyright notice
Description : Hello World in C, Ansi-style
============================================================================
*/
#ifndef libmath
#define libmath

double MatrixDeterminant(double **a, int n);
complex double ComplexDeterminant(complex double **a, int n);
void MatrixCoFactor(double **a, int n, double **b);
void ComplexCoFactor(complex double **a, int n, complex double **b);
void MatrixTranspose(double **a, int n);
void ComplexTranspose(complex double **a, int n);
void MatrixInverse(double a_in[3][3], int n, double b_in[3][3]);
void ComplexInverse(complex double **a, int n, complex double **b);
void MatrixMultiplyVector(double a[][3], int n, double b[], double c[]);
void ComplexMultiplyVector(complex double **a, int n, complex double *b, complex double *c);
void MatrixMultiplyMatrix(double a[3][3], double b[3][3], double c[3][3]);
void SolveLinearSystem(complex double **a, int n, complex double *b, complex double *x);
void VectorCross(double u[], double v[], double x[]);
double VectorSize(double a[]);
void VectorSubtract(double a[], double b[], double c[]);
void MultiplyVectorByConstant(double a[], int n, double val);
void PrintMatrix(double **a, int n);

#endif
