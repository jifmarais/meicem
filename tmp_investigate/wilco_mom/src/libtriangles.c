/*
   ============================================================================
Name        : mom.c
Author      : WJ Strydom
Version     :
Copyright   : Your copyright notice
Description : Hello World in C, Ansi-style
============================================================================
*/

#include <math.h>
#include <complex.h>
#include "libmath.h"

double Distance(double a[], double b[])
{
    double distVec[3];
    VectorSubtract(a, b, distVec);
    return VectorSize(distVec);
}

double TriangleArea(double points[][3])
{
    double AB[3], AC[3], cp[3];
    VectorSubtract(points[0], points[1], AB);
    VectorSubtract(points[0], points[2], AC);
    VectorCross(AB, AC, cp);
    return 0.5*VectorSize(cp);
}

void TriangleCentre (double points[][3], double centre[])
{
    centre[0] = (points[0][0]+points[1][0]+points[2][0])/3.0;
    centre[1] = (points[0][1]+points[1][1]+points[2][1])/3.0;
    centre[2] = (points[0][2]+points[1][2]+points[2][2])/3.0;
}

int TrianglesDontTouch(int *triA, int *triB)
{
    int i,j;
    for (i=0; i<3; i++)
    {
        for (j=0; j<3; j++)
        {
            if (triA[i] == triB[j])
                return 0;
        }
    }
    return 1;
}

double DistanceToTriangle(double TriPoints[3][3], double obsPoint[3])
{
    double dist = 1e9;
    double temp = 0.0;
    int i;
    for (i=0; i<3; i++)
    {
        temp = Distance(TriPoints[i], obsPoint);
        if (temp < dist)
            dist = temp;
    }
    return dist;
}


