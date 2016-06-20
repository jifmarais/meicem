#ifndef libtriangles
#define libtriangles
    
double Distance(double a[], double b[]);
double TriangleArea(double points[][3]);
void TriangleCentre (double points[][3], double centre[]);
int TrianglesDontTouch(int *triA, int *triB);
double DistanceToTriangle(double TriPoints[3][3], double obsPoint[3]);


#endif
