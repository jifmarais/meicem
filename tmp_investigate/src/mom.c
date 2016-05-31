/*
   ============================================================================
Name        : mom.c
Author      : JIF Marais
Version     :
Description : 
Adapted from code written by WJ Strydom
============================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "libmath.h"
#include "libtriangles.h"

/*
   return the integration points and weights based on a gaussian rule
   */
void GaussianQuadrature(double points[][3], int numIntPoints, double **returnPoints)
{
    /*
       Computes quadrature points on a triangle defined by points.
       Returns points and weight.
       */
    int i, j, ii;

    double quadPoints1 [1][4] = {{1.0/3, 1.0/3, 1.0/3, 1}};
    double quadPoints3 [1][4] = {{2.0/3, 1.0/6, 1.0/6, 1.0/3}};
    double quadPoints4 [2][4] = {{1.0/3, 1.0/3, 1.0/3, -0.5625},
        {0.6, 0.2, 0.2, 0.520833333333333}};
    double quadPoints6 [2][4] = {{0.816847572980459, 0.091576213509771, 0.091576213509771, 0.109951743655322},
        {0.108103018168070, 0.445948490915965, 0.445948490915965, 0.223381589678011}};
    double quadPoints7 [3][4] = {{1.0/3, 1.0/3, 1.0/3, 0.225},
        {0.797426985353087, 0.101286507323456, 0.101286507323456, 0.125939180544827},
        {0.059715871789770, 0.470142064105115, 0.470142064105115, 0.132394152788506}};
    double quadPoints12 [4][4] = {{0.873821971016996, 0.063089014491502, 0.063089014491502, 0.050844906370207},
        {0.501426509658179, 0.249286745170910, 0.249286745170910, 0.116786275726379},
        {0.636502499121399, 0.310352451033785, 0.053145049844816, 0.082851075618374},
        {0.636502499121399, 0.053145049844816, 0.310352451033785, 0.082851075618374}};
    double quadPoints13 [5][4] = {{1.0/3, 1.0/3, 1.0/3, -0.149570044467682},
        {0.479308067841920, 0.260345966079040, 0.260345966079040, 0.175615257433208},
        {0.869739794195568, 0.065130102902216, 0.065130102902216, 0.053347235608838},
        {0.048690315425316, 0.312865496004874, 0.638444188569810, 0.077113760890257},
        {0.048690315425316, 0.638444188569810, 0.312865496004874, 0.077113760890257}};
    double quadPoints16 [6][4] = {{1.0/3, 1.0/3, 1.0/3, 0.144315607677787},
        {0.081414823414554, 0.459292588292723, 0.459292588292723, 0.095091634267284},
        {0.658861384496478, 0.170569307751761, 0.170569307751761, 0.103217370534718},
        {0.898905543365938, 0.050547228317031, 0.050547228317031, 0.032458497623198},
        {0.008394777409958, 0.263112829634638, 0.728492392955404, 0.027230314174435},
        {0.008394777409958, 0.728492392955404, 0.263112829634638, 0.027230314174435}};
    double quadPoints25 [9][4] = {{0.333333333333333, 0.333333333333333, 0.333333333333333, 0.090817990382754},
        {0.028844733232685, 0.485577633383657, 0.485577633383657, 0.036725957756467},
        {0.781036849029926, 0.109481575485037, 0.109481575485037, 0.045321059435528},
        {0.141707219414880, 0.307939838764121, 0.550352941820999, 0.072757916845420},
        {0.141707219414880, 0.550352941820999, 0.307939838764121, 0.072757916845420},
        {0.025003534762686, 0.246672560639903, 0.728323904597411, 0.028327242531057},
        {0.025003534762686, 0.728323904597411, 0.246672560639903, 0.028327242531057},
        {0.009540815400299, 0.066803251012200, 0.923655933587500, 0.009421666963733},
        {0.009540815400299, 0.923655933587500, 0.066803251012200, 0.009421666963733}};
    double quadPoints33 [11][4] = {{0.023565220452390, 0.488217389773805, 0.488217389773805, 0.025731066440455},
        {0.120551215411079, 0.439724392294460, 0.439724392294460, 0.043692544538038},
        {0.457579229975768, 0.271210385012116, 0.271210385012116, 0.062858224217885},
        {0.744847708916828, 0.127576145541586, 0.127576145541586, 0.034796112930709},
        {0.957365299093579, 0.021317350453210, 0.021317350453210, 0.006166261051559},
        {0.115343494534698, 0.608943235779788, 0.275713269685514, 0.040371557766381},
        {0.115343494534698, 0.275713269685514, 0.608943235779788, 0.040371557766381},
        {0.022838332222257, 0.695836086787803, 0.281325580989940, 0.022356773202303},
        {0.022838332222257, 0.281325580989940, 0.695836086787803, 0.022356773202303},
        {0.025734050548330, 0.858014033544073, 0.116251915907597, 0.017316231108659},
        {0.025734050548330, 0.116251915907597, 0.858014033544073, 0.017316231108659}};



    double (*qpPointer)[][4] = NULL;

    int qpNum;

    if (numIntPoints > 16 && numIntPoints < 25)
        numIntPoints = 16;

    switch (numIntPoints)
    {
        case 1:
        case 2:
            qpPointer = &quadPoints1;
            qpNum = 1;
            break;
        case 3:
            qpPointer = &quadPoints3;
            qpNum = 1;
            break;
        case 4:
        case 5:
            qpPointer = &quadPoints4;
            qpNum = 2;
            break;
        case 6:
            qpPointer = &quadPoints6;
            qpNum = 2;
            break;
        case 7:
        case 8:
        case 9:
        case 10:
        case 11:
            qpPointer = &quadPoints7;
            qpNum = 3;
            break;
        case 12:
            qpPointer = &quadPoints12;
            qpNum = 4;
            break;
        case 13:
        case 14:
        case 15:
            qpPointer = &quadPoints13;
            qpNum = 5;
            break;
        case 16:
            qpPointer = &quadPoints16;
            qpNum = 6;
            break;
        case 25:
            qpPointer = &quadPoints25;
            qpNum = 9;
            break;
        case 33:
        default:
            qpPointer = &quadPoints33;
            qpNum = 11;
            break;
    }

    double AB[3] = {points[1][0]-points[0][0],
        points[1][1]-points[0][1],
        points[1][2]-points[0][2]};

    double AC[3] = {points[2][0]-points[0][0],
        points[2][1]-points[0][1],
        points[2][2]-points[0][2]};

    double Area = pow(AB[1]*AC[2] - AB[2]*AC[1], 2);
    Area += pow(AB[2]*AC[0] - AB[0]*AC[2], 2);
    Area += pow(AB[0]*AC[1] - AB[1]*AC[0], 2);
    Area = 0.5*sqrt(Area);

    int place = 0;
    int iter;

    for (i=0; i<qpNum; i++)
    {
        if ((*qpPointer)[i][0] != (*qpPointer)[i][1])
            iter = 3;
        else
            iter = 1;

        for (j=0; j<iter; j++)
        {
            for(ii=0; ii<3; ii++)
            {
                returnPoints[place][ii] = (*qpPointer)[i][j]*points[0][ii] + (*qpPointer)[i][(j+1)%3]*points[1][ii] +
                    (*qpPointer)[i][(j+2)%3]*points[2][ii];
            }
            returnPoints[place][3] = Area*(*qpPointer)[i][3];
            place++;
        }
    }

    //free(qpPointer2);
}


/*
   Functions for Radial Angular R1 Sqrt
   */
double wFromXY (double x, double y)
{
    return asinh(x/y);
}
double qFromYWZ (double y, double w, double z)
{
    return sqrt(sqrt(pow(cosh(w),2)*y*y + z*z) - fabs(z));
}
double yFromQWZ (double q, double w, double z)
{
    return sqrt((pow((q*q + fabs(z)),2) - z*z)/pow(cosh(w),2));
}
double xFromYW (double y, double w)
{
    return y*sinh(w);
}

/*
   Calculates the integral on subtriangles.
   The assuption is that the third verice, v3, is the projection of the singularity
   */
void RAR1S_2D(double v1[3], double v2[3], double v3[3], double crux, int numIntPoints, double **returnPoints)
{
    int i,j;

    //Gaussian Quadrature rules in 1D
    double quadPoints1[1][2] = {{2,0}};
    double quadPoints2[2][2] = {{1.0000000000000000,    -0.5773502691896257},
        {1.0000000000000000,    0.5773502691896257}};
    double quadPoints3[3][2] = {{0.8888888888888888,    0.0000000000000000},
        {0.5555555555555556,    -0.7745966692414834},
        {0.5555555555555556,    0.7745966692414834}};
    double quadPoints4[4][2] = {{0.6521451548625461,    -0.3399810435848563},
        {0.6521451548625461,    0.3399810435848563},
        {0.3478548451374538,    -0.8611363115940526},
        {0.3478548451374538,    0.8611363115940526}};
    double quadPoints5[5][2] = {{0.5688888888888889,    0.0000000000000000},
        {0.4786286704993665,    -0.5384693101056831},
        {0.4786286704993665,    0.5384693101056831},
        {0.2369268850561891,    -0.9061798459386640},
        {0.2369268850561891,    0.9061798459386640}};
    double quadPoints6[6][2] = {{0.3607615730481386,    0.6612093864662645},
        {0.3607615730481386,    -0.6612093864662645},
        {0.4679139345726910,    -0.2386191860831969},
        {0.4679139345726910,    0.2386191860831969},
        {0.1713244923791704,    -0.9324695142031521},
        {0.1713244923791704,    0.9324695142031521}};
    double quadPoints7[7][2] = {{0.4179591836734694,    0.0000000000000000},
        {0.3818300505051189,    0.4058451513773972},
        {0.3818300505051189,    -0.4058451513773972},
        {0.2797053914892766,    0.7415311855993945},
        {0.2797053914892766,    -0.7415311855993945},
        {0.1294849661688697,    0.9491079123427585},
        {0.1294849661688697,    -0.9491079123427585}};
    double quadPoints8[8][2] = {{0.10122853629038, -0.96028985649754},
        {0.22238103445338, -0.79666647741362},
        {0.31370664587788, -0.52553240991633},
        {0.36268378337836, -0.18343464249565},
        {0.36268378337836, 0.18343464249565},
        {0.31370664587788, 0.52553240991633},
        {0.22238103445338, 0.79666647741362},
        {0.10122853629038, 0.96028985649754}};
    double quadPoints9[9][2] = {{0.08127438836157, -0.96816023950763},
        {0.18064816069487, -0.83603110732663},
        {0.26061069640291, -0.61337143270059},
        {0.31234707704002, -0.32425342340381},
        {0.33023935500125, 0.0},
        {0.31234707704001, 0.32425342340381},
        {0.26061069640292, 0.61337143270059},
        {0.18064816069487, 0.83603110732663},
        {0.08127438836157, 0.96816023950763}};

    double quadPoints10[10][2] = {{0.06667134430869, -0.97390652851717},
        {0.14945134915058, -0.86506336668899},
        {0.21908636251599, -0.67940956829902},
        {0.26926671930998, -0.43339539412925},
        {0.29552422471476, -0.14887433898163},
        {0.29552422471475, 0.14887433898163},
        {0.26926671930999, 0.43339539412925},
        {0.21908636251599, 0.67940956829902},
        {0.14945134915058, 0.86506336668899},
        {0.06667134430869, 0.97390652851717}};




    //Gaussian Quadrature rules in 1D

    int qpNum;
    double (*qpPointer)[][2];

    switch (numIntPoints)
    {
        case 1:
            qpPointer = &quadPoints1;
            qpNum=1;
            break;
        case 2:
            qpPointer = &quadPoints2;
            qpNum=2;
            break;
        case 3:
            qpPointer = &quadPoints3;
            qpNum=3;
            break;
        case 4:
            qpPointer = &quadPoints4;
            qpNum=4;
            break;
        case 5:
            qpPointer = &quadPoints5;
            qpNum=5;
            break;
        case 6:
            qpPointer = &quadPoints6;
            qpNum=6;
            break;
        case 7:
            qpPointer = &quadPoints7;
            qpNum=7;
            break;
        case 8:
            qpPointer = &quadPoints8;
            qpNum=8;
            break;
        case 9:
            qpPointer = &quadPoints9;
            qpNum=9;
            break;
        case 10:
        default:
            qpNum=10;
            qpPointer = &quadPoints10;
            break;
    }

    double tempV1[3], newV1[3], tempV2[3], newV2[3], pVec[3];

    VectorSubtract(v3, v1, tempV1);        //set origin
    VectorSubtract(v3, v2, tempV2);
    VectorSubtract(tempV1, tempV2, pVec);   //non-origin edge as a vector

    double theta;      //angle to get edge parallel to x-axis
    if (pVec[0]==0)
        theta = M_PI*0.5;
    else
        theta = atan(-1.0*pVec[1]/pVec[0]);

    double rotationMatrix[3][3] = {{cos(theta), -1.0*sin(theta),0},
        {sin(theta), cos(theta),0},
        {0,0,1}};

    MatrixMultiplyVector(rotationMatrix, 3, tempV1, newV1);
    MatrixMultiplyVector(rotationMatrix, 3, tempV2, newV2);         //rotate non-origin points


    int flip = 1.0;
    if (newV1[1] < 0)           //flip triangle if upside down
        flip = -1.0;
    newV1[1] *= flip;
    newV2[1] *= flip;

    rotationMatrix[0][1] *= -1;
    rotationMatrix[1][0] *= -1;                        //rotationMatrix will now undo rotation

    double *leftPoint, *rightPoint;
    if (newV1[0] < newV2[0])
    {
        leftPoint = &newV1[0];
        rightPoint = &newV2[0];
    }
    else
    {
        leftPoint = &newV2[0];              //left most point of edge for integration lower limit
        rightPoint = &newV1[0];
    }

    double w_lower = wFromXY(leftPoint[0], leftPoint[1]);
    double w_upper = wFromXY(rightPoint[0], rightPoint[1]);
    double w_range = w_upper - w_lower;

    for (i=0; i < qpNum; i++)
    {
        double wPoint = w_lower + 0.5*w_range*(1 + (*qpPointer)[i][1]);

        for (j=0; j < qpNum; j++)
        {
            double q_lower = qFromYWZ(0,wPoint,fabs(crux));
            double q_upper = qFromYWZ(leftPoint[1], wPoint, crux);
            double q_range = q_upper - q_lower;

            double qPoint = q_lower + 0.5*q_range*(1 + (*qpPointer)[j][1]);
            double yPoint = yFromQWZ(qPoint, wPoint, crux);
            double xPoint = xFromYW(yPoint, wPoint);

            double tempNewPoint[3] = {xPoint, flip*yPoint, 0.0};
            MatrixMultiplyVector(rotationMatrix, 3, tempNewPoint, returnPoints[i*qpNum + j]);

            returnPoints[i*qpNum + j][0] += v3[0];
            returnPoints[i*qpNum + j][1] += v3[1];
            returnPoints[i*qpNum + j][2] += v3[2];//rotate and translate back to original triangle coordinate system

            returnPoints[i*qpNum + j][3] = (*qpPointer)[j][0]*(*qpPointer)[i][0]*q_range*w_range*0.25*2*qPoint/cosh(wPoint);   //append weighting of quadpoint
        }
    }
}

void RAR1S(double points[][3], double obsPoint[], int numIntPoints, double **returnPoints)
{
    int i, j;
    double triNorm[3], AB[3], AC[3];

    VectorSubtract(points[0], points[1], AB);
    VectorSubtract(points[0], points[2], AC);
    VectorCross(AB, AC, triNorm);                //Find norm of triangle

    if (triNorm[2] < 0)
        MultiplyVectorByConstant(triNorm, 3, -1.0);
    MultiplyVectorByConstant(triNorm, 3, 1/VectorSize(triNorm));

    double thetaNorm = acos(triNorm[2]);
    double phiNorm;

    /*
     *         ThetaNorm and PhiNorm are the required rotation angles not the actual angles in polar coordinates of the norm vector
     */

    if (triNorm[1]== 0)
    {
        if (triNorm[0]>0)
            phiNorm = M_PI/2;
        else
            phiNorm = -M_PI/2;
    }
    else
        phiNorm = atan(triNorm[0]/triNorm[1]);

    if (triNorm[0]<0 && triNorm[1]<0)
        phiNorm += M_PI;
    else if (triNorm[1]<0)
        phiNorm -= M_PI;

    double Rz[3][3] = {{cos(phiNorm), -sin(phiNorm),0},
        {sin(phiNorm), cos(phiNorm), 0},
        {0, 0, 1}};

    double Rx[3][3] = {{1, 0, 0},
        {0, cos(thetaNorm), -sin(thetaNorm)},
        {0, sin(thetaNorm), cos(thetaNorm)}};

    double rotationMatrix[3][3];
    MatrixMultiplyMatrix(Rx, Rz, rotationMatrix);
    double rotationMatrix_inv[3][3];
    MatrixInverse(rotationMatrix, 3, rotationMatrix_inv);

    double newPoints[3][3];

    for (i=0; i< 3; i++)
        MatrixMultiplyVector(rotationMatrix, 3, points[i], newPoints[i]);  //triangle is now rotated

    double zOffset = newPoints[0][2];
    for (i=0; i< 3; i++)
        newPoints[i][2]=0.0;

    double newObsPoint[3];
    MatrixMultiplyVector(rotationMatrix, 3, obsPoint, newObsPoint);       //rotate observation point
    double zObsPoint = newObsPoint[2];
    newObsPoint[2]=0;

    double triCentre[3] = {newPoints[0][0]+newPoints[1][0]+newPoints[2][0],
        newPoints[0][1]+newPoints[1][1]+newPoints[2][1],
        newPoints[0][2]+newPoints[1][2]+newPoints[2][2]};
    MultiplyVectorByConstant(triCentre, 3, 1.0/3);                            //find centre of triangle

    for (i=0; i<3; i++)          //split triangle into three at the projection of the singularity
    {
        double AB[3];
        VectorSubtract(newPoints[i], newPoints[(i+1)%3], AB);

        double AtoCentre[3];
        VectorSubtract(newPoints[i], triCentre, AtoCentre);

        double AtoObs[3];
        VectorSubtract(newPoints[i], newObsPoint, AtoObs);

        double vec1[3], vec2[3];
        VectorCross(AB, AtoCentre, vec1);
        VectorCross(AB, AtoObs, vec2);

        int newNumInt = numIntPoints;
        if (newNumInt>10)
            newNumInt=10;

        double **returnPoints2D = NULL;
        returnPoints2D = malloc(sizeof(double *) * (int)pow(newNumInt,2));
        for (j=0; j < pow(newNumInt, 2); j++)
            returnPoints2D[j] = (double *)malloc(sizeof(double) * 4);

        RAR1S_2D(newPoints[i], newPoints[(i+1)%3], newObsPoint, zObsPoint-zOffset, newNumInt, returnPoints2D);  //get integration points

        int posTri = 1;
        if (vec1[2]*vec2[2] < 0)
            posTri = -1;

        newObsPoint[2] += zObsPoint;
        for (j=0; j < pow(newNumInt,2); j++)
        {
            returnPoints2D[j][2] += zOffset;
            double R = Distance(returnPoints2D[j], newObsPoint);

            MatrixMultiplyVector(rotationMatrix_inv, 3, returnPoints2D[j], returnPoints[(int)pow(newNumInt,2)*i + j]);
            returnPoints[(int)pow(newNumInt,2)*i + j][3] = posTri*R*returnPoints2D[j][3];
        }
        newObsPoint[2] -= zObsPoint;
        for (j=0; j < pow(newNumInt,2); j++)
            free(returnPoints2D[j]);
        free(returnPoints2D);
    }
}

/*
   check if the edge defined by [val1, val2] is a non-boundary edge.
   Return index+1, else return -1
   */
int IsEdgeInMat(int val1, int val2, int **edges, int N)
{
    int i;
    for (i=0; i < N; i++)
    {
        if (edges[i][0] == (int)fmin(val1, val2) && edges[i][1] == (int)fmax(val1,val2))
            return i;
    }
    return -1;
}


void MoM(double freq, int P, double **points, int T, int **triangles, int N, int **edges, complex double **Zmat, complex double *Vvec)
{
    int NumOuterIntPoints = 3;
    int NumInnerIntPoints = 7;
    double prox = 0.33;

    const static double c0 = 299792456.2;
    const static double mu = 1.25663706143592e-06;
    const static double eps = 8.85418781761e-12;
    double k = (2*M_PI*freq)/c0;
    double w = 2*M_PI*freq;

    int i, j, ii, jj, iter;
    int oip, iip;

    int ccount = 0; //used to count occurrences of cancellation quadrature

    double prog = 0.0; //used to monitor progress

    for (i=0; i<T; i++) //Observation triangle
    {
        if ((double)(i+2)/T > prog)
        {
            printf("\r%d%% Completed", (int)(100*prog));
            prog = (double)(i+2)/T;
            fflush(stdout);
        }

        double pPoints[3][3] = {{points[triangles[i][0]][0], points[triangles[i][0]][1], points[triangles[i][0]][2]},
            {points[triangles[i][1]][0], points[triangles[i][1]][1], points[triangles[i][1]][2]},
            {points[triangles[i][2]][0], points[triangles[i][2]][1], points[triangles[i][2]][2]}};

        double pArea = TriangleArea(pPoints);
        double **OuterIntPoints;
        OuterIntPoints = malloc(sizeof(double *) * NumOuterIntPoints);
        for (iter=0; iter < NumOuterIntPoints; iter++)
            OuterIntPoints[iter] = malloc(sizeof(double) * 4);

        GaussianQuadrature(pPoints, NumOuterIntPoints, OuterIntPoints);
        for (j=0; j<T; j++) //Testing triangle
        {
            double qPoints[3][3] = {{points[triangles[j][0]][0], points[triangles[j][0]][1], points[triangles[j][0]][2]},
                {points[triangles[j][1]][0], points[triangles[j][1]][1], points[triangles[j][1]][2]},
                {points[triangles[j][2]][0], points[triangles[j][2]][1], points[triangles[j][2]][2]}};

            double qArea = TriangleArea(qPoints);
            double centrePoint[3];
            TriangleCentre(qPoints, centrePoint);

            double **InnerIntPoints = NULL;
            int newNumInt = 0;
            int doGauss = TrianglesDontTouch(triangles[i], triangles[j]);

            if (doGauss)
            {
                newNumInt = NumInnerIntPoints;
                InnerIntPoints = malloc(sizeof(double *) * newNumInt);
                for (iter=0; iter < newNumInt; iter++)
                    InnerIntPoints[iter] = malloc(sizeof(double) * 4);
                GaussianQuadrature(qPoints, NumInnerIntPoints, InnerIntPoints);
            }
            else
            {
                newNumInt = 3*(int)pow(NumInnerIntPoints,2);
                if (newNumInt > 300)
                    newNumInt = 300;
                InnerIntPoints = malloc(sizeof(double *)*newNumInt);
                for (iter=0; iter < newNumInt; iter++)
                    InnerIntPoints[iter] = malloc(sizeof(double) * 4);
            }


            for (ii=0; ii<3; ii++) //Edges of observation triangle
            {
                int pEdgeIndex = IsEdgeInMat(triangles[i][ii], triangles[i][(ii+1)%3], edges, N);
                if (pEdgeIndex + 1)
                {
                    double lp = Distance(pPoints[ii], pPoints[(ii+1)%3]);

                    //outer integral
                    for (oip=0; oip<NumOuterIntPoints; oip++) //outer integral points
                    {
                        double proximity = DistanceToTriangle(qPoints, OuterIntPoints[oip]);
                        //Distance(OuterIntPoints[oip], centrePoint);
                        if (!doGauss)
                        {
                            if ((proximity > prox*Distance(qPoints[0], qPoints[1])) && (i != j))
                            {
                                GaussianQuadrature(qPoints, 25, InnerIntPoints);
                                newNumInt = 25;
                            }
                            else
                            {
                                RAR1S(qPoints, OuterIntPoints[oip], NumInnerIntPoints, InnerIntPoints);


                                if(pEdgeIndex==581 && IsEdgeInMat(triangles[j][0], triangles[j][1], edges, N)==584)
                                {
                                    printf("%lf %lf, %d %d\n", OuterIntPoints[oip][0], OuterIntPoints[oip][1], pEdgeIndex, IsEdgeInMat(triangles[j][0], triangles[j][1], edges, N));
                                }


                                newNumInt = 3*(int)pow(NumInnerIntPoints,2);
                                if (newNumInt > 300)
                                    newNumInt = 300;
                                ccount ++;
                            }
                        }

                        double pRho[3];
                        VectorSubtract(pPoints[(ii+2)%3], OuterIntPoints[oip], pRho);
                        MultiplyVectorByConstant(pRho, 3, triangles[i][3+(ii+2)%3]);

                        if (j==0)
                            Vvec[pEdgeIndex] += (OuterIntPoints[oip][3]/pArea)*(-lp/2)*(cexp(I*k*OuterIntPoints[oip][2])*pRho[0]);

                        for (jj=0; jj<3; jj++) //Edges of testing triangle
                        {
                            /*if (i==0 && j==0)
                              printf("\n%lf", Distance(qPoints[jj], qPoints[(jj+1)%3]));
                              */

                            int qEdgeIndex = IsEdgeInMat(triangles[j][jj], triangles[j][(jj+1)%3], edges, N);
                            if (qEdgeIndex + 1)
                            {
                                complex double A[3] = {0.0, 0.0, 0.0};
                                complex double Phi = 0+0*I;

                                //inner integral
                                for (iip=0; iip < newNumInt; iip++)
                                {
                                    double qRho[3];
                                    VectorSubtract(qPoints[(jj+2)%3], InnerIntPoints[iip], qRho);
                                    MultiplyVectorByConstant( qRho, 3, triangles[j][3+(jj+2)%3]);

                                    double Rm = Distance(OuterIntPoints[oip], InnerIntPoints[iip]);
                                    double lq = Distance(qPoints[jj], qPoints[(jj+1)%3]);

                                    complex double temp = (InnerIntPoints[iip][3]/(Rm*qArea))*cexp(-I*k*Rm);

                                    A[0] += (qRho[0]*I*lq*mu*w*temp)/(16*M_PI);
                                    A[1] += (qRho[1]*I*lq*mu*w*temp)/(16*M_PI);
                                    A[2] += (qRho[2]*I*lq*mu*w*temp)/(16*M_PI);
                                    Phi += triangles[j][3+(jj+2)%3]*temp*lq/(4*I*M_PI*w*eps);
                                }
                                complex double AdotpRho = A[0]*pRho[0] + A[1]*pRho[1] + A[2]*pRho[2];
                                Zmat[pEdgeIndex][qEdgeIndex] += (OuterIntPoints[oip][3]/pArea)*lp*(AdotpRho + triangles[i][3+(ii+2)%3]*Phi);

                            }
                        }
                    }
                }
            }
            for (iter=0; iter < newNumInt; iter++)
                free(InnerIntPoints[iter]);
            free(InnerIntPoints);
        }
        for (iter=0; iter < NumOuterIntPoints; iter++)
            free(OuterIntPoints[iter]);
        free(OuterIntPoints);
    }
    printf("\nNumber of times RAR1S used: %d\n", ccount/(3*NumOuterIntPoints));
}


int main(void)
{
    /*double points[3][3] = {{0.0,0.0,0.0},
      {10.0,12.0,0.0},
      {20.0,2.0,0.0}};
      double z_tilde[3] = {15.1, 3.1, 1.95};

      int itercancellation

      double **rpg;
      rpg = malloc(sizeof(double *) * 1);
      rpg[0] = malloc(sizeof(double)*4);

      double **rpr;
      rpr = malloc(sizeof(double *)*27);
      for (iter=0; iter < 27; iter++)
      rpr[iter] = malloc(sizeof(double) * 4);


      GaussianQuadrature(points,1,rpg);
      RAR1S(points,z_tilde,3,rpr);

      int i;
      double area = 0.0;
      for (i=0; i<27; i++)
      {
      area += rpr[i][3];
      printf("%lf, %lf, %lf\n",rpr[i][0],rpr[i][1],rpr[i][2]);

      }
      printf("g: %lf, r: %lf\n", rpg[0][3], area);*/





    FILE *myfile;
    int N, P, T;  //N edges (Zmat matrix is NxN) P points and T triangles
    double freq;
    double **points = NULL;
    int **triangles = NULL;
    int **edges = NULL;

    complex double **Zmat = NULL;
    complex double *Vvec = NULL;

    int i=0, j=0;

    myfile = fopen("../inputoutput/input.txt","r");
    if (!myfile)
    {
        printf("No such file (1)");
        return 1;
    }


    fscanf(myfile, "%lf", &freq);

    fscanf(myfile, "%d", &P);
    points = malloc(sizeof(double *) * P);
    for (i=0; i < P; i++)
    {
        points[i] = (double *)malloc(sizeof(double) * 3);
        fscanf(myfile, "%lf %lf %lf", &points[i][0], &points[i][1], &points[i][2]);
    }

    fscanf(myfile, "%d", &T);
    triangles = malloc(sizeof(int *) * T);
    for (i=0; i < T; i++)
    {
        triangles[i] = malloc(sizeof(int) * 6);
        fscanf(myfile, "%d %d %d %d %d %d", &triangles[i][0], &triangles[i][1], &triangles[i][2],
                &triangles[i][3], &triangles[i][4], &triangles[i][5]);
    }

    fscanf(myfile, "%d", &N);
    Zmat = malloc(sizeof(complex double *) * N);
    Vvec = malloc(sizeof(complex double) * N);
    edges = malloc(sizeof(int *) * N);

    for (i=0; i < N; i++)
    {
        Zmat[i] = (complex double *)malloc(sizeof(complex double) * N);
        edges[i] = (int *)malloc(sizeof(int) * 2);
        fscanf(myfile, "%d %d", &edges[i][0], &edges[i][1]);
    }

    fclose(myfile);

    MoM(freq, P, points, T, triangles, N, edges, Zmat, Vvec);


    myfile = fopen("../inputoutput/output.txt","w");
    if (!myfile)
    {
        printf("No such file (2)");
        return 1;
    }


    for (i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            fprintf(myfile, "%d--", i);
            fprintf(myfile, "%d--", j);
            fprintf(myfile, "%.15g--", creal(Zmat[i][j]));
            fprintf(myfile, "%.15g--\n", cimag(Zmat[i][j]));
        }
    }
    for (i=0; i<N; i++)
    {
        fprintf(myfile, "%d--", i);
        fprintf(myfile, "%.15g--", creal(Vvec[i]));
        fprintf(myfile, "%.15g--\n", cimag(Vvec[i]));
    }


    for(i=0; i < P; i++)
        free(points[i]);

    for(i=0; i < T; i++)
        free(triangles[i]);

    for(i = 0; i < N; i++)
    {
        free(Zmat[i]);
        free(edges[i]);
    }
    free(points);
    free(triangles);
    free(edges);
    free(Zmat);
    free(Vvec);

    printf("\nIf it were done when 'tis done, then 'twere well It were done quickly\n");
    return 0;
}
