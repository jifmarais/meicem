#pragma once
#include <vector>
#include "TriangleContainer.hpp"
#include "Node.hpp"

class Quadrature
{
    public:
        struct WeightedPoint{
                Node   node;
                double weight;
        };
        typedef std::vector<WeightedPoint> WeightedPointList_type;
        struct WeightedPoint1D{
                double node;
                double weight;
        };
        typedef std::vector<WeightedPoint1D> WeightedPoint1DList_type;

        Quadrature();
        virtual         ~Quadrature();

        static WeightedPoint1DList_type get1DGaussianQuadraturePoints(unsigned maxNumberOfPoints);
        static WeightedPointList_type 	getTriangleSimplexGaussianQuadraturePoints(unsigned maxNumberOfPoints);
        static WeightedPointList_type   RAR1S(const Triangle& T, Node observationPoint, unsigned maxNumberOfPoints);
        static WeightedPointList_type   RAR1S_2D(Triangle T, double offset, unsigned maxNumberOfPoints);
        static Quadrature::WeightedPointList_type getTriangleGaussianQuadraturePoints(Triangle &T, unsigned maxNumberOfPoints);
protected:
private:
        static double RAR1SuFromXY(double x, double y);
        static double RAR1SvFromYUZ(double y, double u, double z);
        static double RAR1SyFromVUZ(double v, double u, double z);
        static double RAR1SxFromYU(double y, double u);
        static ComplexMatrix getZrotationMatrix(double phi);
        static ComplexMatrix getXrotationMatrix(double theta);
};

