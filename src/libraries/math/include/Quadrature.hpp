#pragma once
#include <vector>
#include "TriangleContainer.hpp"
#include "Node.hpp"
#include <armadillo>

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

        WeightedPoint1DList_type get1DGaussianQuadraturePoints(unsigned maxNumberOfPoints);
        WeightedPointList_type 	getTriangleSimplexGaussianQuadraturePoints(unsigned maxNumberOfPoints);
        WeightedPointList_type   RAR1S(const Triangle& T, Node observationPoint, unsigned maxNumberOfPoints);
        WeightedPointList_type   RAR1S_2D(Triangle T, double offset, unsigned maxNumberOfPoints);
        Quadrature::WeightedPointList_type getTriangleGaussianQuadraturePoints(Triangle &T, unsigned maxNumberOfPoints);
protected:
private:
        static double RAR1SuFromXY(double x, double y);
        static double RAR1SvFromYUZ(double y, double u, double z);
        static double RAR1SyFromVUZ(double v, double u, double z);
        static double RAR1SxFromYU(double y, double u);
        static double RAR1Sdxdy(double u, double v, double R);
        static arma::mat getZrotationMatrix(double phi);
        static arma::mat getXrotationMatrix(double theta);

        unsigned m_cacheNumberOfPointsFor1DGaussianQuadraturePoints{0};
        WeightedPoint1DList_type m_cache1DGaussianQuadraturePoints;
        unsigned m_cacheNumberOfPointsForSimplexTriangleGaussianQuadraturePoints {0};
        WeightedPointList_type m_cacheSimplexTriangleGaussianQuadraturePoints;
};

