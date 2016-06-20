#pragma once
#include <vector>
#include "Node.hpp"

class Quadrature
{
    public:
        struct WeightedPoint{
                Node   node;
                double weight;
        };
        typedef std::vector<WeightedPoint> WeightedPointList_type;

        Quadrature();
        virtual         ~Quadrature();

        static WeightedPointList_type getTriangleSimplexGaussianQuadraturePoints(unsigned maxNumberOfPoints);

    protected:
    private:
};

