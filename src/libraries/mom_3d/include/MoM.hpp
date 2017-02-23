#pragma once
#include "ComplexMatrix.hpp"
#include "TriangleContainer.hpp"
#include "EdgeContainer.hpp"

class MoM
{
    public:
        MoM(TriangleContainer& tContainer);
        virtual         ~MoM();

        void            setFrequency(double freq);
        void            setTriangleContainer(const TriangleContainer& tContainer);
        ComplexMatrix 	fillZmatrixTriangle(double sourceIntegrationPoints, double testIntegrationPoints);
        ComplexMatrix 	calculateRHS(double numberOfIntegrationPoints);
        void		writeCurrentsToOS(std::string fname, ComplexMatrix solutionMatrix) const;

        ComplexMatrix fillZmatrixTriangleEfficient();
        ComplexMatrix fillZmatrixTriangleInefficient1();
        ComplexMatrix fillZmatrixTriangleInefficient(double sourceIntegrationPoints, double testIntegrationPoints);
protected:
    private:
        double 					m_frequency;
        TriangleContainer& 		m_tContainer;
        std::complex<double> 	G0(const double R, const double k) const;
        double 					RWGBasisFunction(const Triangle T) const;
        double 					divRWGBasisFunction(const Triangle T, const double sign) const;
};

