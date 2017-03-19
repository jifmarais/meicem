#pragma once
#include <armadillo>
#include "ComplexMatrix.hpp"
#include "TriangleContainer.hpp"
#include "EdgeContainer.hpp"
#include "PlaneWave.hpp"

class MoM
{
    public:
        MoM(TriangleContainer& tContainer);
        virtual         ~MoM();

        void            setFrequency(double freq);
        void            setTriangleContainer(const TriangleContainer& tContainer);
        arma::cx_mat    fillZmatrixTriangle(double sourceIntegrationPoints, double testIntegrationPoints);
        arma::cx_vec    calculateRHS(double numberOfIntegrationPoints, PlaneWave pw);
        void		    writeCurrentsToOS(std::string fname, arma::cx_vec solutionMatrix) const;

        arma::cx_mat    fillZmatrixTriangleEfficient(double sourceIntegrationPoints, double testIntegrationPoints);
        arma::cx_mat    fillZmatrixTriangleInefficient(double sourceIntegrationPoints, double testIntegrationPoints);
protected:
    private:
        double 					m_frequency;
        TriangleContainer& 		m_tContainer;
        std::complex<double> 	G0(const double R, const double k) const;
        double 					RWGBasisFunction(const Triangle T, const Edge E) const;
        double 					divRWGBasisFunction(const Triangle T, const double sign) const;
        double 					divRWGBasisFunction(const Triangle T, const Edge E, const double sign) const;
        double 					RWGBasisFunctionSign(unsigned baseIndex, unsigned ohterTriangleIndex) const;
};

