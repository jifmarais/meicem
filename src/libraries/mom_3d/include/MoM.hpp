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
        void            setNumberOfSourceIntegrationPoints(unsigned n);
        void            setNumberOfTestIntegrationPoints(unsigned n);
        void            setTriangleContainer(const TriangleContainer& tContainer);
        arma::cx_mat    fillZmatrixTriangle();
        arma::cx_vec    calculateRHS(PlaneWave pw);
        void		    writeCurrentsToOS(std::string fname, const arma::cx_vec& solutionMatrix) const;

        arma::cx_mat    fillZmatrixTriangleEfficient();
        arma::cx_mat    fillZmatrixTriangleInefficient();
protected:
    private:
        double 					m_frequency;
        double 					m_omega;
        double 					m_k;
        TriangleContainer& 		m_tContainer;
        unsigned				m_numberOfSourceIntegrationPoints;
        unsigned				m_numberOfTestIntegrationPoints;

        std::complex<double> 	G0(const double R) const;
        double 					RWGBasisFunction(const Triangle &T, const Edge &E) const;
        double 					divRWGBasisFunction(const Triangle& T, const Edge& E, const double sign) const;
        double 					RWGBasisFunctionSign(unsigned baseIndex, unsigned ohterTriangleIndex) const;
};

