#pragma once
#include "ComplexMatrix.hpp"
#include "TriangleContainer.hpp"

class MoM 
{
    public:
        MoM(TriangleContainer& tContainer);
        virtual         ~MoM();

        void            setFrequency(double freq);
        void            setTriangleContainer(const TriangleContainer& tContainer);
        ComplexMatrix 	fillZmatrixTriangle();
        ComplexMatrix 	calculateRHS();

    protected:
    private:
        double m_frequency;
        TriangleContainer& m_tContainer;

        ComplexMatrix fillZmatrixTriangleInefficient();
};

