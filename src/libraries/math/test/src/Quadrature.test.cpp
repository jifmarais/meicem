#include "Quadrature.hpp"
#include <iostream>

//#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(Quadrature_BasicTests)

double funct1(double x, unsigned numPoints)
{
    std::vector<double> mult {5, -12, 3, 7, -15, 2, 1.1, -20.3, 1.0, 12, 7, -1.0};
    double val = 0.0;

    for (unsigned ii=0; ii < numPoints; ++ii)
    {
        val += std::pow(x, numPoints - 1 - ii) * mult.at(ii);
    }
    return val;
}

double funct2(Node n, unsigned numPoints)
{
    std::vector<double> multx {5, -12, 3, 7, -15, 2, 1.1, -20.3, 1.0, 12, 7, -1.0};
    std::vector<double> multy {15, -2, -3, 1.7, -1.5, 2.8, 1.1, -2.3, 1.9, 1.2, 5.7, -1.2};
    double val = 0.0;

    for (unsigned ii=0; ii < numPoints; ++ii)
    {
        val += std::pow(n.x(), numPoints - 1 - ii) * multx.at(ii);
        val += std::pow(n.y(), numPoints - 1 - ii) * multy.at(ii);
    }
    return val;
}

BOOST_AUTO_TEST_CASE(test1DQuadrature2Up)
{
    double answer = 9.0 + (1.0/3.0);
    unsigned min = 2;
    unsigned max = 12;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = Quadrature::get1DGaussianQuadraturePoints(numPoints);
        double integral = 0;
        for (unsigned ii=0; ii < wps.size(); ++ii)
        {
            double x = wps.at(ii).node;
            integral += wps.at(ii).weight * funct1(x, min+1);
        }
//        std::cout << numPoints<< ":  " << integral << std::endl;
        BOOST_CHECK_MESSAGE(std::fabs(integral - answer) < 1e-6, "1D quadratures failed.");
    }
}

BOOST_AUTO_TEST_CASE(test1DQuadrature3Up)
{
    double answer = 6.0;
    unsigned min = 3;
    unsigned max = 12;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = Quadrature::get1DGaussianQuadraturePoints(numPoints);
        double integral = 0;
        for (unsigned ii=0; ii < wps.size(); ++ii)
        {
            double x = wps.at(ii).node;
            integral += wps.at(ii).weight * funct1(x, min+1);
        }
//        std::cout << numPoints<< ":  " << integral << std::endl;
        BOOST_CHECK_MESSAGE(std::fabs(integral - answer) < 1e-6, "1D quadratures failed.");
    }
}

BOOST_AUTO_TEST_CASE(test1DQuadrature4Up)
{
    double answer = -26;
    unsigned min = 4;
    unsigned max = 12;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = Quadrature::get1DGaussianQuadraturePoints(numPoints);
        double integral = 0;
        for (unsigned ii=0; ii < wps.size(); ++ii)
        {
            double x = wps.at(ii).node;
            integral += wps.at(ii).weight * funct1(x, min+1);
        }
//        std::cout << numPoints<< ":  " << integral << std::endl;
        BOOST_CHECK_MESSAGE(std::fabs(integral - answer) < 1e-6, "1D quadratures failed.");
    }
}

BOOST_AUTO_TEST_CASE(test1DQuadrature5Up)
{
    double answer = 3.8666666666666667;
    unsigned min = 5;
    unsigned max = 12;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = Quadrature::get1DGaussianQuadraturePoints(numPoints);
        double integral = 0;
        for (unsigned ii=0; ii < wps.size(); ++ii)
        {
            double x = wps.at(ii).node;
            integral += wps.at(ii).weight * funct1(x, min+1);
        }
//        std::cout << numPoints<< ":  " << integral << std::endl;
        BOOST_CHECK_MESSAGE(std::fabs(integral - answer) < 1e-6, "1D quadratures failed.");
    }
}

BOOST_AUTO_TEST_CASE(test1DQuadrature6Up)
{
    double answer = -5.17142857142857;
    unsigned min = 6;
    unsigned max = 12;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = Quadrature::get1DGaussianQuadraturePoints(numPoints);
        double integral = 0;
        for (unsigned ii=0; ii < wps.size(); ++ii)
        {
            double x = wps.at(ii).node;
            integral += wps.at(ii).weight * funct1(x, min+1);
        }
//        std::cout << std::setprecision (15)<< numPoints<< ":  " << integral << std::endl;
        BOOST_CHECK_MESSAGE(std::fabs(integral - answer) < 1e-6, "1D quadratures failed.");
    }
}

BOOST_AUTO_TEST_CASE(test1DQuadrature7Up)
{
    double answer = -39.8952380952381;
    unsigned min = 7;
    unsigned max = 12;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = Quadrature::get1DGaussianQuadraturePoints(numPoints);
        double integral = 0;
        for (unsigned ii=0; ii < wps.size(); ++ii)
        {
            double x = wps.at(ii).node;
            integral += wps.at(ii).weight * funct1(x, min+1);
        }
//        std::cout << std::setprecision (15)<< numPoints<< ":  " << integral << std::endl;
        BOOST_CHECK_MESSAGE(std::fabs(integral - answer) < 1e-6, "1D quadratures failed.");
    }
}

BOOST_AUTO_TEST_CASE(test1DQuadrature8Up)
{
    double answer = -1.29841269841274;
    unsigned min = 8;
    unsigned max = 12;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = Quadrature::get1DGaussianQuadraturePoints(numPoints);
        double integral = 0;
        for (unsigned ii=0; ii < wps.size(); ++ii)
        {
            double x = wps.at(ii).node;
            integral += wps.at(ii).weight * funct1(x, min+1);
        }
//        std::cout << std::setprecision (15)<< numPoints<< ":  " << integral << std::endl;
        BOOST_CHECK_MESSAGE(std::fabs(integral - answer) < 1e-6, "1D quadratures failed.");
    }
}

BOOST_AUTO_TEST_CASE(test1DQuadrature9Up)
{
    double answer = 10.5999999999999;
    unsigned min = 9;
    unsigned max = 12;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = Quadrature::get1DGaussianQuadraturePoints(numPoints);
        double integral = 0;
        for (unsigned ii=0; ii < wps.size(); ++ii)
        {
            double x = wps.at(ii).node;
            integral += wps.at(ii).weight * funct1(x, min+1);
        }
//        std::cout << std::setprecision (15)<< numPoints<< ":  " << integral << std::endl;
        BOOST_CHECK_MESSAGE(std::fabs(integral - answer) < 1e-6, "1D quadratures failed.");
    }
}

BOOST_AUTO_TEST_CASE(test1DQuadrature10Up)
{
    double answer = 12.3967099567099;
    unsigned min = 10;
    unsigned max = 12;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = Quadrature::get1DGaussianQuadraturePoints(numPoints);
        double integral = 0;
        for (unsigned ii=0; ii < wps.size(); ++ii)
        {
            double x = wps.at(ii).node;
            integral += wps.at(ii).weight * funct1(x, min+1);
        }
//        std::cout << std::setprecision (15)<< numPoints<< ":  " << integral << std::endl;
        BOOST_CHECK_MESSAGE(std::fabs(integral - answer) < 1e-6, "1D quadratures failed.");
    }
}

BOOST_AUTO_TEST_CASE(test1DQuadrature11Up)
{
    double answer = -2.1748340548341;
    unsigned min = 11;
    unsigned max = 12;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = Quadrature::get1DGaussianQuadraturePoints(numPoints);
        double integral = 0;
        for (unsigned ii=0; ii < wps.size(); ++ii)
        {
            double x = wps.at(ii).node;
            integral += wps.at(ii).weight * funct1(x, min+1);
        }
//        std::cout << std::setprecision (15)<< numPoints<< ":  " << integral << std::endl;
        BOOST_CHECK_MESSAGE(std::fabs(integral - answer) < 1e-6, "1D quadratures failed.");
    }
}

BOOST_AUTO_TEST_CASE(testRAR1S_2D_1)
{
    double answer = 1.0;
    unsigned min = 1;
    unsigned max = 1;
    double offset = 0.0;

    Node n1 {-1.0, 0.0, 0.0};
    Node n2 { 1.0, 0.0, 0.0};
    Node n3 { 0.0, 1.0, 0.0};
    Triangle T {n1, n2, n3};
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPointList_type wps = Quadrature::RAR1S_2D(T, offset, numPoints);
        double integral = 0;
        for (unsigned ii=0; ii < wps.size(); ++ii)
        {
            Node n = wps.at(ii).node;
            integral += wps.at(ii).weight * 1.0; //funct2(n, 1);
        }
        std::cout << std::setprecision (15)<< numPoints<< ":  " << integral << std::endl;
        BOOST_CHECK_MESSAGE(std::fabs(integral - answer) < 1e-6, "RAR1S 2D failed.");
    }
}


BOOST_AUTO_TEST_SUITE_END()

//EOF

