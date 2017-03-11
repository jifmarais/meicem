#include <string>
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

double funct2(double x, unsigned numPoints)
{
    std::vector<double> mult {15, -2, -3, 1.7, -1.5, 2.8, 1.1, -2.3, 1.9, 1.2, 5.7, -1.2};
    double val = 0.0;

    for (unsigned ii=0; ii < numPoints; ++ii)
    {
        val += std::pow(x, numPoints - 1 - ii) * mult.at(ii);
    }
    return val;
}

BOOST_AUTO_TEST_CASE(test1DQuadrature2Up)
{
    double answer = 9.0 + (1.0/3.0);
    unsigned min = 2;
    unsigned max = 12;
    Quadrature quadrature;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = quadrature.get1DGaussianQuadraturePoints(numPoints);
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
    Quadrature quadrature;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = quadrature.get1DGaussianQuadraturePoints(numPoints);
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
    Quadrature quadrature;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = quadrature.get1DGaussianQuadraturePoints(numPoints);
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
    Quadrature quadrature;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = quadrature.get1DGaussianQuadraturePoints(numPoints);
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
    Quadrature quadrature;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = quadrature.get1DGaussianQuadraturePoints(numPoints);
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
    Quadrature quadrature;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = quadrature.get1DGaussianQuadraturePoints(numPoints);
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
    Quadrature quadrature;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = quadrature.get1DGaussianQuadraturePoints(numPoints);
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
    Quadrature quadrature;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = quadrature.get1DGaussianQuadraturePoints(numPoints);
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
    Quadrature quadrature;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = quadrature.get1DGaussianQuadraturePoints(numPoints);
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
    Quadrature quadrature;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPoint1DList_type wps = quadrature.get1DGaussianQuadraturePoints(numPoints);
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

BOOST_AUTO_TEST_CASE(testTriangleGauss01)
{
    double answer = 1.0;
    unsigned min = 1;
    unsigned max = 11;

    Node n1 {-1.0, 0.0, 0.0};
    Node n2 { 1.0, 0.0, 0.0};
    Node n3 { 0.0, 1.0, 0.0};
    Triangle T {n1, n2, n3};
    Quadrature quadrature;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPointList_type wps = quadrature.getTriangleSimplexGaussianQuadraturePoints(numPoints);
        double integral = 0;
        for (unsigned ii=0; ii < wps.size(); ++ii)
        {
//            Node n = wps.at(ii).node;
            integral += wps.at(ii).weight * 1.0;
        }
//        std::cout << std::setprecision (15)<< numPoints<< ":  " << integral << std::endl;
        BOOST_CHECK_MESSAGE(std::fabs(integral - answer) < 1e-6, "Triangle Gauss failed.");
    }
}

BOOST_AUTO_TEST_CASE(testTriangleGauss02)
{
    double answer = 2.0;
    unsigned min = 1;
    unsigned max = 11;

    Node n1 {-1.0, 0.0, 0.0};
    Node n2 { 1.0, 0.0, 0.0};
    Node n3 { 0.0, 2.0, 0.0};
    Triangle T {n1, n2, n3};
    Quadrature quadrature;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPointList_type wps = quadrature.getTriangleGaussianQuadraturePoints(T, numPoints);
        double integral = 0;
        for (unsigned ii=0; ii < wps.size(); ++ii)
        {
            Node n = wps.at(ii).node;
            integral += wps.at(ii).weight * 1.0;
        }
//        std::cout << std::setprecision (15)<< numPoints<< ":  " << integral << std::endl;
        BOOST_CHECK_MESSAGE(std::fabs(integral - answer) < 1e-6, "Triangle Gauss failed.");
    }
}

BOOST_AUTO_TEST_CASE(testTriangleGaussSimplex_loop1)
{
    std::vector<double> answer {0, 75, -36 , -5.5, 10.276188, 8.456582};
    unsigned minmin = 1;
    unsigned minmax = 5;
    unsigned max = 11;
    double tol = 1e-6;

    Quadrature quadrature;
    for (unsigned min = minmin; min <= minmax; ++min)
    {
        Node n1 {-1.0, 0.0, 0.0};
        Node n2 { 1.0, 0.0, 0.0};
        Node n3 { 0.0, 1.0, 0.0};
        Triangle T {n1, n2, n3};
        for (unsigned numPoints = 2*min; numPoints <= max; ++numPoints)
        {
            Quadrature::WeightedPointList_type wps = quadrature.getTriangleSimplexGaussianQuadraturePoints(numPoints);
            double integral = 0;
            for (unsigned ii=0; ii < wps.size(); ++ii)
            {
                Node n = T.fromSimplex(wps.at(ii).node);
                integral += wps.at(ii).weight * funct1(n.x(), min)*funct2(n.y(), min);
            }
            std::string s;
            s.append("Triangle Gauss failed. [");
            s.append(std::to_string(integral));
            s.append(" != ");
            s.append(std::to_string(answer.at(min)));
            s.append("] (tolerance = ");
            s.append(std::to_string(tol));
            s.append(")");
            s.append(" (num points = ");
            s.append(std::to_string(numPoints));
            s.append(")");
            BOOST_CHECK_MESSAGE(std::fabs((integral - answer.at(min))/answer.at(min)) < tol, s);
        }
    }
}

BOOST_AUTO_TEST_CASE(testTriangleGauss_loop2)
{
    std::vector<double> answer {0, 150, -192, 34.333333333, 131.552381, -360.107815}; // Need to check values
    unsigned minmin = 1;
    unsigned minmax = 5;
    unsigned max = 11;
    double tol = 1e-6;

    Quadrature quadrature;
    for (unsigned min = minmin; min <= minmax; ++min)
    {
        Node n1 {-1.0, 0.0, 0.0};
        Node n2 { 1.0, 0.0, 0.0};
        Node n3 { 0.0, 2.0, 0.0};
        Triangle T {n1, n2, n3};
        for (unsigned numPoints = 2*min; numPoints <= max; ++numPoints)
        {
            Quadrature::WeightedPointList_type wps = quadrature.getTriangleGaussianQuadraturePoints(T, numPoints);
            double integral = 0;
            for (unsigned ii=0; ii < wps.size(); ++ii)
            {
                Node n = wps.at(ii).node;
                integral += wps.at(ii).weight * funct1(n.x(), min)*funct2(n.y(), min);
            }
            std::string s;
            s.append("Triangle Gauss failed. [");
            s.append(std::to_string(integral));
            s.append(" != ");
            s.append(std::to_string(answer.at(min)));
            s.append("] (tolerance = ");
            s.append(std::to_string(tol));
            s.append(")");
            s.append(" (num points = ");
            s.append(std::to_string(numPoints));
            s.append(")");
            BOOST_CHECK_MESSAGE(std::fabs((integral - answer.at(min))/answer.at(min)) < tol, s);
        }
    }
}

BOOST_AUTO_TEST_CASE(testRAR1S_2D_0)
{
    double answer = 2.0;
    unsigned min = 2; // Does not work for 1pt - not sure why not.
    unsigned max = 10;
    double offset = 0.0;
    double tol = 1e-2;

//    Node n1 {-1.0, 0.0, 0.0};
//    Node n2 { 0.0, 1.0, 0.0};
//    Node n3 { 0.0, 0.5, 0.0};
//    Node n1 {-1.0, 0.0, 0.0};
//    Node n2 { 1.0, 0.0, 0.0};
//    Node n3 { 0.0, 0.5, 0.0};
    Node n1 {-1.0, 0.0, 0.0};
    Node n2 { 1.0, 0.0, 0.0};
    Node n3 { 0.0, 2.0, 0.0};
    Triangle T {n1, n2, n3};
    Quadrature quadrature;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPointList_type wps = quadrature.RAR1S_2D(T, offset, numPoints);
        double integral = 0;
        for (unsigned ii=0; ii < wps.size(); ++ii)
        {
            Node n = wps.at(ii).node;
            integral += wps.at(ii).weight * 1.0;
        }
//        std::cout << std::setprecision (15)<< numPoints<< ":  " << integral << std::endl;
        std::string s;
        s.append("RAR1S 2D offset failed. [");
        s.append(std::to_string(integral));
        s.append(" != ");
        s.append(std::to_string(answer));
        s.append("] (tolerance = ");
        s.append(std::to_string(tol));
        s.append(")");
        s.append(" (num points = ");
        s.append(std::to_string(numPoints));
        s.append(")");
        BOOST_CHECK_MESSAGE(std::fabs((integral - answer)/answer) < tol, s);
    }
}

BOOST_AUTO_TEST_CASE(testRAR1S_2D_0_offset)
{
    double answer = 2.0;
    unsigned min = 2; // Does not work for 1pt - not sure why not.
    unsigned max = 10;
    double offset = 1.0;  // What is the effect of the offset suppose to be?
    double tol = 1e-2;

    Node n1 {-1.0, 0.0, 0.0};
    Node n2 { 1.0, 0.0, 0.0};
    Node n3 { 0.0, 2.0, 0.0};
    Triangle T {n1, n2, n3};
    Quadrature quadrature;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPointList_type wps = quadrature.RAR1S_2D(T, offset, numPoints);
        double integral = 0;
        for (unsigned ii=0; ii < wps.size(); ++ii)
        {
            Node n = wps.at(ii).node;
            integral += wps.at(ii).weight * 1.0;
        }
        std::string s;
        s.append("RAR1S 2D offset failed. [");
        s.append(std::to_string(integral));
        s.append(" != ");
        s.append(std::to_string(answer));
        s.append("] (tolerance = ");
        s.append(std::to_string(tol));
        s.append(")");
        s.append(" (num points = ");
        s.append(std::to_string(numPoints));
        s.append(")");
        BOOST_CHECK_MESSAGE(std::fabs((integral - answer)/answer) < tol, s);
    }
}

BOOST_AUTO_TEST_CASE(testRAR1S_2D_loop1)
{
    std::vector<double> answer {0, 75, -36 , -5.5, 10.276188, 8.459682};
//    std::vector<double> answer {0, 75, -33.0833333333333, 2.5, 12.4693499622071, 7.54417812022769};
    unsigned minmin = 1;
    unsigned minmax = 5;
    unsigned max = 11;
    double offset = 0.0;
    double tol = 1e-2;

    Quadrature quadrature;
    for (unsigned min = minmin; min <= minmax; ++min)
    {
        Node n1 {-1.0, 0.0, 0.0};
        Node n2 { 1.0, 0.0, 0.0};
        Node n3 { 0.0, 1.0, 0.0};
        Triangle T {n1, n2, n3};
        for (unsigned numPoints = 2*min; numPoints <= max; ++numPoints)
        {
            Quadrature::WeightedPointList_type wps = quadrature.RAR1S_2D(T, offset, numPoints);
            double integral = 0;
            for (unsigned ii=0; ii < wps.size(); ++ii)
            {
                Node n = wps.at(ii).node;
                integral += wps.at(ii).weight * funct1(n.x(), min)*funct2(n.y(), min);
            }
            std::string s;
            s.append("RAR1S 2D failed. [");
            s.append(std::to_string(integral));
            s.append(" != ");
            s.append(std::to_string(answer.at(min)));
            s.append("] (tolerance = ");
            s.append(std::to_string(tol));
            s.append(")");
            s.append(" (num points = ");
            s.append(std::to_string(numPoints));
            s.append(")");
            BOOST_CHECK_MESSAGE(std::fabs((integral - answer.at(min))/answer.at(min)) < tol, s);
        }
    }
}

BOOST_AUTO_TEST_CASE(testRAR1S_1)
{
    double answer = 1.0;
    unsigned min = 2;
    unsigned max = 8;

    Node n1 {-1.0, 0.0, 0.0};
    Node n2 { 1.0, 0.0, 0.0};
    Node n3 { 0.0, 1.0, 0.0};
    Node observationPoint { 9.0, 2.5, 10.0};
//    Node observationPoint { 0.0, 0.5, 0.0};
    Triangle T {n1, n2, n3};
    Quadrature quadrature;
    for (unsigned numPoints = min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPointList_type wps = quadrature.RAR1S(T, observationPoint, numPoints);
        double integral = 0;
        for (unsigned ii=0; ii < wps.size(); ++ii)
        {
            Node n = wps.at(ii).node;
            integral += wps.at(ii).weight * 1.0;
        }
//        std::cout << std::setprecision (15)<< numPoints<< ":  " << integral << std::endl;
        BOOST_CHECK_MESSAGE(std::fabs(integral - answer) < 1e-3, "RAR1S 3D failed.");
    }
}

BOOST_AUTO_TEST_CASE(testRAR1S_2)
{
    double answer = 75;
    unsigned min = 1;
    unsigned max = 8;

    Node n1 {-1.0, 0.0, 0.0};
    Node n2 { 1.0, 0.0, 0.0};
    Node n3 { 0.0, 1.0, 0.0};
//    Node observationPoint { 9.0, 2.5, 0.0};
    Node observationPoint { 0.0, 0.5, 1.0};
    Triangle T {n1, n2, n3};
    Quadrature quadrature;
    for (unsigned numPoints = 2*min; numPoints <= max; ++numPoints)
    {
        Quadrature::WeightedPointList_type wps = quadrature.RAR1S(T, observationPoint, numPoints);
        double integral = 0;
        for (unsigned ii=0; ii < wps.size(); ++ii)
        {
            Node n = wps.at(ii).node;
            integral += wps.at(ii).weight * funct1(n.x(), min)*funct2(n.y(), min);
//            integral += wps.at(ii).weight * 1.0;
        }
//        std::cout << std::setprecision (15)<< numPoints<< ":  " << integral << std::endl;
        BOOST_CHECK_MESSAGE(std::fabs((integral - answer)/answer) < 2e-2, "RAR1S 3D failed.");
    }
}

BOOST_AUTO_TEST_CASE(testRAR1S_loop1)
{
    std::vector<double> answer {0, 75, -36 , -5.5, 10.276188, 8.459682};
//    std::vector<double> answer {0, 75, -33.0833333333333, 2.5, 12.4693499622071, 7.54417812022769};
    unsigned minmin = 1;
    unsigned minmax = 5;
    unsigned max = 11;

    Quadrature quadrature;
    for (unsigned min = minmin; min <= minmax; ++min)
    {
        Node n1 {-1.0, 0.0, 1.0};
        Node n2 { 1.0, 0.0, 1.0};
        Node n3 { 0.0, 1.0, 1.0};
        Node observationPoint { 0.0, 0.5, -10.0};
//		Node observationPoint { 9.0, 2.5, 0.0};
        Triangle T {n1, n2, n3};
        for (unsigned numPoints = 2*min; numPoints <= max; ++numPoints)
        {
            Quadrature::WeightedPointList_type wps = quadrature.RAR1S(T, observationPoint, numPoints);
            double integral = 0;
            for (unsigned ii=0; ii < wps.size(); ++ii)
            {
                Node n = wps.at(ii).node;
                integral += wps.at(ii).weight * funct1(n.x(), min)*funct2(n.y(), min);
            }
            double tol = std::pow(1e-1, min);
            std::string s;
            s.append("RAR1S failed. [");
            s.append(std::to_string(integral));
            s.append(" != ");
            s.append(std::to_string(answer.at(min)));
            s.append("] (tolerance = ");
            s.append(std::to_string(tol));
            s.append(")");
            s.append(" (num points = ");
            s.append(std::to_string(numPoints));
            s.append(")");
            BOOST_CHECK_MESSAGE(std::fabs((integral - answer.at(min))/answer.at(min)) < tol, s);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

//EOF

