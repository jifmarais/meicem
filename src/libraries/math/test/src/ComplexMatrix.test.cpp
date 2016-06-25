#include "ComplexMatrix.hpp"
//#include <complex.h>
#include <cmath>
#include <math.h>
#include <complex>

//#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(ComplexMatrix_BasicTests)


BOOST_AUTO_TEST_CASE(constructor)
{
    ComplexMatrix m1 {1, 1, {1.0, 1.0}};
    std::complex<double> answer {1.0, 1.0};
    BOOST_CHECK_MESSAGE(m1.getRowCount() == 1, "Should be 1 row.");
    BOOST_CHECK_MESSAGE(m1.getColumnCount() == 1, "Should be 1 column.");
    BOOST_CHECK_MESSAGE(m1(0, 0) == answer, "Initialised to incorrect value.");

    ComplexMatrix m2 {100, 100};
    answer  = std::complex<double> {0.0 , 0.0};
    BOOST_CHECK_MESSAGE(m2.getRowCount() == 100, "Should be 100 rows.");
    BOOST_CHECK_MESSAGE(m2.getColumnCount() == 100, "Should be 100 columns.");
    BOOST_CHECK_MESSAGE(m2(0, 0) == answer, "Initialised to incorrect value.");
    BOOST_CHECK_MESSAGE(m2(99, 99) == answer, "Initialised to incorrect value.");
}

BOOST_AUTO_TEST_CASE(indexedAssignments)
{
    ComplexMatrix m1 {10, 10, {0.0, 0.0}};
    std::complex<double> zeros {0.0, 0.0};
    std::complex<double> ones {1.0, 1.0};

    BOOST_CHECK_MESSAGE(m1(0, 0) == zeros, "Should be zero.");
    m1(0, 0)  = ones;
    BOOST_CHECK_MESSAGE(m1(0, 0) == ones, "Should be ones.");

    BOOST_CHECK_MESSAGE(m1(4, 7) == zeros, "Should be zero.");
    m1(4, 7)  = ones;
    BOOST_CHECK_MESSAGE(m1(4, 7) == ones, "Should be ones.");

    BOOST_CHECK_MESSAGE(m1(9, 9) == zeros, "Should be zero.");
    m1(9, 9)  = ones;
    BOOST_CHECK_MESSAGE(m1(9, 9) == ones, "Should be ones.");

}

BOOST_AUTO_TEST_CASE(matrixAssignments)
{
    std::complex<double> zeros {0.0, 0.0};
    std::complex<double> ones {1.0, 1.0};
    std::complex<double> twos {2.0, 2.0};
    ComplexMatrix m1 {10, 10, zeros};
    ComplexMatrix m2 {5, 7, ones};
    ComplexMatrix m3 {11, 18, twos};

    m1 = m2;
    BOOST_CHECK_MESSAGE(m1.getRowCount() == 5, "Should be 5 rows.");
    BOOST_CHECK_MESSAGE(m1.getColumnCount() == 7, "Should be 7 columns.");
    // Checking that all values can be accessed and have correct value
    for (unsigned rr=0 ; rr < m1.getRowCount() ; ++rr )
    {
        for (unsigned cc=0 ; cc < m1.getRowCount() ; ++cc )
        {
            BOOST_CHECK_MESSAGE(m1(rr, cc) == ones, "Initialised to incorrect value.");
        }
    }

    m1 = m3;
    BOOST_CHECK_MESSAGE(m1.getRowCount() == 11, "Should be 11 rows.");
    BOOST_CHECK_MESSAGE(m1.getColumnCount() == 18, "Should be 18 columns.");
    // Checking that all values can be accessed and have correct value
    for (unsigned rr=0 ; rr < m1.getRowCount() ; ++rr )
    {
        for (unsigned cc=0 ; cc < m1.getRowCount() ; ++cc )
        {
            BOOST_CHECK_MESSAGE(m1(rr, cc) == twos, "Initialised to incorrect value.");
        }
    }

    m3(5, 6) = ones;
    BOOST_CHECK_MESSAGE(m1(5, 6) == twos, "Initialised to incorrect value.");
}

BOOST_AUTO_TEST_CASE(matrixComparison)
{
    std::complex<double> zeros {0.0, 0.0};
    std::complex<double> ones {1.0, 1.0};
    std::complex<double> twos {2.0, 2.0};
    ComplexMatrix m1a {10, 10, zeros};
    ComplexMatrix m2a {5, 7, ones};
    ComplexMatrix m3a {11, 18, twos};

    ComplexMatrix m1b {10, 10, zeros};
    ComplexMatrix m2b {5, 7, ones};
    ComplexMatrix m3b {11, 18, twos};

    ComplexMatrix& m3br {m3b};

    BOOST_CHECK_MESSAGE(m1a == m1b, "Should be equal.");
    BOOST_CHECK_MESSAGE(m2a == m2b, "Should be equal.");
    BOOST_CHECK_MESSAGE(m3a == m3b, "Should be equal.");

    BOOST_CHECK_MESSAGE(m3a == m3br, "Should be equal.");
    BOOST_CHECK_MESSAGE(m3b == m3br, "Should be equal.");

    BOOST_CHECK_MESSAGE(m1a != m2a, "Should not be equal.");
    BOOST_CHECK_MESSAGE(m1a != m3a, "Should not be equal.");
    BOOST_CHECK_MESSAGE(m2a != m3a, "Should not be equal.");

    m1b(2, 2) = 10;
    m2b(2, 2) = 10;
    m3b(2, 2) = 10;
    BOOST_CHECK_MESSAGE(m1a != m1b, "Should not be equal.");
    BOOST_CHECK_MESSAGE(m2a != m2b, "Should not be equal.");
    BOOST_CHECK_MESSAGE(m3a != m3b, "Should not be equal.");
    BOOST_CHECK_MESSAGE(m3b == m3br, "Should be equal.");
}

BOOST_AUTO_TEST_CASE(matrixAddition)
{
    std::complex<double> zeros {0.0, 0.0};
    std::complex<double> ones {1.0, 1.0};
    std::complex<double> twos {2.0, 2.0};
    std::complex<double> threes {3.0, 3.0};
    ComplexMatrix m0 {5, 7, zeros};
    ComplexMatrix m1 {5, 7, ones};
    ComplexMatrix m2 {5, 7, twos};
    ComplexMatrix m3 {5, 7, threes};

    BOOST_CHECK_MESSAGE(m1 == m1 + m0, "One plus zero should be one.");
    BOOST_CHECK_MESSAGE(m2 == m1 + m1, "One plus one should be two.");
    BOOST_CHECK_MESSAGE(m3 == m1 + m2, "One plus two should be three.");

    m0 += m1;
    m0 += m2;
    BOOST_CHECK_MESSAGE(m3 == m0, "Zero plus one plus two should be three.");

    m0 = m1 + std::complex<double> {2.0, 2.0};
    BOOST_CHECK_MESSAGE(m3 == m0, "One plus two should be three.");

    m0 = std::complex<double> {2.0, 2.0} + m1;
    BOOST_CHECK_MESSAGE(m3 == m0, "Two plus one should be three.");
}

BOOST_AUTO_TEST_CASE(matrixSubtraction)
{
    std::complex<double> zeros {0.0, 0.0};
    std::complex<double> ones {1.0, 1.0};
    std::complex<double> twos {2.0, 2.0};
    std::complex<double> threes {3.0, 3.0};
    ComplexMatrix m0 {5, 7, zeros};
    ComplexMatrix m1 {5, 7, ones};
    ComplexMatrix m2 {5, 7, twos};
    ComplexMatrix m3 {5, 7, threes};

    BOOST_CHECK_MESSAGE(m1 == m1 - m0, "One minus zero should be one.");
    BOOST_CHECK_MESSAGE(m0 == m1 - m1, "One minus one should be zero.");
    BOOST_CHECK_MESSAGE(m1 == m3 - m2, "Three minus two should be one.");

    m3 -= m1;
    m3 -= m2;
    BOOST_CHECK_MESSAGE(m3 == m0, "Three minus one minus two should be zero.");

    m3 = m2 - std::complex<double> {2.0, 2.0};
    BOOST_CHECK_MESSAGE(m3 == m0, "Two minus two should be zero.");

    m3 = std::complex<double> {2.0, 2.0} - m2;
    BOOST_CHECK_MESSAGE(m3 == m0, "Two minus two should be zero.");
}

BOOST_AUTO_TEST_CASE(matrixMultiplication)
{
    std::complex<double> zeros {0.0, 0.0};
    std::complex<double> ones {1.0, 0.0};
    std::complex<double> twos {2.0, 0.0};
    std::complex<double> fours {4.0, 0.0};
    std::complex<double> thirtytwo {32.0, 0.0};
    ComplexMatrix m0 {3, 2, zeros};
    ComplexMatrix m1 {2, 3, ones};
    ComplexMatrix m2 {3, 4, twos};
    ComplexMatrix m3 {4, 1, fours};

    BOOST_CHECK_MESSAGE((m0*m1).getRowCount() == 3, "Should have 3 rows.");
    BOOST_CHECK_MESSAGE((m0*m1).getColumnCount() == 3, "Should have 3 columns.");
    BOOST_CHECK_MESSAGE((m0*m1)(0 ,0) == zeros, "Values should all be 0.");

    BOOST_CHECK_MESSAGE((m1*m0).getRowCount() == 2, "Should have 2 rows.");
    BOOST_CHECK_MESSAGE((m1*m0).getColumnCount() == 2, "Should have 2 columns.");

    BOOST_CHECK_MESSAGE((m2*m3).getRowCount() == 3, "Should have 3 rows.");
    BOOST_CHECK_MESSAGE((m2*m3).getColumnCount() == 1, "Should have 1 column.");

    BOOST_CHECK_MESSAGE((m2*m3)(0,0) == thirtytwo, "The value should be 32.");
    BOOST_CHECK_MESSAGE((m2*m3)(1,0) == thirtytwo, "The value should be 32.");
    BOOST_CHECK_MESSAGE((m2*m3)(1,0) == thirtytwo, "The value should be 32.");
}

BOOST_AUTO_TEST_CASE(matrixElementWiseMultiplication)
{
    std::complex<double> zeros {0.0, 0.0};
    std::complex<double> ones {1.0, 0.0};
    std::complex<double> twos {2.0, 0.0};
    std::complex<double> fours {4.0, 0.0};
    ComplexMatrix m0 {5, 7, zeros};
    ComplexMatrix m1 {5, 7, ones};
    ComplexMatrix m2 {5, 7, twos};
    ComplexMatrix m3 {5, 7, fours};

    BOOST_CHECK_MESSAGE(m0 == m1.hadamard(m0), "One times zero is zero.");
    BOOST_CHECK_MESSAGE(m1 == m1.hadamard(m1), "One times one is one.");
    BOOST_CHECK_MESSAGE(m2 == m2.hadamard(m1), "Two times one is two.");
    BOOST_CHECK_MESSAGE(m3 == m2.hadamard(m2), "Two times two is four.");

    m0 = m2 * std::complex<double> {2.0, 0.0};
    BOOST_CHECK_MESSAGE(m3 == m0, "Two times two is four.");

    m0 = std::complex<double> {2.0, 0.0} * m2;
    BOOST_CHECK_MESSAGE(m3 == m0, "Two times two is four.");
}

BOOST_AUTO_TEST_CASE(matrixElementWiseDivision)
{
    std::complex<double> zeros {0.0, 0.0};
    std::complex<double> ones {1.0, 0.0};
    std::complex<double> twos {2.0, 0.0};
    std::complex<double> fours {4.0, 0.0};
    ComplexMatrix m0 {5, 7, zeros};
    ComplexMatrix m1 {5, 7, ones};
    ComplexMatrix m2 {5, 7, twos};
    ComplexMatrix m3 {5, 7, fours};

    BOOST_CHECK_MESSAGE(m0 == m0 / m1, "Zero divide one is zero.");
    BOOST_CHECK_MESSAGE(m1 == m2 / m2, "Two divide two is one.");
    BOOST_CHECK_MESSAGE(m2 == m3 / m2, "Four divide two is two.");
    BOOST_CHECK_MESSAGE(m3 == m3 / m1, "Four divide one is four.");

    m0 = m2 / std::complex<double> {2.0, 0.0};
    BOOST_CHECK_MESSAGE(m1 == m0, "Two divide two is one.");

    m0 = std::complex<double> {2.0, 0.0} / m2;
    BOOST_CHECK_MESSAGE(m1 == m0, "Two divide two is one.");
}

BOOST_AUTO_TEST_CASE(matrixSubMatrix)
{
    ComplexMatrix m1 {15, 17};
    for (unsigned rr = 0; rr < m1.getRowCount() ; ++rr)
    {
        for (unsigned cc = 0; cc < m1.getColumnCount() ; ++cc)
        {
            m1(rr, cc) = ComplexMatrix::complex_type {(double)rr, (double)cc};
        }
    }

    ComplexMatrix m0 {5, 7};
    for (unsigned rr = 0; rr < m0.getRowCount() ; ++rr)
    {
        for (unsigned cc = 0; cc < m0.getColumnCount() ; ++cc)
        {
            m0(rr, cc) = ComplexMatrix::complex_type {(double)rr, (double)cc};
        }
    }
    BOOST_CHECK_MESSAGE(m1.subMatrix(0,4, 0, 6) == m0, "Should be equal.");


    unsigned offset_row = 2;
    unsigned offset_col = 3;
    for (unsigned rr = 0; rr < m0.getRowCount() ; ++rr)
    {
        for (unsigned cc = 0; cc < m0.getColumnCount() ; ++cc)
        {
            m0(rr, cc) = ComplexMatrix::complex_type {(double)(rr+offset_row), (double)(cc+offset_col)};
        }
    }
    BOOST_CHECK_MESSAGE(m1.subMatrix(offset_row,4+offset_row, offset_col, 6+offset_col) == m0, "Should be equal.");
}

BOOST_AUTO_TEST_CASE(matrixMinor)
{
    ComplexMatrix ma0 {2, 2, {1.0, 10.0}};
    ComplexMatrix ma1 {1, 1};
    ma0(0, 0) = ComplexMatrix::complex_type {10.0, -15.0};
    ma1(0, 0) = ComplexMatrix::complex_type {10.0, -15.0};
    BOOST_CHECK_MESSAGE(ma0.minorMatrix(1, 1) == ma1, "Should be equal.");

    ComplexMatrix mb0 {3, 5};
    ComplexMatrix mb1 {2, 4};
    unsigned rrRemove = 1;
    unsigned ccRemove = 2;
    for (unsigned rr = 0; rr < mb0.getRowCount() ; ++rr)
    {
        for (unsigned cc = 0; cc < mb0.getColumnCount() ; ++cc)
        {
            mb0(rr, cc) = ComplexMatrix::complex_type {(double)rr, (double)cc};
            if ( (rr != rrRemove)  && (cc != ccRemove) )
            {
                mb1(rr - (rr > rrRemove), cc - (cc > ccRemove)) = ComplexMatrix::complex_type {(double)rr, (double)cc};
            }
        }
    }
    BOOST_CHECK_MESSAGE(mb0.minorMatrix(rrRemove, ccRemove) == mb1, "Should be equal.");
}


BOOST_AUTO_TEST_CASE(matrixTranspose)
{
    ComplexMatrix ma0 {1, 1};
    ComplexMatrix ma1 {1, 1};
    ma0(0, 0) = ComplexMatrix::complex_type {10.0, -15.0};
    ma1(0, 0) = ComplexMatrix::complex_type {10.0, -15.0};
    BOOST_CHECK_MESSAGE(ma1.transpose() == ma0, "Should be equal.");
    BOOST_CHECK_MESSAGE(ma1.transpose().transpose() == ma1, "Should be equal.");

    ComplexMatrix mb0 {3, 5};
    ComplexMatrix mb1 {5, 3};
    for (unsigned rr = 0; rr < mb0.getRowCount() ; ++rr)
    {
        for (unsigned cc = 0; cc < mb0.getColumnCount() ; ++cc)
        {
            mb0(rr, cc) = ComplexMatrix::complex_type {(double)rr, (double)cc};
            mb1(cc, rr) = ComplexMatrix::complex_type {(double)rr, (double)cc};
        }
    }
    BOOST_CHECK_MESSAGE(mb1.transpose() == mb0, "Should be equal.");
    BOOST_CHECK_MESSAGE(mb0.transpose() == mb1, "Should be equal.");
    BOOST_CHECK_MESSAGE(mb1.transpose().transpose() == mb1, "Should be equal.");
}

BOOST_AUTO_TEST_CASE(matrixConjugateTranspose)
{
    ComplexMatrix ma0 {1, 1};
    ComplexMatrix ma1 {1, 1};
    ma0(0, 0) = ComplexMatrix::complex_type {10.0, -15.0};
    ma1(0, 0) = ComplexMatrix::complex_type {10.0, 15.0};
    BOOST_CHECK_MESSAGE(ma1.conjugateTranspose() == ma0, "Should be equal.");
    BOOST_CHECK_MESSAGE(ma1.conjugateTranspose().conjugateTranspose() == ma1, "Should be equal.");

    ComplexMatrix mb0 {3, 5};
    ComplexMatrix mb1 {5, 3};
    for (unsigned rr = 0; rr < mb0.getRowCount() ; ++rr)
    {
        for (unsigned cc = 0; cc < mb0.getColumnCount() ; ++cc)
        {
            mb0(rr, cc) = ComplexMatrix::complex_type {(double)rr, (double)cc};
            mb1(cc, rr) = ComplexMatrix::complex_type {(double)rr, -(double)cc};
        }
    }
    BOOST_CHECK_MESSAGE(mb1.conjugateTranspose() == mb0, "Should be equal.");
    BOOST_CHECK_MESSAGE(mb0.conjugateTranspose() == mb1, "Should be equal.");
    BOOST_CHECK_MESSAGE(mb1.conjugateTranspose().conjugateTranspose() == mb1, "Should be equal.");
}

BOOST_AUTO_TEST_CASE(matrixDeterminant)
{
    ComplexMatrix ma0 {1, 1};
    ma0(0, 0) = ComplexMatrix::complex_type {10.0, -15.0};
    ComplexMatrix::complex_type answer {10.0, -15.0};
    BOOST_CHECK_MESSAGE(ma0.determinant() == answer, "Should be equal.");

    ComplexMatrix ma1 {2, 2};
    ma1(0, 0) = ComplexMatrix::complex_type {5.0, 3.0};
    ma1(0, 1) = ComplexMatrix::complex_type {2.0, -2.0};
    ma1(1, 0) = ComplexMatrix::complex_type {3.0, -1.0};
    ma1(1, 1) = ComplexMatrix::complex_type {4.0, 4.0};
    answer = ComplexMatrix::complex_type {4, 40};
    BOOST_CHECK_MESSAGE(ma1.determinant() == answer, "Should be equal.");

    ComplexMatrix ma2 {3, 3};
    ma2(0, 0) = ComplexMatrix::complex_type {5.0, 3.0};
    ma2(0, 1) = ComplexMatrix::complex_type {2.0, -2.0};
    ma2(0, 2) = ComplexMatrix::complex_type {4.0, -3.0};
    ma2(1, 0) = ComplexMatrix::complex_type {3.0, -1.0};
    ma2(1, 1) = ComplexMatrix::complex_type {4.0, 4.0};
    ma2(1, 2) = ComplexMatrix::complex_type {2.0, -5.0};
    ma2(2, 0) = ComplexMatrix::complex_type {6.0, -3.0};
    ma2(2, 1) = ComplexMatrix::complex_type {4.0, 2.0};
    ma2(2, 2) = ComplexMatrix::complex_type {5.0, 3.0};
    answer = ComplexMatrix::complex_type {-434.0, 198};
    BOOST_CHECK_MESSAGE(ma2.determinant() == answer, "Should be equal.");

}

double getMaxMagnitudeInMatrix(const ComplexMatrix m)
{
    double maxValue = 0.0;
    ComplexMatrix::complex_type cval;
    for ( unsigned rr = 0; rr < m.getRowCount() ; ++rr )
    {
        for ( unsigned cc = 0; cc < m.getColumnCount() ; ++cc )
        {
            cval = m(rr, cc);
            maxValue = std::max(maxValue, sqrt(cval.real()*cval.real() +
                                               cval.imag()*cval.imag()));
        }
    }

    return maxValue;
}

BOOST_AUTO_TEST_CASE(matrixInverse)
{
    ComplexMatrix ma0 {1, 1};
    ma0(0, 0) = ComplexMatrix::complex_type {10.0, -15.0};
    BOOST_CHECK_MESSAGE(ma0.inverse()(0, 0) == 1.0/ma0(0 ,0), "Should be equal.");

    ComplexMatrix ma1 {2, 2};
    ma1(0, 0) = ComplexMatrix::complex_type {5.0, 3.0};
    ma1(0, 1) = ComplexMatrix::complex_type {2.0, -2.0};
    ma1(1, 0) = ComplexMatrix::complex_type {3.0, -1.0};
    ma1(1, 1) = ComplexMatrix::complex_type {4.0, 4.0};
    BOOST_CHECK_MESSAGE(getMaxMagnitudeInMatrix(ma1.inverse()*ma1 - ComplexMatrix::identity(2)) <= 1e-6, "Should be equal.");

    ComplexMatrix ma2 {3, 3};
    ma2(0, 0) = ComplexMatrix::complex_type {5.0, 3.0};
    ma2(0, 1) = ComplexMatrix::complex_type {2.0, -2.0};
    ma2(0, 2) = ComplexMatrix::complex_type {4.0, -3.0};
    ma2(1, 0) = ComplexMatrix::complex_type {3.0, -1.0};
    ma2(1, 1) = ComplexMatrix::complex_type {4.0, 4.0};
    ma2(1, 2) = ComplexMatrix::complex_type {2.0, -5.0};
    ma2(2, 0) = ComplexMatrix::complex_type {6.0, -3.0};
    ma2(2, 1) = ComplexMatrix::complex_type {4.0, 2.0};
    ma2(2, 2) = ComplexMatrix::complex_type {5.0, 3.0};
    BOOST_CHECK_MESSAGE(getMaxMagnitudeInMatrix(ma2.inverse()*ma2 - ComplexMatrix::identity(3)) <= 1e-6, "Should be equal.");

    ComplexMatrix ma3 {4, 4};
    ma3(0, 0) = ComplexMatrix::complex_type {5.0, 3.0};
    ma3(0, 1) = ComplexMatrix::complex_type {2.0, -2.0};
    ma3(0, 2) = ComplexMatrix::complex_type {4.0, -3.0};
    ma3(0, 3) = ComplexMatrix::complex_type {-14.0, 13.0};

    ma3(1, 0) = ComplexMatrix::complex_type {3.0, -1.0};
    ma3(1, 1) = ComplexMatrix::complex_type {4.0, 4.0};
    ma3(1, 2) = ComplexMatrix::complex_type {2.0, -5.0};
    ma3(1, 3) = ComplexMatrix::complex_type {12.0, -15.0};

    ma3(2, 0) = ComplexMatrix::complex_type {6.0, -3.0};
    ma3(2, 1) = ComplexMatrix::complex_type {4.0, 2.0};
    ma3(2, 2) = ComplexMatrix::complex_type {5.0, 3.0};
    ma3(2, 3) = ComplexMatrix::complex_type {8.0, -1.0};

    ma3(3, 0) = ComplexMatrix::complex_type {6.0, -3.0};
    ma3(3, 1) = ComplexMatrix::complex_type {4.0, 2.0};
    ma3(3, 2) = ComplexMatrix::complex_type {5.0, 3.0};
    ma3(3, 3) = ComplexMatrix::complex_type {7.0, -3.0};
    BOOST_CHECK_MESSAGE(getMaxMagnitudeInMatrix(ma3.inverse()*ma3 - ComplexMatrix::identity(4)) <= 1e-6, "Should be equal.");

}

BOOST_AUTO_TEST_SUITE_END()

//EOF

