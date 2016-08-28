#include "NearFieldValue.hpp"

//#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(NearFieldValue_BasicTests)

BOOST_AUTO_TEST_CASE(constructor)
{
    NearFieldValue fval;
    std::complex<double> zero {0.0, 0.0};

    BOOST_CHECK_MESSAGE(fval.getX() == zero, "Default value should be zero.");
    BOOST_CHECK_MESSAGE(fval.getY() == zero, "Default value should be zero.");
    BOOST_CHECK_MESSAGE(fval.getZ() == zero, "Default value should be zero.");
}

BOOST_AUTO_TEST_CASE(getAndSetValue)
{
    std::complex<double> val1 {1.0, 2.0};
    std::complex<double> val2 {-1.0, 3.0};
    std::complex<double> val3 {4.0, -5.0};

    NearFieldValue fval;
    fval.set(val1, val2, val3);
    BOOST_CHECK_MESSAGE(fval.getX() == val1, "Value not correct.");
    BOOST_CHECK_MESSAGE(fval.getY() == val2, "Value not correct.");
    BOOST_CHECK_MESSAGE(fval.getZ() == val3, "Value not correct.");

    fval.setX(val2);
    fval.setY(val3);
    fval.setZ(val1);
    BOOST_CHECK_MESSAGE(fval.getX() == val2, "Value not correct.");
    BOOST_CHECK_MESSAGE(fval.getY() == val3, "Value not correct.");
    BOOST_CHECK_MESSAGE(fval.getZ() == val1, "Value not correct.");
}

BOOST_AUTO_TEST_CASE(assignments)
{
    std::complex<double> val1 {1.0, 2.0};
    std::complex<double> val2 {-1.0, 3.0};
    std::complex<double> val3 {10.0, 15.0};

    NearFieldValue fval1;
    NearFieldValue fval2;
    fval1.set(val1, val2, val3);
    fval2 = fval1;
    BOOST_CHECK_MESSAGE(fval2.getX() == val1, "Value not correct.");
    BOOST_CHECK_MESSAGE(fval2.getY() == val2, "Value not correct.");
    BOOST_CHECK_MESSAGE(fval2.getZ() == val3, "Value not correct.");
}

BOOST_AUTO_TEST_CASE(equality)
{
    std::complex<double> val1 {1.0, 2.0};
    std::complex<double> val2 {-1.0, 3.0};
    std::complex<double> val3 {10.0, 15.0};

    NearFieldValue fval1;
    NearFieldValue fval2;
    NearFieldValue fval3;
    fval1.set(val1, val2, val3);
    fval2.set(val1, val2, val3);
    fval3.set(val3, val1, val2);
    BOOST_CHECK_MESSAGE(fval1 == fval2, "Value not correct.");
    BOOST_CHECK_MESSAGE(fval1 != fval3, "Value not correct.");
}

BOOST_AUTO_TEST_CASE(addition)
{
    std::complex<double> val1 {1.0, 2.0};
    std::complex<double> val2 {-1.0, 3.0};
    std::complex<double> val3 {10.0, 15.0};

    std::complex<double> val11 {2.0, 4.0};
    std::complex<double> val12 {-2.0, 6.0};
    std::complex<double> val13 {20.0, 30.0};

    NearFieldValue fval1;
    NearFieldValue fval2;
    NearFieldValue fval3;
    NearFieldValue fval4;
    fval1.set(val1, val2, val3);
    fval2.set(val1, val2, val3);
    fval3.set(val11, val12, val13);

    fval4 = fval1 + fval2;
    BOOST_CHECK_MESSAGE(fval4 == fval3, "Value not correct.");

    fval4 = fval1;
    fval4 += fval1;
    BOOST_CHECK_MESSAGE(fval4 == fval3, "Value not correct.");
}

BOOST_AUTO_TEST_CASE(subtraction)
{
    std::complex<double> val1 {1.0, 2.0};
    std::complex<double> val2 {-1.0, 3.0};
    std::complex<double> val3 {10.0, 15.0};

    std::complex<double> val11 {2.0, 4.0};
    std::complex<double> val12 {-2.0, 6.0};
    std::complex<double> val13 {20.0, 30.0};

    NearFieldValue fval1;
    NearFieldValue fval2;
    NearFieldValue fval3;
    NearFieldValue fval4;
    fval1.set(val1, val2, val3);
    fval2.set(val1, val2, val3);
    fval3.set(val11, val12, val13);

    fval4 = fval3 - fval2;
    BOOST_CHECK_MESSAGE(fval4 == fval1, "Value not correct.");

    fval4 = fval3;
    fval4 -= fval1;
    BOOST_CHECK_MESSAGE(fval4 == fval2, "Value not correct.");
}

BOOST_AUTO_TEST_SUITE_END()

//EOF

