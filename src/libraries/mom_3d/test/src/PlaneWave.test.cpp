#include "PlaneWave.hpp"
#include <iostream>
//#include "EMcost.hpp"

//#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(PlaneWave_BasicTests)

BOOST_AUTO_TEST_CASE(basic_simple_points)
{
    PlaneWave pw;
    NearFieldValue nfval;
    const double c0 {299792456.2};
    std::complex<double> zero {0.0, 0.0};
    std::complex<double> one  {1.0, 0.0};
    std::complex<double> imgOne  {0.0, 1.0};

    pw.setAmplitude(1.0);
    pw.setFrequency(c0);
    pw.setFieldPoint(0.0, 0.0, 0.0);
    nfval.set(one, zero, zero);

    BOOST_CHECK_MESSAGE(pw.getElectricField().tolerantEqualTo(nfval), "The electric field at the origin should be 1.0 (x-directed).");

    // Testing 1 wavelength away
    pw.setFieldPoint(0.0, 0.0, 1.0);
    nfval.set(one, zero, zero);
    BOOST_CHECK_MESSAGE(pw.getElectricField().tolerantEqualTo(nfval), "The electric field one wavelength away should be 1.0 (x-directed).");

    // Testing 0.5 wavelength away
    pw.setFieldPoint(0.0, 0.0, 0.5);
    nfval.set(-one, zero, zero);
    BOOST_CHECK_MESSAGE(pw.getElectricField().tolerantEqualTo(nfval), "The electric field half a wavelength away should be 180 degrees out of phase (x-directed).");

    // Testing 0.25 wavelength away
    pw.setFieldPoint(0.0, 0.0, 0.25);
    nfval.set(-imgOne, zero, zero);
    BOOST_CHECK_MESSAGE(pw.getElectricField().tolerantEqualTo(nfval), "The electric field a quarter wavelength away should be imaginary (-x-directed).");
}

BOOST_AUTO_TEST_CASE(checkPolarisations)
{
    PlaneWave pw;
    NearFieldValue nfval;
    const double c0 {299792456.2};

    std::complex<double> zero {0.0, 0.0};
    std::complex<double> one  {1.0, 0.0};

    pw.setAmplitude(1.0);
    pw.setFrequency(c0);

    pw.setPolarisationAngle(45);
    pw.setFieldPoint(0.0, 0.0, 0.0);
    nfval.set(sqrt(0.5)*one, sqrt(0.5)*one, zero);
    BOOST_CHECK_MESSAGE(pw.getElectricField().tolerantEqualTo(nfval), "Electric field should be the same in x and y direction.");

    pw.setPolarisationAngle(90);
    pw.setFieldPoint(100.0, 12.20, 0.0);
    nfval.set(zero, one, zero);
    BOOST_CHECK_MESSAGE(pw.getElectricField().tolerantEqualTo(nfval), "Field should be in y direction in entire plane");

    pw.setPolarisationAngle(180);
    pw.setFieldPoint(0.0, 0.0, 0.0);
    nfval.set(-one, zero, zero);
    BOOST_CHECK_MESSAGE(pw.getElectricField().tolerantEqualTo(nfval), "The field should be in negative x direction.");
}

BOOST_AUTO_TEST_CASE(checkAngleOfIncidence)
{
    PlaneWave pw;
    NearFieldValue nfval;
    const double c0 {299792456.2};

    std::complex<double> zero {0.0, 0.0};
    std::complex<double> one  {1.0, 0.0};

    pw.setAmplitude(1.0);
    pw.setFrequency(c0);

    pw.setAngleOfIncidence(-45.0, 0.0);
    pw.setFieldPoint(0.0, 0.0, 0.0);
    nfval.set(sqrt(0.5)*one, zero, sqrt(0.5)*one);
    BOOST_CHECK_MESSAGE(pw.getElectricField().tolerantEqualTo(nfval), "Electric field should be the same in x and z direction.");

    pw.setAngleOfIncidence(-45.0, 45.0);
    pw.setFieldPoint(0.0, 0.0, 0.0);
    nfval.set(0.5*one, 0.5*one, sqrt(0.5)*one);
    BOOST_CHECK_MESSAGE(pw.getElectricField().tolerantEqualTo(nfval), "Electric field should be the same in x, y and z direction.");
}

BOOST_AUTO_TEST_CASE(magnetic_field_tests)
{
    PlaneWave pw;
    NearFieldValue nfval;
    const double c0 {299792456.2};
    const double eps0 {8.85418781761e-12};
    double Zf = 1 / (eps0*c0);
    std::complex<double> zero {0.0, 0.0};
    std::complex<double> one  {1.0, 0.0};
    std::complex<double> imgOne  {0.0, 1.0};

    pw.setAmplitude(1.0);
    pw.setFrequency(c0);
    pw.setFieldPoint(0.0, 0.0, 0.0);
    nfval.set(zero, -one/Zf, zero);

    BOOST_CHECK_MESSAGE(pw.getMagneticField().tolerantEqualTo(nfval), "The electric field at the origin should be 1.0 (x-directed).");
}

BOOST_AUTO_TEST_SUITE_END()

//EOF

