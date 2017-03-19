#include <armadillo>
#include "MoM.hpp"
#include "TriangleContainer.hpp"
#include "NodeContainer.hpp"
#include "reader_nastran.hpp"
#include "iostream"

//#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(MoM_BasicTests)

BOOST_AUTO_TEST_CASE(testTwoTriangles_value)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
    NastranReader reader;
    std::string baseTestFilesDirectory = "../src/libraries/mom_3d/test/models/test1_2triangles/";
    reader.setFile(baseTestFilesDirectory + "test1_2triangles.nas");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();

    double freq = 1e8;

    MoM MoMSetup {tContainer};
    MoMSetup.setFrequency(freq);
    PlaneWave pw;
    pw.setFrequency(freq);
    pw.setAngleOfIncidence(0.0, 0.0);
    MoMSetup.setNumberOfSourceIntegrationPoints(6);
    MoMSetup.setNumberOfTestIntegrationPoints(3);
    arma::cx_mat Zmatrix = MoMSetup.fillZmatrixTriangle();
    arma::cx_vec Vvector = MoMSetup.calculateRHS(pw);
    arma::cx_vec Vsolution((unsigned)Vvector.n_rows);
    Vsolution = arma::solve(Zmatrix, Vvector);

//    arma::cout.precision(11);
//    Vsolution.raw_print();

    arma::cx_vec ans(1);
    ans(0) = {-1.1331559565e-07, -0.00037046587346};

    BOOST_CHECK_MESSAGE(Vsolution.n_rows == 1, "Should have only one basis function.");
    BOOST_CHECK_MESSAGE(approx_equal(Vsolution, ans, "reldiff", 1e-3), "Answer is not correct.");
}

BOOST_AUTO_TEST_CASE(testTwoTriangles_rotation)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
    NastranReader reader;
    std::string baseTestFilesDirectory = "../src/libraries/mom_3d/test/models/test1_2triangles/";
    reader.setFile(baseTestFilesDirectory + "test1_2triangles.nas");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();
    double freq = 1e8;
    MoM MoMSetup {tContainer};
    MoMSetup.setFrequency(freq);
    PlaneWave pw;
    pw.setFrequency(freq);
    pw.setAngleOfIncidence(0.0, 0.0);
    MoMSetup.setNumberOfSourceIntegrationPoints(6);
    MoMSetup.setNumberOfTestIntegrationPoints(3);
    arma::cx_mat Zmatrix = MoMSetup.fillZmatrixTriangle();
    arma::cx_vec Vvector = MoMSetup.calculateRHS(pw);
    arma::cx_vec Vsolution1((unsigned)Vvector.n_rows);
    Vsolution1 = arma::solve(Zmatrix, Vvector);

    reader.setFile(baseTestFilesDirectory + "test1_2triangles_rotate.nas");
    TriangleContainer tContainer2(pContainer);
    reader.setTriangleContainer(&tContainer2);
    reader.importModel();
    MoM MoMSetup2 {tContainer2};
    MoMSetup2.setFrequency(freq);
    pw.setFrequency(freq);
    pw.setAngleOfIncidence(30.0, 20.0);
    Zmatrix = MoMSetup2.fillZmatrixTriangle();
    Vvector = MoMSetup2.calculateRHS(pw);
    arma::cx_vec Vsolution2((unsigned)Vvector.n_rows);
    Vsolution2 = arma::solve(Zmatrix, Vvector);

    BOOST_CHECK_MESSAGE(approx_equal(Vsolution1, Vsolution2, "reldiff", 1e-3), "Rotation should have no effect - everything rotated.");
}

BOOST_AUTO_TEST_CASE(testFourTriangles_value1)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
    NastranReader reader;
    std::string baseTestFilesDirectory = "../src/libraries/mom_3d/test/models/test2_4triangles/";
    reader.setFile(baseTestFilesDirectory + "test2_4triangles.nas");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();

    double freq = 1e8;

    MoM MoMSetup {tContainer};
    MoMSetup.setFrequency(freq);
    PlaneWave pw;
    pw.setFrequency(freq);
    pw.setAngleOfIncidence(0.0, 0.0);
    MoMSetup.setNumberOfSourceIntegrationPoints(6);
    MoMSetup.setNumberOfTestIntegrationPoints(3);
    arma::cx_mat Zmatrix = MoMSetup.fillZmatrixTriangle();
    arma::cx_vec Vvector = MoMSetup.calculateRHS(pw);
    arma::cx_vec Vsolution((unsigned)Vvector.n_rows);
    Vsolution = arma::solve(Zmatrix, Vvector);

    arma::cx_vec ans(4);
    ans(0) = {-1.2094701946e-07, -0.00038272310298};
    ans(1) = {-1.2094701946e-07, -0.00038272310298};
    ans(2) = {-1.2094701946e-07, -0.00038272310298};
    ans(3) = {-1.2094701946e-07, -0.00038272310298};

    BOOST_CHECK_MESSAGE(Vsolution.n_rows == 4, "Should have only one basis function.");
    BOOST_CHECK_MESSAGE(approx_equal(Vsolution, ans, "reldiff", 1e-3), "Answer is not correct.");

    pw.setAngleOfIncidence(30.0, 20.0);
    Zmatrix = MoMSetup.fillZmatrixTriangle();
    Vvector = MoMSetup.calculateRHS(pw);
    Vsolution = arma::solve(Zmatrix, Vvector);

    ans(0) = {-2.7053098135e-06, -0.00042470442645};
    ans(1) = {-2.6336804069e-06, -0.00019804292658};
    ans(2) = { 2.5085105306e-06, -0.00019804292709};
    ans(3) = { 2.4368821478e-06, -0.00042470442696};

    BOOST_CHECK_MESSAGE(approx_equal(Vsolution, ans, "reldiff", 1e-3), "Answer is not correct.");
}

BOOST_AUTO_TEST_CASE(testFourTriangles_rotation)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
    NastranReader reader;
    std::string baseTestFilesDirectory = "../src/libraries/mom_3d/test/models/test2_4triangles/";
    reader.setFile(baseTestFilesDirectory + "test2_4triangles.nas");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();
    double freq = 1e8;
    MoM MoMSetup {tContainer};
    MoMSetup.setFrequency(freq);
    PlaneWave pw;
    pw.setFrequency(freq);
    pw.setAngleOfIncidence(0.0, 0.0);
    MoMSetup.setNumberOfSourceIntegrationPoints(6);
    MoMSetup.setNumberOfTestIntegrationPoints(3);
    arma::cx_mat Zmatrix = MoMSetup.fillZmatrixTriangle();
    arma::cx_vec Vvector = MoMSetup.calculateRHS(pw);
    arma::cx_vec Vsolution1((unsigned)Vvector.n_rows);
    Vsolution1 = arma::solve(Zmatrix, Vvector);

    reader.setFile(baseTestFilesDirectory + "test2_4triangles_rotate.nas");
    TriangleContainer tContainer2(pContainer);
    reader.setTriangleContainer(&tContainer2);
    reader.importModel();
    MoM MoMSetup2 {tContainer2};
    MoMSetup2.setFrequency(freq);
    pw.setFrequency(freq);
    pw.setAngleOfIncidence(30.0, 20.0);
    Zmatrix = MoMSetup2.fillZmatrixTriangle();
    Vvector = MoMSetup2.calculateRHS(pw);
    arma::cx_vec Vsolution2((unsigned)Vvector.n_rows);
    Vsolution2 = arma::solve(Zmatrix, Vvector);

    BOOST_CHECK_MESSAGE(approx_equal(Vsolution1, Vsolution2, "reldiff", 1e-3), "Rotation should have no effect - everything rotated.");
}

BOOST_AUTO_TEST_CASE(testPlate_value1)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
    NastranReader reader;
    std::string baseTestFilesDirectory = "../src/libraries/mom_3d/test/models/test3_plate/";
    reader.setFile(baseTestFilesDirectory + "test3_plate.nas");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();

    double freq = 1e8;

    MoM MoMSetup {tContainer};
    MoMSetup.setFrequency(freq);
    PlaneWave pw;
    pw.setFrequency(freq);
    pw.setAngleOfIncidence(0.0, 0.0);
    MoMSetup.setNumberOfSourceIntegrationPoints(6);
    MoMSetup.setNumberOfTestIntegrationPoints(3);
    arma::cx_mat Zmatrix = MoMSetup.fillZmatrixTriangle();
    arma::cx_vec Vvector = MoMSetup.calculateRHS(pw);
    arma::cx_vec Vsolution((unsigned)Vvector.n_rows);
    Vsolution = arma::solve(Zmatrix, Vvector);

//    Vsolution.save(baseTestFilesDirectory + "ans");
    arma::cx_vec ans;
    ans.load(baseTestFilesDirectory + "ans");

    BOOST_CHECK_MESSAGE(Vsolution.n_rows == ans.n_rows, "Check for  number of basis functions. Assumes the mesh does not change.");
    BOOST_CHECK_MESSAGE(approx_equal(Vsolution, ans, "reldiff", 1e-3), "Answer is not correct.");
}

BOOST_AUTO_TEST_CASE(testPlate_rotation)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
    NastranReader reader;
    std::string baseTestFilesDirectory = "../src/libraries/mom_3d/test/models/test3_plate/";
    reader.setFile(baseTestFilesDirectory + "test3_plate.nas");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();
    double freq = 1e8;
    MoM MoMSetup {tContainer};
    MoMSetup.setFrequency(freq);
    PlaneWave pw;
    pw.setFrequency(freq);
    pw.setAngleOfIncidence(0.0, 0.0);
    MoMSetup.setNumberOfSourceIntegrationPoints(6);
    MoMSetup.setNumberOfTestIntegrationPoints(3);
    arma::cx_mat Zmatrix = MoMSetup.fillZmatrixTriangle();
    arma::cx_vec Vvector = MoMSetup.calculateRHS(pw);
    arma::cx_vec Vsolution1((unsigned)Vvector.n_rows);
    Vsolution1 = arma::solve(Zmatrix, Vvector);

    reader.setFile(baseTestFilesDirectory + "test3_plate_rotate.nas");
    TriangleContainer tContainer2(pContainer);
    reader.setTriangleContainer(&tContainer2);
    reader.importModel();
    MoM MoMSetup2 {tContainer2};
    MoMSetup2.setFrequency(freq);
    pw.setFrequency(freq);
    pw.setAngleOfIncidence(30.0, 20.0);
    Zmatrix = MoMSetup2.fillZmatrixTriangle();
    Vvector = MoMSetup2.calculateRHS(pw);
    arma::cx_vec Vsolution2((unsigned)Vvector.n_rows);
    Vsolution2 = arma::solve(Zmatrix, Vvector);

    BOOST_CHECK_MESSAGE(approx_equal(Vsolution1, Vsolution2, "reldiff", 1e-3), "Rotation should have no effect - everything rotated.");
}

BOOST_AUTO_TEST_CASE(testMixed_value1)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
    NastranReader reader;
    std::string baseTestFilesDirectory = "../src/libraries/mom_3d/test/models/test4_mixed/";
    reader.setFile(baseTestFilesDirectory + "test4_mixed.nas");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();

    double freq = 1e8;

    MoM MoMSetup {tContainer};
    MoMSetup.setFrequency(freq);
    PlaneWave pw;
    pw.setFrequency(freq);
    pw.setAngleOfIncidence(0.0, 0.0);
    MoMSetup.setNumberOfSourceIntegrationPoints(6);
    MoMSetup.setNumberOfTestIntegrationPoints(3);
    arma::cx_mat Zmatrix = MoMSetup.fillZmatrixTriangle();
    arma::cx_vec Vvector = MoMSetup.calculateRHS(pw);
    arma::cx_vec Vsolution((unsigned)Vvector.n_rows);
    Vsolution = arma::solve(Zmatrix, Vvector);

//    Vsolution.save(baseTestFilesDirectory + "ans");
    arma::cx_vec ans;
    ans.load(baseTestFilesDirectory + "ans");

    BOOST_CHECK_MESSAGE(Vsolution.n_rows == ans.n_rows, "Check for  number of basis functions. Assumes the mesh does not change.");
    BOOST_CHECK_MESSAGE(approx_equal(Vsolution, ans, "reldiff", 1e-3), "Answer is not correct.");
}

BOOST_AUTO_TEST_CASE(testMixed_rotation)
{
    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
    NastranReader reader;
    std::string baseTestFilesDirectory = "../src/libraries/mom_3d/test/models/test4_mixed/";
    reader.setFile(baseTestFilesDirectory + "test4_mixed.nas");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();
    double freq = 1e8;
    MoM MoMSetup {tContainer};
    MoMSetup.setFrequency(freq);
    PlaneWave pw;
    pw.setFrequency(freq);
    pw.setAngleOfIncidence(0.0, 0.0);
    MoMSetup.setNumberOfSourceIntegrationPoints(6);
    MoMSetup.setNumberOfTestIntegrationPoints(3);
    arma::cx_mat Zmatrix = MoMSetup.fillZmatrixTriangle();
    arma::cx_vec Vvector = MoMSetup.calculateRHS(pw);
    arma::cx_vec Vsolution1((unsigned)Vvector.n_rows);
    Vsolution1 = arma::solve(Zmatrix, Vvector);

    reader.setFile(baseTestFilesDirectory + "test4_mixed_rotate.nas");
    TriangleContainer tContainer2(pContainer);
    reader.setTriangleContainer(&tContainer2);
    reader.importModel();
    MoM MoMSetup2 {tContainer2};
    MoMSetup2.setFrequency(freq);
    pw.setFrequency(freq);
    pw.setAngleOfIncidence(30.0, 20.0);
    Zmatrix = MoMSetup2.fillZmatrixTriangle();
    Vvector = MoMSetup2.calculateRHS(pw);
    arma::cx_vec Vsolution2((unsigned)Vvector.n_rows);
    Vsolution2 = arma::solve(Zmatrix, Vvector);

    BOOST_CHECK_MESSAGE(approx_equal(Vsolution1, Vsolution2, "reldiff", 1e-3), "Rotation should have no effect - everything rotated.");
}

BOOST_AUTO_TEST_SUITE_END()

//EOF

