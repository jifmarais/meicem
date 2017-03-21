#include "Node.hpp"
#include "NodeContainer.hpp"
#include "ComplexMatrix.hpp"
#include "TriangleContainer.hpp"
#include "reader_nastran.hpp"
#include "MoM.hpp"
#include "assert.h"
#include <iostream>
#include <chrono>
#include "NearFieldContainer.hpp"
#include "PlaneWave.hpp"
#include <armadillo>

int main()
{
    std::cout << "Let the MEICEM mayhem begin..." << std::endl;

    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
    NastranReader reader;
//    std::string baseTestFilesDirectory = "../models/meicem_simple_plate_test/";
//    std::string baseTestFilesDirectory = "../models/basic_MoM_test_models/";
//    std::string baseTestFilesDirectory = "../src/libraries/mom_3d/test/models/test3_plate/";

    std::string baseTestFilesDirectory = "../src/libraries/mom_3d/test/models/test4_mixed/";
    reader.setFile(baseTestFilesDirectory + "test4_mixed.nas");

//    reader.setFile(baseTestFilesDirectory + "input_cononical_2basisfunction.txt");
//    reader.setFile(baseTestFilesDirectory + "simple_plate_test.nas");
//    reader.setFile(baseTestFilesDirectory + "test2_4triangles.nas");
//    reader.setFile(baseTestFilesDirectory + "test1_2triangles.nas");
//    reader.setFile(baseTestFilesDirectory + "test3_plate.nas");
//    reader.setFile(baseTestFilesDirectory + "test3_plate_rotate.nas");
//    reader.setFile(baseTestFilesDirectory + "mini_plate_test.nas");
//    reader.setFile(baseTestFilesDirectory + "mini_plate_test_rot.nas");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();

    double freq = 1e8;

    MoM MoMSetup {tContainer};
    MoMSetup.setFrequency(freq);
    MoMSetup.setNumberOfSourceIntegrationPoints(6);
    MoMSetup.setNumberOfTestIntegrationPoints(3);
    PlaneWave pw;
    pw.setFrequency(freq);
    pw.setAngleOfIncidence(0.0, 0.0);
//    pw.setAngleOfIncidence(45.0, 0.0);
//    pw.setAngleOfIncidence(30, 20);
//    pw.setPolarisationAngle(0.0);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    arma::cx_mat Zmatrix = MoMSetup.fillZmatrixTriangle();
    std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds> (end - begin).count() <<std::endl;

    arma::cx_vec Vvector = MoMSetup.calculateRHS(pw);

    std::cout << std::endl << "Zimatrix" << std::endl;
//    Zmatrix.print();

    std::cout << std::endl << "Vvector" << std::endl;
//    Vvector.print();

    std::cout << std::endl << "Coefficients" << std::endl;
    arma::cx_vec Vsolution((unsigned)Vvector.n_rows);
    Vsolution = arma::solve(Zmatrix, Vvector);
//    Vsolution.print();

    // Calculate currents (and export them)
    MoMSetup.writeCurrentsToOS("test1", Vsolution);

    // Calculate near fields and export them
    NearFieldContainer nfContainer;
    NearFieldContainer nfContainer2;
    nfContainer.setFieldType(NearFieldContainer::ELECTICFIELD);
    nfContainer2.setFieldType(NearFieldContainer::MAGNETICFIELD);
    Node startCorner {-1.0, -1.0, -1.0};
    Node endCorner   {+1.0, +1.0, +1.0};
    nfContainer.setPoints(startCorner, endCorner, 11, 11, 51);
    nfContainer2.setPoints(startCorner, endCorner, 11, 11, 51);
    for ( unsigned pIndex = 0; pIndex < nfContainer.size(); ++pIndex )
    {
            pw.setFieldPoint(nfContainer.getPointAt(pIndex));
            nfContainer.setValueAt(pIndex, pw.getElectricField());

            pw.setFieldPoint(nfContainer2.getPointAt(pIndex));
            nfContainer2.setValueAt(pIndex, pw.getMagneticField());
    }

//    std::ofstream file ('test.efe');
    nfContainer.writeToEFEHFE("test1");
    nfContainer2.writeToEFEHFE("test1");


    return EXIT_SUCCESS;
}

