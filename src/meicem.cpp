#include "Node.hpp"

#include "NodeContainer.hpp"
#include "ComplexMatrix.hpp"
#include "TriangleContainer.hpp"
#include "reader_nastran.hpp"
//#include "reader_wilco_input.hpp"
#include "MoM.hpp"
//#include "math.h"
#include "assert.h"
#include <iostream>
#include "NearFieldContainer.hpp"
#include "PlaneWave.hpp"
#include <armadillo>
//#include <iterator>
//#include <complex>
//#include <cmath>

//#include <boost/program_options.hpp>
//namespace po = boost::program_options;


//int main( int argc, char *argv[] )
int main()
{
    std::cout << "Let the MEICEM mayhem begin..." << std::endl;

    NodeContainer pContainer;
    TriangleContainer tContainer(pContainer);
//    WilcoInputReader reader;
    NastranReader reader;
    //std::string baseTestFilesDirectory = "../src/libraries/reader_wilco_input/test/test_files/";
//    std::string baseTestFilesDirectory = "../models/meicem_simple_plate_test/";
    std::string baseTestFilesDirectory = "../models/basic_MoM_test_models/";
//    reader.setFile(baseTestFilesDirectory + "input_cononical_2basisfunction.txt");
//    reader.setFile(baseTestFilesDirectory + "simple_plate_test.nas");
//    reader.setFile(baseTestFilesDirectory + "test2_4triangles.nas");
    reader.setFile(baseTestFilesDirectory + "test4_mixed.nas");
//    reader.setFile(baseTestFilesDirectory + "mini_plate_test.nas");
//    reader.setFile(baseTestFilesDirectory + "mini_plate_test_rot.nas");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();

//    double freq = reader.getFrequency();
    double freq = 1e8;

    MoM MoMSetup {tContainer};
    MoMSetup.setFrequency(freq);
    PlaneWave pw;
    pw.setFrequency(freq);
    pw.setAngleOfIncidence(0.0, 0.0);
//    pw.setAngleOfIncidence(45.0, 0.0);
//    pw.setAngleOfIncidence(30, 20);
//    pw.setPolarisationAngle(0.0);
    arma::cx_mat Zmatrix = MoMSetup.fillZmatrixTriangle(6, 3);
    arma::cx_vec Vvector = MoMSetup.calculateRHS(6, pw);

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

