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
//    std::string baseTestFilesDirectory = "../src/libraries/reader_wilco_input/test/test_files/";
    std::string baseTestFilesDirectory = "../models/meicem_simple_plate_test/";
//    reader.setFile(baseTestFilesDirectory + "input_cononical_2basisfunction.txt");
    reader.setFile(baseTestFilesDirectory + "simple_plate_test.nas");
    reader.setTriangleContainer(&tContainer);
    reader.importModel();

//    double freq = reader.getFrequency();
    double freq = 1e8;

    MoM MoMSetup {tContainer};
    MoMSetup.setFrequency(freq);
    ComplexMatrix Zmatrix = MoMSetup.fillZmatrixTriangle();
    ComplexMatrix Vvector = MoMSetup.calculateRHS();

    std::cout << std::endl << "Zmatrix" << std::endl;
//    Zmatrix.print();

    std::cout << std::endl << "Vvector" << std::endl;
//    Vvector.print();

//    (Zmatrix.inverse()*Zmatrix).print();
    std::cout << std::endl << "Coefficients" << std::endl;
    ComplexMatrix Vsolution {(unsigned)Vvector.getRowCount(), 1};
    Vsolution = Zmatrix.inverse()*Vvector;
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
    PlaneWave pw;
    pw.setAngleOfIncidence(0.0, 0.0);
    pw.setPolarisationAngle(0.0);
    pw.setFrequency(1.5e8);
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

